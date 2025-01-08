// Standard Library Imports
use std::collections::{hash_map::DefaultHasher, HashMap, HashSet};
use std::env;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc, RwLock,
};
use std::thread;

// External Crate Imports
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

// File Type
#[derive(Debug, Clone, Copy)]
enum FileType {
    FASTA,
    FASTQ,
}

impl FileType {
    fn from_filename(filename: &str) -> io::Result<Self> {
        let fastq_exts = [".fastq.gz", ".fq.gz", ".fastq", ".fq"];
        let fasta_exts = [
            ".fasta",
            ".fa",
            ".faa",
            ".fna",
            ".fasta.gz",
            ".fa.gz",
            ".faa.gz",
            ".fna.gz",
        ];
        if fastq_exts.iter().any(|ext| filename.ends_with(ext)) {
            Ok(FileType::FASTQ)
        } else if fasta_exts.iter().any(|ext| filename.ends_with(ext)) {
            Ok(FileType::FASTA)
        } else {
            Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Unrecognized file extension for '{}'.", filename),
            ))
        }
    }
}

#[derive(Debug, Clone)]
struct Seq {
    id: String,
    sequence: String,
    quality: Option<String>,
}

fn hash_string(s: &str) -> u64 {
    let mut hasher = DefaultHasher::new();
    s.hash(&mut hasher);
    hasher.finish()
}

fn read_sequences_to_raw(file_path: &str, file_type: FileType) -> io::Result<Vec<String>> {
    let file = File::open(file_path)?;
    let reader: Box<dyn BufRead> = if file_path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut unique_seqs: HashMap<u64, Seq> = HashMap::new();
    match file_type {
        FileType::FASTQ => {
            let mut lines = reader.lines();
            while let (Some(header), Some(seq), Some(_), Some(qual)) =
                (lines.next(), lines.next(), lines.next(), lines.next())
            {
                let header = header?;
                let seq = seq?;
                let qual = qual?;
                let seq_hash = hash_string(&seq);
                unique_seqs.entry(seq_hash).or_insert(Seq {
                    id: header[1..].to_string(),
                    sequence: seq,
                    quality: Some(qual),
                });
            }
        }
        FileType::FASTA => {
            let mut lines = reader.lines().peekable();
            while let Some(header) = lines.next() {
                let header = header?;
                if !header.starts_with('>') {
                    continue; // Skip invalid headers
                }
                let seqid = header[1..].to_string();
                let mut seq = String::new();
                while let Some(Ok(line)) = lines.peek() {
                    if line.starts_with('>') {
                        break;
                    }
                    seq.push_str(lines.next().unwrap()?.trim());
                }
                let seq_hash = hash_string(&seq);
                unique_seqs.entry(seq_hash).or_insert(Seq {
                    id: seqid,
                    sequence: seq,
                    quality: None,
                });
            }
        }
    }
    Ok(unique_seqs
        .into_iter()
        .map(|(_, seq)| seq.sequence)
        .collect())
}

fn write_contigs_to_file(sequences: &[String], file_path: &str) -> io::Result<()> {
    let file = File::create(file_path)?;
    let mut writer: Box<dyn Write> = if file_path.ends_with(".gz") {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(file)
    };
    for (i, seq) in sequences.iter().enumerate() {
        writeln!(writer, ">{}\n{}", format!("{}|{}", i, seq.len()), seq)?;
    }
    Ok(())
}

fn build_local_graph(sequences: &[String], k: usize) -> Vec<String> {
    let mut graph = Vec::new();
    for sequence in sequences {
        for i in 0..=sequence.len() - k {
            let kmer = &sequence[i..i + k];
            let annotated_kmer = if i == 0 && i + k == sequence.len() {
                format!("^{}$", kmer) // Start and end
            } else if i == 0 {
                format!("^{}", kmer) // Start node
            } else if i + k == sequence.len() {
                format!("{}$", kmer) // End node
            } else {
                kmer.to_string() // Intermediate node
            };
            if !graph.contains(&annotated_kmer) {
                graph.push(annotated_kmer);
            }
        }
    }
    graph
}

fn build_global_graph(
    sequences: Vec<String>,
    k: usize,
    num_threads: usize,
) -> (Vec<String>, Vec<String>, Vec<String>) {
    let chunk_size = (sequences.len() + num_threads - 1) / num_threads;
    let sequence_chunks: Vec<_> = sequences
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    let global_graph = Arc::new(RwLock::new(HashSet::<String>::new()));
    let progress = Arc::new(AtomicUsize::new(0));
    let total_sequences = sequences.len();
    let mut handles = Vec::new();
    // Parallel processing
    for chunk in sequence_chunks {
        let global_graph_clone = Arc::clone(&global_graph);
        let progress = Arc::clone(&progress);
        let handle = thread::spawn(move || {
            let local_graph = build_local_graph(&chunk, k); // Generate local kmers
            let mut global_graph_lock = global_graph_clone.write().unwrap();
            for local_kmer in local_graph {
                global_graph_lock.insert(local_kmer); // Add directly
            }
            // Update progress
            let processed_sequences = chunk.len();
            let current_progress =
                progress.fetch_add(processed_sequences, Ordering::Relaxed) + processed_sequences;
            let approximate_progress = (current_progress * 100) / total_sequences;
            println!(
                "Progress: {}% ({}/{}) sequences processed",
                approximate_progress, current_progress, total_sequences
            );
        });
        handles.push(handle);
    }
    // Wait for all threads to finish
    for handle in handles {
        handle.join().unwrap();
    }
    // Separate the nodes
    let global_graph_lock = global_graph.read().unwrap();
    let mut start_nodes = Vec::new();
    let mut end_nodes = Vec::new();
    let mut intermediate_nodes = HashSet::new();
    // First pass: Collect intermediate nodes
    for kmer in global_graph_lock.iter() {
        if !kmer.starts_with('^') && !kmer.ends_with('$') {
            intermediate_nodes.insert(kmer.trim_matches(&['^', '$']).to_string());
        }
    }
    // Second pass: Filter start and end nodes based on intermediates
    for kmer in global_graph_lock.iter() {
        let raw_kmer = kmer.trim_matches(&['^', '$']).to_string();
        if kmer.starts_with('^') && !intermediate_nodes.contains(&raw_kmer) {
            start_nodes.push(kmer.clone());
        } else if kmer.ends_with('$') && !intermediate_nodes.contains(&raw_kmer) {
            end_nodes.push(kmer.clone());
        }
    }
    (
        start_nodes,
        intermediate_nodes.into_iter().collect(),
        end_nodes,
    )
}

fn construct_contigs(
    start_nodes: Vec<String>,
    intermediate_nodes: Vec<String>,
    end_nodes: Vec<String>,
    num_threads: usize,
) -> Vec<String> {
    // Shared data structures for threads
    let all_contigs = Arc::new(RwLock::new(Vec::new()));
    let progress = Arc::new(AtomicUsize::new(0)); // Progress tracker
    let total_start_nodes = start_nodes.len();
    let start_node_chunks: Vec<_> = start_nodes
        .chunks((total_start_nodes + num_threads - 1) / num_threads)
        .map(|chunk| chunk.to_vec())
        .collect();
    let end_nodes = Arc::new(end_nodes);
    let mut handles = Vec::new();
    for chunk in start_node_chunks {
        let all_contigs = Arc::clone(&all_contigs);
        let end_nodes = Arc::clone(&end_nodes);
        let local_intermediate_nodes = intermediate_nodes.clone();
        let progress = Arc::clone(&progress);
        let handle = thread::spawn(move || {
            let mut local_contigs = Vec::new();
            for (i, start_node) in chunk.iter().enumerate() {
                let mut current_paths = vec![start_node.clone()]; // Initialize with the current start node
                let mut remaining_nodes = local_intermediate_nodes.clone(); // Local Vec for thread
                while !current_paths.is_empty() {
                    let mut new_paths = Vec::new(); // For branching paths
                    for path in current_paths {
                        // Get the last kmer in the current path (strip '^')
                        let last_kmer = path.trim_start_matches('^');
                        // Find matching nodes based on overlap
                        let mut matches = Vec::new();
                        for node in remaining_nodes.iter().chain(end_nodes.iter()) {
                            let overlap_length = if end_nodes.contains(node) {
                                node.len() - 2 // End nodes have a shorter overlap
                            } else {
                                node.len() - 1
                            };
                            if last_kmer.ends_with(&node[..overlap_length]) {
                                matches.push(node.clone());
                            }
                        }
                        // Process matches
                        for matched_node in matches {
                            if end_nodes.contains(&matched_node) {
                                // Finalize the path if it's an end node
                                let mut new_contig = path.clone();
                                new_contig.push_str(&matched_node[matched_node.len() - 2..]);
                                local_contigs.push(new_contig);
                            } else {
                                // Extend the path with intermediate nodes
                                let mut extended_path = path.clone();
                                extended_path.push_str(&matched_node[matched_node.len() - 1..]);
                                new_paths.push(extended_path);
                                // Remove the used intermediate node
                                remaining_nodes.retain(|n| n != &matched_node);
                            }
                        }
                    }
                    // Update the current paths with the new paths
                    current_paths = new_paths;
                }
                // Update progress
                let current_progress = progress.fetch_add(1, Ordering::Relaxed) + 1;
                let approximate_progress = ((current_progress * 100) / total_start_nodes).min(100); // Prevent over 100%
                if current_progress % 10 == 0 || current_progress == total_start_nodes {
                    println!(
                        "Progress: {}% ({}/{}) contigs processed",
                        approximate_progress, current_progress, total_start_nodes
                    );
                }
            }
            // Write local results to the shared contigs
            let mut global_contigs = all_contigs.write().unwrap();
            global_contigs.extend(local_contigs);
        });

        handles.push(handle);
    }
    // Wait for threads to finish
    for handle in handles {
        handle.join().unwrap();
    }
    // Filter valid contigs (must start with "^" and end with "$")
    let global_contigs = all_contigs.read().unwrap();
    global_contigs
        .iter()
        .cloned()
        .filter(|c| c.starts_with('^') && c.ends_with('$'))
        .map(|c| c.trim_matches(&['^', '$']).to_string())
        .collect()
}

fn main() -> io::Result<()> {
    // Command-line argument handling
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "Usage: {} <input_file> <output_file> [k] [num_threads]",
            args[0]
        );
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Insufficient arguments",
        ));
    }
    let input_file = &args[1];
    let output_file = &args[2];
    if !Path::new(input_file).exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "Input file not found",
        ));
    }
    let k: usize = args
        .get(3)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(55);
    let available_threads = thread::available_parallelism().map_or(1, |n| n.get());
    let num_threads = std::cmp::min(
        available_threads,
        args.get(4)
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(4),
    );
    let raw_sequences = read_sequences_to_raw(input_file, FileType::from_filename(&input_file)?)?;
    let (start_nodes, intermediate_nodes, end_nodes) =
        build_global_graph(raw_sequences, k, num_threads);
    let contigs = construct_contigs(start_nodes, intermediate_nodes, end_nodes, num_threads);
    write_contigs_to_file(&contigs, output_file)?;
    Ok(())
}
