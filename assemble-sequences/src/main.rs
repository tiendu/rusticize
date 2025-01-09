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

fn build_local_graph(sequences: &[String], k: usize) -> HashSet<String> {
    let mut graph = HashSet::new();
    for sequence in sequences {
        if sequence.len() < k {
            continue; // Skip sequences shorter than k
        }
        for i in 0..=sequence.len() - k {
            let kmer = &sequence[i..i + k];
            let annotated_kmer = if i == 0 {
                format!("^{}", kmer) // Start node
            } else if i + k == sequence.len() {
                format!("{}$", kmer) // End node
            } else {
                kmer.to_string() // Intermediate node
            };
            graph.insert(annotated_kmer);
        }
    }
    graph
}

fn build_global_graph(sequences: &Vec<String>, k: usize, num_threads: usize) -> HashSet<String> {
    let chunk_size = (sequences.len() + num_threads - 1) / num_threads;
    let sequence_chunks: Vec<_> = sequences
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    let global_graph = Arc::new(RwLock::new(HashSet::<String>::new()));
    let progress = Arc::new(AtomicUsize::new(0));
    let total_sequences = sequences.len();
    let mut handles = Vec::new();
    for chunk in sequence_chunks {
        let global_graph_clone = Arc::clone(&global_graph);
        let progress = Arc::clone(&progress);
        let handle = thread::spawn(move || {
            let local_graph = build_local_graph(&chunk, k); // Generate local kmers
            let mut global_graph_lock = global_graph_clone.write().unwrap();
            for local_kmer in local_graph {
                global_graph_lock.insert(local_kmer); // Add directly
            }
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
    let global_graph_lock = global_graph.read().unwrap();
    global_graph_lock.clone()
}

fn separate_nodes(
    global_graph: &HashSet<String>,
) -> (HashSet<String>, HashSet<String>, HashSet<String>) {
    // Step 1: Split into start, intermediate, and end nodes
    let mut start_nodes = HashSet::new();
    let mut end_nodes = HashSet::new();
    let mut intermediate_nodes = HashSet::new();
    for kmer in global_graph {
        if kmer.starts_with('^') {
            start_nodes.insert(kmer.clone());
        } else if kmer.ends_with('$') {
            end_nodes.insert(kmer.clone());
        } else {
            intermediate_nodes.insert(kmer.clone());
        }
    }
    // Step 2: Validate start nodes
    let mut valid_start_nodes = HashSet::new();
    let mut intermediate_nodes_to_add = HashSet::new();
    for start_node in &start_nodes {
        let raw_kmer = start_node.trim_start_matches('^').to_string();
        if intermediate_nodes.contains(&raw_kmer)
            || end_nodes
                .iter()
                .any(|n| n.trim_end_matches('$') == raw_kmer)
        {
            // If the start node has an equivalent in intermediate_nodes or end_nodes, mark it as invalid
            // If it also exists in end_nodes, add it to intermediate_nodes
            if end_nodes
                .iter()
                .any(|n| n.trim_end_matches('$') == raw_kmer)
            {
                intermediate_nodes_to_add.insert(raw_kmer);
            }
        } else {
            // Otherwise, it's valid
            valid_start_nodes.insert(start_node.clone());
        }
    }
    // Add new intermediates
    intermediate_nodes.extend(intermediate_nodes_to_add.into_iter());
    // Step 3: Validate end nodes
    let valid_end_nodes: HashSet<String> = end_nodes
        .into_iter()
        .filter(|end_node| {
            let raw_kmer = end_node.trim_end_matches('$').to_string();
            !intermediate_nodes.contains(&raw_kmer)
        })
        .collect();
    (valid_start_nodes, intermediate_nodes, valid_end_nodes)
}

fn construct_contigs(
    start_nodes: &HashSet<String>,
    intermediate_nodes: &HashSet<String>,
    end_nodes: &HashSet<String>,
    threshold: usize,
    num_threads: usize,
) -> Vec<String> {
    let all_contigs = Arc::new(RwLock::new(Vec::new()));
    let progress = Arc::new(AtomicUsize::new(0));
    let total_start_nodes = start_nodes.len();
    let start_node_chunks: Vec<_> = start_nodes
        .iter()
        .cloned()
        .collect::<Vec<String>>()
        .chunks((total_start_nodes + num_threads - 1) / num_threads)
        .map(|chunk| chunk.to_vec())
        .collect();
    let end_nodes = Arc::new(end_nodes.clone());
    let mut handles = Vec::new();
    for chunk in start_node_chunks {
        let all_contigs = Arc::clone(&all_contigs);
        let end_nodes = Arc::clone(&end_nodes);
        let local_intermediate_nodes = intermediate_nodes.clone(); // Local copy of intermediate nodes
        let progress = Arc::clone(&progress);
        let handle = thread::spawn(move || {
            let mut local_contigs = Vec::new(); // Contigs generated by this thread
            for start_node in chunk {
                let mut current_paths = vec![start_node.clone()]; // Initialize paths with the start node
                let mut remaining_nodes = local_intermediate_nodes.clone(); // Local intermediate nodes
                                                                            // Iteratively extend paths until no new paths can be created
                while !current_paths.is_empty() {
                    let mut new_paths = Vec::new(); // Store new paths generated in this iteration
                    for path in current_paths {
                        let last_kmer = path.trim_start_matches('^'); // Get the suffix of the current path
                        let mut matches = Vec::new(); // Nodes that can extend this path
                                                      // Find nodes that overlap with the last kmer in the path
                        for node in remaining_nodes.iter().chain(end_nodes.iter()) {
                            let overlap_length = if end_nodes.contains(node) {
                                node.len() - 2 // End nodes have shorter overlap
                            } else {
                                node.len() - 1 // Intermediate nodes have full overlap
                            };
                            if last_kmer.ends_with(&node[..overlap_length]) {
                                matches.push(node.clone()); // Add matching node
                            }
                        }
                        // Extend the path for each matching node
                        for matched_node in matches {
                            if end_nodes.contains(&matched_node) {
                                // If the matched node is an end node, finalize the contig
                                let mut new_contig = path.clone();
                                new_contig.push_str(&matched_node[matched_node.len() - 2..]);
                                local_contigs.push(new_contig);
                            } else {
                                // If the matched node is intermediate, extend the path
                                let mut extended_path = path.clone();
                                extended_path.push_str(&matched_node[matched_node.len() - 1..]);
                                new_paths.push(extended_path);
                                // Remove the matched node from the remaining intermediate nodes
                                remaining_nodes.remove(&matched_node);
                            }
                        }
                    }
                    // Update current paths with newly generated paths
                    current_paths = new_paths;
                }
                let current_progress = progress.fetch_add(1, Ordering::Relaxed) + 1;
                if current_progress % 10 == 0 || current_progress == total_start_nodes {
                    println!(
                        "Progress: {}% ({}/{}) contigs processed",
                        (current_progress * 100) / total_start_nodes,
                        current_progress,
                        total_start_nodes
                    );
                }
            }
            let mut global_contigs = all_contigs.write().unwrap();
            global_contigs.extend(local_contigs);
        });
        handles.push(handle);
    }
    for handle in handles {
        handle.join().unwrap();
    }
    // Filter and clean up contigs: must start with "^" and end with "$"
    let global_contigs = all_contigs.read().unwrap();
    global_contigs
        .iter()
        .cloned()
        .filter(|c| c.starts_with('^') && c.ends_with('$') && c.len() > threshold + 2) // Keep only valid contigs, accounted for two extra chars '^', '$'
        .map(|c| c.trim_matches(&['^', '$']).to_string()) // Remove annotations from contigs
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
    let mut min_length = usize::MAX; // Initialize with the maximum possible usize value
    let _: Vec<_> = raw_sequences
        .clone()
        .into_iter()
        .map(|seq| {
            let len = seq.len();
            if len < min_length {
                min_length = len; // Update the minimum length
            }
        })
        .collect();
    let graph = build_global_graph(&raw_sequences, k, num_threads);
    let (start_nodes, intermediate_nodes, end_nodes) = separate_nodes(&graph);
    let contigs = construct_contigs(
        &start_nodes,
        &intermediate_nodes,
        &end_nodes,
        min_length,
        num_threads,
    );
    write_contigs_to_file(&contigs, output_file)?;
    Ok(())
}
