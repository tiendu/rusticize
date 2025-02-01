// Standard Library Imports
use std::collections::{hash_map::DefaultHasher, HashMap, HashSet};
use std::env;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc, Mutex, RwLock,
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

fn write_contigs_to_file(sequences: Vec<String>, file_path: &str) -> io::Result<()> {
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
    let mut local_graph = HashMap::<u64, String>::new(); // Map hashes to k-mers
    for sequence in sequences {
        if sequence.len() <= k {
            continue; // Skip sequences shorter than k
        }
        for i in 0..=sequence.len() - k {
            let kmer = &sequence[i..i + k];
            let hash = hash_string(kmer); // Compute the hash of the raw k-mer
            let annotated_kmer = if i == 0 {
                format!("^{}", kmer) // Start node
            } else if i + k == sequence.len() {
                format!("{}$", kmer) // End node
            } else {
                kmer.to_string() // Intermediate node
            };
            // Resolve conflicts locally
            match local_graph.get(&hash) {
                Some(existing_kmer) => {
                    if existing_kmer.contains('^') || existing_kmer.contains('$') {
                        // Prefer unannotated k-mer if present
                        local_graph.insert(hash, kmer.to_string());
                    }
                }
                None => {
                    // No existing k-mer; insert the annotated one
                    local_graph.insert(hash, annotated_kmer);
                }
            }
        }
    }
    // Extract the values from the map as a HashSet
    local_graph.into_values().collect::<HashSet<String>>()
}

fn build_global_graph(sequences: Vec<String>, k: usize, num_threads: usize) -> HashSet<String> {
    let chunk_size = 10_000;
    let sequence_chunks: Vec<_> = sequences
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    let global_graph = Arc::new(RwLock::new(HashMap::<u64, String>::new()));
    let progress = Arc::new(AtomicUsize::new(0));
    let total_sequences = sequences.len();
    for chunk in sequence_chunks {
        let chunk_len = chunk.len();
        let sub_chunks: Vec<_> = if chunk_len < num_threads {
            vec![chunk.to_vec()] // Keep the original chunk as is
        } else {
            chunk
                .chunks((chunk_len + num_threads - 1) / num_threads)
                .map(|sub_chunk| sub_chunk.to_vec())
                .collect()
        };
        let mut handles = Vec::new();
        for sub_chunk in sub_chunks {
            let global_graph_clone = Arc::clone(&global_graph);
            let progress = Arc::clone(&progress);
            let handle = thread::spawn(move || {
                let local_graph = build_local_graph(&sub_chunk, k); // Generate local k-mers
                let mut global_graph_lock = global_graph_clone.write().unwrap();
                for local_kmer in local_graph {
                    let raw_local_kmer = local_kmer.trim_matches(&['^', '$']);
                    let hash = hash_string(raw_local_kmer);
                    match global_graph_lock.get(&hash) {
                        Some(existing_kmer) => {
                            // Resolve conflict if annotations differ
                            if existing_kmer.contains('^') || existing_kmer.contains('$') {
                                global_graph_lock.insert(hash, raw_local_kmer.to_string());
                            }
                        }
                        None => {
                            // Insert the k-mer directly to the global graph
                            global_graph_lock.insert(hash, local_kmer);
                        }
                    }
                }
                // Update progress
                let processed_sequences = sub_chunk.len();
                let current_progress = progress.fetch_add(processed_sequences, Ordering::Relaxed)
                    + processed_sequences;
                let approximate_progress = (current_progress * 100) / total_sequences;
                println!(
                    "Progress: {}% ({}/{}) sequences processed",
                    approximate_progress,
                    format_with_delimiters(current_progress),
                    format_with_delimiters(total_sequences)
                );
            });
            handles.push(handle);
        }
        for handle in handles {
            handle.join().unwrap();
        }
    }
    // Extract the values from the map as a HashSet
    let global_graph_lock = global_graph.read().unwrap();
    global_graph_lock
        .values()
        .cloned()
        .collect::<HashSet<String>>()
}

fn construct_contigs(kmers: HashSet<String>, min_length: usize, num_threads: usize) -> Vec<String> {
    let (start_nodes, end_nodes, intermediate_nodes): (Vec<_>, Vec<_>, HashSet<_>) = {
        let mut start_nodes = Vec::new();
        let mut end_nodes = Vec::new();
        let mut intermediate_nodes = HashSet::new();
        for kmer in kmers {
            if kmer.starts_with('^') {
                start_nodes.push(kmer.clone());
            } else if kmer.ends_with('$') {
                end_nodes.push(kmer.clone());
            } else {
                intermediate_nodes.insert(kmer.clone());
            }
        }
        (start_nodes, end_nodes, intermediate_nodes)
    };
    let all_contigs = Arc::new(RwLock::new(Vec::new()));
    let progress = Arc::new(AtomicUsize::new(0));
    let total_start_nodes = start_nodes.len();
    let suffix_map: HashMap<String, Vec<String>> = {
        let mut map = HashMap::<String, Vec<String>>::new();
        for node in intermediate_nodes.iter().chain(end_nodes.iter()) {
            let overlap_length = if node.ends_with('$') {
                node.len() - 2
            } else {
                node.len() - 1
            };
            let suffix = node[..overlap_length].to_string();
            map.entry(suffix).or_default().push(node.clone());
        }
        map
    };
    let chunk_size = (total_start_nodes + num_threads - 1) / num_threads;
    let start_node_chunks: Vec<_> = start_nodes
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    let mut handles = Vec::new();
    for chunk in start_node_chunks {
        let all_contigs = Arc::clone(&all_contigs);
        let suffix_map = Arc::new(Mutex::new(suffix_map.clone()));
        let progress = Arc::clone(&progress);
        let handle = thread::spawn(move || {
            let mut local_contigs = Vec::new();
            for start_node in chunk {
                let mut current_paths = vec![start_node.clone()];
                while !current_paths.is_empty() {
                    let mut new_paths = Vec::new();
                    for path in current_paths {
                        let last_kmer = path.trim_start_matches('^');
                        let suffix_key = &last_kmer[path.len() - start_node.len() + 1..];
                        if let Some(matches) = suffix_map.lock().unwrap().get_mut(suffix_key) {
                            for matched_node in matches.drain(..) {
                                if matched_node.ends_with('$') {
                                    let mut new_contig = path.clone();
                                    new_contig.push_str(&matched_node[matched_node.len() - 2..]);
                                    local_contigs.push(new_contig);
                                } else {
                                    let mut extended_path = path.clone();
                                    extended_path.push_str(&matched_node[matched_node.len() - 1..]);
                                    new_paths.push(extended_path);
                                }
                            }
                        }
                    }
                    current_paths = new_paths;
                }
                let current_progress = progress.fetch_add(1, Ordering::Relaxed) + 1;
                if current_progress % 10 == 0 || current_progress == total_start_nodes {
                    println!(
                        "Progress: {}% ({}/{}) contigs processed",
                        (current_progress * 100) / total_start_nodes,
                        format_with_delimiters(current_progress),
                        format_with_delimiters(total_start_nodes)
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
    let global_contigs = all_contigs.read().unwrap();
    global_contigs
        .iter()
        .filter(|c| c.starts_with('^') && c.ends_with('$') && c.len() > min_length + 2)
        .map(|c| c.trim_matches(&['^', '$']).to_string())
        .collect()
}

fn format_with_delimiters(number: usize) -> String {
    let num_str = number.to_string();
    let mut formatted = String::new();
    let mut count = 0;
    for ch in num_str.chars().rev() {
        if count > 0 && count % 3 == 0 {
            formatted.push(',');
        }
        formatted.push(ch);
        count += 1;
    }
    formatted.chars().rev().collect()
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
    let min_length = raw_sequences.iter().map(|seq| seq.len()).min().unwrap();
    let graph = build_global_graph(raw_sequences, k, num_threads);
    let contigs = construct_contigs(graph, min_length, num_threads);
    write_contigs_to_file(contigs, output_file)?;
    Ok(())
}
