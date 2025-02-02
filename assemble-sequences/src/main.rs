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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct KmerCompact {
    bits: u128, // Bit-packed k-mer.
    k: u8, // k-mer length.
    is_start: bool,
    is_end: bool,
}

impl KmerCompact {
    /// Convert the compact k-mer back into a String, adding annotations if needed.
    fn to_string(&self) -> String {
        let mut s = String::with_capacity(self.k as usize);
        // Decode each nucleotide (2 bits per nucleotide).
        for i in 0..(self.k as usize) {
            let shift = 2 * (self.k as usize - i - 1);
            let code = (self.bits >> shift) & 0b11;
            let ch = match code {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => 'N',
            };
            s.push(ch);
        }
        // Add annotations if this k-mer is at a sequence boundary.
        if self.is_start {
            s = format!("^{}", s);
        }
        if self.is_end {
            s.push('$');
        }
        s
    }
}

/// Encodes a DNA k-mer string into a bit-packed u128.
/// Returns None if an invalid character is encountered.
fn encode_kmer(kmer: &str) -> Option<u128> {
    let mut bits: u128 = 0;
    for ch in kmer.chars() {
        bits <<= 2;
        bits |= match ch {
            'A' | 'a' => 0,
            'C' | 'c' => 1,
            'G' | 'g' => 2,
            'T' | 't' => 3,
            _ => return None,
        };
    }
    Some(bits)
}

/// Estimates the total number of k-mers that will be extracted from the sequences.
fn estimate_kmer_count(sequences: &[String], k: usize) -> usize {
    sequences
        .iter()
        .map(|s| if s.len() > k { s.len() - k + 1 } else { 0 })
        .sum()
}

fn build_local_graph(sequences: &[String], k: usize) -> HashSet<KmerCompact> {
    // Estimate the capacity to reduce reallocation.
    let estimated_capacity = estimate_kmer_count(sequences, k);
    let mut local_graph = HashMap::<u128, KmerCompact>::with_capacity(estimated_capacity);
    // Iterate over each sequence.
    for sequence in sequences {
        if sequence.len() <= k {
            continue; // Skip sequences that are too short.
        }
        // Slide a window of length k over the sequence.
        for i in 0..=sequence.len() - k {
            let kmer_slice = &sequence[i..i + k];
            if let Some(raw_bits) = encode_kmer(kmer_slice) {
                let is_start = i == 0;
                let is_end = i + k == sequence.len();
                let kmer_compact = KmerCompact {
                    bits: raw_bits,
                    k: k as u8,
                    is_start,
                    is_end,
                };
                // If a k-mer is already present but annotated, replace it with the unannotated version.
                match local_graph.get(&raw_bits) {
                    Some(existing) => {
                        if existing.is_start || existing.is_end {
                            let unannotated = KmerCompact {
                                bits: raw_bits,
                                k: k as u8,
                                is_start: false,
                                is_end: false,
                            };
                            local_graph.insert(raw_bits, unannotated);
                        }
                    }
                    None => {
                        local_graph.insert(raw_bits, kmer_compact);
                    }
                }
            }
        }
    }
    local_graph.into_values().collect()
}

fn build_global_graph(
    sequences: Vec<String>,
    k: usize,
    num_threads: usize,
) -> HashSet<KmerCompact> {
    let chunk_size = 10_000;
    let sequences_arc = Arc::new(sequences);
    let total_sequences = sequences_arc.len();
    // Create index ranges for chunks.
    let mut chunk_ranges = Vec::new();
    let mut start = 0;
    while start < total_sequences {
        let end = std::cmp::min(start + chunk_size, total_sequences);
        chunk_ranges.push((start, end));
        start = end;
    }
    // Global graph shared across threads.
    let global_graph = Arc::new(RwLock::new(HashMap::<u128, KmerCompact>::new()));
    let progress = Arc::new(AtomicUsize::new(0));
    // Process each chunk.
    for (chunk_start, chunk_end) in chunk_ranges {
        let chunk_len = chunk_end - chunk_start;
        // Partition the current chunk into sub-chunks.
        let mut subchunk_ranges = Vec::new();
        if chunk_len < num_threads {
            subchunk_ranges.push((chunk_start, chunk_end));
        } else {
            let sub_chunk_size = (chunk_len + num_threads - 1) / num_threads;
            let mut sub_start = chunk_start;
            while sub_start < chunk_end {
                let sub_end = std::cmp::min(sub_start + sub_chunk_size, chunk_end);
                subchunk_ranges.push((sub_start, sub_end));
                sub_start = sub_end;
            }
        }
        let mut handles = Vec::new();
        // Spawn threads for each sub-chunk.
        for (sub_start, sub_end) in subchunk_ranges {
            let sequences_clone = Arc::clone(&sequences_arc);
            let global_graph_clone = Arc::clone(&global_graph);
            let progress_clone = Arc::clone(&progress);
            let handle = thread::spawn(move || {
                // Create a slice for the sub-chunk.
                let sub_chunk = &sequences_clone[sub_start..sub_end];
                // Build the local graph for this sub-chunk.
                let local_graph = build_local_graph(sub_chunk, k);
                {
                    // Merge local k-mers into the global graph.
                    let mut global_lock = global_graph_clone.write().unwrap();
                    for kmer in local_graph {
                        let raw_bits = kmer.bits;
                        match global_lock.get(&raw_bits) {
                            Some(existing) => {
                                if existing.is_start || existing.is_end {
                                    // Replace with unannotated version.
                                    let unannotated = KmerCompact {
                                        bits: raw_bits,
                                        k: kmer.k,
                                        is_start: false,
                                        is_end: false,
                                    };
                                    global_lock.insert(raw_bits, unannotated);
                                }
                            }
                            None => {
                                global_lock.insert(raw_bits, kmer);
                            }
                        }
                    }
                }
                // Update and print progress.
                let sub_chunk_len = sub_end - sub_start;
                let current_progress =
                    progress_clone.fetch_add(sub_chunk_len, Ordering::Relaxed) + sub_chunk_len;
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
        // Wait for all threads in this chunk.
        for handle in handles {
            handle.join().unwrap();
        }
    }
    // Return the global graph as a HashSet.
    let global_graph_lock = global_graph.read().unwrap();
    global_graph_lock.values().cloned().collect()
}

fn construct_contigs(
    kmers: HashSet<KmerCompact>,
    min_length: usize,
    num_threads: usize,
) -> Vec<String> {
    // Convert the compact k-mers to Strings.
    let kmers_str: HashSet<String> = kmers.into_iter().map(|k| k.to_string()).collect();
    // Separate nodes into start, end, and intermediate.
    let (start_nodes, end_nodes, intermediate_nodes): (Vec<String>, Vec<String>, HashSet<String>) = {
        let mut start_nodes = Vec::new();
        let mut end_nodes = Vec::new();
        let mut intermediate_nodes = HashSet::new();
        for kmer in kmers_str {
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
    // Shared vector to collect contigs.
    let all_contigs = Arc::new(RwLock::new(Vec::new()));
    let progress = Arc::new(AtomicUsize::new(0));
    let total_start_nodes = start_nodes.len();
    // Build a suffix map: map from an overlap to nodes having that suffix.
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
    // Partition start nodes for parallel processing.
    let chunk_size = (total_start_nodes + num_threads - 1) / num_threads;
    let start_node_chunks: Vec<Vec<String>> = start_nodes
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    let mut handles = Vec::new();
    for chunk in start_node_chunks {
        let all_contigs_clone = Arc::clone(&all_contigs);
        let suffix_map = Arc::new(Mutex::new(suffix_map.clone()));
        let progress_clone = Arc::clone(&progress);
        let handle = thread::spawn(move || {
            let mut local_contigs = Vec::new();
            for start_node in chunk {
                let mut current_paths = vec![start_node.clone()];
                // Extend paths until no further extension is possible.
                while !current_paths.is_empty() {
                    let mut new_paths = Vec::new();
                    for path in current_paths {
                        let last_kmer = path.trim_start_matches('^');
                        // Calculate the overlap key.
                        let suffix_key = &last_kmer[path.len() - start_node.len() + 1..];
                        if let Some(matches) = suffix_map.lock().unwrap().get_mut(suffix_key) {
                            for matched_node in matches.drain(..) {
                                if matched_node.ends_with('$') {
                                    // End node: complete contig.
                                    let mut new_contig = path.clone();
                                    new_contig.push_str(&matched_node[matched_node.len() - 2..]);
                                    local_contigs.push(new_contig);
                                } else {
                                    // Intermediate: extend the current path.
                                    let mut extended_path = path.clone();
                                    extended_path.push_str(&matched_node[matched_node.len() - 1..]);
                                    new_paths.push(extended_path);
                                }
                            }
                        }
                    }
                    current_paths = new_paths;
                }
                // Update progress for contig construction.
                let current_progress = progress_clone.fetch_add(1, Ordering::Relaxed) + 1;
                if current_progress % 10 == 0 || current_progress == total_start_nodes {
                    println!(
                        "Progress: {}% ({}/{}) contigs processed",
                        (current_progress * 100) / total_start_nodes,
                        format_with_delimiters(current_progress),
                        format_with_delimiters(total_start_nodes)
                    );
                }
            }
            // Merge local contigs into the global vector.
            let mut global_contigs = all_contigs_clone.write().unwrap();
            global_contigs.extend(local_contigs);
        });
        handles.push(handle);
    }
    // Wait for all contig construction threads.
    for handle in handles {
        handle.join().unwrap();
    }
    // Filter and return only valid contigs.
    let global_contigs = all_contigs.read().unwrap();
    global_contigs
        .iter()
        .filter(|c| c.starts_with('^') && c.ends_with('$') && c.len() > min_length + 2)
        .map(|c| c.trim_matches(&['^', '$'][..]).to_string())
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
