use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc, RwLock,
};
use std::thread;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

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

fn read_sequences(file_path: &str, file_type: FileType) -> io::Result<Vec<Seq>> {
    let file = File::open(file_path)?;
    let reader: Box<dyn BufRead> = if file_path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut seqs: Vec<Seq> = Vec::new();
    match file_type {
        FileType::FASTQ => {
            let mut lines = reader.lines();
            while let (Some(header), Some(seq), Some(_), Some(qual)) =
                (lines.next(), lines.next(), lines.next(), lines.next())
            {
                let header = header?;
                let seq = seq?;
                let qual = qual?;
                seqs.push(Seq {
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
                seqs.push(Seq {
                    id: seqid,
                    sequence: seq,
                    quality: None,
                });
            }
        }
    }
    Ok(seqs)
}

fn write_sequences_to_file(sequences: &[String], file_path: &str) -> io::Result<()> {
    let file = File::create(file_path)?;
    let mut writer: Box<dyn Write> = if file_path.ends_with(".gz") {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(file)
    };
    for (i, seq) in sequences.iter().enumerate() {
        writeln!(writer, ">{}\n{}", format!("Contig_{}", i), seq)?;
    }
    Ok(())
}

fn build_de_bruijn_graph(reads: &[String], k: usize) -> HashMap<String, Vec<String>> {
    let mut graph: HashMap<String, Vec<String>> = HashMap::new();
    for read in reads {
        for i in 0..=read.len() - k {
            let kmer = &read[i..i + k];
            let prefix = &kmer[..k - 1];
            let suffix = &kmer[1..];
            graph
                .entry(prefix.to_string())
                .and_modify(|suffixes| {
                    if !suffixes.contains(&suffix.to_string()) {
                        suffixes.push(suffix.to_string());
                    }
                })
                .or_insert_with(|| vec![suffix.to_string()]);
        }
    }
    graph
}

fn generate_contigs_from_graph(
    graph: &RwLock<HashMap<String, Vec<String>>>,
    k: usize,
    start_nodes: Vec<String>,
) -> Vec<String> {
    let mut contigs = Vec::new();
    let mut visited_edges: HashSet<(String, String)> = HashSet::new();
    // Helper to build a contig from a path
    let build_contig = |path: &[String]| -> String {
        path.iter()
            .skip(1) // Skip the first node since it is already included as the base
            .fold(path[0].clone(), |mut contig, node| {
                // For each subsequent node in the path, append its unique suffix (k-1 overlapping portion)
                contig.push_str(&node[k - 1..]);
                contig
            })
    };
    // Traverse the graph starting from each node in the chunk
    for start_node in start_nodes {
        let mut stack = vec![start_node.clone()];
        let mut path = Vec::new();
        while let Some(current_node) = stack.pop() {
            path.push(current_node.clone());
            let graph_read = graph.read().unwrap(); // Acquire read lock
            if let Some(suffixes) = graph_read.get(&current_node) {
                let mut sorted_suffixes = suffixes.clone();
                sorted_suffixes.sort();
                let mut extended = false;
                for suffix in sorted_suffixes {
                    let edge = (current_node.clone(), suffix.clone());
                    if !visited_edges.contains(&edge) {
                        visited_edges.insert(edge);
                        stack.push(suffix.clone());
                        extended = true;
                    }
                }
                if !extended {
                    contigs.push(build_contig(&path));
                    path.pop(); // Backtrack
                }
            } else {
                contigs.push(build_contig(&path));
                path.pop(); // Backtrack
            }
        }
    }
    contigs
}

fn process_reads_to_contigs(sequences: Vec<String>, k: usize, num_threads: usize) -> Vec<String> {
    let chunk_size = (sequences.len() + num_threads - 1) / (num_threads * 10);
    let sequence_chunks: Vec<_> = sequences
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    let graph = Arc::new(RwLock::new(HashMap::new()));
    let start_nodes = Arc::new(RwLock::new(HashSet::new()));
    let total_chunks = sequence_chunks.len();
    let progress = Arc::new(AtomicUsize::new(0));
    let mut handles = Vec::new();
    for chunk in sequence_chunks {
        let graph_clone = Arc::clone(&graph);
        let start_nodes_clone = Arc::clone(&start_nodes);
        let progress_clone = Arc::clone(&progress);
        let handle = thread::spawn(move || {
            let local_graph = build_de_bruijn_graph(&chunk, k);
            let mut local_start_nodes = HashSet::new();
            let mut global_graph = graph_clone.write().unwrap(); // Acquire write lock
            for (prefix, suffixes) in local_graph {
                global_graph
                    .entry(prefix.to_owned())
                    .and_modify(|existing_suffixes: &mut Vec<_>| {
                        for suffix in &suffixes {
                            if !existing_suffixes.contains(suffix) {
                                existing_suffixes.push(suffix.to_owned());
                            }
                        }
                    })
                    .or_insert_with(|| suffixes.clone());
                if !global_graph
                    .values()
                    .any(|suffixes| suffixes.contains(&prefix))
                {
                    local_start_nodes.insert(prefix);
                }
                for suffix in &suffixes {
                    if !global_graph.contains_key(suffix) {
                        local_start_nodes.insert(suffix.clone());
                    }
                }
            }
            let mut start_nodes_write = start_nodes_clone.write().unwrap();
            start_nodes_write.extend(local_start_nodes);
            let completed = progress_clone.fetch_add(1, Ordering::Relaxed) + 1;
            println!(
                "Graph construction progress: {}/{}",
                completed, total_chunks
            );
        });
        handles.push(handle);
    }
    for handle in handles {
        handle.join().unwrap();
    }
    let start_nodes_read = start_nodes.read().unwrap();
    let start_nodes: Vec<String> = start_nodes_read.clone().into_iter().collect();
    let chunk_size = if start_nodes.len() > 0 {
        (start_nodes.len() + num_threads - 1) / (num_threads * 10)
    } else {
        1 // Default to 1 if there are no start nodes
    };
    let contigs = Arc::new(RwLock::new(Vec::new()));
    let start_chunks: Vec<_> = start_nodes
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    let total_start_chunks = start_chunks.len();
    let progress = Arc::new(AtomicUsize::new(0));
    let mut handles = Vec::new();
    for chunk in start_chunks {
        let graph_clone = Arc::clone(&graph);
        let contigs_clone = Arc::clone(&contigs);
        let progress_clone = Arc::clone(&progress);
        let handle = thread::spawn(move || {
            let local_contigs = generate_contigs_from_graph(&graph_clone, k - 1, chunk);
            contigs_clone.write().unwrap().extend(local_contigs); // Acquire write lock
            let completed = progress_clone.fetch_add(1, Ordering::Relaxed) + 1;
            println!(
                "Contig generation progress: {}/{}",
                completed, total_start_chunks
            );
        });
        handles.push(handle);
    }
    for handle in handles {
        handle.join().unwrap();
    }
    let contigs_read = contigs.read().unwrap(); // Acquire read lock
    (*contigs_read).clone()
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
    let sequences = read_sequences(input_file, FileType::from_filename(&input_file)?)?;
    let unique_sequences: HashSet<String> = sequences.into_iter().map(|s| s.sequence).collect();
    let contigs =
        process_reads_to_contigs(Vec::from_iter(unique_sequences.into_iter()), k, num_threads);
    write_sequences_to_file(&contigs, output_file)?;
    Ok(())
}
