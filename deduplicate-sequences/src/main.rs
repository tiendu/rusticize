// Standard Library Imports
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
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
        match filename {
            // Detect the file type from a given filename
            f if f.ends_with(".fastq.gz")
                || f.ends_with(".fq.gz")
                || f.ends_with(".fastq")
                || f.ends_with(".fq") =>
            {
                Ok(FileType::FASTQ)
            }
            f if f.ends_with(".fasta")
                || f.ends_with(".fa")
                || f.ends_with(".faa")
                || f.ends_with(".fna") =>
            {
                Ok(FileType::FASTA)
            }
            _ => {
                eprintln!("Unrecognized file extension for '{}'.", filename);
                Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Unrecognized file extension",
                ))
            }
        }
    }
}

#[derive(Debug, Clone)]
struct Seq {
    id: String,
    sequence: String,
    quality: Option<String>,
}

impl Seq {
    fn length(&self) -> usize {
        self.sequence.len()
    }
}

fn hash_string(s: &str) -> u64 {
    let mut hasher = DefaultHasher::new();
    s.hash(&mut hasher);
    hasher.finish()
}

fn read_sequences(file_path: &str, file_type: FileType) -> io::Result<Vec<Seq>> {
    let file = File::open(file_path)?;
    let reader: Box<dyn BufRead> = if file_path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut count = 0;
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
                count += 1;
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
                count += 1;
            }
        }
    }
    let mut sorted_seqs: Vec<Seq> = unique_seqs.into_values().collect();
    sorted_seqs.sort_by(|a, b| a.length().cmp(&b.length()));
    println!("Start: {}", count);
    Ok(sorted_seqs)
}

fn write_sequences_to_file(sequences: &[Seq], file_path: &str) -> io::Result<()> {
    let file = File::create(file_path)?;
    let mut writer: Box<dyn Write> = if file_path.ends_with(".gz") {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(file)
    };
    for seq in sequences {
        writeln!(
            writer,
            "{}{}\n{}",
            if seq.quality.is_some() { "@" } else { ">" },
            seq.id,
            seq.sequence
        )?;
        if let Some(quality) = &seq.quality {
            writeln!(writer, "+\n{}", quality)?;
        }
    }
    Ok(())
}

fn deduplicate_sequences(
    sequences: Vec<Seq>,
    num_threads: usize,
    similarity: f64,
) -> io::Result<Vec<Seq>> {
    let length_change_indices: Vec<usize> = sequences
        .windows(2)
        .enumerate()
        .filter_map(|(i, window)| {
            if window[0].sequence.len() != window[1].sequence.len() {
                Some(i + 1)
            } else {
                None
            }
        })
        .collect();
    let total_sequences = sequences.len();
    let sequences = Arc::new(sequences); // Shared immutable data
    let mut handles = Vec::new();
    let progress = Arc::new(AtomicUsize::new(0)); // Thread-safe counter for progress
    for thread_id in 0..num_threads {
        let sequences = Arc::clone(&sequences); // Each thread gets its own Arc reference
        let length_change_indices = length_change_indices.clone();
        let progress = Arc::clone(&progress);
        let handle = thread::spawn(move || {
            let mut local_results = Vec::new();
            for i in (thread_id..total_sequences).step_by(num_threads) {
                let next_index = match length_change_indices.binary_search(&i) {
                    Ok(idx) => length_change_indices[idx],
                    Err(idx) => length_change_indices
                        .get(idx)
                        .copied()
                        .unwrap_or(total_sequences),
                };
                let current_seq = &sequences[i];
                if similarity == 100.0 {
                    if !sequences[next_index + 1..]
                        .iter()
                        .any(|s| s.sequence.contains(&current_seq.sequence))
                    {
                        local_results.push(i);
                    }
                } else {
                    if !sequences[i + 1..].iter().any(|s| {
                        let current_len = current_seq.sequence.len();
                        let candidate_len = s.sequence.len();
                        let max_mismatches = ((1.0 - similarity / 100.0)
                            * usize::min(current_len, candidate_len) as f64)
                            .ceil() as usize;
                        (0..=candidate_len.saturating_sub(current_len)).any(|start| {
                            let window = &s.sequence[start..start + current_len];
                            current_seq
                                .sequence
                                .chars()
                                .zip(window.chars())
                                .filter(|(a, b)| a != b) // Compare each char
                                .count()
                                <= max_mismatches
                        })
                    }) {
                        local_results.push(i);
                    }
                }
                let current_progress = progress.fetch_add(1, Ordering::SeqCst) + 1;
                if current_progress % 100 == 0 || current_progress == total_sequences {
                    println!("Progress: {}/{}", current_progress, total_sequences);
                }
            }
            local_results
        });
        handles.push(handle);
    }
    // Collect and combine results from threads
    let deduplicated: Vec<Seq> = handles
        .into_iter()
        .flat_map(|handle| handle.join().unwrap())
        .map(|idx| sequences[idx].clone())
        .collect();
    Ok(deduplicated)
}

fn main() -> io::Result<()> {
    // Command-line argument handling
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <input_file> <output_file> <similarity_threshold> [num_threads]",
            args[0]
        );
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Insufficient arguments",
        ));
    }
    let input_file = &args[1];
    let output_file = &args[2];
    let similarity_threshold: f64 = args[3]
        .parse::<f64>()
        .expect("Please provide a valid similarity value as a floating-point number.");
    let max_cpus = thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1); // Fallback to 1 if unavailable
    let num_threads: usize = args.get(6)
        .and_then(|s| s.parse().ok())
        .unwrap_or_else(|| {
            let threads = 4;
            eprintln!("Invalid number of threads or not provided. Using default value: {}", threads);
            std::cmp::min(threads, max_cpus)
        });
    if !Path::new(input_file).exists() {
        eprintln!("The input file '{}' does not exist.", input_file);
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "Input file not found",
        ));
    }
    if let Err(e) = File::create(output_file) {
        eprintln!("Error creating output file '{}': {}", output_file, e);
        return Err(e);
    }
    // Detect file type
    let file_type = FileType::from_filename(input_file)?;
    let sequences = read_sequences(input_file, file_type)?;
    if sequences.is_empty() {
        eprintln!("No sequences detected!");
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "No sequences detected",
        ));
    }
    let deduplicated_sequences =
        deduplicate_sequences(sequences, num_threads, similarity_threshold)?;
    println!("End: {}", deduplicated_sequences.len());
    write_sequences_to_file(&deduplicated_sequences, output_file)?;
    Ok(())
}
