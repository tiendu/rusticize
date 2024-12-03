// Standard Library Imports
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::env;
use std::fs::{File};
use std::hash::{Hash, Hasher};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{Arc};
use std::thread;

// External Crate Imports
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

// Sequence Struct Definition
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

// Hash Function for Sequence Deduplication
fn hash_string(s: &str) -> u64 {
    let mut hasher = DefaultHasher::new();
    s.hash(&mut hasher);
    hasher.finish()
}

// Function to Read Sequences from a File
fn read_sequences(file_path: &str, file_type: &str) -> io::Result<Vec<Seq>> {
    let file = File::open(file_path)?;
    let reader: Box<dyn BufRead> = if file_path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    
    let mut unique_seqs: HashMap<u64, Seq> = HashMap::new();
    let mut count = 0;

    match file_type {
        "FASTQ" => {
            let mut lines = reader.lines();
            while let (Some(header), Some(seq), Some(_), Some(qual)) = (
                lines.next(),
                lines.next(),
                lines.next(),
                lines.next(),
            ) {
                let header = header?;
                let seq = seq?;
                let qual = qual?;
                count += 1;
                let seq_hash = hash_string(&seq);
                unique_seqs.entry(seq_hash).or_insert(Seq {
                    id: header[1..].to_string(),
                    sequence: seq,
                    quality: Some(qual),
                });
            }
        }
        "FASTA" => {
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
                count += 1;
                let seq_hash = hash_string(&seq);
                unique_seqs.entry(seq_hash).or_insert(Seq {
                    id: seqid,
                    sequence: seq,
                    quality: None,
                });
            }
        }
        _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Unknown file type")),
    }

    println!("Start: {}", count);

    let mut sorted_seqs: Vec<Seq> = unique_seqs.into_values().collect();
    sorted_seqs.sort_by(|a, b| a.length().cmp(&b.length()));

    Ok(sorted_seqs)
}

// Function to Write Sequences to a File
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

// Function to Deduplicate Sequences
fn deduplicate_sequences(sequences: Vec<Seq>, num_threads: usize) -> io::Result<Vec<Seq>> {
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
    
    for thread_id in 0..num_threads {
        let sequences = Arc::clone(&sequences); // Each thread gets its own Arc reference
        let length_change_indices = length_change_indices.clone();
        
        let handle = thread::spawn(move || {
            let mut local_results = Vec::new();

            for i in (thread_id..total_sequences).step_by(num_threads) {
                let next_index = match length_change_indices.binary_search(&i) {
                    Ok(idx) => length_change_indices[idx],
                    Err(idx) => length_change_indices.get(idx).copied().unwrap_or(total_sequences),
                };

                let current_seq = &sequences[i];
                if !sequences[next_index..]
                    .iter()
                    .any(|s| s.sequence.contains(&current_seq.sequence))
                {
                    local_results.push(i);
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

    println!("End: {}", deduplicated.len());
    Ok(deduplicated)
}

// Main Function
fn main() -> io::Result<()> {
    // Command-line argument handling
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <input_file> <output_file> [num_threads]", args[0]);
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Insufficient arguments"));
    }

    let input_file = &args[1];
    let output_file = &args[2];
    let num_threads: usize = if args.len() > 3 {
        args[3].parse().unwrap_or_else(|_| {
            eprintln!("Invalid number of threads. Using default value: 4");
            4
        })
    } else {
        4
    };

    if !Path::new(input_file).exists() {
        eprintln!("The input file '{}' does not exist.", input_file);
        return Err(io::Error::new(io::ErrorKind::NotFound, "Input file not found"));
    }

    if let Err(e) = File::create(output_file) {
        eprintln!("Error creating output file '{}': {}", output_file, e);
        return Err(e);
    }

    let file_type = if input_file.ends_with(".fastq.gz") || input_file.ends_with(".fq.gz") || input_file.ends_with(".fastq") || input_file.ends_with(".fq") {
        "FASTQ"
    } else if input_file.ends_with(".fasta") || input_file.ends_with(".fa") || input_file.ends_with(".faa") || input_file.ends_with(".fna") {
        "FASTA"
    } else {
        eprintln!("Unrecognized file extension for '{}'.", input_file);
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Unrecognized file extension",
        ));
    };

    let sequences = read_sequences(input_file, file_type)?;
    if sequences.is_empty() {
        eprintln!("No sequences detected!");
        return Err(io::Error::new(io::ErrorKind::InvalidData, "No sequences detected"));
    }

    let adjusted_threads = std::cmp::max(1, std::cmp::min(num_threads, sequences.len() / 10));
    println!("Using {} threads for deduplication.", adjusted_threads);

    let deduplicated_sequences = deduplicate_sequences(sequences, adjusted_threads)?;
    write_sequences_to_file(&deduplicated_sequences, output_file)?;

    println!("Deduplication completed. Output written to '{}'.", output_file);

    Ok(())
}
