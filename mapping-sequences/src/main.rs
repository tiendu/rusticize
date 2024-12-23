// Standard Library Imports
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;

// External Crate Imports
use flate2::read::GzDecoder;

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

fn dna_to_bits(dna: &str) -> Vec<u8> {
    let nucl_to_bits = [
        ('A', 0b0001),
        ('T', 0b0010),
        ('G', 0b0100),
        ('C', 0b1000),
        ('W', 0b0011),
        ('R', 0b0101),
        ('Y', 0b1010),
        ('S', 0b1100),
        ('K', 0b0110),
        ('M', 0b1001),
        ('B', 0b1110),
        ('D', 0b0111),
        ('H', 0b1011),
        ('V', 0b1101),
        ('N', 0b1111),
    ]
    .iter()
    .cloned()
    .collect::<HashMap<_, _>>();

    dna.chars()
        .map(|c| {
            *nucl_to_bits.get(&c).unwrap_or_else(|| {
                eprintln!("Invalid base '{}'. Setting to N.", c);
                &0b1111 // Default to 'N'
            })
        })
        .collect()
}

// Reverse complement of a DNA sequence
fn reverse_complement(dna: &str) -> String {
    let complement = [
        ('A', 'T'),
        ('T', 'A'),
        ('G', 'C'),
        ('C', 'G'),
        ('R', 'Y'),
        ('Y', 'R'),
        ('S', 'S'),
        ('W', 'W'),
        ('K', 'M'),
        ('M', 'K'),
        ('B', 'V'),
        ('D', 'H'),
        ('H', 'D'),
        ('V', 'B'),
        ('N', 'N'),
    ]
    .iter()
    .cloned()
    .collect::<HashMap<_, _>>();

    dna.chars()
        .rev()
        .map(|c| {
            *complement.get(&c).unwrap_or_else(|| {
                eprintln!("Invalid base: '{}'. Setting to N.", c);
                &'N' // Default to 'N'
            })
        })
        .collect()
}

fn search(
    text: &str,
    pattern: &str,
    is_nucl: bool,
    is_circ: bool,
    tolerance: usize,
) -> Vec<(String, f64, usize, String, String)> {
    let mut matches = Vec::new();
    if is_nucl {
        // Prepare bit representations
        let extended_text = if is_circ {
            let concatenated = text.to_owned() + text;
            concatenated
        } else {
            text.to_string()
        };
        let text_bits = dna_to_bits(&extended_text);
        let pattern_bits = dna_to_bits(pattern);
        let rev_pattern_bits = dna_to_bits(&reverse_complement(pattern));
        // Sliding window search
        let len_limit = if is_circ {
            text.len()
        } else {
            text_bits.len() - pattern_bits.len() + 1
        };
        for i in 0..len_limit {
            // Forward match
            let mismatch = text_bits[i..i + pattern_bits.len()]
                .iter()
                .zip(&pattern_bits)
                .filter(|(&t, &p)| (t & p) == 0)
                .count();
            if mismatch <= tolerance {
                let score = (pattern_bits.len() - mismatch) as f64 / pattern_bits.len() as f64;
                let end_pos = i + pattern.len();
                let topology = if end_pos > text.len() {
                    "circular"
                } else {
                    "linear"
                };
                matches.push((
                    format!(
                        "{}..{}",
                        i + 1,
                        if end_pos > text.len() {
                            end_pos - text.len()
                        } else {
                            end_pos
                        }
                    ),
                    score,
                    mismatch,
                    "forward".to_string(),
                    topology.to_string(),
                ));
            }
            // Reverse match
            let rev_mismatch = text_bits[i..i + rev_pattern_bits.len()]
                .iter()
                .zip(&rev_pattern_bits)
                .filter(|(&t, &p)| (t & p) == 0)
                .count();
            if rev_mismatch <= tolerance && rev_mismatch <= mismatch {
                let score =
                    (rev_pattern_bits.len() - rev_mismatch) as f64 / rev_pattern_bits.len() as f64;
                let end_pos = i + pattern.len();
                let topology = if end_pos > text.len() {
                    "circular"
                } else {
                    "linear"
                };
                matches.push((
                    format!(
                        "{}..{}",
                        i + 1,
                        if end_pos > text.len() {
                            end_pos - text.len()
                        } else {
                            end_pos
                        }
                    ),
                    score,
                    rev_mismatch,
                    "reverse".to_string(),
                    topology.to_string(),
                ));
            }
        }
    } else {
        // Amino acid search
        for i in 0..=text.len() - pattern.len() {
            let mismatch = text[i..i + pattern.len()]
                .chars()
                .zip(pattern.chars())
                .filter(|(t, p)| t != p)
                .count();
            if mismatch <= tolerance {
                let score = (pattern.len() - mismatch) as f64 / pattern.len() as f64;
                matches.push((
                    format!("{}..{}", i + 1, i + pattern.len()),
                    score,
                    mismatch,
                    ".".to_string(), // No reverse complement
                    ".".to_string(), // No circularity
                ));
            }
        }
    }
    println!("{:#?}", matches);
    matches
}

fn compare_sequences(
    seqs1: Vec<Seq>,
    seqs2: Vec<Seq>,
    sim_thres: f64,
    cov_thres: f64,
    is_nucl: bool,
    is_circ: bool,
    num_threads: usize,
) -> Vec<(String, String, String, f64, usize, String, String)> {
    let seqs1 = Arc::new(seqs1);
    let seqs2 = Arc::new(seqs2);
    let results = Arc::new(Mutex::new(Vec::with_capacity(seqs1.len() * seqs2.len())));
    let mut handles = Vec::with_capacity(num_threads);

    // We will compare each unique pair of seq1 and seq2 (i < j) to avoid redundant comparisons
    let total_combinations = seqs1.len() * seqs2.len();
    let chunk_size = total_combinations / num_threads + 1;

    for i in 0..num_threads {
        let results = Arc::clone(&results);
        let seqs1 = Arc::clone(&seqs1);
        let seqs2 = Arc::clone(&seqs2);
        let start_index = i * chunk_size;
        let end_index = std::cmp::min((i + 1) * chunk_size, total_combinations);

        let handle = thread::spawn(move || {
            let mut local_results = Vec::with_capacity(chunk_size);

            for index in start_index..end_index {
                let seq1 = &seqs1[index / seqs2.len()]; // seq1 comes from seqs1
                let seq2 = &seqs2[index % seqs2.len()]; // seq2 comes from seqs2
                println!("{:#?} {:#?}", seq1, seq2);
                // Determine which sequence is the pattern and which is the text
                let (pattern, text) = if seq1.sequence.len() < seq2.sequence.len() {
                    (seq1, seq2)
                } else {
                    (seq2, seq1)
                };

                // Calculate the tolerance based on similarity threshold
                let tolerance = (sim_thres * pattern.sequence.len() as f64).round() as usize;

                let search_results = search(
                    &text.sequence,
                    &pattern.sequence,
                    is_nucl,
                    is_circ,
                    tolerance,
                );

                for (position, score, mismatch, orientation, topology) in search_results {
                    // Push the result with the required format including the IDs
                    if pattern.sequence.len() - mismatch
                        >= (text.sequence.len() as f64 * cov_thres).round() as usize
                    {
                        local_results.push((
                            pattern.id.clone(), // Query ID (from pattern)
                            text.id.clone(),    // Reference ID (from text)
                            position,           // Position (start..end)
                            score,              // Match score
                            mismatch,           // Number of mismatches
                            orientation,        // Orientation (forward/reverse)
                            topology,           // Topology (circular/linear)
                        ))
                    };
                }
            }
            let mut results_lock = results.lock().unwrap();
            results_lock.extend(local_results);
        });
        handles.push(handle);
    }
    // Wait for all threads to complete
    for handle in handles {
        handle.join().unwrap();
    }
    let results_lock = results.lock().unwrap();
    results_lock.clone()
}

fn write_results_to_tsv(
    results: &[(String, String, String, f64, usize, String, String)],
    output_file: &str,
) -> io::Result<()> {
    let mut file = File::create(output_file)?;
    // Write headers
    writeln!(
        file,
        "Query\tReference\tPosition\tScore\tMismatch\tOrientation\tTopology"
    )?;
    // Write the results
    for (query_id, ref_id, position, score, mismatch, orientation, topology) in results.iter() {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            query_id, ref_id, position, score, mismatch, orientation, topology
        )?;
    }
    Ok(())
}

fn main() -> io::Result<()> {
    // Command-line argument handling
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <query_file> <reference_file> <output_file> [similarity_threshold] [coverage_threshold] [num_threads] [nu/aa] [circular/linear]", args[0]);
    }
    let input_file1 = &args[1];
    let input_file2 = &args[2];
    let output_file = &args[3];
    let similarity_threshold: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(0.95);
    let coverage_threshold: f64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(0.5);
    let max_cpus = thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1); // Fallback to 1 if unavailable
    let num_threads: usize = args.get(6).and_then(|s| s.parse().ok()).unwrap_or_else(|| {
        let threads = 4;
        eprintln!(
            "Invalid number of threads or not provided. Using default value: {}",
            threads
        );
        std::cmp::min(threads, max_cpus)
    });
    if !Path::new(input_file1).exists() {
        eprintln!("The input file '{}' does not exist.", input_file1);
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "Input file not found",
        ));
    }
    if !Path::new(input_file2).exists() {
        eprintln!("The input file '{}' does not exist.", input_file2);
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "Input file not found",
        ));
    }
    if let Err(e) = File::create(output_file) {
        eprintln!("Error creating output file '{}': {}", output_file, e);
        return Err(e);
    }
    let file1_type = FileType::from_filename(input_file1)?;
    let file2_type = FileType::from_filename(input_file2)?;
    let sequences1 = read_sequences(input_file1, file1_type)?;
    let sequences2 = read_sequences(input_file2, file2_type)?;
    if sequences1.is_empty() {
        eprintln!("No sequences detected!");
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "No sequences detected",
        ));
    }
    if sequences2.is_empty() {
        eprintln!("No sequences detected!");
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "No sequences detected",
        ));
    }
    let mode = args.get(7).map(String::as_str).unwrap_or("nu");
    let topology = args.get(8).map(String::as_str).unwrap_or("linear");
    if mode != "nu" && mode != "aa" {
        panic!("Invalid mode '{}'. Expected 'nu' or 'aa'.", mode);
    }
    if topology != "circular" && topology != "linear" {
        panic!(
            "Invalid topology '{}'. Expected 'circular' or 'linear'",
            topology
        );
    }
    let is_circ = if topology == "circular" { true } else { false };
    let is_nucl = if mode == "nu" { true } else { false };
    if !is_nucl && is_circ {
        panic!("Circular topology is not available for amino acids");
    }
    let results = compare_sequences(
        sequences1,
        sequences2,
        similarity_threshold,
        coverage_threshold,
        is_nucl,
        is_circ,
        num_threads,
    );
    write_results_to_tsv(&results, output_file)?;
    Ok(())
}