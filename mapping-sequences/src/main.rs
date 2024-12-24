// Standard Library Imports
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

// External Crate Imports
use flate2::read::GzDecoder;

#[derive(Debug)]
struct ArgsStruct {
    query_file: String,
    reference_file: String,
    output_file: String,
    similarity_threshold: Option<f64>,
    coverage_threshold: Option<f64>,
    num_threads: Option<usize>,
    is_nucleotide: bool,
    is_circular: bool,
}

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

fn dna_to_bits(dna: &str) -> Vec<u8> {
    dna.bytes()
        .map(|c| match c {
            b'A' => 0b0001,
            b'T' => 0b0010,
            b'G' => 0b0100,
            b'C' => 0b1000,
            b'W' => 0b0011,
            b'R' => 0b0101,
            b'Y' => 0b1010,
            b'S' => 0b1100,
            b'K' => 0b0110,
            b'M' => 0b1001,
            b'B' => 0b1110,
            b'D' => 0b0111,
            b'H' => 0b1011,
            b'V' => 0b1101,
            b'N' => 0b1111,
            _ => {
                eprintln!("Invalid base '{}'. Setting to N.", c as char);
                0b1111
            }
        })
        .collect()
}

fn reverse_complement(dna: &str) -> String {
    dna.bytes()
        .rev()
        .map(|c| match c {
            b'A' => 'T',
            b'T' => 'A',
            b'G' => 'C',
            b'C' => 'G',
            b'R' => 'Y',
            b'Y' => 'R',
            b'S' => 'S',
            b'W' => 'W',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'D' => 'H',
            b'H' => 'D',
            b'V' => 'B',
            b'N' => 'N',
            _ if (c as char).is_ascii() => {
                eprintln!("Invalid base '{}'. Setting to N.", c as char);
                'N'
            }
            _ => {
                eprintln!("Non-ASCII character detected. Setting to N.");
                'N'
            }
        })
        .collect()
}

fn compare(
    text: &str,
    pattern: &str,
    is_nucl: bool,
    is_circ: bool,
    max_mismatches: usize,
    min_coverage: usize,
) -> Vec<(String, f64, f64, usize, String, String)> {
    let mut matches = Vec::new();
    if is_nucl {
        // Prepare bit representations for nucleotide sequences
        let extended_text = if is_circ {
            text.to_owned() + text
        } else {
            text.to_string()
        };
        let text_bits = dna_to_bits(&extended_text);
        let pattern_bits = dna_to_bits(pattern);
        let rev_pattern_bits = dna_to_bits(&reverse_complement(pattern));
        let len_limit = if is_circ {
            text.len()
        } else {
            text_bits.len() - pattern_bits.len() + 1
        };
        for i in 0..len_limit {
            // Forward match calculation
            let forward_mismatch = text_bits[i..i + pattern_bits.len()]
                .iter()
                .zip(&pattern_bits)
                .filter(|(&t, &p)| (t & p) == 0)
                .count();
            if forward_mismatch <= max_mismatches
                && pattern.len() - forward_mismatch >= min_coverage
            {
                let score = (pattern.len() - forward_mismatch) as f64 / pattern.len() as f64;
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
                    (pattern.len() - forward_mismatch) as f64 / text.len() as f64,
                    forward_mismatch,
                    "forward".to_string(),
                    topology.to_string(),
                ));
            }
            // Reverse match calculation
            let reverse_mismatch = text_bits[i..i + rev_pattern_bits.len()]
                .iter()
                .zip(&rev_pattern_bits)
                .filter(|(&t, &p)| (t & p) == 0)
                .count();
            if reverse_mismatch <= max_mismatches
                && pattern.len() - reverse_mismatch >= min_coverage
            {
                let score = (pattern.len() - reverse_mismatch) as f64 / pattern.len() as f64;
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
                    (pattern.len() - reverse_mismatch) as f64 / text.len() as f64,
                    reverse_mismatch,
                    "reverse".to_string(),
                    topology.to_string(),
                ));
            }
        }
    } else {
        // Amino acid sequence search
        for i in 0..=text.len() - pattern.len() {
            let mismatch_count = text[i..i + pattern.len()]
                .chars()
                .zip(pattern.chars())
                .filter(|(t, p)| t != p)
                .count();
            if mismatch_count <= max_mismatches && pattern.len() - mismatch_count >= min_coverage {
                let score = (pattern.len() - mismatch_count) as f64 / pattern.len() as f64;
                matches.push((
                    format!("{}..{}", i + 1, i + pattern.len()),
                    score,
                    (pattern.len() - mismatch_count) as f64 / text.len() as f64,
                    mismatch_count,
                    ".".to_string(),
                    ".".to_string(),
                ));
            }
        }
    }
    matches
}

fn search(
    query_seqs: Vec<Seq>,
    target_seqs: Vec<Seq>,
    similarity_threshold: f64,
    coverage_threshold: f64,
    is_nucl: bool,
    is_circ: bool,
    num_threads: usize,
) -> Vec<(String, String, String, f64, usize, f64, String, String)> {
    let query_seqs = Arc::new(query_seqs);
    let target_seqs = Arc::new(target_seqs);
    let results = Arc::new(Mutex::new(Vec::new()));
    let progress = Arc::new(AtomicUsize::new(0));
    let total_combinations = query_seqs.len() * target_seqs.len();
    let chunk_size = total_combinations / num_threads + 1;
    let mut handles = Vec::with_capacity(num_threads);
    for thread_id in 0..num_threads {
        let results = Arc::clone(&results);
        let progress = Arc::clone(&progress);
        let query_seqs = Arc::clone(&query_seqs);
        let target_seqs = Arc::clone(&target_seqs);
        let start_index = thread_id * chunk_size;
        let end_index = std::cmp::min((thread_id + 1) * chunk_size, total_combinations);
        let handle = thread::spawn(move || {
            let mut local_results = Vec::new();
            for index in start_index..end_index {
                let query = &query_seqs[index / target_seqs.len()];
                let target = &target_seqs[index % target_seqs.len()];
                let (shorter_seq, longer_seq) = if query.sequence.len() < target.sequence.len() {
                    (query, target)
                } else {
                    (target, query)
                };
                let max_mismatches = (shorter_seq.sequence.len() as f64
                    * (1.0 - similarity_threshold))
                    .round() as usize;
                let min_coverage =
                    (longer_seq.sequence.len() as f64 * coverage_threshold).round() as usize;
                let search_results = compare(
                    &longer_seq.sequence,
                    &shorter_seq.sequence,
                    is_nucl,
                    is_circ,
                    max_mismatches,
                    min_coverage,
                );
                for (position, score, coverage, mismatch_count, orientation, topology) in
                    search_results
                {
                    local_results.push((
                        shorter_seq.id.clone(),
                        longer_seq.id.clone(),
                        position,
                        score,
                        mismatch_count,
                        coverage,
                        orientation,
                        topology,
                    ));
                }
                let current_progress = progress.fetch_add(1, Ordering::Relaxed) + 1;
                if current_progress % 100 == 0 || current_progress == total_combinations {
                    println!("Progress: {}/{}", current_progress, total_combinations);
                }
            }
            let mut results_lock = results.lock().unwrap();
            results_lock.extend(local_results);
        });
        handles.push(handle);
    }
    for handle in handles {
        handle.join().unwrap();
    }
    let results_lock = results.lock().unwrap();
    results_lock.clone()
}

fn write_results_to_tsv(
    results: &[(String, String, String, f64, usize, f64, String, String)],
    output_file: &str,
) -> io::Result<()> {
    let mut file = File::create(output_file)?;
    // Write headers
    writeln!(
        file,
        "Query\tReference\tPosition\tScore\tMismatch\tCoverage\tOrientation\tTopology"
    )?;
    // Write the results
    for (query_id, ref_id, position, score, mismatch, coverage, orientation, topology) in
        results.iter()
    {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            query_id, ref_id, position, score, mismatch, coverage, orientation, topology
        )?;
    }
    Ok(())
}

fn parse_args() -> io::Result<ArgsStruct> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!(
            "Usage: {} <query_file> <reference_file> <output_file> [similarity_threshold] [coverage_threshold] [num_threads] [nu/aa] [circular/linear]",
            args[0]
        );
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Insufficient arguments",
        ));
    }
    let similarity_threshold = args.get(4).and_then(|s| s.parse().ok());
    let coverage_threshold = args.get(5).and_then(|s| s.parse().ok());
    let num_threads = args.get(6).and_then(|s| s.parse().ok());
    let sequence_type = args.get(7).map(String::as_str).unwrap_or("nu");
    let sequence_topology = args.get(8).map(String::as_str).unwrap_or("linear");
    Ok(ArgsStruct {
        query_file: args[1].clone(),
        reference_file: args[2].clone(),
        output_file: args[3].clone(),
        similarity_threshold,
        coverage_threshold,
        num_threads,
        is_nucleotide: sequence_type == "nu",
        is_circular: sequence_topology == "circular",
    })
}

fn validate_file_exists(filename: &str) -> io::Result<()> {
    if !Path::new(filename).exists() {
        Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("File '{}' not found", filename),
        ))
    } else {
        Ok(())
    }
}

fn main() -> io::Result<()> {
    // Command-line argument handling
    let args = parse_args()?;
    validate_file_exists(&args.query_file)?;
    validate_file_exists(&args.reference_file)?;
    let query_file_type = FileType::from_filename(&args.query_file)?;
    let reference_file_type = FileType::from_filename(&args.reference_file)?;
    let query_sequences = read_sequences(&args.query_file, query_file_type)?;
    let reference_sequences = read_sequences(&args.reference_file, reference_file_type)?;
    let output_file = &args.output_file;
    let available_threads = thread::available_parallelism().map_or(1, |n| n.get());
    let num_threads = std::cmp::min(args.num_threads.unwrap_or(4), available_threads);
    if let Err(e) = File::create(output_file) {
        eprintln!("Error creating output file '{}': {}", output_file, e);
        return Err(e);
    }
    if !args.is_nucleotide && args.is_circular {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Circular topology is not supported for amino acids.",
        ));
    }
    let results = search(
        query_sequences,
        reference_sequences,
        args.similarity_threshold.unwrap_or(0.95),
        args.coverage_threshold.unwrap_or(0.5),
        args.is_nucleotide,
        args.is_circular,
        num_threads,
    );
    write_results_to_tsv(&results, &args.output_file)?;
    Ok(())
}
