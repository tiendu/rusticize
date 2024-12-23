// Standard Library Imports
use std::env;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use std::sync::Arc;
use std::thread;

// Enum to represent the sequence comparison mode
#[derive(Debug, PartialEq)]
enum SearchMode {
    Nucleotide,
    AminoAcid,
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
        .map(|c| *nucl_to_bits.get(&c).unwrap_or_else(|| {
            eprintln!("Invalid base '{}'. Setting to N.", c);
            &0b1111 // Default to 'N'
        }))
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
        .map(|c| *complement.get(&c).unwrap_or_else(|| {
            eprintln!("Invalid base: '{}'. Setting to N.", c);
            &'N' // Default to 'N'
        }))
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
        let len_limit = if is_circ { text.len() } else { text_bits.len() - pattern_bits.len() + 1 };
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
                let topology = if end_pos > text.len() { "circular" } else { "linear" };
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
                let score = (rev_pattern_bits.len() - rev_mismatch) as f64 / rev_pattern_bits.len() as f64;
                let end_pos = i + pattern.len();
                let topology = if end_pos > text.len() { "circular" } else { "linear" };
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
    matches
}

fn main() -> io::Result<()> {
    // Command-line argument handling
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <query_file> <reference_file> <output_file> [similarity_threshold] [coverage_threshold] [num_threads] [nu/aa] [circular/linear]", args[0]);
    }

    let query_file = &args[1];
    let reference_file = &args[2];
    let output_file = &args[3];
    let similarity_threshold: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(95.0);
    let coverage_threshold: f64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(50.0);
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
    let mode = args.get(7).map(String::as_str).unwrap_or("nu");
    let is_circular = args.get(8).map(String::as_str).unwrap_or("linear");

    if mode != "nu" && mode != "aa" {
        panic!("Invalid mode '{}'. Expected 'nu' or 'aa'.", mode);
    }
    let search_mode = if mode == "nu" {
        SearchMode::Nucleotide
    } else {
        SearchMode::AminoAcid
    };

    if search_mode == SearchMode::AminoAcid && is_circular {
        panic!("Circular topology is not available for amino acids");
    }

    if !Path::new(query_file).exists() {
        eprintln!("The query file '{}' does not exists.", query_file);
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "Query file not found",
        ));
    }

    if !Path::new(reference_file).exists() {
        eprintln!("The reference file '{}' does not exists.", reference_file);
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "Reference file not found",
        ));
    }

    // Read query and reference sequences
    let query_sequences = read_sequences(query_file)?;
    let reference_sequences = read_sequences(reference_file)?;

    write_results_to_csv(&results, output_path)?;

    Ok(())
}
