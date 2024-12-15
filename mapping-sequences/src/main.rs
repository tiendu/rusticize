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

// Convert DNA sequence to bits representation
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
    .collect::<std::collections::HashMap<_, _>>();
    dna.chars()
        .map(|c| *nucl_to_bits.get(&c).unwrap())
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
    .collect::<std::collections::HashMap<_, _>>();
    dna.chars()
        .rev()
        .map(|c| *complement.get(&c).unwrap())
        .collect()
}

fn search(
    text: SearchMode,
    pattern: SearchMode,
    tolerance: usize,
) -> HashSet<(String, f64, usize)> {
    let mut matches = HashSet::new();

    match (text, pattern) {
        (SearchMode::Nucleotide(text), SearchMode::Nucleotide(pattern)) => {
            // Convert DNA strings to bit representations
            let text_bits = dna_to_bits(text);
            let pattern_bits = dna_to_bits(pattern);

            for i in 0..=text_bits.len() - pattern_bits.len() {
                let mismatch = text_bits[i..i + pattern_bits.len()]
                    .iter()
                    .zip(&pattern_bits)
                    .filter(|(&t, &p)| (t & p) == 0)
                    .count();
                if mismatch <= tolerance {
                    let score = (pattern_bits.len() - mismatch) as f64 / pattern_bits.len() as f64;
                    matches.insert((
                        format!("{}..{}", i + 1, i + pattern_bits.len()),
                        score,
                        mismatch,
                    ));
                }
            }
        }
        (SearchMode::AminoAcid(text), SearchMode::AminoAcid(pattern)) => {
            for i in 0..=text.len() - pattern.len() {
                let mismatch = text[i..i + pattern.len()]
                    .chars()
                    .zip(pattern_str.chars())
                    .filter(|(t, p)| t != p)
                    .count();
                if mismatch <= tolerance {
                    let score = (pattern.len() - mismatch) as f64 / pattern.len() as f64;
                    matches.insert((
                        format!("{}..{}", i + 1, i + pattern_str.len()),
                        score,
                        mismatch,
                    ));
                }
            }
        }
        _ => panic!("Mismatched search modes: both must be Nucleic acid or Amino acid."),
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
