use std::collections::{BTreeMap, HashSet};
use std::sync::{Arc, Mutex};
use std::thread;

fn build_de_bruijn_graph_for_read(read: &String, k: usize) -> BTreeMap<String, Vec<String>> {
    let mut graph: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for i in 0..=read.len() - k {
        let kmer = &read[i..i + k];
        let prefix = &kmer[..k - 1];
        let suffix = &kmer[1..];
        graph
            .entry(prefix.to_string())
            .or_insert_with(Vec::new)
            .push(suffix.to_string());
    }
    graph
}

fn generate_contigs(graph: &BTreeMap<String, Vec<String>>, k: usize) -> Vec<String> {
    let mut contigs = Vec::new();
    let mut seen_contigs = HashSet::new();

    for prefix in graph.keys() {
        let mut path = vec![prefix];
        let mut current_node = prefix;

        loop {
            if let Some(suffixes) = graph.get(current_node) {
                if suffixes.is_empty() {
                    break; // No outgoing edges, dead-end, break the loop
                }
                let mut extended = false;
                for suffix in suffixes {
                    path.push(suffix);
                    current_node = suffix;
                    extended = true;
                    break;
                }
                if !extended {
                    break;
                }
            } else {
                break;
            }
        }

        if path.len() > 1 {
            let mut contig = path[0].to_string();
            for node in &path[1..] {
                contig.push_str(&node[k - 1..]);
            }

            if contig.len() > k + 1 && !is_substring_of_existing(&contig, &seen_contigs) {
                contigs.push(contig.clone());
                seen_contigs.insert(contig);
            }
        }
    }
    println!("{:#?}", seen_contigs);
    contigs
}

fn is_substring_of_existing(contig: &str, seen_contigs: &HashSet<String>) -> bool {
    for existing in seen_contigs {
        if existing.contains(contig) {
            return true;
        }
    }
    false
}

fn main() {
    let reads = vec![
        "TTAGGG".to_string(),
        "ATGCT".to_string(),
        "TGCTTA".to_string(),
        "CTTACC".to_string(),
    ];

    let mut sorted_reads = reads.clone();
    sorted_reads.sort_by(|a, b| {
        let len_cmp = b.len().cmp(&a.len());
        if len_cmp == std::cmp::Ordering::Equal {
            a.cmp(b) // Lexicographic sorting for equal length
        } else {
            len_cmp
        }
    });

    let k = 4; // k-mer size
    let num_threads = 2; // Fixed number of threads
    let chunk_size = (sorted_reads.len() + num_threads - 1) / num_threads; // Calculate chunk size

    // Split the sorted reads into chunks
    let chunks: Vec<Vec<String>> = sorted_reads
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();

    // Initialize the shared graph
    let graph = Arc::new(Mutex::new(BTreeMap::new()));
    let mut handles = vec![];

    // Spawn a thread for each chunk of reads
    for chunk in chunks {
        let graph = Arc::clone(&graph);
        let handle = thread::spawn(move || {
            for read in chunk {
                let local_graph = build_de_bruijn_graph_for_read(&read, k);
                let mut graph = graph.lock().unwrap();
                for (prefix, suffixes) in local_graph {
                    graph
                        .entry(prefix)
                        .or_insert_with(Vec::new)
                        .extend(suffixes);
                }
            }
        });
        handles.push(handle);
    }

    // Wait for all threads to finish building the graph
    for handle in handles {
        handle.join().unwrap();
    }

    // Once all threads finish, generate contigs using the graph
    let graph = graph.lock().unwrap();
    let contigs = generate_contigs(&graph, k - 1);

    println!("Generated Contigs:");
    for contig in contigs {
        println!("{}", contig);
    }
}
