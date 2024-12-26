use std::collections::{HashMap, HashSet};
use std::sync::{Arc, RwLock};
use std::thread;

fn build_de_bruijn_graph(reads: &[String], k: usize) -> HashMap<String, Vec<String>> {
    let mut graph: HashMap<String, Vec<String>> = HashMap::new();
    for read in reads {
        for i in 0..=read.len() - k {
            let kmer = &read[i..i + k];
            let prefix = &kmer[..k - 1];
            let suffix = &kmer[1..];
            graph
                .entry(prefix.to_string())
                .or_insert_with(Vec::new)
                .push(suffix.to_string());
        }
    }
    graph
}

fn generate_contigs(
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

fn main() {
    let reads = vec![
        "TTAGGG".to_string(),
        "ATGCT".to_string(),
        "TGCTTA".to_string(),
        "CTTACC".to_string(),
        "AAATAT".to_string(),
        "AAATAC".to_string(),
    ];
    let k = 4; // k-mer size
    let max_cpus = thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1); // Fallback to 1 if unavailable
    let num_threads = 4;
    let chunk_size = (reads.len() + num_threads - 1) / num_threads;
    let reads_chunks: Vec<_> = reads
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    let graph = Arc::new(RwLock::new(HashMap::new()));
    let mut handles = Vec::new();
    for chunk in reads_chunks {
        let graph_clone = Arc::clone(&graph);
        let handle = thread::spawn(move || {
            let local_graph = build_de_bruijn_graph(&chunk, k);
            let mut global_graph = graph_clone.write().unwrap(); // Acquire write lock
            for (prefix, suffixes) in local_graph {
                global_graph
                    .entry(prefix)
                    .or_insert_with(Vec::new)
                    .extend(suffixes);
            }
        });
        handles.push(handle);
    }
    for handle in handles {
        handle.join().unwrap();
    }

    let graph_read = graph.read().unwrap(); // Acquire read lock
    println!("De Bruijn Graph:");
    for (node, edges) in &*graph_read {
        println!("{} -> {:?}", node, edges);
    }

    let contigs = Arc::new(RwLock::new(Vec::new()));
    let start_nodes: Vec<String> = graph_read
        .keys()
        .filter(|node| {
            // For each node, check if it has no incoming edges (i.e., it's not a suffix of any other node)
            !graph_read
                .values()
                .flat_map(|suffixes| suffixes.iter())
                .any(|suffix| suffix == *node)
        })
        .cloned()
        .collect();
    let chunk_size = (start_nodes.len() + num_threads - 1) / num_threads;
    let start_chunks: Vec<_> = start_nodes
        .chunks(chunk_size)
        .map(|chunk| chunk.to_vec())
        .collect();
    let mut handles = Vec::new();
    for chunk in start_chunks {
        let graph_clone = Arc::clone(&graph);
        let contigs_clone = Arc::clone(&contigs);
        let handle = thread::spawn(move || {
            let local_contigs = generate_contigs(&graph_clone, k - 1, chunk);
            contigs_clone.write().unwrap().extend(local_contigs); // Acquire write lock
        });
        handles.push(handle);
    }
    for handle in handles {
        handle.join().unwrap();
    }
    let contigs_read = contigs.read().unwrap(); // Acquire read lock
    println!("Contigs:");
    for contig in &*contigs_read {
        println!("{}", contig);
    }
}
