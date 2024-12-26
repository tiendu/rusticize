use std::collections::{HashMap, HashSet};

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

fn generate_contigs(graph: &HashMap<String, Vec<String>>, k: usize) -> Vec<String> {
    let mut contigs = Vec::new();
    let mut visited_edges: HashSet<(String, String)> = HashSet::new();
    // Determine incoming edges for each node
    let incoming_edges: HashMap<&String, usize> = graph
        .values()
        .flat_map(|suffixes| suffixes.iter())
        .fold(HashMap::new(), |mut acc, suffix| {
            *acc.entry(suffix).or_insert(0) += 1;
            acc
        });
    // Identify start nodes
    let mut start_nodes: Vec<String> = graph
        .keys()
        .filter(|node| !incoming_edges.contains_key(*node))
        .cloned()
        .collect();
    // Sort start nodes for consistent traversal
    start_nodes.sort();
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
    // Traverse the graph starting from each node
    for start_node in start_nodes {
        let mut stack = vec![start_node.clone()];
        let mut path = Vec::new();
        while let Some(current_node) = stack.pop() {
            path.push(current_node.clone());
            if let Some(suffixes) = graph.get(&current_node) {
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
    ];
    println!("Reads:");
    for read in &reads {
        println!("{}", read);
    }

    let k = 4; // k-mer size
    let dbg = build_de_bruijn_graph(&reads, k);
    println!("De Bruijn Graph:");
    for (node, edges) in &dbg {
        println!("{} -> {:?}", node, edges);
    }

    let contigs = generate_contigs(&dbg, k - 1);
    for contig in contigs {
        println!("{}", contig);
    }
}
