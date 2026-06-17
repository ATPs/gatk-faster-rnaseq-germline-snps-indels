use std::collections::HashMap;

/// Represents a K-mer in the assembly graph.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Kmer {
    pub sequence: Vec<u8>,
}

impl Kmer {
    pub fn new(sequence: &[u8]) -> Self {
        Self {
            sequence: sequence.to_vec(),
        }
    }
}

/// A vertex in the read threading graph, representing a unique K-mer.
#[derive(Debug, Clone)]
pub struct Vertex {
    pub sequence: Vec<u8>,
}

impl Vertex {
    pub fn new(sequence: &[u8]) -> Self {
        Self {
            sequence: sequence.to_vec(),
        }
    }
}

/// An edge connecting two K-mer vertices.
#[derive(Debug, Clone)]
pub struct Edge {
    pub is_ref: bool,
    pub multiplicity: usize,
}

impl Edge {
    pub fn new(is_ref: bool, multiplicity: usize) -> Self {
        Self {
            is_ref,
            multiplicity,
        }
    }
}

/// The De Bruijn Graph used for local assembly.
pub struct ReadThreadingGraph {
    pub kmer_size: usize,
    pub vertices: Vec<Vertex>,
    pub kmer_to_vertex: HashMap<Kmer, usize>,
    // Adjacency lists for out-edges and in-edges.
    pub out_edges: HashMap<usize, HashMap<usize, Edge>>,
    pub in_edges: HashMap<usize, HashMap<usize, Edge>>,
}

impl ReadThreadingGraph {
    pub fn new(kmer_size: usize) -> Self {
        Self {
            kmer_size,
            vertices: Vec::new(),
            kmer_to_vertex: HashMap::new(),
            out_edges: HashMap::new(),
            in_edges: HashMap::new(),
        }
    }

    /// Adds a vertex for the given k-mer sequence if it doesn't exist, and returns its index.
    pub fn get_or_create_vertex(&mut self, kmer_seq: &[u8]) -> usize {
        let kmer = Kmer::new(kmer_seq);
        if let Some(&v_idx) = self.kmer_to_vertex.get(&kmer) {
            v_idx
        } else {
            let v_idx = self.vertices.len();
            self.vertices.push(Vertex::new(kmer_seq));
            self.kmer_to_vertex.insert(kmer, v_idx);
            self.out_edges.insert(v_idx, HashMap::new());
            self.in_edges.insert(v_idx, HashMap::new());
            v_idx
        }
    }

    /// Adds an edge between two vertices.
    pub fn add_edge(&mut self, source: usize, target: usize, is_ref: bool) {
        // Update out-edge from source to target
        let source_out = self.out_edges.get_mut(&source).unwrap();
        let edge = source_out.entry(target).or_insert_with(|| Edge::new(is_ref, 0));
        edge.multiplicity += 1;
        if is_ref {
            edge.is_ref = true;
        }

        // Update in-edge from target to source
        let target_in = self.in_edges.get_mut(&target).unwrap();
        let edge = target_in.entry(source).or_insert_with(|| Edge::new(is_ref, 0));
        edge.multiplicity += 1;
        if is_ref {
            edge.is_ref = true;
        }
    }

    /// Extends a chain in the graph by one base.
    /// Given the previous vertex and the sequence containing the next k-mer,
    /// gets or creates the target vertex and adds an edge.
    pub fn extend_chain_by_one(
        &mut self,
        prev_vertex: usize,
        sequence: &[u8],
        next_kmer_start: usize,
        is_ref: bool,
    ) -> usize {
        let next_kmer_seq = &sequence[next_kmer_start..next_kmer_start + self.kmer_size];
        let next_vertex = self.get_or_create_vertex(next_kmer_seq);
        self.add_edge(prev_vertex, next_vertex, is_ref);
        next_vertex
    }

    /// Adds a sequence to the graph, connecting its kmers.
    pub fn add_sequence(&mut self, sequence: &[u8], is_ref: bool) {
        if sequence.len() < self.kmer_size {
            return;
        }
        
        let first_kmer = &sequence[0..self.kmer_size];
        let mut prev_vertex = self.get_or_create_vertex(first_kmer);

        for start_pos in 1..=sequence.len() - self.kmer_size {
            prev_vertex = self.extend_chain_by_one(prev_vertex, sequence, start_pos, is_ref);
        }
    }
}

use std::collections::BinaryHeap;
use std::cmp::Ordering;

#[derive(Clone, Debug)]
pub struct KBestPath {
    pub edges: Vec<(usize, usize)>, // (source, target)
    pub score: f64,
    pub is_reference: bool,
    pub last_vertex: usize,
}

impl PartialEq for KBestPath {
    fn eq(&self, other: &Self) -> bool {
        self.score == other.score
    }
}

impl Eq for KBestPath {}

impl PartialOrd for KBestPath {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.score.partial_cmp(&other.score)
    }
}

impl Ord for KBestPath {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}

impl ReadThreadingGraph {
    pub fn find_best_haplotypes(&self, source: usize, sink: usize, max_haplotypes: usize) -> Vec<KBestPath> {
        let mut result = Vec::new();
        let mut queue = BinaryHeap::new();
        
        queue.push(KBestPath {
            edges: Vec::new(),
            score: 0.0,
            is_reference: true,
            last_vertex: source,
        });

        let mut vertex_counts: HashMap<usize, usize> = HashMap::new();
        let max_expansions = 10_000_usize; // Hard limit to avoid blowup
        let mut expansions = 0_usize;

        while let Some(path) = queue.pop() {
            expansions += 1;
            if expansions > max_expansions {
                break;
            }

            if path.last_vertex == sink {
                result.push(path);
                if result.len() >= max_haplotypes {
                    break;
                }
            } else {
                let count = vertex_counts.entry(path.last_vertex).or_insert(0);
                if *count < max_haplotypes {
                    *count += 1;

                    // Skip paths that are unreasonably long (likely cycling)
                    if path.edges.len() > self.vertices.len() * 2 {
                        continue;
                    }

                    if let Some(out_edges) = self.out_edges.get(&path.last_vertex) {
                        let total_multiplicity: usize = out_edges.values().map(|e| e.multiplicity).sum();
                        let total_f64 = total_multiplicity as f64;

                        for (&target, edge) in out_edges {
                            if edge.multiplicity == 0 {
                                continue;
                            }
                            
                            let mut new_edges = path.edges.clone();
                            new_edges.push((path.last_vertex, target));
                            
                            let edge_prob = (edge.multiplicity as f64).log10() - total_f64.log10();
                            let new_score = path.score + edge_prob;

                            queue.push(KBestPath {
                                edges: new_edges,
                                score: new_score,
                                is_reference: path.is_reference && edge.is_ref,
                                last_vertex: target,
                            });
                        }
                    }
                }
            }
        }
        
        result
    }
    
    pub fn reconstruct_sequence(&self, path: &KBestPath) -> Vec<u8> {
        if path.edges.is_empty() {
            return self.vertices[path.last_vertex].sequence.clone();
        }
        
        let mut seq = self.vertices[path.edges[0].0].sequence.clone();
        
        for &(_, target) in &path.edges {
            let target_seq = &self.vertices[target].sequence;
            seq.push(*target_seq.last().unwrap());
        }
        
        seq
    }
}
