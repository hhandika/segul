pub enum SeqFormat {
    Phylip,
    Fasta,
}

// Alphabeth for dna. Include IUPAC characters plus missing (?)
pub fn valid_dna() -> String {
    String::from("ACGTRYSWKMBDHVNacgtryswkmbdhvn.-?")
}
