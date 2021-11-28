use ahash::AHashMap as HashMap;

pub struct NcbiTables {
    translation: Vec<(String, String)>,
}

impl NcbiTables {
    pub fn new() -> Self {
        Self {
            translation: vec![
                (String::from("TTT"), String::from("F")),
                (String::from("TTC"), String::from("F")),
                (String::from("TTA"), String::from("L")),
                (String::from("TTG"), String::from("L")),
                (String::from("CTT"), String::from("L")),
                (String::from("CTC"), String::from("L")),
                (String::from("CTA"), String::from("L")),
                (String::from("CTG"), String::from("L")),
                (String::from("ATT"), String::from("I")),
                (String::from("ATC"), String::from("I")),
                (String::from("ATA"), String::from("I")),
                (String::from("ATG"), String::from("M")),
                (String::from("GTT"), String::from("V")),
                (String::from("GTC"), String::from("V")),
                (String::from("GTA"), String::from("V")),
                (String::from("GTG"), String::from("V")),
                (String::from("TCT"), String::from("S")),
                (String::from("TCC"), String::from("S")),
                (String::from("TCA"), String::from("S")),
                (String::from("TCG"), String::from("S")),
                (String::from("CCT"), String::from("P")),
                (String::from("CCC"), String::from("P")),
                (String::from("CCA"), String::from("P")),
                (String::from("CCG"), String::from("P")),
                (String::from("ACT"), String::from("T")),
                (String::from("ACC"), String::from("T")),
                (String::from("ACA"), String::from("T")),
                (String::from("ACG"), String::from("T")),
                (String::from("GCT"), String::from("A")),
                (String::from("GCC"), String::from("A")),
                (String::from("GCA"), String::from("A")),
                (String::from("GCG"), String::from("A")),
                (String::from("TAT"), String::from("Y")),
                (String::from("TAC"), String::from("Y")),
                (String::from("TAA"), String::from("*")),
                (String::from("TAG"), String::from("*")),
                (String::from("CAT"), String::from("H")),
                (String::from("CAC"), String::from("H")),
                (String::from("CAA"), String::from("Q")),
                (String::from("CAG"), String::from("Q")),
                (String::from("AAT"), String::from("N")),
                (String::from("AAC"), String::from("N")),
                (String::from("AAA"), String::from("K")),
                (String::from("AAG"), String::from("K")),
                (String::from("GAT"), String::from("D")),
                (String::from("GAC"), String::from("D")),
                (String::from("GAA"), String::from("E")),
                (String::from("GAG"), String::from("E")),
                (String::from("TGT"), String::from("C")),
                (String::from("TGC"), String::from("C")),
                (String::from("TGA"), String::from("*")),
                (String::from("TGG"), String::from("W")),
                (String::from("CGT"), String::from("R")),
                (String::from("CGC"), String::from("R")),
                (String::from("CGA"), String::from("R")),
                (String::from("CGG"), String::from("R")),
                (String::from("AGT"), String::from("S")),
                (String::from("AGC"), String::from("S")),
                (String::from("AGA"), String::from("R")),
                (String::from("AGG"), String::from("R")),
                (String::from("GGT"), String::from("G")),
                (String::from("GGC"), String::from("G")),
                (String::from("GGA"), String::from("G")),
                (String::from("GGG"), String::from("G")),
                (String::from("NNN"), String::from("X")),
                (String::from("???"), String::from("?")),
                (String::from("---"), String::from("-")),
            ],
        }
    }

    // NCBI Table 1. NCBI Table 2 uses the same translation table as NCBI Table 1.
    pub fn standard_code(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            code.insert(codon.to_string(), protein.to_string());
        });
        code
    }

    // NCBI Table 2.
    pub fn vert_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "AGA" => code.insert(codon.to_string(), String::from("*")),
                "AGG" => code.insert(codon.to_string(), String::from("*")),
                "ATA" => code.insert(codon.to_string(), String::from("M")),
                "TGA" => code.insert(codon.to_string(), String::from("W")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    // NCBI Table 3.
    pub fn yeast_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "ATA" => code.insert(codon.to_string(), String::from("M")),
                "CTT" => code.insert(codon.to_string(), String::from("T")),
                "CTC" => code.insert(codon.to_string(), String::from("T")),
                "CTA" => code.insert(codon.to_string(), String::from("T")),
                "CTG" => code.insert(codon.to_string(), String::from("T")),
                "TGA" => code.insert(codon.to_string(), String::from("W")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    // NCBI Table 4.
    pub fn moldprotocoe_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TGA" => code.insert(codon.to_string(), String::from("W")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    // NCBI Table 5.
    pub fn invert_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "AGA" => code.insert(codon.to_string(), String::from("S")),
                "AGG" => code.insert(codon.to_string(), String::from("S")),
                "ATA" => code.insert(codon.to_string(), String::from("M")),
                "TGA" => code.insert(codon.to_string(), String::from("W")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    // NCBI Table 6.
    pub fn cildashex_nudna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TAA" => code.insert(codon.to_string(), String::from("Q")),
                "TAG" => code.insert(codon.to_string(), String::from("Q")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    // NCBI Table 9.
    pub fn echiflatworm_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "AAA" => code.insert(codon.to_string(), String::from("N")),
                "AGA" => code.insert(codon.to_string(), String::from("S")),
                "AGG" => code.insert(codon.to_string(), String::from("S")),
                "TGA" => code.insert(codon.to_string(), String::from("W")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    // NCBI Table 10.
    pub fn euplotid_nudna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TGA" => code.insert(codon.to_string(), String::from("C")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }
}
