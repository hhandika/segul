//! Translation tables for the standard genetic code and the mitochondrial genetic code.
/// The Translation struct contains the translation tables for the standard genetic code and the
/// mitochondrial genetic code.
/// The translation tables are based on the NCBI translation tables.
/// See: <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>.
/// The translation tables are stored in a vector of tuples.
/// The first element of the tuple is the codon and the second element is the amino acid.
/// This implementation may change in the future.
///
/// ## Implemented Tables:
/// 1. The Standard Code (Table 1)
/// 2. The Vertebrate Mitochondrial Code (Table 2)
/// 3. The Yeast Mitochondrial Code (Table 3)
/// 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (Table 4)
/// 5. The Invertebrate Mitochondrial Code (Table 5)
/// 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code (Table 6)
/// 7. The Echinoderm and Flatworm Mitochondrial Code (Table 9)
/// 8. The Euplotid Nuclear Code (Table 10)
/// 9. The Bacterial and Plant Plastid Code (Table 11)
/// 10. The Bacterial, Archaeal and Plant Plastid Code (Table 12)
/// 11. The Alternative Yeast Nuclear Code (Table 13)
/// 12. The Ascidian Mitochondrial Code (Table 14)
/// 13. Chlorophycean Mitochondrial Code (Table 16)
/// 14. Trematode Mitochondrial Code (Table 21)
/// 15. Mesodinium Nuclear Code (Table 29)
/// 16. Peritrich Nuclear Code (Table 30)
/// 17. Blastocrithidia Nuclear Code (Table 31)
/// 18. Cephalodiscidae Mitochondrial UAA-Tyr Code (Table 33)
///
/// ## Example
/// ```
/// use segul::helper::translation::NcbiTables;
///
/// let translation = NcbiTables::new();
/// let standard_code = translation.standard_code();
///
/// let codon = "ATG";
/// let amino_acid = standard_code.get(codon);
/// assert_eq!(amino_acid, Some(&String::from("M")));
///```
use ahash::AHashMap as HashMap;
use lazy_static::lazy_static;

lazy_static! {
    /// A version of lazyly initialized of standard translation code.
    /// It may provide better performance than the `NcbiTables::standard_code()` method.
    pub static ref STANDARD_CODE: HashMap<&'static str, &'static str> = {
        let mut m = HashMap::new();
        for (codon, amino_acid) in TRANSLATION {
            m.insert(codon.as_ref(), amino_acid.as_ref());
        }
        m
    };
}

const TRANSLATION: &[(&str, &str)] = &[
    ("TTT", "F"),
    ("TTC", "F"),
    ("TTA", "L"),
    ("TTG", "L"),
    ("TCT", "S"),
    ("TCC", "S"),
    ("TCA", "S"),
    ("TCG", "S"),
    ("TAT", "Y"),
    ("TAC", "Y"),
    ("TAA", "*"),
    ("TAG", "*"),
    ("TGT", "C"),
    ("TGC", "C"),
    ("TGA", "*"),
    ("TGG", "W"),
    ("CTT", "L"),
    ("CTC", "L"),
    ("CTA", "L"),
    ("CTG", "L"),
    ("CCT", "P"),
    ("CCC", "P"),
    ("CCA", "P"),
    ("CCG", "P"),
    ("CAT", "H"),
    ("CAC", "H"),
    ("CAA", "Q"),
    ("CAG", "Q"),
    ("CGT", "R"),
    ("CGC", "R"),
    ("CGA", "R"),
    ("CGG", "R"),
    ("ATT", "I"),
    ("ATC", "I"),
    ("ATA", "I"),
    ("ATG", "M"),
    ("ACT", "T"),
    ("ACC", "T"),
    ("ACA", "T"),
    ("ACG", "T"),
    ("AAT", "N"),
    ("AAC", "N"),
    ("AAA", "K"),
    ("AAG", "K"),
    ("AGT", "S"),
    ("AGC", "S"),
    ("AGA", "R"),
    ("AGG", "R"),
    ("GTT", "V"),
    ("GTC", "V"),
    ("GTA", "V"),
    ("GTG", "V"),
    ("GCT", "A"),
    ("GCC", "A"),
    ("GCA", "A"),
    ("GCG", "A"),
    ("GAT", "D"),
    ("GAC", "D"),
    ("GAA", "E"),
    ("GAG", "E"),
    ("GGT", "G"),
    ("GGC", "G"),
    ("GGA", "G"),
    ("GGG", "G"),
    ("NNN", "X"),
    ("---", "-"),
    ("???", "?"),
];

pub struct NcbiTables<'a> {
    translation: &'a [(&'a str, &'a str)],
}

impl<'a> NcbiTables<'a> {
    /// The constructor for the NCBI Table struct.
    /// The constructor initializes the translation tables.
    pub fn new() -> Self {
        Self {
            translation: TRANSLATION,
        }
    }

    /// Returns a HashMap of the standard genetic code (NCBI Table 1)
    pub fn standard_code(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            code.insert(codon.to_string(), protein.to_string());
        });
        code
    }

    /// Returns a HashMap of the vertebrate mitochondrial genetic code (NCBI Table 2)
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

    /// Returns a HashMap of the yeast mitochondrial genetic code (NCBI Table 3)
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

    /// Returns a HashMap of the mold, protozoan, and coelenterate mitochondrial genetic code (NCBI Table 4)
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

    /// Returns a HashMap of the invertebrate mitochondrial genetic code (NCBI Table 5)
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

    /// Returns a HashMap of the ciliate, dasycladacean and hexamita nuclear genetic code (NCBI Table 6)
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

    /// Returns a HashMap of the echinoderm and flatworm mitochondrial genetic code (NCBI Table 9)
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

    /// Returns a HashMap of the euplotid nuclear genetic code (NCBI Table 10)
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

    /// Returns a HashMap of the alternative yeast nuclear genetic code (NCBI Table 11 and 12)
    /// NCBI Table 12. NCBI Table 11 is the same as the Table 1.
    /// Uses by Candida albicans, Candida cylindracea, Candida melibiosica,
    /// Candida parapsilosis, and Candida rugosa (Ohama et al., 1993).
    /// Other yeast, uses the standard code.
    pub fn alt_yeast_nu(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "CTG" => code.insert(codon.to_string(), String::from("S")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Returns a HashMap of the ascidian mitochondrial genetic code (NCBI Table 13)
    pub fn ascidian_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "AGA" => code.insert(codon.to_string(), String::from("G")),
                "AGG" => code.insert(codon.to_string(), String::from("G")),
                "ATA" => code.insert(codon.to_string(), String::from("M")),
                "TGA" => code.insert(codon.to_string(), String::from("W")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Returns a HashMap of the alternative flatworm mitochondrial genetic code (NCBI Table 14)
    pub fn alt_flatworm_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "AAA" => code.insert(codon.to_string(), String::from("N")),
                "AGA" => code.insert(codon.to_string(), String::from("S")),
                "AGG" => code.insert(codon.to_string(), String::from("S")),
                "TAA" => code.insert(codon.to_string(), String::from("Y")),
                "TGA" => code.insert(codon.to_string(), String::from("W")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Returns a HashMap of the blepharisma nuclear genetic code (NCBI Table 15)
    pub fn chlorophycean_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TAG" => code.insert(codon.to_string(), String::from("L")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Returns a HashMap of the alternative blepharisma mitochondrial genetic code (NCBI Table 16)
    pub fn trematode_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TGA" => code.insert(codon.to_string(), String::from("W")),
                "ATA" => code.insert(codon.to_string(), String::from("M")),
                "AGA" => code.insert(codon.to_string(), String::from("S")),
                "AGG" => code.insert(codon.to_string(), String::from("S")),
                "AAA" => code.insert(codon.to_string(), String::from("N")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Returns a HashMap of the alternative bryophyte mitochondrial genetic code (NCBI Table 17)
    pub fn scenedesmus_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TCA" => code.insert(codon.to_string(), String::from("*")),
                "TAG" => code.insert(codon.to_string(), String::from("L")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Returns a HashMap of the alternative pterobranchia mitochondrial genetic code (NCBI Table 18)
    pub fn thraustochytrium_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TTA" => code.insert(codon.to_string(), String::from("*")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Returns a HashMap of the alternative echinoderm mitochondrial genetic code (NCBI Table 19)
    pub fn rhabdopleuridae_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "AGA" => code.insert(codon.to_string(), String::from("S")),
                "AGG" => code.insert(codon.to_string(), String::from("K")),
                "TGA" => code.insert(codon.to_string(), String::from("W")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Returns a HashMap of the alternative euplotid mitochondrial genetic code (NCBI Table 20)
    pub fn candid_div_sr1_gracil(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TGA" => code.insert(codon.to_string(), String::from("G")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Returns a HashMap of the alternative euplotid mitochondrial genetic code (NCBI Table 21)
    pub fn pachysolen_tanno_nu(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "CTG" => code.insert(codon.to_string(), String::from("A")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Return a HashMap of the alternative flatworm mitochondrial genetic code (NCBI Table 29)
    pub fn mesodinium_nu(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TAA" => code.insert(codon.to_string(), String::from("Y")),
                "TAG" => code.insert(codon.to_string(), String::from("Y")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Return a HashMap of the alternative flatworm mitochondrial genetic code (NCBI Table 30)
    pub fn peritrich_nu(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TAA" => code.insert(codon.to_string(), String::from("E")),
                "TAG" => code.insert(codon.to_string(), String::from("E")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }

    /// Return a HashMap of the alternative flatworm mitochondrial genetic code (NCBI Table 33)
    pub fn cephalodiscidae_mtdna(&self) -> HashMap<String, String> {
        let mut code = HashMap::new();

        self.translation.iter().for_each(|(codon, protein)| {
            match codon.as_ref() {
                "TAA" => code.insert(codon.to_string(), String::from("Y")),
                "TGA" => code.insert(codon.to_string(), String::from("W")),
                "AGA" => code.insert(codon.to_string(), String::from("S")),
                "AGG" => code.insert(codon.to_string(), String::from("K")),
                _ => code.insert(codon.to_string(), protein.to_string()),
            };
        });
        code
    }
}
