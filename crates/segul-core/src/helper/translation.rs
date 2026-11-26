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
/// use segul_core::helper::translation::NcbiTables;
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

use crate::helper::types::GeneticCodes;

/// Direct 64-entry codon lookup table indexed by `(b0 << 4) | (b1 << 2) | b2`,
/// where `A/a = 0`, `C/c = 1`, `G/g = 2`, `T/t = 3`.
pub type CodonTable = [u8; 64];

const INVALID_BASE: u8 = 0xFF;

const BASE_2BIT: [u8; 256] = {
    let mut table = [INVALID_BASE; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

const fn encode_codon(c: &[u8]) -> Option<usize> {
    if c.len() != 3 {
        return None;
    }
    let b0 = BASE_2BIT[c[0] as usize];
    let b1 = BASE_2BIT[c[1] as usize];
    let b2 = BASE_2BIT[c[2] as usize];
    if b0 < 4 && b1 < 4 && b2 < 4 {
        Some(((b0 as usize) << 4) | ((b1 as usize) << 2) | (b2 as usize))
    } else {
        None
    }
}

const fn build_codon_table(overrides: &[(&[u8; 3], u8)]) -> CodonTable {
    let mut table = [b'X'; 64];
    let mut i = 0;
    while i < STD_CODE.len() {
        let (codon, aa) = STD_CODE[i];
        if let Some(idx) = encode_codon(codon.as_bytes()) {
            table[idx] = aa.as_bytes()[0];
        }
        i += 1;
    }
    let mut j = 0;
    while j < overrides.len() {
        let (codon, aa) = overrides[j];
        if let Some(idx) = encode_codon(codon) {
            table[idx] = aa;
        }
        j += 1;
    }
    table
}

/// Return the compile-time 64-entry byte lookup table for `code`.
pub const fn codon_table(code: GeneticCodes) -> CodonTable {
    match code {
        GeneticCodes::StandardCode | GeneticCodes::BacArchPlantPlast => build_codon_table(&[]),
        GeneticCodes::VertMtDna => build_codon_table(&[
            (b"AGA", b'*'),
            (b"AGG", b'*'),
            (b"ATA", b'M'),
            (b"TGA", b'W'),
        ]),
        GeneticCodes::YeastMtDna => build_codon_table(&[
            (b"ATA", b'M'),
            (b"CTT", b'T'),
            (b"CTC", b'T'),
            (b"CTA", b'T'),
            (b"CTG", b'T'),
            (b"TGA", b'W'),
        ]),
        GeneticCodes::MoldProtCoelMtDna => build_codon_table(&[(b"TGA", b'W')]),
        GeneticCodes::InvertMtDna => build_codon_table(&[
            (b"AGA", b'S'),
            (b"AGG", b'S'),
            (b"ATA", b'M'),
            (b"TGA", b'W'),
        ]),
        GeneticCodes::CilDasHexNu => build_codon_table(&[(b"TAA", b'Q'), (b"TAG", b'Q')]),
        GeneticCodes::EchiFlatwormMtDna => build_codon_table(&[
            (b"AAA", b'N'),
            (b"AGA", b'S'),
            (b"AGG", b'S'),
            (b"TGA", b'W'),
        ]),
        GeneticCodes::EuplotidNu => build_codon_table(&[(b"TGA", b'C')]),
        GeneticCodes::AltYeastNu => build_codon_table(&[(b"CTG", b'S')]),
        GeneticCodes::AsciMtDna => build_codon_table(&[
            (b"AGA", b'G'),
            (b"AGG", b'G'),
            (b"ATA", b'M'),
            (b"TGA", b'W'),
        ]),
        GeneticCodes::AltFlatwormMtDna => build_codon_table(&[
            (b"AAA", b'N'),
            (b"AGA", b'S'),
            (b"AGG", b'S'),
            (b"TAA", b'Y'),
            (b"TGA", b'W'),
        ]),
        GeneticCodes::ChlorMtDna => build_codon_table(&[(b"TAG", b'L')]),
        GeneticCodes::TrematodeMtDna => build_codon_table(&[
            (b"TGA", b'W'),
            (b"ATA", b'M'),
            (b"AGA", b'S'),
            (b"AGG", b'S'),
            (b"AAA", b'N'),
        ]),
        GeneticCodes::ScenedesmusMtDna => build_codon_table(&[(b"TCA", b'*'), (b"TAG", b'L')]),
        GeneticCodes::ThrausMtDna => build_codon_table(&[(b"TTA", b'*')]),
        GeneticCodes::RhabdopMtDna => {
            build_codon_table(&[(b"AGA", b'S'), (b"AGG", b'K'), (b"TGA", b'W')])
        }
        GeneticCodes::CaDivSR1GraciBac => build_codon_table(&[(b"TGA", b'G')]),
        GeneticCodes::PachyNu => build_codon_table(&[(b"CTG", b'A')]),
        GeneticCodes::MesodiniumNu => build_codon_table(&[(b"TAA", b'Y'), (b"TAG", b'Y')]),
        GeneticCodes::PeritrichNu => build_codon_table(&[(b"TAA", b'E'), (b"TAG", b'E')]),
        GeneticCodes::CephalodiscidaeMtDna => build_codon_table(&[
            (b"TAA", b'Y'),
            (b"TGA", b'W'),
            (b"AGA", b'S'),
            (b"AGG", b'K'),
        ]),
    }
}

/// Translate nucleotide sequence bytes into amino acid bytes using `table`,
/// starting at 1-based reading `frame` and appending into `out`.
#[inline]
pub fn translate_bytes(seq: &[u8], frame: usize, table: &CodonTable, out: &mut Vec<u8>) {
    let offset = frame.saturating_sub(1).min(seq.len());
    let slice = &seq[offset..];
    out.reserve(slice.len().div_ceil(3));
    for chunk in slice.chunks(3) {
        if chunk.len() < 3 {
            out.push(b'X');
            continue;
        }
        let aa = match chunk {
            [b'-', b'-', b'-'] => b'-',
            [b'?', b'?', b'?'] => b'?',
            [b0, b1, b2] => {
                let i0 = BASE_2BIT[*b0 as usize];
                let i1 = BASE_2BIT[*b1 as usize];
                let i2 = BASE_2BIT[*b2 as usize];
                if (i0 | i1 | i2) < 4 {
                    table[((i0 as usize) << 4) | ((i1 as usize) << 2) | (i2 as usize)]
                } else {
                    b'X'
                }
            }
            _ => unreachable!(),
        };
        out.push(aa);
    }
}

lazy_static! {
    /// A version of lazyly initialized of standard translation code.
    /// It may provide better performance than the `NcbiTables::standard_code()` method.
    pub static ref STD_CODE_MAP: HashMap<&'static str, &'static str> = {
        let mut m = HashMap::new();
        for (codon, amino_acid) in STD_CODE {
            m.insert(*codon, amino_acid.as_ref());
        }
        m
    };
}

const STD_CODE: &[(&str, &str)] = &[
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

impl Default for NcbiTables<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl NcbiTables<'_> {
    /// The constructor for the NCBI Table struct.
    /// The constructor initializes the translation tables.
    pub fn new() -> Self {
        Self {
            translation: STD_CODE,
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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
            match *codon {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codon_tables_match_ncbi_tables_exhaustively() {
        let tables = NcbiTables::new();
        let variants = [
            (GeneticCodes::StandardCode, tables.standard_code()),
            (GeneticCodes::VertMtDna, tables.vert_mtdna()),
            (GeneticCodes::YeastMtDna, tables.yeast_mtdna()),
            (GeneticCodes::MoldProtCoelMtDna, tables.moldprotocoe_mtdna()),
            (GeneticCodes::InvertMtDna, tables.invert_mtdna()),
            (GeneticCodes::CilDasHexNu, tables.cildashex_nudna()),
            (GeneticCodes::EchiFlatwormMtDna, tables.echiflatworm_mtdna()),
            (GeneticCodes::EuplotidNu, tables.euplotid_nudna()),
            (GeneticCodes::BacArchPlantPlast, tables.standard_code()),
            (GeneticCodes::AltYeastNu, tables.alt_yeast_nu()),
            (GeneticCodes::AsciMtDna, tables.ascidian_mtdna()),
            (GeneticCodes::AltFlatwormMtDna, tables.alt_flatworm_mtdna()),
            (GeneticCodes::ChlorMtDna, tables.chlorophycean_mtdna()),
            (GeneticCodes::TrematodeMtDna, tables.trematode_mtdna()),
            (GeneticCodes::ScenedesmusMtDna, tables.scenedesmus_mtdna()),
            (GeneticCodes::ThrausMtDna, tables.thraustochytrium_mtdna()),
            (GeneticCodes::RhabdopMtDna, tables.rhabdopleuridae_mtdna()),
            (
                GeneticCodes::CaDivSR1GraciBac,
                tables.candid_div_sr1_gracil(),
            ),
            (GeneticCodes::PachyNu, tables.pachysolen_tanno_nu()),
            (GeneticCodes::MesodiniumNu, tables.mesodinium_nu()),
            (GeneticCodes::PeritrichNu, tables.peritrich_nu()),
            (
                GeneticCodes::CephalodiscidaeMtDna,
                tables.cephalodiscidae_mtdna(),
            ),
        ];

        let bases = *b"ACGT";
        for (code, legacy) in variants {
            let table = codon_table(code);
            let mut out = Vec::new();
            for &b0 in &bases {
                for &b1 in &bases {
                    for &b2 in &bases {
                        let codon = [b0, b1, b2];
                        let codon_str = std::str::from_utf8(&codon).unwrap();
                        let expected = legacy.get(codon_str).unwrap().as_bytes()[0];
                        out.clear();
                        translate_bytes(&codon, 1, &table, &mut out);
                        assert_eq!(out, vec![expected], "{code:?} codon {codon_str}");

                        let lower = [
                            b0.to_ascii_lowercase(),
                            b1.to_ascii_lowercase(),
                            b2.to_ascii_lowercase(),
                        ];
                        out.clear();
                        translate_bytes(&lower, 1, &table, &mut out);
                        assert_eq!(out, vec![expected], "{code:?} lower {codon_str}");
                    }
                }
            }
            for (special, expected) in [(b"---", b'-'), (b"???", b'?'), (b"NNN", b'X')] {
                out.clear();
                translate_bytes(special, 1, &table, &mut out);
                assert_eq!(out, vec![expected], "{code:?} special");
            }
        }
    }
}
