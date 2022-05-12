use anyhow::{anyhow, Result};

use strum_macros::{Display, EnumIter, EnumString, EnumVariantNames};

use crate::parser::lex::{lex, Token};

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, PartialEq, Eq, Hash, EnumString, EnumVariantNames, EnumIter, Display)]
pub enum Term {
    #[strum(serialize = ".")]
    ID,
    CDD,
    Coils,
    Gene3D,
    MobiDBLite,
    PANTHER,
    Pfam,
    PIRSF,
    PIRSR,
    PRINTS,
    ProSitePatterns,
    ProSiteProfiles,
    SFLD,
    SMART,
    SUPERFAMILY,
    TIGRFAM,
    GoTerm,
    Reactome,
    MetaCyc,
    InterPro,
}

impl Term {
    pub fn try_infer(name: &str) -> Result<Term> {
        let term = match name {
            "mobidb-lite" => Term::MobiDBLite,
            // C
            n if n.starts_with("cd") => Term::CDD,
            // G
            n if n.starts_with("G3DSA") => Term::Gene3D,
            n if n.starts_with("GO") => Term::GoTerm,
            // I
            n if n.starts_with("IRR") => Term::InterPro,
            // P
            n if n.starts_with("PTHR") => Term::PANTHER,
            n if n.starts_with("PWY") => Term::MetaCyc,
            n if n.starts_with("PS") => Term::ProSiteProfiles,
            n if n.starts_with("PF") => Term::Pfam,
            // S
            n if n.starts_with("SSF") => Term::SUPERFAMILY,
            n if n.starts_with("SM") => Term::SMART,
            // T
            n if n.starts_with("TIGR") => Term::TIGRFAM,
            // R
            n if n.starts_with("R-") => Term::Reactome,

            _ => return Err(anyhow!(format!("No term is match: {}", name))),
        };

        Ok(term)
    }

    pub fn try_from_tokens(tokens: &[Token]) -> Vec<Self> {
        tokens
            .into_iter()
            .filter_map(|token| {
                let term = match token {
                    Token::Name(name) => Term::try_infer(&name).ok(),
                    _ => None,
                };
                term
            })
            .collect()
    }

    pub fn try_from_expr(expr: &str) -> Result<Vec<Self>> {
        let tokens = lex(&expr)?;
        Ok(Self::try_from_tokens(&tokens))
    }
}
