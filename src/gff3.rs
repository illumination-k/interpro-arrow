use std::{
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
    str::FromStr,
};

use anyhow::{anyhow, Result};
use flate2::read::MultiGzDecoder;

use crate::records::{DomainRecords, Term};

fn is_compressed<P: AsRef<Path>>(p: &P) -> bool {
    let ext = p.as_ref().extension();

    ext == Some(OsStr::new("gz"))
}

pub fn read_with_gz<P: AsRef<Path>>(p: &P) -> Result<Box<dyn BufRead>> {
    let file = File::open(p)?;
    let reader: Box<dyn BufRead> = if is_compressed(p) {
        let gz = MultiGzDecoder::new(file);
        Box::new(BufReader::new(gz))
    } else {
        Box::new(BufReader::new(file))
    };

    Ok(reader)
}

pub fn parse_gffrecord_line(line: &str, domain_records: &mut DomainRecords) -> Result<()> {
    let line = line.trim();

    let records: Vec<&str> = line.split('\t').collect();
    if records.len() != 9 {
        return Err(anyhow!(format!("Invalid line: {}", line)));
    }

    let id = records[0];

    let source = Term::from_str(records[1])?;

    if source == Term::ID {
        return Ok(());
    }

    let start: i16 = records[3].parse()?;
    let end: i16 = records[4].parse()?;

    let mut domain_name = None;
    let mut domain_desc = None;
    for attr in records[8].split(';') {
        let attr_records: Vec<&str> = attr.split('=').collect();

        if attr_records.len() != 2 {
            continue;
        }

        let (key, values) = (attr_records[0], attr_records[1]);

        if key == "Name" {
            domain_name = Some(values);
        } else if key == "signature_desc" {
            domain_desc = Some(values)
        } else if key == "Ontology_term" {
            // Ontology_term="GO:0008654","GO:0016020","GO:0016780"
            for term_name in values.split(",") {
                let term_name = term_name.trim_matches('\"');
                domain_records.push(
                    Term::GoTerm,
                    start,
                    end,
                    term_name.to_string(),
                    None,
                    id.to_string(),
                )?;
            }
        } else if attr_records[0] == "Dbxref" {
            for type_term in values.split(",") {
                let type_term: Vec<&str> = type_term.trim_matches('\"').split(":").collect();
                let (type_, term_name) = (type_term[0], type_term[1]);

                match type_ {
                    "InterPro" => domain_records.push(
                        Term::InterPro,
                        start,
                        end,
                        term_name.to_string(),
                        None,
                        id.to_string(),
                    )?,
                    "MetaCyc" => domain_records.push(
                        Term::MetaCyc,
                        start,
                        end,
                        term_name.to_string(),
                        None,
                        id.to_string(),
                    )?,

                    "Reactome" => domain_records.push(
                        Term::Reactome,
                        start,
                        end,
                        term_name.to_string(),
                        None,
                        id.to_string(),
                    )?,
                    _ => unreachable!(),
                }
            }
        }
    }

    if let Some(domain_name) = domain_name {
        domain_records.push(
            source,
            start,
            end,
            domain_name.to_string(),
            domain_desc.map(|s| s.to_string()),
            id.to_string(),
        )?;

        Ok(())
    } else {
        Err(anyhow!("domain name is required"))
    }
}
pub struct Reader {
    reader: Box<dyn BufRead>,
}

impl Reader {
    pub fn from_path<P: AsRef<Path>>(p: &P) -> Result<Self> {
        Ok(Self {
            reader: read_with_gz(p)?,
        })
    }

    pub fn finish(self) -> Result<()> {
        let comment = '#';
        let fasta_line = "##FASTA";
        
        let mut domain_records = DomainRecords::new(5000);
        
        let mut is_fasta = false;
        for line in self.reader.lines() {
            let line = line?;
            if line.starts_with(fasta_line) {
                is_fasta = true;
                continue;
            }

            if line.starts_with(comment) {
                continue;
            }

            if line.len() <= 1 {
                continue;
            }

            if !is_fasta {
                parse_gffrecord_line(&line, &mut domain_records)?;
            }
        }

        dbg!(&domain_records);


        Ok(())
    }
}
