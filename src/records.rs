use anyhow::{anyhow, Result};
use arrow2::array::ArrayRef;
use polars_core::POOL;
use strum::IntoEnumIterator;

use std::{
    collections::HashMap,
    fmt::Debug,
    fs::{self, File},
    io::BufWriter,
    path::{Path, PathBuf},
    sync::Arc,
};

use arrow2::{
    array::{Int16Array, Utf8Array},
    chunk::Chunk,
    datatypes::{DataType, Field, Schema},
    io::ipc::write::{self, Compression, FileWriter},
};
use strum_macros::{Display, EnumIter, EnumString, EnumVariantNames};

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
}

pub struct GeneRecords {
    gene_ids: Vec<String>,
    seqs: Vec<String>,
    desc: Vec<Option<String>>,
    organism: Vec<String>,
    schema: Schema,
    chunk_size: u32,
    chunks: Vec<Chunk<ArrayRef>>,
}

pub fn gene_records_schema() -> Schema {
    Schema::from(vec![
        Field::new("gene_id", DataType::Utf8, false),
        Field::new("seq", DataType::Utf8, false),
        Field::new("desc", DataType::Utf8, true),
        Field::new("organism", DataType::Utf8, false),
    ])
}

impl GeneRecords {
    pub fn new(chunk_size: u32) -> Self {
        let schema = gene_records_schema();
        Self {
            gene_ids: Vec::new(),
            seqs: Vec::new(),
            desc: Vec::new(),
            organism: Vec::new(),
            schema,
            chunk_size,
            chunks: Vec::new(),
        }
    }

    pub fn push(
        &mut self,
        gene_id: String,
        seq: String,
        desc: Option<String>,
        organism: String,
    ) -> Result<()> {
        self.rechunk()?;
        self.gene_ids.push(gene_id);
        self.seqs.push(seq);
        self.desc.push(desc);
        self.organism.push(organism);
        Ok(())
    }
    pub fn len(&self) -> usize {
        self.gene_ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn to_chunk(&self) -> Result<Chunk<ArrayRef>> {
        Ok(Chunk::try_new(vec![
            Arc::new(Utf8Array::<i32>::from_slice(&self.gene_ids)) as ArrayRef,
            Arc::new(Utf8Array::<i32>::from_slice(&self.seqs)) as ArrayRef,
            Arc::new(Utf8Array::<i32>::from(&self.desc)) as ArrayRef,
            Arc::new(Utf8Array::<i32>::from_slice(&self.organism)) as ArrayRef,
        ])?)
    }

    fn init(&mut self) {
        self.gene_ids = Vec::new();
        self.seqs = Vec::new();
        self.desc = Vec::new();
        self.organism = Vec::new()
    }

    fn finish(&mut self) -> Result<()> {
        if !self.is_empty() {
            let chunk = self.to_chunk()?;
            self.chunks.push(chunk);
            self.init();
        }
        Ok(())
    }

    fn rechunk(&mut self) -> Result<()> {
        if self.chunk_size as usize != self.len() {
            return Ok(());
        }

        let chunk = self.to_chunk()?;

        self.chunks.push(chunk);

        self.init();
        Ok(())
    }

    pub fn write(mut self, path: &Path) -> Result<()> {
        self.finish()?;
        let file = File::create(path)?;
        let options = write::WriteOptions {
            compression: Some(Compression::LZ4),
        };

        let mut writer = FileWriter::try_new(BufWriter::new(file), &self.schema, None, options)?;
        // writer.start()?;

        for chunk in self.chunks.iter() {
            writer.write(chunk, None)?;
        }

        writer.finish()?;

        Ok(())
    }
}

pub fn domain_record_schema() -> Schema {
    Schema::from(vec![
        Field::new("start", DataType::Int16, false),
        Field::new("end", DataType::Int16, false),
        Field::new("domain_name", DataType::Utf8, false),
        Field::new("domain_desc", DataType::Utf8, true),
        Field::new("gene_id", DataType::Utf8, false),
    ])
}

struct DomainRecord {
    starts: Vec<i16>,
    ends: Vec<i16>,
    domain_names: Vec<String>,
    domain_descs: Vec<Option<String>>,
    gene_ids: Vec<String>,
}

impl DomainRecord {
    fn new() -> Self {
        Self {
            starts: Vec::new(),
            ends: Vec::new(),
            domain_names: Vec::new(),
            domain_descs: Vec::new(),
            gene_ids: Vec::new(),
        }
    }

    fn push(
        &mut self,
        start: i16,
        end: i16,
        domain_name: String,
        domain_desc: Option<String>,
        gene_id: String,
    ) {
        self.starts.push(start);
        self.ends.push(end);
        self.domain_names.push(domain_name);
        self.domain_descs.push(domain_desc);
        self.gene_ids.push(gene_id);
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn len(&self) -> usize {
        self.gene_ids.len()
    }

    fn to_chunk(&self) -> Result<Chunk<ArrayRef>> {
        Ok(Chunk::try_new(vec![
            Arc::new(Int16Array::from_slice(&self.starts)) as ArrayRef,
            Arc::new(Int16Array::from_slice(&self.ends)) as ArrayRef,
            Arc::new(Utf8Array::<i32>::from_slice(&self.domain_names)) as ArrayRef,
            Arc::new(Utf8Array::<i32>::from(&self.domain_descs)) as ArrayRef,
            Arc::new(Utf8Array::<i32>::from_slice(&self.gene_ids)) as ArrayRef,
        ])?)
    }
}

impl Default for DomainRecord {
    fn default() -> Self {
        Self::new()
    }
}

pub struct DomainRecords {
    sources: HashMap<Term, DomainRecord>,
    schema: Schema,
    chunks: HashMap<Term, Vec<Chunk<ArrayRef>>>,
    chunk_size: u32,
}

impl DomainRecords {
    pub fn new(chunk_size: u32) -> Self {
        Self {
            sources: Term::iter().map(|t| (t, DomainRecord::default())).collect(),
            schema: domain_record_schema(),
            chunks: Term::iter().map(|t| (t, Vec::new())).collect(),
            chunk_size,
        }
    }

    fn rechunk(&mut self) -> Result<()> {
        for (k, v) in self.sources.iter_mut() {
            if self.chunk_size as usize != v.len() {
                continue;
            }
            let chunk = v.to_chunk()?;

            self.chunks.get_mut(k).unwrap().push(chunk);
            *v = DomainRecord::default();
        }
        Ok(())
    }

    pub fn push(
        &mut self,
        source: Term,
        start: i16,
        end: i16,
        domain_name: String,
        domain_desc: Option<String>,
        gene_id: String,
    ) -> Result<()> {
        if let Some(x) = self.sources.get_mut(&source) {
            x.push(start, end, domain_name, domain_desc, gene_id)
        };

        self.rechunk()?;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        for (k, v) in self.sources.iter_mut() {
            if !v.is_empty() {
                let chunk = v.to_chunk()?;

                self.chunks.get_mut(k).unwrap().push(chunk);
                *v = DomainRecord::default();
            }
        }
        Ok(())
    }

    pub fn write(mut self, dir: PathBuf) -> Result<()> {
        self.finish()?;
        println!("------ {} ------", dir.display());
        POOL.install(|| {
            self.chunks
                .into_iter()
                .map(|(source, batches)| {
                    let mut path = PathBuf::from(&dir);
                    path.push(format!("source={}", source));

                    fs::create_dir_all(&path)?;

                    path.push(format!("{}.ipc", uuid::Uuid::new_v4()));
                    let file = File::create(path)?;
                    let options = write::WriteOptions {
                        compression: Some(Compression::LZ4),
                    };
                    let mut writer =
                        FileWriter::try_new(BufWriter::new(file), &self.schema, None, options)?;

                    for chunk in batches.iter() {
                        writer.write(chunk, None)?;
                    }
                    writer.finish()?;
                    Ok(())
                })
                .collect::<Result<Vec<_>>>()
        })?;

        Ok(())
    }
}
