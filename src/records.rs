use anyhow::Result;
use std::{fmt::Debug, fs::File, sync::Arc};

use arrow2::{
    array::{Array, Int16Array, Utf8Array},
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

pub struct GeneRecords {
    gene_ids: Vec<String>,
    seqs: Vec<String>,
    desc: Vec<Option<String>>,
    organism: Vec<String>,
    schema: Schema,
    chunk_size: u32,
    chunks: Vec<Chunk<Arc<dyn Array>>>,
}

impl Debug for GeneRecords {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let table = arrow2::io::print::write(
            &self.chunks,
            &self
                .schema
                .fields
                .iter()
                .map(|f| f.name.to_string())
                .collect::<Vec<String>>(),
        );
        f.write_str(&table)?;
        Ok(())
    }
}


impl GeneRecords {
    pub fn new(chunk_size: u32) -> Self {
        let schema = Schema::from(vec![
            Field::new("gene_id", DataType::Utf8, false),
            Field::new("seq", DataType::Utf8, false),
            Field::new("desc", DataType::Utf8, true),
            Field::new("organism", DataType::Utf8, false),
        ]);

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

    fn to_chunk(&self) -> Result<Chunk<Arc<dyn Array>>> {
        Ok(Chunk::try_new(vec![
            Arc::new(Utf8Array::<i32>::from_slice(&self.gene_ids)) as Arc<dyn Array>,
            Arc::new(Utf8Array::<i32>::from_slice(&self.seqs)) as Arc<dyn Array>,
            Arc::new(Utf8Array::<i32>::from(&self.desc)) as Arc<dyn Array>,
            Arc::new(Utf8Array::<i32>::from_slice(&self.organism)) as Arc<dyn Array>,
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
        if !self.chunk_size as usize == self.len() {
            return Ok(());
        }

        let chunk = self.to_chunk()?;

        self.chunks.push(chunk);

        self.init();
        Ok(())
    }
}

pub struct DomainRecords {
    sources: Vec<String>,
    starts: Vec<i16>,
    ends: Vec<i16>,
    domain_names: Vec<String>,
    domain_descs: Vec<Option<String>>,
    gene_ids: Vec<String>,
    schema: Schema,
    chunks: Vec<Chunk<Arc<dyn Array>>>,
    chunk_size: u32,
}

impl Debug for DomainRecords {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let table = arrow2::io::print::write(
            &self.chunks,
            &self
                .schema
                .fields
                .iter()
                .map(|f| f.name.to_string())
                .collect::<Vec<String>>(),
        );
        f.write_str(&table)?;
        Ok(())
    }
}

pub fn domain_records_schema() -> Schema {
    Schema::from(vec![
        Field::new("source", DataType::Utf8, false),
        Field::new("start", DataType::Int16, false),
        Field::new("end", DataType::Int16, false),
        Field::new("domain_name", DataType::Utf8, false),
        Field::new("domain_descs", DataType::Utf8, true),
        Field::new("gene_id", DataType::Utf8, false),
    ])
}

impl DomainRecords {
    pub fn new(chunk_size: u32) -> Self {
        Self {
            sources: Vec::new(),
            starts: Vec::new(),
            ends: Vec::new(),
            domain_names: Vec::new(),
            domain_descs: Vec::new(),
            gene_ids: Vec::new(),
            schema: domain_records_schema(),
            chunks: Vec::new(),
            chunk_size,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn len(&self) -> usize {
        self.sources.len()
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
        self.rechunk()?;

        self.sources.push(source.to_string());
        self.starts.push(start);
        self.ends.push(end);
        self.domain_names.push(domain_name);
        self.domain_descs.push(domain_desc);
        self.gene_ids.push(gene_id);

        Ok(())
    }

    fn init(&mut self) {
        self.sources = Vec::new();
        self.starts = Vec::new();
        self.ends = Vec::new();
        self.domain_names = Vec::new();
        self.domain_descs = Vec::new();
        self.gene_ids = Vec::new();
    }

    fn to_chunk(&self) -> Result<Chunk<Arc<dyn Array>>> {
        Ok(Chunk::try_new(vec![
            Arc::new(Utf8Array::<i32>::from_slice(&self.sources)) as Arc<dyn Array>,
            Arc::new(Int16Array::from_slice(&self.starts)) as Arc<dyn Array>,
            Arc::new(Int16Array::from_slice(&self.ends)) as Arc<dyn Array>,
            Arc::new(Utf8Array::<i32>::from_slice(&self.domain_names)) as Arc<dyn Array>,
            Arc::new(Utf8Array::<i32>::from(&self.domain_descs)) as Arc<dyn Array>,
            Arc::new(Utf8Array::<i32>::from_slice(&self.gene_ids)) as Arc<dyn Array>,
        ])?)
    }

    fn rechunk(&mut self) -> Result<()> {
        if !self.chunk_size as usize == self.len() {
            return Ok(());
        }

        let chunk = self.to_chunk()?;

        self.chunks.push(chunk);

        self.init();

        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        if !self.is_empty() {
            let chunk = self.to_chunk()?;
            self.chunks.push(chunk);
            self.init();
        }
        Ok(())
    }

    pub fn write(mut self, path: &str) -> Result<()> {
        self.finish()?;
        let file = File::create(path)?;
        let options = write::WriteOptions {
            compression: Some(Compression::LZ4),
        };

        let mut writer = FileWriter::new(file, self.schema, None, options);
        writer.start()?;

        for chunk in self.chunks.iter() {
            writer.write(chunk, None)?;
        }

        writer.finish()?;

        Ok(())
    }
}
