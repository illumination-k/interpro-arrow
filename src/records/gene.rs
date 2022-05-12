use anyhow::Result;
use arrow2::array::ArrayRef;

use std::{fs::File, io::BufWriter, path::Path, sync::Arc};

use arrow2::{
    array::{Int16Array, Utf8Array},
    chunk::Chunk,
    datatypes::{DataType, Field, Schema},
    io::ipc::write::{self, Compression, FileWriter},
};

pub struct GeneRecords {
    gene_ids: Vec<String>,
    seqs: Vec<String>,
    lengths: Vec<i16>,
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
        Field::new("length", DataType::Int16, false),
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
            lengths: Vec::new(),
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
        let length = seq.len() as i16;
        self.gene_ids.push(gene_id);
        self.seqs.push(seq);
        self.lengths.push(length);
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
            Arc::new(Int16Array::from_slice(&self.lengths)) as ArrayRef,
            Arc::new(Utf8Array::<i32>::from(&self.desc)) as ArrayRef,
            Arc::new(Utf8Array::<i32>::from_slice(&self.organism)) as ArrayRef,
        ])?)
    }

    fn init(&mut self) {
        self.gene_ids = Vec::new();
        self.seqs = Vec::new();
        self.lengths = Vec::new();
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
