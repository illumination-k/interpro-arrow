mod args;
mod gff3;
mod parser;
mod partition;
mod records;

use std::path::PathBuf;

use anyhow::Result;
use polars::{
    chunked_array::ChunkedArray,
    datatypes::BooleanType,
    io::SerReader,
    prelude::{ChunkApply, IpcReader},
};

use crate::{parser::Expr, partition::PartitionedIpcReader, records::Term};

fn main() -> Result<()> {
    println!("Hello, world!");
    // let reader = gff3::Reader::from_path(&"tests/Olucimarinus_231_v2.0.protein.fa.gff3.gz")?;
    // let (dr, gr) = reader.finish()?;
    // dr.write(PathBuf::from("tmp/domain/org=Olucimarinus"))?;
    // gr.write("tests/genes.ipc")?;

    let df = PartitionedIpcReader::new(PathBuf::from("tmp/domain"))
        .with_org(Some(vec!["Olucimarinus".to_string()]))
        .with_source(Some(vec![Term::GoTerm, Term::Pfam]))
        .finish()?
        .groupby(["gene_id"])?
        .agg_list()?;

    dbg!(&df);

    let expr = Expr::from_string("GO:0004497").unwrap();

    let mask: ChunkedArray<BooleanType> = df["domain_name_agg_list"]
        .list()?
        .into_iter()
        .map(|l| {
            if let Some(s) = l {
                let s: Vec<&str> = s.0.utf8().unwrap().into_iter().flatten().collect();
                let bool = expr.matches(&s).unwrap();
                Ok(bool)
            } else {
                Ok(false)
            }
        })
        .collect::<Result<_>>()?;

    dbg!(&df.filter(&mask));
    Ok(())
}
