mod args;
mod gff3;
mod parser;
mod partition;
mod records;

use std::fs;

use anyhow::Result;
use args::OutFormat;
use polars::{chunked_array::ChunkedArray, datatypes::BooleanType};
use structopt::StructOpt;

use crate::{
    args::{Opt, SubCommands},
    parser::{
        lex::{lex, Token},
        Expr,
    },
    partition::PartitionedIpcReader,
    records::Term,
};

fn main() -> Result<()> {
    let opt = Opt::from_args();

    match &opt.subcommands {
        SubCommands::Register { input, org, dir } => {
            let (dr, gr) = gff3::Reader::from_path(input)?.finish()?;
            let orgname = format!("org={}", org);
            let mut domain_dir = dir.clone();
            domain_dir.push("domain");
            domain_dir.push(&orgname);
            dr.write(domain_dir)?;

            let mut gene_path = dir.clone();
            gene_path.push("gene");
            gene_path.push(&orgname);
            fs::create_dir_all(&gene_path)?;
            gene_path.push(format!("{}.ipc", uuid::Uuid::new_v4()));
            gr.write(&gene_path)?;
        }
        SubCommands::Find {
            dir,
            expr,
            org,
            format,
        } => {
            let format = format.as_ref().unwrap_or(&args::OutFormat::Id);
            let tokens = lex(expr)?;
            let mut source = vec![];
            for token in tokens.iter() {
                match token {
                    Token::Name(s) => source.push(Term::try_infer(s)?),
                    _ => {}
                }
            }

            let mut domain_dir = dir.clone();
            domain_dir.push("domain");
            let domain_df = PartitionedIpcReader::new(domain_dir)
                .with_org(org.to_owned())
                .with_source(Some(source))
                .finish()?
                .groupby(["gene_id"])?
                .agg_list()?;
            let expr = Expr::from_string(expr).unwrap();
            let mask: ChunkedArray<BooleanType> = domain_df["domain_name_agg_list"]
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

            match format {
                OutFormat::Id => {
                    let df = domain_df.filter(&mask)?;
                    println!(
                        "{}",
                        df["gene_id"]
                            .0
                            .utf8()
                            .unwrap()
                            .into_iter()
                            .flatten()
                            .collect::<Vec<&str>>()
                            .join("\n")
                    )
                }
                OutFormat::Fasta => {
                    let mut gene_dir = dir.clone();
                    gene_dir.push("gene");
                    let gene_df = PartitionedIpcReader::new(gene_dir).with_org(org.to_owned()).finish()?;
                    dbg!(&gene_df);
                },
            };
        }
    };

    Ok(())
}
