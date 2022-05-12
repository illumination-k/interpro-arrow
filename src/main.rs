mod args;
mod gff3;
mod parser;
mod partition;
mod records;

use std::fs;

use anyhow::{anyhow, Result};
use args::OutFormat;
use polars::{
    chunked_array::ChunkedArray,
    datatypes::{BooleanType, Utf8Chunked},
};
use structopt::StructOpt;

use crate::{
    args::{Opt, SubCommands},
    parser::Expr,
    partition::PartitionedIpcReader,
    records::Term,
};

fn main() -> Result<()> {
    let opt = Opt::from_args();

    match &opt.subcommands {
        SubCommands::Register { input, org, dir } => {
            let orgname = format!("org={}", org);
            let domain_dir = dir.join("domain").join(&orgname);

            if domain_dir.exists() {
                return Err(anyhow!(format!("{} is already registered", org)));
            }
            let (dr, gr) = gff3::Reader::from_path(input)?.finish()?;
            let orgname = format!("org={}", org);
            let domain_dir = dir.join("domain").join(&orgname);

            dr.write(domain_dir)?;

            let mut gene_path = dir.join("gene").join(&orgname);

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
            let sources = Term::try_from_expr(expr)?;

            let mut domain_dir = dir.clone();
            domain_dir.push("domain");
            let domain_df = PartitionedIpcReader::new(domain_dir)
                .with_org(org.to_owned())
                .with_source(Some(sources))
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

            let df = domain_df.filter(&mask)?;

            match format {
                OutFormat::Id => {
                    println!(
                        "{}",
                        df["gene_id"]
                            .0
                            .utf8()?
                            .into_iter()
                            .flatten()
                            .collect::<Vec<&str>>()
                            .join("\n")
                    )
                }
                OutFormat::Fasta => {
                    let gene_dir = dir.join("gene");

                    let gene_df = PartitionedIpcReader::new(gene_dir)
                        .with_org(org.to_owned())
                        .finish()?;
                    let gene_df = df
                        .join(
                            &gene_df,
                            ["gene_id"],
                            ["gene_id"],
                            polars::prelude::JoinType::Inner,
                            None,
                        )?
                        .select(["gene_id", "seq"])?;
                    let len = gene_df.height();
                    let header = &Utf8Chunked::from_iter(std::iter::repeat(">").take(len))
                        + gene_df["gene_id"].utf8()?;
                    let header_with_n =
                        header + Utf8Chunked::from_iter(std::iter::repeat("\n").take(len));
                    let fasta = &header_with_n + gene_df["seq"].utf8()?;

                    println!(
                        "{}",
                        fasta
                            .into_iter()
                            .flatten()
                            .collect::<Vec<&str>>()
                            .join("\n")
                    );
                }
            };
        }
    };

    Ok(())
}
