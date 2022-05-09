/*
register
- gff3
- write ipc dir
- org
find
- roots dir
- orgs
- domain expr
- format
    - fasta
    - ids

infer source (for faster find)
*/

use std::path::PathBuf;
use structopt::{clap, clap::arg_enum, StructOpt};

#[derive(Debug, StructOpt)]
#[structopt(name = "ipsr")]
#[structopt(long_version(option_env!("LONG_VERSION").unwrap_or(env!("CARGO_PKG_VERSION"))))]
#[structopt(setting(clap::AppSettings::ColoredHelp))]
pub struct Opt {
    #[structopt(subcommand)]
    pub subcommands: SubCommands,
}

#[derive(Debug, StructOpt)]
pub enum SubCommands {
    #[structopt(name = "register", about = "register gff")]
    #[structopt(setting(clap::AppSettings::ColoredHelp))]
    Register {
        #[structopt(short = "-i", long = "input", about = "input GFF3 (can gzipped)")]
        input: PathBuf,
        #[structopt(short = "-o", long = "org", about = "Organism name")]
        org: String,
        #[structopt(short = "-d", long = "dir", about = "output dir")]
        dir: PathBuf,
    },
    #[structopt(name = "find", about = "find gene(s) which has the specific domain")]
    #[structopt(setting(clap::AppSettings::ColoredHelp))]
    Find {
        #[structopt(short = "-d", long = "dir", about = "output dir")]
        dir: PathBuf,
        #[structopt(short = "-e", long = "domain-expr", about = "domain expr")]
        expr: String,
        #[structopt(short = "-o", long = "org", about = "organism name")]
        org: Option<Vec<String>>,
        #[structopt(short = "-f", long = "fmt", about = "output format")]
        format: Option<OutFormat>,
    },
}

arg_enum! {
    #[derive(Debug)]
    pub enum OutFormat {
        Id,
        Fasta,
    }
}
