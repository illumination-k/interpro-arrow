use anyhow::Result;
use polars::{
    io::SerReader,
    prelude::{DataFrame, IpcReader},
};
use polars_core::{utils::accumulate_dataframes_vertical, POOL};
use std::{
    collections::HashMap,
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use rayon::prelude::*;

use crate::records::Term;

pub struct PartitionedIpcReader {
    dir: PathBuf,
    org: Option<Vec<String>>,
    source: Option<Vec<String>>,
}

fn check_path(path: &Path, map: &HashMap<String, Vec<String>>) -> bool {
    let mut flags = Vec::new();

    for p in path.iter() {
        if let Some(p) = p.to_str() {
            let kv: Vec<&str> = p.split('=').collect();
            if kv.len() != 2 {
                continue;
            }

            let (k, v) = (kv[0], kv[1]);

            if let Some(expected) = map.get(k) {
                flags.push(expected.contains(&v.to_string()))
            } else {
                flags.push(true);
            }
        } else {
            continue;
        }
    }

    flags.into_iter().all(|b| b)
}

pub fn select_paths(glob: &str, map: &HashMap<String, Vec<String>>) -> Result<Vec<PathBuf>> {
    let mut ret = vec![];

    for path in glob::glob(glob)? {
        let path = path?;
        if check_path(&path, map) {
            dbg!(&path);
            ret.push(path)
        }
    }

    Ok(ret)
}

impl PartitionedIpcReader {
    pub fn new(dir: PathBuf) -> Self {
        Self {
            dir,
            org: None,
            source: None,
        }
    }

    pub fn with_org(mut self, org: Option<Vec<String>>) -> Self {
        self.org = org;
        self
    }

    pub fn with_source(mut self, source: Option<Vec<Term>>) -> Self {
        self.source = source.map(|s| s.into_iter().map(|x| x.to_string()).collect());
        self
    }

    fn map(&self) -> HashMap<String, Vec<String>> {
        let mut map = HashMap::new();

        if let Some(org) = self.org.clone() {
            map.entry("org".to_string()).or_insert(org);
        }

        if let Some(source) = self.source.clone() {
            map.entry("source".to_string()).or_insert(source);
        }

        map
    }

    pub fn finish(&self) -> Result<DataFrame> {
        let map = self.map();

        let paths = select_paths(&format!("{}/**/*.ipc", self.dir.display()), &map)?;

        let parsed_dfs = POOL.install(|| {
            paths
                .into_par_iter()
                .map(|path| {
                    let df = IpcReader::new(BufReader::new(File::open(path)?)).finish()?;
                    Ok(df)
                })
                .collect::<Result<Vec<_>>>()
        })?;

        let df = accumulate_dataframes_vertical(parsed_dfs)?;
        Ok(df)
    }
}
