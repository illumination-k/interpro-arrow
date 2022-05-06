mod gff3;
mod records;
mod args;
mod parser;

use anyhow::Result;

fn main() -> Result<()> {
    println!("Hello, world!");
    let reader = gff3::Reader::from_path(&"tests/Olucimarinus_231_v2.0.protein.fa.gff3.gz")?;
    let (dr, gr) = reader.finish()?;
    dr.write("tests/Olucimarinus_231_v2.0.protein.ipc")?;
    gr.write( "tests/genes.ipc")?;
    Ok(())
}
