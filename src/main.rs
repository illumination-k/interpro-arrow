mod gff3;
mod records;

use anyhow::Result;

fn main() -> Result<()> {
    println!("Hello, world!");
    let reader = gff3::Reader::from_path(&"tests/Olucimarinus_231_v2.0.protein.fa.gff3.gz")?;
    reader.finish()?;
    Ok(())
}
