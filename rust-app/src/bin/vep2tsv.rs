use vcf::{VCFReader, U8Vec, VCFHeaderFilterAlt, VCFError, VCFRecord};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::BufReader;


fn main() -> Result<(), VCFError> {
    let input_vcf: &str = "../../resources/examples/nf.vep.vcf";

    let open_file = File::open(input_vcf)?;
    let gzdecoder = MultiGzDecoder::new(open_file);
    let reader = VCFReader::new(BufReader::new(gzdecoder))?;

    println!(reaer.header());
    
    // let bufreader = BufReader::new(inner)

    Ok(())
}

