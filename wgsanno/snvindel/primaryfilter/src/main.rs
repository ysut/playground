use csv::ReaderBuilder;
use std::collections::HashMap;
use std::env;
use serde::{Serialize, Deserialize};
use std::error::Error;
use std::fs::File;


#[derive(Debug, Serialize, Deserialize)]
struct Thresholds {
    gnomad40_genome_af: f64,
    gnomad40_genome_af_afr: f64,
    gnomad40_genome_af_ami: f64,
    gnomad40_genome_af_amr: f64,
    gnomad40_genome_af_asj: f64,
    gnomad40_genome_af_eas: f64,
    gnomad40_genome_af_fin: f64,
    gnomad40_genome_af_mid: f64,
    gnomad40_genome_af_nfe: f64,
    gnomad40_genome_af_remaining: f64,
    gnomad40_genome_af_sas: f64,
    all_sites_2015_08: f64,
    hrc_af: f64,
    hrc_non1000g_af: f64,
    kaviar_af: f64,
    gme_af: f64,
    gme_nwa: f64,
    gme_nea: f64,
    gme_ap: f64,
    gme_israel: f64,
    gme_sd: f64,
    gme_tp: f64,
    gme_ca: f64
}


fn read_thresholds(yaml_file: &str) -> Result<Thresholds, Box<dyn Error>> {
    let file = File::open(yaml_file)?;
    let thresholds: Thresholds = serde_yaml::from_reader(file)?;
    Ok(thresholds)
}


fn parse_value_or_min(value: Option<&String>) -> f64 {
    value
    .and_then(|v| v.parse::<f64>().ok())
    .unwrap_or(f64::MIN)
}


fn is_lower_threshold_or_min(value: f64, threshold: f64) -> bool {
    value < threshold || value == f64::MIN
}


fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <inital_genome_summary.txt> <config.yaml>", args[0]); 
        std::process::exit(1);
    }

    let input_file = &args[1];
    let yaml_file = &args[2];
    
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_reader(File::open(input_file)?);

    // Get column names (header) from input genome summary.
    let headers = rdr.headers()?.clone();

    // Read thresholds from config.yaml and store them in Thresholds.
    let thresholds = read_thresholds(yaml_file)?;
    
    for result in rdr.records() {
        // Pick up each record and store it in HashMap.
        let record = result?;
        let mut map: HashMap<String, String> = HashMap::new();

        for (i, header) in headers.iter().enumerate() {
            map.insert(
                header.to_string(),
                record.get(i).unwrap_or(".").to_string(),
            );
        }
        
        // Intialize all_conditions_met as true.
        let mut all_conditions_met: bool = true;

        // The set of column names and thresholds are iterated over.
        for (key, threshold) in [
            ("gnomad40_genome_AF", thresholds.gnomad40_genome_af),
            ("gnomad40_genome_AF_afr", thresholds.gnomad40_genome_af_afr),
            ("gnomad40_genome_AF_ami", thresholds.gnomad40_genome_af_ami),
            ("gnomad40_genome_AF_amr", thresholds.gnomad40_genome_af_amr),
            ("gnomad40_genome_AF_asj", thresholds.gnomad40_genome_af_asj),
            ("gnomad40_genome_AF_eas", thresholds.gnomad40_genome_af_eas),
            ("gnomad40_genome_AF_fin", thresholds.gnomad40_genome_af_fin),
            ("gnomad40_genome_AF_mid", thresholds.gnomad40_genome_af_mid),
            ("gnomad40_genome_AF_nfe", thresholds.gnomad40_genome_af_nfe),
            ("gnomad40_genome_AF_remaining", thresholds.gnomad40_genome_af_remaining),
            ("gnomad40_genome_AF_sas", thresholds.gnomad40_genome_af_sas),
            ("ALL.sites.2015_08", thresholds.all_sites_2015_08),
            ("HRC_AF", thresholds.hrc_af),
            ("HRC_non1000G_AF", thresholds.hrc_non1000g_af),
            ("Kaviar_AF", thresholds.kaviar_af),
            ("GME_AF", thresholds.gme_af),
            ("gme_NWA", thresholds.gme_nwa),
            ("gme_NEA", thresholds.gme_nea),
            ("gme_AP", thresholds.gme_ap),
            ("GME_Israel", thresholds.gme_israel),
            ("GME_SD", thresholds.gme_sd),
            ("GME_TP", thresholds.gme_tp),
            ("GME_CA", thresholds.gme_ca)
            ] {
                // Get the value from HashMap.
                let value: f64 = parse_value_or_min(map.get(key));
                
                /*
                1. The above column name and threshold pairs are enter 
                   the function, named "is_lower_threshold_or_min()".
                   Then, check the value whether it is lower than 
                   the threshold or not.
                2. When the value is f64::MIN, it is also lower 
                   than the threshold.
                3. If the all values are lower than the threshold, 
                   "all_conditions_met" will not be changed (remain true).
                */

                if !is_lower_threshold_or_min(value, threshold) {
                    /*
                    If any value is higher than the threshold,
                    "all_conditions_met" will be changed to false.
                    Then, the record will not be printed to stdout
                    and the next record will be picked up
                    (This loop will be broken.).
                    */
                    all_conditions_met = false;
                    break;
                }
            }

        if all_conditions_met {
            let line: String= record.iter().collect::<Vec<&str>>().join("\t");
            println!("{}", line);
        }
    }
    Ok(())
}

