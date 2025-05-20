use tsv_utils::bed::{BedRecord, BedRecordVecEx};
use tsv_utils::peak::{PeakRecord};
use tsv_utils::tsv::{TsvRecord, read_tsv_file, read_tsv_gz_file};
use std::fs::File;

fn main() {

//    let bed_file = "./data/sample/sample.bed";
//    let bed = read_tsv_file::<BedRecord>(bed_file).unwrap();


    // get all .bed.gz file from ./data/peak
    let bed_gz_files = glob::glob("./data/peak/*.bed.gz").unwrap();
    for bed_gz_file in bed_gz_files {
        let bed_gz_file = bed_gz_file.unwrap();
        println!("{:?}", bed_gz_file);

        let bed = read_tsv_gz_file::<PeakRecord>(bed_gz_file.to_str().unwrap()).unwrap();
        bed.iter().for_each(|record| {
            println!("{:?}", record);
        });
    }

}
