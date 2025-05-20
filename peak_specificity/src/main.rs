use serde::Deserialize;
use std::fs::File;
use std::io::{BufReader, Write, BufRead};
use std::process::Command;
use std::error::Error;
use std::path::{PathBuf, Path};
use tsv_utils::bed::{BedRecord, BedRecordVecEx};
use tsv_utils::peak::{PeakRecord};
use tsv_utils::tsv::{TsvRecord, read_tsv_file, read_tsv_gz_file, write_tsv_file};

#[derive(Debug, Deserialize)]
struct Setting {
    target_peak_directory: PathBuf,
    target_peak_file: PathBuf,
    source_bigwig_directory: PathBuf,
    output_dir_base: PathBuf,
    log_dir: PathBuf,
}

fn read_setting(path: &str) -> Result<Setting, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let setting: Setting = serde_yaml::from_reader(reader)?;
    Ok(setting)
}

fn main() -> Result<(), Box<dyn Error>> {
    let setting = read_setting("setting.yaml").unwrap();
    println!("{:#?}", setting.target_peak_file);

    // eg. Follicular
    let target_str : String = setting.target_peak_file
        .file_stem()  // 1回目: "Follicular.bed"
        .and_then(|s| Path::new(s).file_stem()) // 2回目: "Follicular"
        .unwrap()
        .to_string_lossy()
        .to_string();
    let output_dir : PathBuf = setting.output_dir_base.join(target_str);

    if !output_dir.exists() {
        std::fs::create_dir_all(&output_dir).unwrap();
    }
    let output_csv_path = output_dir.join("merged_depth.csv");
    let output_zscore_path = output_dir.join("depth_zscore.csv");
    let region_bed_path = output_dir.join("regions.bed");


    let bed_gz_files : Vec<PeakRecord> = read_tsv_gz_file::<PeakRecord>(&setting.target_peak_directory.join(&setting.target_peak_file)).unwrap();

    let bed_records: Vec<BedRecord> = bed_gz_files.iter()
    .filter(|record| {
        record.chrom.starts_with("chr")
    })
    .filter_map(|record| {
        BedRecord::from_fields(record.to_fields())
    }).collect();
    write_tsv_file(&region_bed_path, &bed_records).unwrap();


    for cell_type_file in glob::glob(&setting.source_bigwig_directory.join("*.bw").to_string_lossy()).unwrap() {
        let cell_type_file = cell_type_file.unwrap();
        let cell_type_file_stem = cell_type_file.file_stem().unwrap().to_string_lossy().to_string();
        let output_file = output_dir.join(format!("{}_depth.tsv", &cell_type_file_stem));
        if output_file.exists() {
            println!("File already exists: {:?}", output_file);
            continue;
        }

        println!("Megadepth start! {:?}", cell_type_file);

        let output = Command::new("./megadepth")
            .arg(&cell_type_file)
            .arg("--annotation")
            .arg(output_dir.join("regions.bed"))
            .output()
            .expect("failed to execute megadepth");

        if !output.status.success() {
            eprintln!("Megadepth failed:\n{}", String::from_utf8_lossy(&output.stderr));
            return Ok(());
        }

        // 出力をファイルに保存
        let mut outputfile = File::create(output_file)?;
        outputfile.write_all(&output.stdout)?;
        println!("{:?} Finished", cell_type_file);
    }

    println!("Megadepth finished successfully. {}", String::from_utf8_lossy("Follicular_depth.tsv".as_bytes()));


    // *_depth.tsv ファイル群を読み込み
    let depth_pattern = output_dir.join("*_depth.tsv").to_string_lossy().to_string();
    let mut readers: Vec<_> = glob::glob(&depth_pattern)?
        .filter_map(Result::ok)
        .map(|path| {
            let stem = path.file_stem().unwrap().to_string_lossy().to_string();
            let lines = BufReader::new(File::open(path).unwrap()).lines();
            (stem.strip_suffix("_depth").unwrap().to_string(), lines)
        })
        .collect();

    let mut csv_writer = File::create(&output_csv_path)?;
    let mut zscore_writer = File::create(&output_zscore_path)?;

    let mut header = vec!["chr", "start", "end"].into_iter().map(String::from).collect::<Vec<_>>();
    
    header.extend(readers.iter().map(|(name, _)| name.clone()));
    writeln!(csv_writer, "{}", header.join(","))?;

    // データ行処理（行ごとに全ファイルを横断）
    loop {
        let mut row_values = Vec::new();
        let mut chrom = String::new();
        let mut start = String::new();
        let mut end = String::new();
        let mut all_done = true;

        for (idx, (_, reader)) in readers.iter_mut().enumerate() {
            if let Some(Ok(line)) = reader.next() {
                all_done = false;
                let fields: Vec<&str> = line.split('\t').collect();
                if idx == 0 {
                    chrom = fields[0].to_string();
                    start = fields[1].to_string();
                    end = fields[2].to_string();
                }

                let val: f64 = fields[3].trim().parse().unwrap();  // 保証されている
                row_values.push(format!("{:.4}", val));
            }
        }

        if all_done {
            break;
        }

        let mut output_row = vec![chrom, start, end];
        output_row.extend(row_values);
        writeln!(csv_writer, "{}", output_row.join(","))?;
    }


    Ok(())
    // chunk_iter(bed_gz_file.into_iter(), 1000).for_each(|chunk| {
    //     println!("{:?}", chunk.len());
    // });


}




fn chunk_iter<T, I>(mut iter: I, size: usize) -> impl Iterator<Item = Vec<T>>
where
    I: Iterator<Item = T>,
{
    std::iter::from_fn(move || {
        let mut chunk = Vec::with_capacity(size);
        for item in iter.by_ref().take(size) {
            chunk.push(item);
        }
        if chunk.is_empty() {
            None
        } else {
            Some(chunk)
        }
    })
}