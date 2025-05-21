use std::collections::HashMap;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::{BufRead, BufReader, Write};

fn main() {
    println!("GTFファイルの読み取りサンプル");

    let file_path = "data/gencode.v48.annotation.gtf";
    let file = File::open(file_path).expect("GTFファイルを開けませんでした");
    let reader = BufReader::new(file);

    // gene_id -> gene_name, chr, start, end, strand, gene_type, level, source, その他
    let mut gene_info: HashMap<String, (String, String, String, String, String, String, String, String, Vec<String>)> = HashMap::new();

    for line in reader.lines() {
        let line = line.expect("行の読み込みに失敗");
        if line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 { continue; }
        let chr = fields[0].to_string();
        let source = fields[1].to_string();
        let feature_type = fields[2];
        let start = fields[3].to_string();
        let end = fields[4].to_string();
        let strand = fields[6].to_string();
        let attrs = fields[8];
        if feature_type != "gene" { continue; }

        let mut gene_id = String::new();
        let mut gene_name = String::new();
        let mut gene_type = String::new();
        let mut level = String::new();
        let mut other_info = Vec::new();
        for attr in attrs.split(';') {
            let attr = attr.trim();
            if attr.starts_with("gene_id ") {
                gene_id = attr[8..].trim_matches('"').to_string();
            } else if attr.starts_with("gene_name ") {
                gene_name = attr[10..].trim_matches('"').to_string();
            } else if attr.starts_with("gene_type ") {
                gene_type = attr[10..].trim_matches('"').to_string();
            } else if attr.starts_with("level ") {
                level = attr[6..].trim().to_string();
            } else if !attr.is_empty() {
                other_info.push(attr.to_string());
            }
        }
        if !gene_id.is_empty() {
            gene_info.insert(
                gene_id,
                (gene_name, chr, start, end, strand, gene_type, level, source, other_info),
            );
        }
    }

    // 出力ファイル作成
    let output_path = "output/genes.tsv";
    let mut output = OpenOptions::new().create(true).write(true).truncate(true).open(output_path)
        .expect("出力ファイルを作成できませんでした");
    // ヘッダー
    writeln!(output, "gene_id\tgene_name\tchr\tstart\tend\tstrand\tgene_type\tlevel\tsource\tother").unwrap();
    for (gene_id, (gene_name, chr, start, end, strand, gene_type, level, source, other_info)) in &gene_info {
        let other = other_info.join("; ");
        writeln!(output, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", gene_id, gene_name, chr, start, end, strand, gene_type, level, source, other).unwrap();
    }
    println!("genes.tsvを{}に出力しました", output_path);
}
