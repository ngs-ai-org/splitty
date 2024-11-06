//! ## vcf_bnd_2fasta ##
//!---------------------
//!Extracts from paired BND VCF files the sequence in the proper direction based on the events.
//!It reconstructs based on biological limitations in silico the sequence which might have given rise to the event.

use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
use std::env;
// our library which is within the same project
extern crate genefusion;
use genefusion::lib::hts_lib_based::{*};
//use genefusion::lib::common::{*};
use rust_htslib::bcf::Read as BcfRead;
use std::str;

use bio::io::fasta::Writer;
//use bio::alphabets::dna;

fn main() {
    pretty_env_logger::init();

    // now the next is not really for any argument 
    // parsing but simply to get the command which 
    // was used to execute as I cant get this from clap
    let args: Vec<String> = env::args().collect();
    let _args_string = args.join(" ");
    let matches = app_from_crate!()
            .about("This tool takes BND events and extracts annotated FASTA files with the in silico combined event.\
                    It expects pairs of BND events and needs the matching reference indexed fasta sequence. \
                    It cant currently work with events other than BND and will ignore these. \
                    By default the header of a FASTA entry will contain the name of the BND events. \
                    The comment field in the header contains additionally information about the position of the fusion point (FP) \
                    and the range chosen for the in silico edges.\
                    The sequence is split where the first half is upper case and the second part of the sequence is lower-case. \
                    This faciliates the identification of the event point \
                    Entries need to be sorted by the identifier of the BND events and pairs must follow each other, e.g.:
                    
                    \n\n
        E.g.: awk \'/^#/{print $0 ; next}{ print $0 | \"sort -V -k3,3\" }\' File.vcf
            ")
            .arg(Arg::with_name("VCF")
                .short("i")
                .long("vcf")
                .value_name("file")
                .help("VCF file with BND events")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("REF")
                .short("r")
                .long("reference")
                .value_name("file")
                .help("the already indexed FASTA sequence of the VCF matching reference")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("OUT")
                .short("o")
                .long("out")
                .value_name("file")
                .help("the in silico sequence of the combined BND events")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("SIMPLE")
                .short("s")
                .long("simple-id")
                .value_name("bool")
                .help("simplify ID and description to InSilico_n")
                .takes_value(false)
                .required(false))
            // .arg(Arg::with_name("OUT2")
            //     .short("b")
            //     .long("bed")
            //     .value_name("file")
            //     .help("a matching BED file for the FASTA output which annotates the name and position of the FP event")
            //     .takes_value(true)
            //     .required(false))    
            // .arg(Arg::with_name("THREADS")
            //     .short("t")
            //     .long("threads")
            //     .value_name("int")
            //     .help("vcf reading threads (files need to be huge to make an impact)")
            //     .takes_value(true)
            //     .default_value("2")
            //     .required(false))
            .arg(Arg::with_name("RANGE")
                .short("e")
                .long("extent")
                .value_name("int")
                .help("the range for both BND position till which the in silico sequence will extend")
                .takes_value(true)
                .default_value("500")
                .required(false))
            .arg(Arg::with_name("PREFIX")
                .short("p")
                .long("prefix")
                .value_name("string")
                .help("a prefix for the fasta header, can be e.g. sample name - will default to file name without prefix. Simple argument will set default to empty")
                .takes_value(true)
                .required(false))
            // .arg(Arg::with_name("STRAND")
            //     .short("s")
            //     .long("strand-ignore")
            //     .value_name("bool")
            //     .help("if activated ignores strand and generates 2 entries for each with both directions")
            //     .takes_value(false)
            //     .required(false))    
            .get_matches();

    let file      = matches.value_of("VCF").unwrap();
    let assembly  = matches.value_of("REF").unwrap();
    let range      = matches.value_of("RANGE").unwrap().parse::<u32>().expect("ERROR: could not parse range option correctly!");
    let _stranded = matches.is_present("STRAND");
    let output    = matches.value_of("OUT").unwrap();
    let simple_name= matches.is_present("SIMPLE");
    let mut n_records : i32 = 0;
    

    // lets take now the input and read the VCF correctly:
    let mut vcf_file     = rust_htslib::bcf::Reader::from_path(file).expect("ERROR: could not open vcf file!");
    let mut writer = Writer::to_file(output).expect("ERROR: could not open output file");
    let vcf_header   = vcf_file.header().samples().to_vec();
    if vcf_header.len() > 1 {
        panic!("ERROR: Currently only single sample vcf accepted!");
    };
    
    // we need to initialize an empty record first
    let mut prev_record = vcf_file.empty_record();
    // now we need the FASTA writer ready and  iterate
    // over VCF entries. Important: must be sorted by BND pairs
    for (n,entry) in vcf_file.records().enumerate() {
        // Turns out first entry is 0 not 1,
        // so we need to count  
        if let 0=n%2 { 
            prev_record = entry.expect("ERROR: could not read 1st record of vcf file!"); 
            continue
        }else{
            n_records+=1;
            let cur_record = entry.expect("ERROR: could not read 1st record of vcf file!");
            let id_1     = &prev_record.id() ;
            let id_2     = &cur_record.id() ;

            let id1  = str::from_utf8(id_1).expect("ERROR: cant extract record id!");
            if id1 == "." {
                panic!("ERROR: all entries must have an annotated entry ID!");
            };
            let id2  = str::from_utf8(id_2).expect("ERROR: cant extract record id!");
            if id2 == "." {
                panic!("ERROR: all entries must have an annotated entry ID!");
            };

            let mate_id       = cur_record.info(b"MATEID").string().unwrap().unwrap();
            let mate_id2      = mate_id[0];
            if (id_1 != mate_id2) || (mate_id2.is_empty()) {
                eprintln!("Non-matching entries pair entries {:?} and {:?}",str::from_utf8(id_1),str::from_utf8(mate_id2));
                panic!("ERROR: VCF entries are either mal-formated or not sorted by matching IDs!");
            }
            // Now we would only return 1 entry logically speaking,
            // but we want to provide the possibility to get both senses if needed.
            // E.g. RNAseq has no real sense
            // Therefore lets make here a list 
            let record : bio::io::fasta::Record = vcfpair2fa_record(
                &prev_record,
                &cur_record,
                assembly,
                range,
            );
            if simple_name {
                let simple_record = bio::io::fasta::Record::with_attrs(
                    &format!("InSilico_{}",n_records), 
                    record.desc(), 
                    record.seq());
                    writer.write_record(&simple_record).expect("ERROR: could not write fasta record!");
            }else{
                writer.write_record(&record).expect("ERROR: could not write fasta record!");
            }
            // if stranded {
            //     let tmp = bio::alphabets::dna::revcomp(record.seq());
            //     let new_record = bio::io::fasta::Record::with_attrs(record.id(), Some(&format!("{:?}/revComp",record.desc())),record.seq());
            //     writer.write_record(&new_record).expect("ERROR: could not write fasta record!");
            // }
        };
    }
    println!("INFO: {} records written",n_records);

}
    