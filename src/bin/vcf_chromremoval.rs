//! ## vcf_chromremoval ##
//! -----------------------
//! A tool which allows to remove/keep VCF entries which match a list of chromosome entries.
//! Used to remove e.g. decoy chromosome events. Can treat correctly as well BND entries, meaning to remove
//! both entries if one of them is in the list.
use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
use std::env;
//use rustc_hash::FxHashMap;
//use bio::data_structures::interval_tree;
// our library which is within the same project
extern crate genefusion;
//use genefusion::lib::hts_lib_based::{*};
use genefusion::lib::common::{*};
use rust_htslib::bcf::header::Header;
use rust_htslib::bcf::{Read as BcfRead};
use rust_htslib::bcf::{Format};
use std::str::from_utf8;
use std::str;
use rustc_hash::FxHashMap;
extern crate pretty_env_logger;

// essentially we will do for the 1st vcf a first run population
// a hashmap were the keys are concatenated iD:mateID + the link to the vcf entry 
// together with a structure containging very simplified intormation 
// In a second step we do the same for the 2nd vcf but instead of the 2nd structure
// we compare it directly to the 1st strucure and add it accordingly 

fn main() {
    pretty_env_logger::init();

    // now the next is not really for any argument 
    // parsing but simply to get the command which 
    // was used to execute as I cant get this from clap
    let args: Vec<String> = env::args().collect();
    let args_string = args.join(" ");
    let matches = app_from_crate!()
    .about("This tools will analyze a VCF and only keep events provided in a list of chromosomes.\n\
    Importantly this will work both for all type of events (e.g. INS,DEL,SNV,DUPL,..) but as well paired BND events.\n\
    For the latter it will remove both pairs if one of them is on the black list. \
    The chromosome section of the VCF header will **not be altered**.
    

    ")
    .arg(Arg::with_name("VCF")
            .short("a")
            .long("vcf")
            .value_name("FILE")
            .help("the VCF file to subset")
            .takes_value(true)
            .required(true))
    .arg(Arg::with_name("CHROMS")
            .short("c")
            .long("chroms")
            .value_name("FILE")
            .help("list of chromosomes to keep, one entry per line")
            .takes_value(true)
            .required(true))
    .arg(Arg::with_name("OUT")
            .short("o")
            .long("out")
            .value_name("FILE")
            .help("outfile name for the VCF")
            .takes_value(true)
            .required(false))
    .arg(Arg::with_name("INVERSE")
            .short("i")
            .long("inverse")
            .value_name("BOOL")
            .help("will do inverse selection, removing all entries matching the provided list ")
            .takes_value(false))
    .arg(Arg::with_name("SIMPLE")
            .short("s")
            .long("simple")
            .value_name("BOOL")
            .help("will not do a 2-pass but filter directly. Allows missing VCF-ID. Will fail on paired BND events")
            .takes_value(false))
    .get_matches();

    let vcf_file       = matches.value_of("VCF").unwrap();
    let out_vcf= matches.value_of("OUT");
    let chroms         = matches.value_of("CHROMS").unwrap();
    let inverse        = matches.is_present("INVERSE");
    let simple               = matches.is_present("SIMPLE");
    // here we just get the chromosomes as a simple hash
    let selected_chroms = parse_chroms_txt(chroms);
    eprintln!("INFO: chromosome list contained {} unique entries",selected_chroms.len());

    // some information which we inject in the header
    // to document the VCF being altered
    let infos =  &VersionInfo {
        program: &String::from("vcf_chromremoval"),
        version: crate_version!(),
        author : crate_authors!(),
        command: &args_string,
    };
    let header_source_line= format!("{}:{}",r#"##source=splitty"# ,infos.version);    
    let header_author_line= format!("{}{}",r#"##source_author="# ,infos.author);
    let header_cmd_line   = format!("{}{}",r#"##command="# ,infos.command);

    // So next we want now to parse the VCF file get it again in pairs.
    // Then we check if the pair is in one of the annotated regions pairs
    // and if this is the case we write the VCF entry in a new file
    // lets take now the input and read the VCF correctly:
    let mut vcf_reader1 = rust_htslib::bcf::Reader::from_path(vcf_file).expect("ERROR: could not open vcf file!");
    let mut vcf_header : rust_htslib::bcf::Header = Header::from_template(vcf_reader1.header());
    vcf_header.push_record(header_source_line.as_bytes());
    vcf_header.push_record(header_author_line.as_bytes());
    vcf_header.push_record(header_cmd_line.as_bytes());

    // either we write to file or STDOUT
    let mut writer : rust_htslib::bcf::Writer = match out_vcf {
        Some(x) => rust_htslib::bcf::Writer::from_path(
            x,
            &vcf_header,
            true,
            Format::Vcf).expect("ERROR: cant write VCF!"),
        None    => rust_htslib::bcf::Writer::from_stdout(
            &vcf_header,
            true,
            Format::Vcf).expect("ERROR: cant write VCF!"),
    };

    // this will be our list with IDs that should be removed    
    let mut black_list_id : FxHashMap<String,u8> = FxHashMap::default();
    // now we will iterate twice over it.
    // once we will find all IDs and MATE-IDs which are associated
    // with unwanted chromosomes and then a second time to read and write directly

    if !simple {
    // 1st path defining entries to remove
        for entry in vcf_reader1.records() {
            let record = entry.expect("ERROR: could not read entry!");
            // Turns out first entry is 0 not 1,
            let id1 = &record.id() ;
            let id  = from_utf8(id1).expect("ERROR: cant extract record id!");
            if id == "." {
                panic!("ERROR: all entries must have an annotated entry ID!");
            };

            // I wanted to do that instead with a result and define
            // none but somehow it really butchered the resulting string
            // occasionally - no idea how
            let mate_id = if record.info(b"MATEID").string().is_err() || record.info(b"MATEID").string().unwrap().is_none() {
                String::from(".")
            }else{
                let mate_id1       = record.info(b"MATEID").string().unwrap().unwrap();
                str::from_utf8(mate_id1[0]).expect("ERROR: could not extract MateID string!").to_string()
            };

            
            // next we get the chromosome names via the IDs from the header
            let rid1 = record.rid().unwrap();
            let chr1 : &str  = str::from_utf8(record.header().rid2name(rid1).expect("ERROR: could not convert rid to chromosome name!")).unwrap();
            //eprintln!("Found ID {} and mate {} for chromosome {}",&id,&mate_id,&chr1);
            // now comes the decision if the ID goes into
            // our black-list or not
            if (!selected_chroms.contains_key(chr1) & !inverse) || ( selected_chroms.contains_key(chr1) & inverse ){
                //eprintln!("Pushing now ID {} and mate {} for chromosome {} into bad bunch",&id,&mate_id,&chr1);
                black_list_id.insert(id.to_string(), 1);
                // now if we have as well the mate, then 
                // this one should be on the black-list, too
                if mate_id != "." {
                    black_list_id.insert(mate_id, 1);
                }
            };
            //writer.write(&prev_record).expect("ERROR: could not write VCF records!");
        } ;
        eprintln!("INFO: IDs to remove contained {} unique entries",black_list_id.len());
    }
    // now we go simply again over all entries and 
    // remove the ones which are on the black-list
    let mut vcf_reader2 = rust_htslib::bcf::Reader::from_path(vcf_file).expect("ERROR: could not open vcf file!");
    let mut simple_count = 0_u32;
    for entry2 in vcf_reader2.records() {
        let record = entry2.expect("ERROR: could not read entry!");
        // Turns out first entry is 0 not 1,
        let id = &record.id() ;
        let id1  = from_utf8(id).expect("ERROR: cant extract record id!");
        let rid = record.rid().unwrap();
        let chr : &str  = str::from_utf8(record.header().rid2name(rid).expect("ERROR: could not convert rid to chromosome name!")).unwrap();
        if !simple {
            if id1 == "." {
                panic!("ERROR: all entries must have an annotated entry ID!");
            };
            if black_list_id.contains_key(id1) {
                continue
            }else {
                writer.write(&record).expect("ERROR: could not write VCF records!");
            }
        }else{
            if (selected_chroms.contains_key(chr) & !inverse) || ( !selected_chroms.contains_key(chr) & inverse ){
                writer.write(&record).expect("ERROR: could not write VCF records!");
            }else{
                simple_count += 1;
            };
        }
        
    } ;
    if simple {
        eprintln!("INFO: Simple approach removed {} unique entries",simple_count);
    }
}