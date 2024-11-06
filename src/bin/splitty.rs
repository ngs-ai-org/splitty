//! 
//! ## Fragment-based ##
//! -------------
//! Splitty uses mapping of cDNA or transcript fragments to a genome in order to identify fusion points.
//! Currently exclusively for longRead/cDNA-on-genome-mapping, with the abilitiy to annotate as well the
//! genes if it falls into the region of one. GTF annotation from either CHESS or Gencode is supported
//! 
//! ## Paired-read based ##
//! --------------------
//! Splitty uses paired-read mapping onto a genome in order to identify fusion points.
//! GTF annotation from either CHESS or Gencode is supported.
//! The output is organized by position in constrast to the read-based result of splitty.
//! 
//! ## id_fusion_reads ##
//!---------------------
//!This tool is designed for paired-read mapping data and identifies reads
//!which are indicative for fusions. Read must be mapped onto a transcriptome to make this
//!work effectively.
//!It does not suggest fusion points but purely returns the subset of reads.

use rust_htslib::{bam, bam::Read};
use std::path::Path;
use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,SubCommand,Arg};

use std::io::{BufRead, BufReader};
use std::str;
use std::env;
use std::process;
use rustc_hash::FxHashMap;
use std::fs::File;
extern crate pretty_env_logger;
#[macro_use] extern crate log;


// our library which is within the same project
extern crate genefusion;
use genefusion::lib::common::{*};
use genefusion::lib::hts_lib_based::{*};
use genefusion::lib::noodles_based::{*};


fn file_per_line_hashmap(
    my_file: &str
) -> Option<FxHashMap<String, i32>> {
    // this function takes a file and returns a hashmap
    // of all entries in the form entry=>1
    // this is used later to serve as list of entries
    // which are to be ignored
    // it does not accept space between entries
    // and will exit in that case

    let mut tmp_map = FxHashMap::default();
    eprintln!("INFO: blacklist {} provided, reading entries...", my_file);
    assert!(
        Path::new(my_file).exists(),
        "ERROR: blacklist file {:?} does not exist!",
        my_file
    );
    let input  = File::open(my_file).expect("ERROR: Could not open blacklist file!");
    let reader = BufReader::new(input);
    for line in reader.lines() {
        let l = line.expect("ERROR: could not read line!");
        // now we split by empty space
        let e: Vec<&str> = l.split(' ').collect();
        if e.len() != 1 {
            panic!("ERROR: your blacklist contains space delimited entries!");
        }
        tmp_map.insert(e[0].to_string(), 1);
    }
    Some(tmp_map)
}

fn run_fragment_base(
    matches: &clap::ArgMatches ,
    arg_string: &str 
){
        ////////////////////////
        ////  prep options  ////
        ////////////////////////
        let stranded        = matches.is_present("STRAND");
        let bam_threads     = matches.value_of("THREAD").unwrap().parse::<u32>().unwrap();
        let filter          = matches.value_of("BLACKLIST").unwrap_or("NA");
        let annot           = matches.value_of("GTF").unwrap_or("NA");
        let annot_filter    = matches.value_of("MATCH").unwrap_or("transcript");
        let mut keep_all    = false;
        let ref_file        = matches.value_of("REF").unwrap_or("NONE");
        let vcf_file        = matches.value_of("VCF").unwrap_or("NONE");
        let multi_mapper    = matches.is_present("MULTI");
        let ratio           = match matches.value_of("RATIO").unwrap().parse::<f32>().expect("ERROR: could not parse ratio value"){
            0.0 => None,
            x => Some(x)
        };

        if matches.is_present("KEEP") {
            keep_all = true;
        };

        let identity        = match matches.value_of("SIM").unwrap_or("98.0").parse::<f32>() {
            Ok(v) => v,
            _ => panic!("ERROR: could not parse \"SIM\"-field correctly"),
        };

        let intra_dist: u64  = match matches.value_of("DIST").unwrap_or("10000").parse::<u64>() {
            Ok(v) => v,
            _ => panic!("ERROR: could not parse \"DIST\"-filed correctly"),
        };

        if matches.is_present("BAMBAM") {
            bambam::bam_bam_inda_house();
        }

        // we can just unwrap as mandatory argument
        let bam_file = matches.value_of("BAM").unwrap();
        eprintln!("Input BAM file is {}", &bam_file);
        // if file not exist quit
        assert!(
            Path::new(&bam_file).exists(),
            "ERROR: input file {:?} does not exit!",
            &bam_file
        );

        ////////////////////////
        ////  bam checking  ////
        ////////////////////////
        let mut bam     = bam::Reader::from_path(bam_file)
            .expect("ERROR: Could not open BAM file ");
        
        let header = bam::Header::from_template(bam.header());
        let header_view = header.to_hashmap();
    
    
        // we need this later to look up the ID of a target-id
        // as these are otherwise only i32
        
        //////////////////////////
        ////  prep blacklist  ////
        //////////////////////////
        let blacklist_map: Option<FxHashMap<String, i32>> ;
        if filter == "NA" {
            blacklist_map = None;
            eprintln!("INFO: no blacklist provided");
        } else {
            blacklist_map = file_per_line_hashmap(filter);
        }

        /////////////////////
        //// bam reading ////
        /////////////////////
        
        // generate a thread pool for the multi-threaded reading and writing of BAM files
        let pool = rust_htslib::tpool::ThreadPool::new(bam_threads).unwrap();
        bam.set_thread_pool(&pool).unwrap();


        // here we read now the BAM and populate our read-based
        // structure which documents all reads with SA alignment
        let (read_based_organization,encountered_chroms) = analyze_alignments_from_bam(
            bam,
            blacklist_map,
            identity,
            false,
            multi_mapper,
            ratio
        );

        debug!("Read based organization: {:?}",&read_based_organization);
        // now we can start with "pruning"
        // essentially we check for each read if there
        // are 2 entries supporting a fusion
        // if that is not the case we drop the read
        // E.g. one of both alignments did not match
        // our similarity cut-off
        let pruned_reads = prune_read_based_organization(read_based_organization,keep_all);


        ///////////////////////////
        ////  prep annotation  ////
        ///////////////////////////
        let mut annot_map: Option<FxHashMap<String, Vec<Gtf>>> = None;
        if annot == "NA" {
            eprintln!("INFO: no annotation provided");
        } else {
            if !encountered_chroms.is_empty() {
                let mut tmp : Vec<String> = Vec::new();
                for (my_key,_) in  encountered_chroms.iter(){
                    tmp.push(my_key.to_owned());
                };
                annot_map = Some(
                    parse_gtf_annotation(
                        annot,
                        Some(annot_filter),
                        Some(tmp)
                    )
                );
                let annot_keys = annot_map.as_ref().unwrap().keys().len();
                if annot_keys == 0 {
                    panic!("ERROR: your annotation did not result in any considered feature!");
                }
                eprintln!("INFO: your provided annotation contained {} considered chromosomes",annot_keys);
            }else{
                eprintln!("WARNING: did not encounter any entry with a valid chromosome !");
            }
           
        }
        
        // next we format correctly our output and annotate if a 
        // annotation was provided
        debug!("read_collection: {:?}, gtf: {:?}",pruned_reads,annot_map);
        let final_collection: Vec<PrincipleReadBasedFusionOutput>  = format_and_annotate_read_based_fusions(
            pruned_reads,
            annot_map, 
            stranded,
            intra_dist
        );
        debug!("final_collection: {:?}",final_collection);
        // last but not least we write our output
        match vcf_file {
            "NONE" => { 
                if let Err(err) = write_read_based_tsv_stdout(&final_collection,crate_version!(),crate_authors!(),arg_string) {
                    println!("{}", err);
                    process::exit(1);
                }
            },
            _ => {
                if let Err(err) = write_read_based_vcf_file(&final_collection,
                    &VersionInfo {
                        program: &String::from("Splitty"),
                        version: crate_version!(),
                        author: crate_authors!(),
                        command: arg_string,
                    },
                    header_view,
                    ref_file,
                    vcf_file,
                ){
                    println!("{}", err);
                    process::exit(1);
                }
                if let Err(err) = write_read_based_tsv_stdout(&final_collection,crate_version!(),crate_authors!(),arg_string) {
                    println!("{}", err);
                    process::exit(1);
                }
                /*if let Err(err) = write_vcf_file(&final_collection,x,crate_version!(),crate_authors!(),&args_string) {
                    println!("{}", err);
                    process::exit(1);
                }
                */   
                },
        }
}

fn run_paired_base(
    matches: &clap::ArgMatches ,
    arg_string: &str
){
      ////////////////////////
        ////  prep options  ////
        //////////////////////// 
        let stranded        = matches.is_present("STRAND");
        let bam_threads      = matches.value_of("THREAD").unwrap().parse::<usize>().unwrap();
        let annot           = matches.value_of("GTF").unwrap_or("NA");
        let annot_filter    = matches.value_of("MATCH").unwrap_or("transcript");
        let ref_file        = matches.value_of("REF").unwrap_or("NONE");
        let vcf_file        = matches.value_of("VCF").unwrap_or("NONE");
        let multi_mapper    = matches.is_present("MULTI");

            // the range which we allow as deviation from a reported fusion
        // point to be collapsed with
        let range_sr : u32  = matches.value_of("PRECISION").unwrap().parse::<u32>().unwrap();
        let range_sp : u32  = matches.value_of("RANGE").unwrap().parse::<u32>().unwrap();
        let intra_dist      = match matches.value_of("DIST").unwrap_or("10000").parse::<u64>() {
            Ok(v) => v,
            _ => panic!("ERROR: could not parse \"DIST\"-filed correctly"),
        };
        let identity        = match matches.value_of("SIM").unwrap_or("98.0").parse::<f32>() {
            Ok(v) => v,
            _ => panic!("ERROR: could not parse \"SIM\"-field correctly"),
        };

        if matches.is_present("BAMBAM") {
            bambam::bam_bam_inda_house();
        }

        // we can just unwrap as mandatory argument
        let bam_file = matches.value_of("BAM").unwrap();
        eprintln!("Input BAM file is {}", &bam_file);
        // if file not exist quit
        assert!(
            Path::new(&bam_file).exists(),
            "ERROR: input file {:?} does not exit!",
            &bam_file
        );

        ////////////////////////
        ////  bam checking  ////
        ////////////////////////
        // okay here I am only going for an indexed reader because we need
        // an index for the later step of verification and we want to 
        // avoid that the program runs for some minutes before then
        // complaining that the 2nd stage has no index

        let mut stage_1_bam = bam::Reader::from_path(bam_file).expect("ERROR: Could not open BAM file ");
        let mut stage_2_bam = bam::IndexedReader::from_path(bam_file).expect("ERROR: Could not open BAM file ");

        let header = bam::Header::from_template(stage_1_bam.header());
        let header_view = header.to_hashmap();


        /////////////////////
        //// bam reading ////
        /////////////////////

        // generate a thread pool for the multi-threaded reading and writing of BAM files
        //let pool = rust_htslib::tpool::ThreadPool::new(bam_threads).unwrap();
        stage_1_bam.set_threads(bam_threads).expect("ERROR: could not set correctly read threads");
        stage_2_bam.set_threads(bam_threads).expect("ERROR: could not set correctly read threads");
        // stage_1_bam.set_thread_pool(&pool).unwrap();
        // stage_2_bam.set_thread_pool(&pool).unwrap();

        // here we read now the BAM and populate our read-based
        // structure which documents all reads with SA alignment
        // Essentially I do not report any fusion if there is not
        // at least 1 read which gives a precise position
        let (read_based_organization,encountered_chroms) = analyze_alignments_from_bam(
            stage_1_bam,
            None, // blacklist
            identity, // match similarity
            true, // if paired
            multi_mapper, // if accepting multi-mapper
            None // ratio of 5-and 3-clipping accepted
        );
        
        debug!("Read based organization: {:?}",&read_based_organization);
        // now we can start with "pruning"
        // essentially we check for each read if there
        // are 2 entries supporting a fusion
        // if that is not the case we drop the read
        // E.g. one of both alignments did not match
        // our similarity cut-off
        let pruned_reads      = prune_read_based_organization(read_based_organization,false);

        ///////////////////////////
        ////  prep annotation  ////
        ///////////////////////////
        debug!("Getting the annotation ");
        let mut annot_map: Option<FxHashMap<String, Vec<Gtf>>> = None;
        if annot == "NA" {
            eprintln!("INFO: no annotation provided");
        } else {
            if !encountered_chroms.is_empty() {
                let mut tmp : Vec<String> = Vec::new();
                for (my_key,_) in  encountered_chroms.iter(){
                    tmp.push(my_key.to_owned());
                };
                annot_map = Some(
                    parse_gtf_annotation(
                        annot,
                        Some(annot_filter),
                        Some(tmp)
                    )
                );
            }else{
                panic!("ERROR: did not encounter any entry with a valid chromosome !");
            }
            let annot_keys = annot_map.as_ref().unwrap().keys().len();
            if annot_keys == 0 {
                panic!("ERROR: your annotation did not result in any considered feature!");
            }
            eprintln!("INFO: your provided annotation contained {} considered chromosomes",annot_keys);
        }
        debug!("Following annotation: {:?}",annot_map);
        
        // we take now the results which are still on a per-read-basis
        // and sum them up into a position based system, collapsing identical
        // information and keeping track of the additional support
        let pos_based_results = read_based_2_pos_based_fusions(
            pruned_reads,
            range_sr,
            stranded
        );
        debug!("Read-based 2 postion based result: {:?}",pos_based_results);
        debug!("Now adding split pairs to split reads");
        // As we have now a position based results 
        // which contains for a given fusion a chromosome A + B 
        // together with a position A + B we can now get for each
        // the support of split-pairs as well.
        // Since we now have positions this can be done with 
        // querying directly said region instead of re-analyzing all 
        // reads 1 by 1 again
        let results_with_sp = add_sp_to_sr_fusions(
            & mut stage_2_bam,
            pos_based_results,
            range_sp
        );
        debug!("Now formatting everything accordingly");
        let mut final_collection: Vec<PrinciplePosBasedFusionOutput>  = format_and_annotate_pos_based_fusions(results_with_sp, annot_map, stranded, intra_dist);
        debug!("Final colection: {:?}", final_collection);
        // generate our principle tsv table
        debug!("Now writing TSV output");
        //write_pos_based_tsv_stdout(&final_collection,crate_version!(),crate_authors!(),arg_string).expect("ERROR: failed to write results!");
        // if requested as well the vcf
        if ref_file !="NONE" && vcf_file !="NONE"{
            debug!("Finally as well VCF generation");
            write_pos_based_vcf_file(&mut final_collection,
                &VersionInfo {
                    program: &String::from("Splitty"),
                    version: crate_version!(),
                    author: crate_authors!(),
                    command: arg_string,
                },
                header_view,
                ref_file,
                vcf_file,
            ).expect("ERROR: failed to write results!");
        }
}

fn run_id_reads(
    matches: &clap::ArgMatches
){
        ////////////////////////
        ////  prep options  ////
        ////////////////////////
        let identity        = match matches.value_of("SIM").unwrap_or("98.0").parse::<f32>() {
            Ok(v) => v,
            _ => panic!("ERROR: could not parse \"SIM\"-field correctly"),
        };
        let threads : usize = match matches.value_of("THREAD").unwrap().parse::<usize>() {
            Ok(v) => v,
            _ => panic!("ERROR: could not parse provided thread number!"),
        };
        let source = matches.value_of("REF").expect("ERROR: could not get the source field successfully!");
        if matches.is_present("BAMBAM") {
            bambam::bam_bam_inda_house();
        }
        // we can just unwrap as mandatory argument
        let bam_file = matches.value_of("BAM").unwrap();
        eprintln!("Input BAM file is {}", &bam_file);
        // if file not exist quit
        assert!(
            Path::new(&bam_file).exists(),
            "ERROR: input file {:?} does not exit!",
            &bam_file
        );

        /////////////////////
        //// bam reading ////
        /////////////////////
        // here we read now the BAM and populate our read-based
        // structure which documents all reads with SA alignment
        let fusion_reads = analyze_paired_alignments_from_trans_bam_noodles(bam_file,identity,source,threads);      
        for i in fusion_reads {
            println!("{}",i);
        }
}

fn main() {
    pretty_env_logger::init();
    // define counters for stats
    let read_threads = "1";

    // now the next is not really for any argument 
    // parsing but simply to get the command which 
    // was used to execute as I cant get this from clap
    let args: Vec<String> = env::args().collect();
    let args_string = args.join(" ");
    let matches = app_from_crate!()
        .subcommand(SubCommand::with_name("fragment")
            .about("classifying cDNA/DNA mappings into gene-fusion events, not for paired-end reads")
            .arg(Arg::with_name("BAM")
                .short("b")
                .long("bam")
                .value_name("FILE")
                .help("long read/cDNA alignment vs reference")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("GTF")
                .short("g")
                .long("gtf")
                .value_name("FILE")
                .help("Provide annotation of the genome. \
                        Use option m to specify selected feature \
                        The field \"gene_id\" or, if available \"gene_name\" will be used for reporting \
                        Comply either to CHESS or Gencode standard")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("MATCH")
                .short("m")
                .long("match")
                .value_name("feature-name")
                .help("specifies the feature field in a GTF to be used, \
                    [defaults: \"transcript\" ] \
                    \nactivates automatically option \"gtf\" too")
                .takes_value(true)
                .requires("GTF")
                )
            .arg(Arg::with_name("MULTI")
                .short("u")
                .long("unique-off")
                .help("accepts reads with non-unique mapping mapq=0")
                .takes_value(false)
                )
            .arg(Arg::with_name("BLACKLIST")
                .short("l")
                .long("blacklist")
                .value_name("FILE")
                .help("black-list of queries to ignore (e.g. assembled transcripts identified \
                        as mouse derived), one ID per line")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("SIM")
                .short("i")
                .long("similarity")
                .value_name("float")
                .help("minimum seq similarity of matches [default:98.0]")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("REF")
                .short("r")
                .long("reference")
                .value_name("int")
                .help("reference in fasta format provided, required for vcf output")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("VCF")
                .short("f")
                .long("vcf")
                .value_name("int")
                .help("generate vcf encoded fusion call, requires reference")
                .takes_value(true)
                .requires_all(&["REF"])
                .required(false))
            .arg(Arg::with_name("DIST")
                .short("w")
                .long("distance")
                .help("minimum allowed distance for intra-chromosomal fusions [default=10000]")
                .value_name("int")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("STRAND")
                .short("s")
                .long("stranded")
                .help("implies that used cDNA is stranded \
                    this provides better annotation precision")
                .takes_value(false)
                .required(false))
            .arg(Arg::with_name("BAMBAM")
                .short("x")
                .long("bambam")
                .value_name("IDH")
                .takes_value(false)
                .required(false)
                .hidden(true))	
            .arg(Arg::with_name("THREAD")
                .short("t")
                .long("threads")
                .value_name("int")
                .help("number of threads for reading BAM")
                .takes_value(true)
                .required(false)
                .default_value(read_threads))
            .arg(Arg::with_name("RATIO")
                .short("d")
                .long("ratio")
                .value_name("float")
                .help("ratio of accepted 5' to 3' clipping (and other way around), 0.0 will deactivate check")
                .takes_value(true)
                .required(false)
                .default_value("0.1"))
            .arg(Arg::with_name("KEEP")
                .short("k")
                .long("keep")
                .value_name("bool")
                .help("keep all suggested entries if multiple SAs available, otherwise only most precise\
                    careful though this might increase the number of FPs and precision might become wrong for these entries\
                    in particular for entries which are |----A-->|---B-->|---C--> architecture")
                .takes_value(false)
                .required(false)))
        .subcommand(SubCommand::with_name("paired")
            .about("detecting integrations and genefusions from mapped paired-reads")
                .arg(Arg::with_name("BAM")
                .short("b")
                .long("bam")
                .value_name("FILE")
                .help("long read/cDNA alignment vs reference")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("GTF")
                .short("g")
                .long("gtf")
                .value_name("FILE")
                .help("Provide annotation of the genome. \
                        Use option m to specify selected feature \
                        The field \"gene_id\" or, if available \"gene_name\" will be used for reporting \
                        Comply either to CHESS or Gencode standard")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("MATCH")
                .short("m")
                .long("match")
                .value_name("feature-name")
                .help("specifies the feature field in a GTF to be used, \
                     [defaults: \"transcript\" ] \
                     \nactivates automatically option \"gtf\" too")
                .takes_value(true)
                .requires("GTF")
                )
            .arg(Arg::with_name("MULTI")
                .short("u")
                .long("unique-off")
                .help("accepts reads with non-unique mapping mapq=0")
                .takes_value(false)
                )
            .arg(Arg::with_name("DIST")
                .short("w")
                .long("distance")
                .help("minimum allowed distance for intra-chromosomal fusions [default=10000]")
                .value_name("int")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("STRAND")
                .short("s")
                .long("stranded")
                .help("implies that used cDNA is stranded \
                    this provides better annotation precision")
                .takes_value(false)
                .required(false))
            .arg(Arg::with_name("SIM")
                .short("i")
                .long("similarity")
                .value_name("float")
                .help("minimum seq similarity of matches [default:98.0]")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("BAMBAM")
                .short("x")
                .long("bambam")
                .value_name("IDH")
                .takes_value(false)
                .required(false)
                .hidden(true))	
            .arg(Arg::with_name("THREAD")
                .short("t")
                .long("threads")
                .value_name("int")
                .help("number of threads for reading BAM")
                .takes_value(true)
                .required(false)
                .default_value(read_threads))
            .arg(Arg::with_name("PRECISION")
                .short("p")
                .long("precision")
                .value_name("int")
                .takes_value(true)
                .default_value("5")
                .help("the required precision in bp for split-read derived fusion/integration points \
                Note: increasing might help to pool events in close proximity. \
                This might fuse though independent events and the lowest possible value might be advisable"))
            .arg(Arg::with_name("REF")
                .short("r")
                .long("reference")
                .value_name("int")
                .help("reference in fasta format provided, required for vcf output")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("VCF")
                .short("f")
                .long("vcf")
                .value_name("int")
                .help("generate vcf encoded fusion call, requires reference")
                .takes_value(true)
                .requires_all(&["REF"])
                .required(false))
            .arg(Arg::with_name("RANGE")
                .short("z")
                .long("range")
                .value_name("int")
                .takes_value(true)
                .default_value("500")
                .help("the allowed range in bp within which a split pair is considered support for a fusion point. Applies
                for both sides of the fusion")))
        .subcommand(SubCommand::with_name("id_fusion_reads")
                .about("classifying potential fusion derived reads")
                .arg(Arg::with_name("BAM")
                .short("b")
                .long("bam")
                .value_name("FILE")
                .help("long read/cDNA alignment vs reference")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("SIM")
                .short("i")
                .long("similarity")
                .value_name("float")
                .help("minimum seq similarity of matches [default:98.0]")
                .takes_value(true)
                .required(false))
            .arg(Arg::with_name("BAMBAM")
                .short("x")
                .long("bambam")
                .value_name("IDH")
                .takes_value(false)
                .required(false)
                .hidden(true))	
            .arg(Arg::with_name("THREAD")
                .short("t")
                .long("threads")
                .value_name("int")
                .help("number of threads for multi-threading")
                .takes_value(true)
                .required(false)
                .default_value("10"))
            .arg(Arg::with_name("REF")
                .short("r")
                .long("source")
                .value_name("STRING")
                .help("define the source of your transcriptome data. \
                This is key to identify if 2 FASTA entries are different transcript but from same gene or \
                if they are indeed 2 different genes and are indicative of a fusion")
                .takes_value(true)
                .required(true)
                .possible_values(&["CHESS"])
                .default_value("CHESS")))
		.get_matches();

    if let Some(matches) = matches.subcommand_matches("fragment") {
        run_fragment_base(matches,&args_string);
    }else if let Some(matches) = matches.subcommand_matches("paired") {
        run_paired_base(matches,&args_string);
    }else if let Some(matches) = matches.subcommand_matches("id_fusion_reads") {
        run_id_reads(matches);
    }else{
        eprintln!("Please choose one of the sub-commands, specify --help for more information");
    }
}



#[cfg(test)]
mod tests {
    //use std::collections::HashMap;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use genefusion::lib::common::{*};
    use genefusion::lib::hts_lib_based::{*};
    use rustc_hash::FxHashMap;
   // use crate::splitty;
    use crate::StrandDirection::Rev;
    use crate::StrandDirection::Fwd;

    // convenience function as the data is used
    // in multiple tests
    // we have the fullFusion evidence, princpe fusion evidence non annotated, the GTF evidence and the annotated fusion evidence and one with annotation but nothing matching
    fn build_basic_evidence () -> (FullFusionEvidence,PrincipleReadBasedFusionOutput,FxHashMap<String, Vec<Gtf>>,PrincipleReadBasedFusionOutput,PrincipleReadBasedFusionOutput){
        let n1 = FullFusionEvidence { 
            cigar_pr: String::from("15=1I169=14123N61=6233N103=50300N45=1I78=1I17=13610N95=1I5=33472N123=204S"), 
            cigar_sa: String::from("715S137=1I66="), 
            chrom_pr: String::from("chr3"), 
            chrom_sa: String::from("chr3"), 
            query_name: String::from("m64143_220818_104322/55640273/ccs;full_length_coverage=0;length=951;NpolyAremoved=32;FL"), 
            query_length: 919, 
            target_break_pr: 115520086, 
            target_break_sa: 115792300, 
            query_break_pr: 204, 
            query_break_sa: 203, 
            strand_pr: Rev, 
            strand_sa: Rev, 
            orient_pr_upstream: false, 
            orient_sa_upstream: true, 
            aln_start_pr: 115401638, 
            aln_start_sa: 115792502, 
            precision: Some(0), 
            support: 1 
        };
        let n2 = PrincipleReadBasedFusionOutput { 
            query_name: String::from("m64143_220818_104322/55640273/ccs;full_length_coverage=0;length=951;NpolyAremoved=32;FL"),
            length: 919,
            fusion_point: 203, 
            fusion_genes: vec![String::from("chr3--chr3")], 
            fp_fuzzy: Some(0), 
            fusion_distance: Some(272214), 
            a_chrom: String::from("chr3"), 
            a_start: 115792300, 
            a_end: 115792502, 
            a_fp: 115792300,
            a_strand: Rev,
            b_chrom: String::from("chr3"),
            b_start: 115401638,
            b_end: 115520086,
            b_strand: Rev,
            b_fp: 115520086
        };
        // just the needed info, which is by default first gene_name
        let mut n3 : FxHashMap<String, Vec<Gtf>> = FxHashMap::default();
        let mut att : FxHashMap<String,String> =  FxHashMap::default();
        att.insert(String::from("gene_name"), String::from("AC026341.1"));
        let gtf_entry = Gtf { 
            chromosome: String::from("chr3"), 
            source: String::from("HAVANA"), 
            feature: String::from("transcript"), 
            start: 115147641,
            end: 115564389,
            score: String::from("."),
            strand: Fwd,
            frame: String::from("."),
            attribute: att,
        };
        debug!("{:?}\n",n3);
        n3.insert(String::from("chr3"), vec![gtf_entry]);
        let n4 = PrincipleReadBasedFusionOutput { 
            query_name: String::from("m64143_220818_104322/55640273/ccs;full_length_coverage=0;length=951;NpolyAremoved=32;FL"),
            length: 919,
            fusion_point: 203, 
            fusion_genes: vec![String::from("NA--AC026341.1")], 
            fp_fuzzy: Some(0), 
            fusion_distance: Some(272214), 
            a_chrom: String::from("chr3"), 
            a_start: 115792300, 
            a_end: 115792502, 
            a_fp: 115792300,
            a_strand: Rev,
            b_chrom: String::from("chr3"),
            b_start: 115401638,
            b_end: 115520086,
            b_strand: Rev,
            b_fp: 115520086
        };
        let n5 = PrincipleReadBasedFusionOutput { 
            query_name: String::from("m64143_220818_104322/55640273/ccs;full_length_coverage=0;length=951;NpolyAremoved=32;FL"),
            length: 919,
            fusion_point: 203, 
            fusion_genes: vec![String::from("NA--NA")], 
            fp_fuzzy: Some(0), 
            fusion_distance: Some(272214), 
            a_chrom: String::from("chr3"), 
            a_start: 115792300, 
            a_end: 115792502, 
            a_fp: 115792300,
            a_strand: Rev,
            b_chrom: String::from("chr3"),
            b_start: 115401638,
            b_end: 115520086,
            b_strand: Rev,
            b_fp: 115520086
        };
        (n1, n2, n3,n4,n5)
    }

    /// just a simple test to verify that we can get the gene-name from our function accordingly
    #[test]
    fn get_gene_name_simple1(){
        // getting all the shared values from our convenience function
        // no annotation here needed though
        let (_fusion_evidence,_entry,gtf,_entr_annotated, _entry_nana) = build_basic_evidence();
        let entry:FxHashMap<String,String>  = gtf.get("chr3").expect("ERROR: couldnt find chrom!")[0].attribute.clone();
        let name= get_gene_name(&entry);
        let truth= String::from("AC026341.1");
        assert_eq!(name,truth);
    }

    /// just a simple test to verify that we can get the gene-name from our function accordingly
    /// if no match is found
    #[test]
    fn get_gene_name_simple2(){
        let entry:FxHashMap<String,String>  = FxHashMap::default();
        let name= get_gene_name(&entry);
        let truth= String::from("NA");
        assert_eq!(name,truth);
    }
    //////////////////////////////////////////////////////////////////////////
    /// test that the produced results from fragment annotations are good ////
    //////////////////////////////////////////////////////////////////////////
    #[test]
    fn format_and_annotate_read_based_fusions_noannot() {
        let mut read_collection: FxHashMap<String, Vec<FullFusionEvidence>> = FxHashMap::default();

        // getting all the shared values from our convenience function
        // no annotation here needed though
        let (fusion_evidence,entry,_gtf,_entr_annotated, _entry_nana) = build_basic_evidence();
        let mut truth_collection: Vec<PrincipleReadBasedFusionOutput> = Vec::new();
        truth_collection.push(entry);
        read_collection.insert(String::from("m64143_220818_104322/55640273/ccs;full_length_coverage=0;length=951;NpolyAremoved=32;FL") ,vec![fusion_evidence]);
        
        // prep finished, run function
        let result = format_and_annotate_read_based_fusions(read_collection,None,false,10000);
        assert_eq!(truth_collection,result);
    }

    #[test]
    fn format_and_annotate_read_based_fusions_wannot() {
        let mut read_collection: FxHashMap<String, Vec<FullFusionEvidence>> = FxHashMap::default();

        // getting all the shared values from our convenience function
        // no annotation here needed though
        let (fusion_evidence,_entry,gtf,entr_annotated,_entry_nana) = build_basic_evidence();
        let mut truth_collection: Vec<PrincipleReadBasedFusionOutput> = Vec::new();
        truth_collection.push(entr_annotated);
        read_collection.insert(String::from("m64143_220818_104322/55640273/ccs;full_length_coverage=0;length=951;NpolyAremoved=32;FL") ,vec![fusion_evidence]);
        

        //, gtf: Some({ }    
        let result = format_and_annotate_read_based_fusions(read_collection,Some(gtf),false,10000);
    
        assert_eq!(truth_collection,result);
    }

    // this one assumes that there is actually annotation but it does not
    // cover the region of interest, so the annotation in that case should be "NA--NA"
    #[test]
    fn format_and_annotate_read_based_fusions_lackingannot() {
        let mut read_collection: FxHashMap<String, Vec<FullFusionEvidence>> = FxHashMap::default();

        // getting all the shared values from our convenience function
        // no annotation here needed though
        let (fusion_evidence,_entry,_gtf,_entr_annotated,entry_nana) = build_basic_evidence();
        let mut truth_collection: Vec<PrincipleReadBasedFusionOutput> = Vec::new();
        truth_collection.push(entry_nana);
        read_collection.insert(String::from("m64143_220818_104322/55640273/ccs;full_length_coverage=0;length=951;NpolyAremoved=32;FL") ,vec![fusion_evidence]);
        
        let gtf : FxHashMap<String, Vec<Gtf>> = FxHashMap::default();
        let result = format_and_annotate_read_based_fusions(read_collection,Some(gtf),false,10000);
        
        assert_eq!(truth_collection,result);
    }

    // This one tests if the annotation works properly and adds a gene names intead
    // of the chr3 annotations to the name
    #[test]
    fn add_fusion_annotation_read_based_1(){
        // getting all the shared values from our convenience function
        // no annotation here needed though
        let (_fusion_evidence,entry,gtf,entr_annotated,_entry_nana) = build_basic_evidence();
        println!("Whats the gtf {:?}",&gtf);
        let annotated = add_fusion_annotation_read_based(entry,&Some(gtf),true); 
        assert_eq!(vec![entr_annotated],annotated);
    }
}

 