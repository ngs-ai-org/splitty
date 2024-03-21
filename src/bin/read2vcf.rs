//! ## read2vcf ##
//! ---------------------
//! Takes BED annotated reads with one line describing the reference and the second part describing the alien region.
//! Uses that information and describes the genomic location as a VCF standard formated single BND event.

use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
use std::env;
//use rustc_hash::FxHashMap;

//use bio::data_structures::interval_tree;
// our library which is within the same project
extern crate genefusion;
use genefusion::lib::hts_lib_based::{*};
use genefusion::lib::common::{*};
extern crate pretty_env_logger;




fn main() {
    pretty_env_logger::init();

    // now the next is not really for any argument 
    // parsing but simply to get the command which 
    // was used to execute as I cant get this from clap
    //let args: Vec<String> = env::args().collect();
    //let args_string = args.join(" ");
    let args: Vec<String> = env::args().collect();
    let args_string = args.join(" ");
    let matches = app_from_crate!()
    .about("
This tool is tailored to translated read/contig annotation in a BED file for an alien DNA integration into a valid VCF. \
It expects for any read/contig exactly 2 lines which provide the annotation of alien and human region. \
Using these lines it will determine which part of the read represents the alien and which the human part. \
Next, it will based on the mapping information provide the annotation where in the genome the integration took place. \
This happens in a single BND entry VCF format.

        \n\n
Example of a valid format: \n\n
m64143_211228_145337/63635503/ccs       0       6563    mouse_merged,muERV_merged       0       .
m64143_211228_145337/63635503/ccs       6562    8881    human_merged    +,-     chr10@46668689-46671026,chr10@47671023-47673341
m64143_211228_145337/65471531/ccs       9       1908    human_merged    -       chrX@38858884-38860780
m64143_211228_145337/65471531/ccs       1908    7207    mouse_merged,muERV_merged       0       .

If a line has as shown above has 2 entry this will reported as 2 separated VCF events.

    ")
    .arg(Arg::with_name("VCF")
            .short("o")
            .long("vcf-out")
            .value_name("FILE")
            .help("the VCF formated output")
            .takes_value(true)
            .required(true))
    .arg(Arg::with_name("BED")
            .short("i")
            .long("bed-in")
            .value_name("FILE")
            .help("BED formated input file")
            .takes_value(true)
            .required(true))
    .arg(Arg::with_name("REF")
            .short("r")
            .long("reference")
            .value_name("FILE")
            .help("the reference FASTA genome sequence")
            .takes_value(true)
            .required(true))
    .arg(Arg::with_name("HUMAN")
            .short("m")
            .long("human-match")
            .value_name("STR")
            .help("string or sub-string which identifies in the BED file the human entry")
            .takes_value(true)
            .default_value("human")
            .required(false))
    .arg(Arg::with_name("MOUSE")
            .short("n")
            .long("mouse-match")
            .value_name("STR")
            .help("string or sub-string which identifies in the BED file the mouse/muERV entry")
            .takes_value(true)
            .default_value("muERV")
            .required(false))
    .arg(Arg::with_name("LAP")
            .short("w")
            .long("overlap")
            .value_name("INT")
            .help("Number of bp allowed for features to overlap")
            .takes_value(true)
            .default_value("10")
            .required(false))
    .get_matches();

    let bed_file        = matches.value_of("BED").unwrap();
    let reference       = matches.value_of("REF").unwrap();
    let out_vcf         = matches.value_of("VCF").unwrap();
    let regex_human     = matches.value_of("HUMAN").unwrap();
    let regex_alien     = matches.value_of("MOUSE").unwrap();
    let overlap          = matches.value_of("LAP").unwrap().parse::<u64>().unwrap();

    // get the annotated reads now
    let collected_reads = read_2lines_bed(bed_file,regex_human,regex_alien,&overlap);
    
    //Prepare the VCF header
    let header_fragid_line  = r#"##INFO=<ID=Contig_ID,Number=1,Type=String, Description="Name of the queried fragment">"#;   
    let header_frag1_line   = r#"##INFO=<ID=ReadFragRef,Number=2,Type=Integer, Description="Region of read which is reference derived">"#;    
    let header_frag2_line   = r#"##INFO=<ID=ReadFragAlien,Number=2,Type=Integer, Description="Region of read which is alien">"#; 

    // now we get a VCF BND writer. It does though not contain
    // some fields in the header which we want to add afterwards
    let mut vcf_writer      = prep_bnd_vcf(
        &VersionInfo {
            program: &String::from("read2vcf"),
            version: crate_version!(),
            author : crate_authors!(),
            command: &args_string,
        },
        Some(vec![header_fragid_line,header_frag1_line,header_frag2_line]),
        reference,
        out_vcf
    );

    // Finally write single BND entries
    write_single_bnd_vcf(reference,& mut vcf_writer,collected_reads,"muERV-Integration")
    
}