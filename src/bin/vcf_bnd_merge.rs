//! ## vcf_bnd_merge ## 
//! -------------------
//! This tool is designed for splitty based fusion calling and combines fragment- and 
//! position-based entries into one vcf file. Importantly one needs to supply the fragment 
//! based one as FileA and the position based one as FileB. This will be checked and fails
//! otherwise.
//! 
//! 
use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
use std::env;
use genefusion::lib::common::VersionInfo;
// our library which is within the same project
extern crate genefusion;
use genefusion::lib::hts_lib_based::{*};
// essentially we will do for the 1st vcf a first run population
// a hashmap were the keys are concatenated iD:mateID + the link to the vcf entry 
// together with a structure containging very simplified intormation 
// In a second step we do the same for the 2nd vcf but instead of the 2nd structure
// we compare it directly to the 1st strucure and add it accordingly 
extern crate pretty_env_logger;


fn main() {
    pretty_env_logger::init();

    // now the next is not really for any argument 
    // parsing but simply to get the command which 
    // was used to execute as I cant get this from clap
    let args: Vec<String> = env::args().collect();
    let args_string = args.join(" ");
    let matches = app_from_crate!()
            .about("This tool combines 2 vcf files describing BND. SVTYPE needs to be BND and ID + MATEID has to be set. \
            Important: if you have Fragment and read derived results, please pass them in this order ignoring the \"--simple\" and \"--flag\". \
            To simply check if entries in FileA were supported in FileB use the \"--simple\" flag. \
            In this case the \"--flag\" field will be added as INFO field and be set to 0 if no support was found and 1 if support was detected. \
            One can provide further INFO fields which will be propagated accordingly in a comma separated list. Currently this requires though \
            elements which have a single entrie. E.g. CIPOS which provides 2 numbers typically are not supported.

            Scoring: if only entryA/B but not both are present, low confidence events will get the label \"SPLITTY_FAILED\" . \
            If both evidence are combined the filter label is set to \"PASS\".
            Otherwise no filter field is applied.
            IMPORTANT: currently the strand is not verified of the events!")
            .arg(Arg::with_name("A")
                .short("a")
                .long("FileA")
                .value_name("vcf")
                .help("the first file used to compare BND events")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("B")
                .short("b")
                .long("FileB")
                .value_name("vcf")
                .help("the 2nd file used to compare BND events")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("OUT")
                .short("o")
                .long("out")
                .value_name("vcf")
                .help("the combined out vcf")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("THREADS")
                .short("t")
                .long("threads")
                .value_name("int")
                .help("vcf reading threads (files need to be huge to make an impact)")
                .takes_value(true)
                .default_value("2")
                .required(false))
            .arg(Arg::with_name("CIPOS")
                .short("c")
                .long("cipos")
                .value_name("int")
                .help("the distance within which events are combined")
                .takes_value(true)
                .required(false)
                .default_value("500"))
            .arg(Arg::with_name("INFO_FIELD")
                .short("d")
                .long("info-fields")
                .value_name("string")
                .help("a comma separated list of INFO fields to propagate, fomrat of field must be Integer")
                .takes_value(true)
                .required(false)
                .default_value("SR"))
            .arg(Arg::with_name("SIMPLE")
                .short("s")
                .long("simple")
                .value_name("bool")
                .help("simple true/false comparison if fields in A are supported in B \
                    one can set a flag name which is then added as TRUE if existent in B")
                .takes_value(false)
                .required(false))
            .arg(Arg::with_name("FLAG")
                .short("f")
                .long("flag")
                .value_name("string")
                .help("name of the flag value in case of \"simple\"")
                .takes_value(true)
                .required(false)
                .default_value_if("SIMPLE",None,"SUPPORTED"))    
            .get_matches();

let cipos: u32   = matches.value_of("CIPOS").unwrap().parse::<u32>().unwrap();
let file_a  = matches.value_of("A").unwrap();
let file_b  = matches.value_of("B").unwrap();
let file_o  = matches.value_of("OUT").unwrap();
let threads = matches.value_of("THREADS").unwrap().parse::<usize>().unwrap();
let simple: bool = matches.is_present("SIMPLE");
let simple_flag: Option<&str> = matches.value_of("FLAG");
let info_field_raw: &str= matches.value_of("INFO_FIELD").unwrap();
//let info_frmt_raw: &str = matches.value_of("INFO_FRMT").unwrap();

let info_field: Vec<&str> = info_field_raw.split(",").collect();
//let info_frmt: Vec<&str> = info_frmt_raw.split(",").collect();

// first we parse A and populate 2 structures
let (mut full_vcf_a, simpl_vcf_a) = vcf_a_parsing(
    file_a,
    threads,
    &simple
);
// next we parse B, populate 1 structure and already compare with 2nd A structure
let (mut full_vcf_b, simpl_vcf_a_modified) = vcf_b_parsing(
    file_b,
    cipos,
    simpl_vcf_a,
    threads,
    &simple
);
// before populating a vcf we need to prepare the new one
// with all the header options which we plan to have
let mut vcf = prepare_out_vcf(
    file_o,
    file_a,
    file_b,
    &VersionInfo {
        program: &String::from("vcf_bnd_merge"),
        version: crate_version!(),
        author:  crate_authors!(),
        command: &args_string,
    },
    &simple,
    &simple_flag,
    &info_field,
    //&info_frmt,
);
println!("{:?}",full_vcf_a);
vcf_combine_a_and_b(
    simpl_vcf_a_modified, 
    & mut full_vcf_a, 
    & mut full_vcf_b, 
    & mut vcf,
    &simple,
    &simple_flag,
    &info_field,
    //&info_frmt
);
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn merge_a_and_b_function(){

        use genefusion::lib::common::is_same_file;
        use std::path::Path;
        //let tmp_dir = tempdir().unwrap();
        //let tmp_path = tmp_dir.into_path();
        
        let cipos:u32 = 500;   
        let file_a  = "test/vcf_merge/mini_IntAssemblyResult.vcf";
        let file_b  = "test/vcf_merge/mini_inSilicoIntegrations.vcf";
        let file_verify =Path::new("test/vcf_merge/mini_verify.vcf");
        let file_o  = "test/vcf_merge/tmp.vcf";
        let threads = 1;
        let simple: bool = false;
        let simple_flag: Option<&str> = None;
        let info_field_raw: &str= "INSILICO_SP,INSILICO_SR";       
        let info_field: Vec<&str> = info_field_raw.split(",").collect();
        
        // first we parse A and populate 2 structures
        let (mut full_vcf_a, simpl_vcf_a) = vcf_a_parsing(
            file_a,
            threads,
            &simple
        );
        // next we parse B, populate 1 structure and already compare with 2nd A structure
        let (mut full_vcf_b, simpl_vcf_a_modified) = vcf_b_parsing(
            file_b,
            cipos,
            simpl_vcf_a,
            threads,
            &simple
        );
        // before populating a vcf we need to prepare the new one
        // with all the header options which we plan to have
        let mut vcf = prepare_out_vcf(
            &file_o,
            file_a,
            file_b,
            &VersionInfo {
                program: &String::from("vcf_bnd_merge"),
                version: "none",
                author:  "nobody",
                command: "unittest",
            },
            &simple,
            &simple_flag,
            &info_field,
            //&info_frmt,
        );
        println!("{:?}",full_vcf_a);
        vcf_combine_a_and_b(
            simpl_vcf_a_modified, 
            & mut full_vcf_a, 
            & mut full_vcf_b, 
            & mut vcf,
            &simple,
            &simple_flag,
            &info_field,
            //&info_frmt
        );
        let comparison = is_same_file(&file_verify,Path::new(file_o));
        assert_eq!(comparison.unwrap(), true);
        std::fs::remove_file(file_o).unwrap();
    }


}