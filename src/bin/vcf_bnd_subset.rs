//! ## vcf_bnd_subset ##
//! ---------------------
//! Will filter a paired BND based VCF based on provided paired ranges in a BEDPE file.
//! This allows to chain the clustering tools and extracting common events based on it's results for the original VCF files.

use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
use rustc_hash::FxHashMap;
use std::env;
//use rustc_hash::FxHashMap;
//use bio::data_structures::interval_tree;
// our library which is within the same project
extern crate genefusion;
//use genefusion::lib::hts_lib_based::{*};
use genefusion::lib::common::{*};
use genefusion::lib::hts_lib_based::get_header;
use bio::data_structures::interval_tree::{*};
use rust_htslib::bcf::Read as BcfRead;
use rust_htslib::bcf::Format;
use std::str::from_utf8;
use std::str;
use std::convert::TryFrom;
use std::error::Error;
use itertools::Itertools; // 0.8.2
use human_sort::compare;

extern crate pretty_env_logger;
#[macro_use] extern crate log;

// essentially we will do for the 1st vcf a first run population
// a hashmap were the keys are concatenated iD:mateID + the link to the vcf entry 
// together with a structure containging very simplified intormation 
// In a second step we do the same for the 2nd vcf but instead of the 2nd structure
// we compare it directly to the 1st strucure and add it accordingly 
fn subset_vcf_bnd (
    vcf_file : &str,
    out_vcf : Option<&str>, 
    bedpe_tree: BedpeTrees,
    inverse : bool,
    args: &VersionInfo
) ->   Result<(), Box<dyn Error>> {
    // So next we want now to parse the VCF file get it again in pairs.
    // Then we check if the pair is in one of the annotated regions pairs
    // and if this is the case we write the VCF entry in a new file
    // lets take now the input and read the VCF correctly:
    let mut vcf_reader       = rust_htslib::bcf::Reader::from_path(vcf_file).expect("ERROR: could not open vcf file!");
    let vcf_input_header     = get_header(&vcf_reader);
    let mut vcf_header       = rust_htslib::bcf::Header::from_template(&vcf_input_header);

    //let vcf_header : rust_htslib::bcf::Header = Header::from_template(vcf_reader.header());
    let header_cmd_line   = format!("{}{}",r#"##command_subset="# ,args.command);
    vcf_header.push_record(header_cmd_line.as_bytes());
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
    let samples     = vcf_reader.header().samples().to_vec();
    if samples.len() > 1 {
        panic!("ERROR: Currently only single sample vcf accepted!");
    };

    // now we need the FASTA writer ready and  iterate
    // over VCF entries. Important: must be sorted by BND pairs
    let mut bnd_pair_collection   : FxHashMap<String,Vec<rust_htslib::bcf::Record>> = FxHashMap::default();
    let mut bnd_single_collection : FxHashMap<String,Vec<rust_htslib::bcf::Record>> = FxHashMap::default();
    let mut bnd_indel_collection  : FxHashMap<String,Vec<rust_htslib::bcf::Record>> = FxHashMap::default();
    
    for entry in vcf_reader.records() {
        let mut cur_record = entry.expect("ERROR: could not read record of vcf file!");
        // here we have an ugly work around to fix a problem
        // form Niko's pipeline where entries contain ; in the field
        if cur_record.info(b"full_length_coverage").string().is_ok() {
            cur_record.clear_info_string(b"full_length_coverage").expect("ERROR: could not remove full_length_coverage");
        };
        if cur_record.info(b"length").string().is_ok() {
            cur_record.clear_info_string(b"length").expect("ERROR: could not remove length");
        };
        if cur_record.info(b"NpolyAremoved").string().is_ok(){
            cur_record.clear_info_string(b"NpolyAremoved").expect("ERROR: could not remove NpolyAremoved");
        };
        if cur_record.info(b"FL").string().is_ok() {
            cur_record.clear_info_string(b"FL").expect("ERROR: could not remove FL");
        };
        
        // if cur_record.info(b"length").string().unwrap().is_some()
        // if cur_record.info(b"NpolyAremoved").string().unwrap().is_some()
        // if cur_record.info(b"FL").string().unwrap().is_some()
        // define if a potential pair, single BND or INDEL 
        let sv_full       = cur_record.info(b"SVTYPE").string().unwrap().unwrap();
        let sv_rec        = str::from_utf8(sv_full[0]).expect("ERROR: could not extract SVTYPE string!");
        let id = &cur_record.id() ;
        let id1  = from_utf8(id).expect("ERROR: cant extract record id!");
        if id1 == "." {
            panic!("ERROR: all entries must have an annotated entry ID!");
        };
        let sv_type = match sv_rec{
            "INS" => SVType::INS,
            "DEL" => SVType::DEL,
            "BND" => {
                    if cur_record.info(b"MATEID").string().unwrap().is_some() {
                        SVType::BndPair
                    }else{
                        SVType::BndSingle
                    }
                },
            _ => panic!("ERROR: Currently only these types are supported: BND_PAIR, INS, DEL !"),
        };
        debug!("SV type: {:?}",sv_type);
        let key = id1.trim_end_matches(&['1', '2']);
        debug!("SV key: {:?}",key);
        match sv_type {
            SVType::INS       => bnd_indel_collection.entry(key.to_string()).or_default().push(cur_record),
            SVType::BndPair   => bnd_pair_collection.entry(key.to_string()).or_default().push(cur_record),
            SVType::BndSingle => bnd_single_collection.entry(key.to_string()).or_default().push(cur_record),
            SVType::DEL       => bnd_indel_collection.entry(key.to_string()).or_default().push(cur_record),
            _ => panic!("ERROR: Currently only these types are supported: BND_PAIR, INS, DEL !"),
        }

    }
    // currently we treat only the BND pair events
    debug!("BND PAIR Collection {:?}", bnd_pair_collection);
    debug!("BND SINGLE Collection {:?}", bnd_single_collection);
    debug!("BND INDEL Collection {:?}", bnd_indel_collection);
    eprintln!("INFO: VCF parsed entries: paired BND {}; single BND {}; INDEL {}", bnd_pair_collection.len(),bnd_single_collection.len(), bnd_indel_collection.len());
    // Now we sort the paired entries
    for (_,element) in bnd_pair_collection.iter_mut(){
        element.sort_by(|a,b| a.id().cmp(&b.id()));
    };
    let mut warned : bool = false;
    // now here we deploy natural sorting as we have otherwise very counterintuive
    // number sorting instead, then we loop over all entries and verify that the
    // pairs are present and matching
    let mut kept: usize = 0;
    for element in bnd_pair_collection.keys().sorted_by(|a,b| compare(a, b)) {
        let pair = bnd_pair_collection.get(element).expect("ERROR: could not fetch key in hashmap!");
        if pair.len() != 2 {
            panic!("ERROR: Not pairs found for all BND pair entryies!");
        }
        let cur_record = &pair[1];
        let prev_record = &pair[0];
        // Turns out first entry is 0 not 1,
        // so we need to count  

        let id_1 = &prev_record.id() ;
        let id_2 = &cur_record.id() ;
        let id1  = from_utf8(id_1).expect("ERROR: cant extract record id!");
        if id1 == "." {
            panic!("ERROR: all entries must have an annotated entry ID!");
        };
        let id2  = from_utf8(id_2).expect("ERROR: cant extract record id!");
        if id2 == "." {
            panic!("ERROR: all entries must have an annotated entry ID!");
        };
        let mate_id       = cur_record.info(b"MATEID").string().unwrap().unwrap();
        let mate_id2      = mate_id[0];
        if (id_1 != mate_id2) || (mate_id2.is_empty()) {
            eprintln!("Non-matching entries pair entries {:?} and {:?}",str::from_utf8(id_1),str::from_utf8(mate_id2));
            panic!("ERROR: VCF entries are either mal-formated or not sorted by matching IDs!");
        }
        // now we need to verify that we are indeed comparing proper pairs
        // we expect and ID .1 and ID .2 following each other
        // Otherwise we have either a different ordered file
        // which we cant treat here or we have a missing entry

        // next we get the chromosome names via the IDs from the header
        let rid1 = prev_record.rid().unwrap();
        let chr1 : &str  = str::from_utf8(prev_record.header().rid2name(rid1).expect("ERROR: could not convert rid to chromosome name!")).unwrap();
        let rid2 = cur_record.rid().unwrap();
        let chr2 : &str  = str::from_utf8(cur_record.header().rid2name(rid2).expect("ERROR: could not convert rid to chromosome name!")).unwrap();
        
        // here we get botht the a and b nucleotide
        // for both of the alleles provided
        let alleles1  = prev_record.alleles();
        let alleles2  = cur_record.alleles();
        // here we currently default these but will
        // add the posibility to have as well event types correctly 
        // classified
        if alleles1.len() != 2 || alleles2.len() != 2 {
            panic!("ERROR: encountered more than 2 alleles, which is not supported!")
        };
        // now we need again to get the string instead
        let allele_1 = str::from_utf8(alleles1[1]).expect("ERROR: could not convert allele!");
        let allele_2 = str::from_utf8(alleles2[1]).expect("ERROR: could not convert allele!");

        // now we have a function which  desconstructs the
        // 2nd allele info field in a BND event and returns 
        // already a orgaized structure instead
        // That means if we want the principle information of the primary allele we need
        // to take the complementary pair and get it's second allele
        // We need though later to use then the secondary direction!
        // Careful, coordinates are 1 based in the second allele
        let allele1_desconstr : Bnd2ndAllele = deconstruct_2nd_allele(allele_2,true);
        let allele2_desconstr : Bnd2ndAllele = deconstruct_2nd_allele(allele_1,true);
        // now we need to check if they are within a know pairs of entries
        // and if one matches if the associated match is as well having a proper similar pair
        // Any query will return a result iterator for the 1st of the 2 VCF records
        let tmp_tree = bedpe_tree.paired.get(chr1).expect("ERROR: could not fetch chromosome in tree!");
        //eprintln!("{:?}",tmp_tree);
        let result1 : Vec<Entry<'_, u64, ClusterPairs>> = tmp_tree.find(
            u64::try_from(prev_record.pos()).expect("ERROR: got non positive genomic coordinates!")
            ..
            u64::try_from(prev_record.pos()+1).expect("ERROR: got non positive genomic coordinates!")
        ).collect();
        //eprintln!("start {} end {}, result {:?}",prev_record.pos(),prev_record.pos()+1,result1);
        // now we iterate over them and verify if we find an
        // entry for the 1st VCF record for which we have then as well
        // a match for the 2nd VCF record
        // If one region matches we write directly and break, that way
        // we avoid duplicated entries
        let mut seen   : bool = false;
        
        for entry in result1 {
            // if already seen we just write and quit
            if seen {
                //eprintln!("Entry already encountered");
                break
            }
            let data = entry.data();
            // we check here the first entry but the order could be
            // the other way around. Therefore we have the else if
            // loop
            if (entry.interval().start == data.range1.start) && (entry.interval().end == data.range1.end) {
                // eprintln!("Found a matching entry");
                // now we first need to veify if the sense
                // is actually identical for the 1st pair, otherwise we break
                // We might encounter though in the BEDPE "unknown"
                // sense. I think in this case we want to warn the user but continue
                if allele1_desconstr.prim_direction != data.forward1 {
                    if  data.forward1 == StrandDirection::Unknown {
                        if !warned {
                            eprintln!("WARNING: your BEDPE contained unknown stranded and selection will be not verified for sense!");
                            warned = true ;
                        }
                    } else {
                        continue
                    }
                }
                //eprintln!("Direction is good");
                // here we compare now with the 2nd pair
                if chr2 == data.chr2 &&  data.range2.contains(&u64::try_from(cur_record.pos()).expect("ERROR: got non positive genomic coordinates!")) && 
                    ( allele2_desconstr.prim_direction == data.forward2 || data.forward2 == StrandDirection::Unknown ){
                        //eprintln!("Seen now");
                        seen = true ;
                        continue
                }
            }else if (entry.interval().start == data.range2.start) && (entry.interval().end == data.range2.end) {
                debug!("Checking out second pair");
                // here we compare now with the 1st pair
                // is actually identical for the 2st pair, otherwise we break
                if allele2_desconstr.prim_direction != data.forward2 {
                    if data.forward2 == StrandDirection::Unknown {
                        if !warned {
                            eprintln!("WARNING: your BEDPE contained unknown stranded and selection will be not verified for sense!");
                            warned = true ;
                        }
                    }else {
                        continue
                    }
                }
                if chr2 == data.chr1 && 
                    data.range1.contains(
                        &u64::try_from(cur_record.pos()).expect("ERROR: got non positive genomic coordinates!")
                    ) && 
                    (allele1_desconstr.prim_direction == data.forward1 || data.forward1 == StrandDirection::Unknown )
                    {
                        seen = true ;
                        continue
                }
            }else{
                panic!("ERROR: found incompatible ranges and ClusterPairs!");
            }
        }
        // now we check if both have been seen arccordingly, then we can
        // write our VCF entry accordingly.
        if (seen && !inverse) | (!seen && inverse ) {
            writer.write(prev_record).expect("ERROR: could not write VCF records!");
            writer.write(cur_record).expect("ERROR: could not write VCF records!");
            kept+=1;
        }
    
    }   
    eprintln!("INFO: VCF records kept and written: {}", kept);
    Ok(())
}

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
        .about("This tool allows to subset a VCF file based on ranges in a BEDPE file.\n\
        It's primary application is paired BND events for which the SVTYPE field must be properly labeled and an ID + MATEID must be set.\
        These pairs must have a '.1' and '.2' ending pair BND ID. Single BND events and INDELS are not yet fully implemented.\
        Resulting VCF files will be sorted by BND PAIR entry IDs, not genomic positions.
        ")
        .arg(Arg::with_name("VCF")
                .short("a")
                .long("vcf")
                .value_name("FILE")
                .help("the VCF file to subset")
                .takes_value(true)
                .required(true))
        .arg(Arg::with_name("BEDPE")
                .short("b")
                .long("bedpe")
                .value_name("FILE")
                .help("BEDPE format file indicating 2 regions per line")
                .takes_value(true)
                .required(true))
        .arg(Arg::with_name("CHROM")
                .short("c")
                .long("chroms")
                .value_name("FILE")
                .help("tab-separated file with name in 1st column and chrom length in second column")
                .takes_value(true)
                .required(true))
        .arg(Arg::with_name("OUT")
                .short("o")
                .long("out")
                .value_name("FILE")
                .help("outfile name for the VCF")
                .takes_value(true)
                .required(false))
        .arg(Arg::with_name("SCORE")
                .short("s")
                .long("score")
                .value_name("INT")
                .help("if set, allows to filter BEDPE regions based on scores, will keep only entries higher than this (e.g. number of observations)")
                .takes_value(true)
                .required(false))
        .arg(Arg::with_name("INVERSE")
                .short("i")
                .long("inverse")
                .value_name("BOOL")
                .help("will do inverse selection, reporting all events except regions in BEDPE file")
                .takes_value(false))
        .get_matches();

    let vcf_file     = matches.value_of("VCF").unwrap();
    let chroms       = matches.value_of("CHROM").unwrap();
    let contigs      = parse_chrom_file(chroms);
    let score_filter = matches.value_of("SCORE");
    let out_vcf      = matches.value_of("OUT");
    let bedpe        = matches.value_of("BEDPE").unwrap();
    let inverse      = matches.is_present("INVERSE");
    
    let infos = &VersionInfo {
        program: &String::from("bionano_translate"),
        version: crate_version!(),
        author: crate_authors!(),
        command: &args_string,
    };

    // now get the BEDPE file and parse it . We convert it into
    // pairs of coordinate entries and allow filtering based on a score.
    // If the score comes from our cluster tool it actually indicates
    // number of samples with that region and it can be used to filter
    // directly afterwards based on that
    // This returns then a hasmap with chromosomes as key and values being 
    // interval trees with all observed ranges and the pair information as 
    // value.
    let bedpe_tree   = svtree_from_bedpe(bedpe,score_filter,&contigs);
    subset_vcf_bnd(vcf_file, out_vcf, bedpe_tree, inverse, infos).expect("ERROR: failed to subset the VCF file!");
    
}


#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use genefusion::lib::common::{*};
    //use genefusion::lib::hts_lib_based::{*};
    //use bio::data_structures::interval_tree::{*};
    use crate::subset_vcf_bnd;


    // Note: the file test/vcf_subset/all_bundled_unidir.bedpe
    // and the directional was generated with 100kbp extend

    //////////////////////////////////////////////////////////////////////////
    ///  Test1: tests if we have a proper expected subsetting,  //////////////
    ///   both for uni-directional and directional with all     //////////////
    ///   that should accordingly be reported and verified      //////////////
    ///   against manual curated entries                        //////////////
    //////////////////////////////////////////////////////////////////////////
    /// 
    /////////////////////////////////////////
    ///       UNIDIRECTIONAL     ////////////
    /////////////////////////////////////////
    /// 
    #[test]
    fn test1_subset_unidir_rr() {
        pretty_env_logger::init();

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/all_bundled_unidir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_rr_2.vcf");
        let out_vcf  = String::from("test/vcf_subset/check_rr_uni_2.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test1/manual_rr_uni_2.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }

    #[test]
    fn test1_subset_unidir_ff() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author : &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/all_bundled_unidir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_ff_1.vcf");
        let out_vcf  = String::from("test/vcf_subset/check_ff_uni_1.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test1/manual_ff_uni_2.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }

    #[test]
    fn test1_subset_unidir_rf() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/all_bundled_unidir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_rf_4.vcf");
        let out_vcf  = String::from("test/vcf_subset/check_rf_uni_4.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test1/manual_rf_uni_4.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        //std::fs::remove_file(out_vcf).unwrap();

    }

    #[test]
    fn test1_subset_unidir_fr() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/all_bundled_unidir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_fr_3.vcf");
        let out_vcf  = String::from("test/vcf_subset/check_fr_uni_3.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test1/manual_fr_uni_3.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }

    #[test]
    fn test_subset_unidir_cornercase1() {
        // cornercase discovered by Nikos, if the cluster is supposed to be happening
        // at the very last basepair of the contig
        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test//vcf_subset/cornercase_1.bedpe");
        let chrom_in = String::from("test/vcf_cluster/GRCh38.primary_assembly.genome.chr.txt");
        let vcf_in   = String::from("test/vcf_cluster/test_file.vcf");
        let out_vcf  = String::from("test/vcf_subset/cornercase1.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/manual_corncercase1.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }

    //////////////////////////////////////
    ///       DIRECTIONAL     ////////////
    //////////////////////////////////////
    ///
    #[test]
    fn test1_subset_dir_rr() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/all_bundled_dir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_rr_2.vcf");
        let out_vcf  = String::from("test/vcf_subset//test1/check_rr_dir_2.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test1/manual_rr_dir_2.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }

    #[test]
    fn test1_subset_dir_ff() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/all_bundled_dir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_ff_1.vcf");
        let out_vcf  = String::from("test/vcf_subset/test1/check_ff_dir_1.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test1/manual_ff_dir_2.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }

    #[test]
    fn test1_subset_dir_rf() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/all_bundled_dir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_rf_4.vcf");
        let out_vcf  = String::from("test/vcf_subset//test1/check_rf_dir_4.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test1/manual_rf_dir_4.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }

    #[test]
    fn test1_subset_dir_fr() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/all_bundled_dir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_fr_3.vcf");
        let out_vcf  = String::from("test/vcf_subset//test1/check_fr_dir_3.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test1/manual_fr_dir_3.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }
    /////////////////////////////////////////////////////////////////////////////
    ///  Test2: tests that for a directional setting the opposite  //////////////
    ///   event wont be reported                                   //////////////
    ///  We expect for all tests here an empty VCF, only headers   //////////////
    /////////////////////////////////////////////////////////////////////////////
    /// 
    //////////////////////////////////////
    ///       DIRECTIONAL     ////////////
    //////////////////////////////////////
    #[test]
    fn test2_subset_dir_rr() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/ff_dir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_rr_2.vcf");
        let out_vcf  = String::from("test/vcf_subset//test2/check_rr_dir_2.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test2/manual_rr_dir_2.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }

    #[test]
    fn test2_subset_dir_ff() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/rr_dir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_ff_1.vcf");
        let out_vcf  = String::from("test/vcf_subset//test2/check_ff_dir_1.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test2/manual_ff_dir_1.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }

    #[test]
    fn test2_subset_dir_rf() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/rr_dir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_rf_4.vcf");
        let out_vcf  = String::from("test/vcf_subset//test2/check_rf_dir_4.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test2/manual_rf_dir_4.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }

    #[test]
    fn test2_subset_dir_fr() {

        let infos = &VersionInfo {
            program: &String::from("bionano_translate"),
            version: &String::from("0.0.0"),
            author: &String::from("nobody"),
            command: &String::from("unittest"),
        };        
        let bedpe_in = String::from("test/vcf_subset/rr_dir.bedpe");
        let chrom_in = String::from("test/vcf_subset/chroms.txt");
        let vcf_in   = String::from("test/data/manual_fr_3.vcf");
        let out_vcf  = String::from("test/vcf_subset//test2/check_fr_dir_3.vcf");
        let contigs      = parse_chrom_file(&chrom_in);

        let tree =    svtree_from_bedpe(&bedpe_in,None,&contigs);
        subset_vcf_bnd(&vcf_in, Some(&out_vcf), tree, false, infos).expect("ERROR: failed to subset the VCF file!");
        let vcf_test = String::from("test/vcf_subset/test2/manual_fr_dir_3.vcf");
        let comparison = is_same_file(&vcf_test, &out_vcf);
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(out_vcf).unwrap();

    }
}