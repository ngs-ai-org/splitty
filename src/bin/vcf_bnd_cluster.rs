//! ## vcf_sv_cluster ##
//! ---------------------
//! This tool expects a list of vcf files and will cluster events of similar BND TRA events
//! together and report all samples which contain said events.
use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
use std::env;
use rustc_hash::FxHashMap;
//use bio::data_structures::interval_tree;
// our library which is within the same project
extern crate genefusion;
use genefusion::lib::hts_lib_based::{*};
use genefusion::lib::common::{*};
use bio::data_structures::interval_tree::{*};

extern crate pretty_env_logger;
#[macro_use] extern crate log;

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
    //let args: Vec<String> = env::args().collect();
    //let args_string = args.join(" ");
    let matches = app_from_crate!()
    .about("This tool compares and extends BND based events for single or multiple samples. \
        It expects a list of vcf files as fofn or a single vcf which is parsed and organized in an inverval tree. \n\
        SVTYPE can be BND, INS and DEL. For BND ID + MATEID has to be set. \n\
        Events of same type are clustered together based on a defined maximal distance and potentially thereby extend the original event range. \n\
        The program produces a 0-based BEDPE file with the ranges of the clusters and members that have been observed within.")
    .arg(Arg::with_name("FILES")
            .short("f")
            .long("fofn")
            .value_name("FILE")
            .help("a file with files to parse, one entry per line")
            .takes_value(true)
            .required_unless("VCF"))
    .arg(Arg::with_name("VCF")
            .short("v")
            .long("vcf")
            .value_name("FILE")
            .help("a single vcf file to collapse and cluster entries")
            .takes_value(true)
            .required_unless("FILES"))
    .arg(Arg::with_name("CHROM")
            .short("c")
            .long("chroms")
            .value_name("FILE")
            .help("tab-separated file with name in 1st column and chrom length in second column")
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
    .arg(Arg::with_name("OUTFILE")
            .short("o")
            .long("bedpe")
            .value_name("FILE")
            .help("write to outfile instead stdout")
            .takes_value(true)
            .required(false)
            .hidden(true))
    .arg(Arg::with_name("RANGE")
            .short("r")
            .long("range")
            .value_name("int")
            .help("the distance within which events are combined")
            .takes_value(true)
            .required(false)
            .default_value("500"))
    .arg(Arg::with_name("MULTI")
            .short("m")
            .long("multi-only")
            .help("will only report multi event regions, no singletons")
            .takes_value(false)
            .required(false))
    .arg(Arg::with_name("PURGE")
            .short("p")
            .long("purge")
            .help("will only report one event per sample. Can be used with single files to remove redundant events, too")
            .takes_value(false)
            .required(false))
    .arg(Arg::with_name("UNI")
            .short("u")
            .long("uni-dir")
            .help("uni-directional, ignoring direction of events")
            .takes_value(false)
            .required(false))
    .get_matches();

    let range: u64 = matches.value_of("RANGE").unwrap().parse::<u64>().unwrap();
    let fofn      = matches.value_of("FILES");
    let vcf_in    = matches.value_of("VCF");
    let chroms    = matches.value_of("CHROM").unwrap();
    let threads  = matches.value_of("THREADS").unwrap().parse::<usize>().unwrap();
    let multi     = matches.is_present("MULTI");
    let purge     = matches.is_present("PURGE");
    let direction = matches.is_present("UNI");
    let contigs   = parse_chrom_file(chroms);
    let file_out        = matches.value_of("OUTFILE");

    // first we go through all VCF files and gnerate an interval tree which keeps for all the essential
    // information. It is organized within a hashmap where the keys are the chromosomes
    let sv_tree   : FullSvTrees = match (fofn,vcf_in){
        (Some(x),None) =>  vcf_parsing_tree(
                x,
                threads,
                &contigs,
                direction,
                true
            ),
        (None,Some(x)) =>  vcf_parsing_tree(
            x,
            threads,
            &contigs,
            direction,
            false
        ),
        _ => panic!("ERROR: provide a fofn or single vcf file as input!"),
    };
    debug!("FullTree: {:?}", sv_tree);
    
    ///////////////////////
    // PAIRED BND EVENTS //
    ///////////////////////
    
    // now we can start and get based on our allowed
    // distance between clusters the ranges which define 
    // a cluster of events:
    let bndpair_break_points : FxHashMap<String,Vec<std::ops::Range<u64>>> = chr_dist_based_groups(&sv_tree.bnd_pair,range,&contigs);
    debug!("BND_TRbndpair_break_pointsEE: {:?}",bndpair_break_points);
    let mut bnd_break_tree   : FxHashMap<String,IntervalTree<u64,String>>  = FxHashMap::default();
    
    for (key,value) in bndpair_break_points {
            let local_tree = cluster_to_tree(value);
            bnd_break_tree.insert(key,local_tree);
    };
    debug!("bnd_break_tree {:?}",bnd_break_tree);

    // the next step now defines pairs of grouped events
    let bnd_cluster_pairs   = bnd_pairs_clustered(
        &sv_tree.bnd_pair,
        &bnd_break_tree,
        &contigs
    );
    debug!("bnd_cluster_pairs {:?}",bnd_cluster_pairs);
    // here we find potential clusters which have events that do not match the
    // pair, so we refine the event from the cluster and refine again the boundaries
    let bnd_refined_pairs   = bnd_filter_pairs(bnd_cluster_pairs,range);
    debug!("refined {:?}",bnd_refined_pairs);
    // finally we collapse now the cluster and should have
    // for each event a pair that is finally returned
    let bnd_collapsed_pairs = bnd_collapse_pairs(bnd_refined_pairs);
    // now we can write the pairs in a BEDPE format of paired entries
    bedpe_collapsed_pairs(bnd_collapsed_pairs,multi,purge,file_out).expect("ERROR: could not write results!");

    //////////////////
    // INDEL EVENTS //
    //////////////////
    
    let ins_break_points   : FxHashMap<String,Vec<std::ops::Range<u64>>> = chr_dist_based_groups_indel(
        &sv_tree.indel,
        range,&contigs,
        Some(SVType::INS)
    );
    let del_break_points   : FxHashMap<String,Vec<std::ops::Range<u64>>> = chr_dist_based_groups_indel(
        &sv_tree.indel,
        range,&contigs,
        Some(SVType::DEL)
    );
    let mut ins_break_tree : FxHashMap<String,IntervalTree<u64,String>>  = FxHashMap::default();
    let mut del_break_tree : FxHashMap<String,IntervalTree<u64,String>>  = FxHashMap::default();
    
    for (key,value) in ins_break_points {
        let local_tree = cluster_to_tree(value);
        ins_break_tree.insert(key,local_tree);
    };
    for (key,value) in del_break_points {
        let local_tree = cluster_to_tree(value);
        del_break_tree.insert(key,local_tree);
    };

    // the next step now defines pairs of grouped events
    let clustered_ins = indels_clustered(
        &sv_tree.indel,
        &ins_break_tree,
        &contigs,
        &Some(SVType::INS)
    );
    let clustered_del = indels_clustered(
        &sv_tree.indel,
        &del_break_tree,
        &contigs,
        &Some(SVType::DEL)
    );
    // now we can write the pairs in a BEDPE format of paired entries
    bedpe_sv_clusters(clustered_ins,multi,purge,&String::from("INS"),file_out);
    bedpe_sv_clusters(clustered_del,multi,purge,&String::from("DEL"),file_out);
}




#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    //use bio::io::fasta::IndexedReader;  
    //use std::str::from_utf8;
    use std::iter::FromIterator;
    //use std::str::from_utf8;
    use tempfile::NamedTempFile;

    #[test]
    fn full_test(){
        pretty_env_logger::init();
        let range : u64 = 50;
        // this below is the kind of "unwanted" situation
        // which we later need to accound for
        //let test_tree : IntervalTree<i64, BNDentry> = IntervalTree::new();
        let test_tree_chr1 = IntervalTree::from_iter(
        vec![
            (1000..1300,BNDentry{
                    sample: String::from("S1"), 
                    id : String::from("ID.1.1"), 
                    mid: String::from("ID.1.2"),
                    chr1: String::from("chr1"),
                    fp1: 1200,
                    st1: 1000,
                    end1: 1300,
                    forward1: StrandDirection::Fwd,
                    chr2: String::from("chr7"),
                    fp2: 11275,
                    st2: 10050,
                    end2: 12500,
                    forward2: StrandDirection::Fwd,
                    sv_type : SVType::default(),
                    name: Some(String::from("Z--Y")),
            }),
            (900..1200, BNDentry{
                    sample: String::from("S2"), 
                    id : String::from("ID.2.1"), 
                    mid: String::from("ID.2.2"),
                    chr1: String::from("chr1"),
                    fp1: 1050,
                    st1: 900,
                    end1: 1200,
                    forward1: StrandDirection::Fwd,
                    chr2: String::from("chr7"),
                    fp2: 10050,
                    st2: 10000,
                    end2: 10100,
                    forward2: StrandDirection::Fwd,
                    sv_type : SVType::default(),
                    name: Some(String::from("Z--Y")),
                    }),
            (1000..1100, BNDentry{
                    sample: String::from("S3"), 
                    id : String::from("ID.3.1"), 
                    mid: String::from("ID.3.2"),
                    chr1: String::from("chr1"),
                    fp1: 1050,
                    st1: 1000,
                    end1: 1100,
                    forward1: StrandDirection::Fwd,
                    chr2: String::from("chr7"),
                    fp2: 10200,
                    st2: 10000,
                    end2: 10400,
                    forward2: StrandDirection::Fwd,
                    sv_type : SVType::default(),
                    name: Some(String::from("Z--Y")),
                    }),
            (1400..1600, BNDentry{
                    sample: String::from("S4"), 
                    id : String::from("ID.5.1"), 
                    mid: String::from("ID.5.2"),
                    chr1: String::from("chr1"),
                    fp1: 1500,
                    st1: 1400,
                    end1: 1600,
                    forward1: StrandDirection::Fwd,
                    chr2: String::from("chr7"),
                    fp2: 200,
                    st2: 100,
                    end2: 300,
                    forward2: StrandDirection::Fwd,
                    sv_type : SVType::default(),
                    name: Some(String::from("A--B")),
                    }),
            (1500..1600, BNDentry{
                    sample: String::from("S5"), 
                    id : String::from("ID.6.1"), 
                    mid: String::from("ID.6.2"),
                    chr1: String::from("chr1"),
                    fp1: 1550,
                    st1: 1500,
                    end1: 1600,
                    forward1: StrandDirection::Fwd,
                    chr2: String::from("chr7"),
                    fp2: 200,
                    st2: 100,
                    end2: 300,
                    forward2: StrandDirection::Fwd,
                    sv_type : SVType::default(),
                    name: Some(String::from("A--B")),
                    }),
            (1500..1700, BNDentry{
                    sample: String::from("S6"), 
                    id : String::from("ID.7.1"), 
                    mid: String::from("ID.7.2"),
                    chr1: String::from("chr1"),
                    fp1: 1600,
                    st1: 1500,
                    end1: 1700,
                    forward1: StrandDirection::Fwd,
                    chr2: String::from("chr7"),
                    fp2: 200,
                    st2: 100,
                    end2: 300,
                    forward2: StrandDirection::Fwd,
                    sv_type : SVType::default(),
                    name: Some(String::from("A--B")),
                    }),
            (5000..5500, BNDentry{
                    sample: String::from("S7"), 
                    id : String::from("ID.8.1"), 
                    mid: String::from("ID.8.2"),
                    chr1: String::from("chr1"),
                    fp1: 5250,
                    st1: 5000,
                    end1: 5500,
                    forward1: StrandDirection::Fwd,
                    chr2: String::from("chr7"),
                    fp2: 200,
                    st2: 100,
                    end2: 300,
                    forward2: StrandDirection::Fwd,
                    sv_type : SVType::default(),
                    name: None,
                    }),
            (1000..1250, BNDentry{
                    sample: String::from("S2"), 
                    id : String::from("ID.9.1"), 
                    mid: String::from("ID.9.2"),
                    chr1: String::from("chr1"),
                    fp1: 1100,
                    st1: 1000,
                    end1: 1250,
                    forward1: StrandDirection::Fwd,
                    chr2: String::from("chr7"),
                    fp2: 300,
                    st2: 120,
                    end2: 500,
                    forward2: StrandDirection::Fwd,
                    sv_type : SVType::default(),
                    name: None,
                    }),  
            (4800..10400, BNDentry{
                    sample: String::from("S3"), 
                    id : String::from("ID.10.1"), 
                    mid: String::from("ID.10.2"),
                    chr1: String::from("chr1"),
                    fp1: 5200,
                    st1: 4800,
                    end1: 6000,
                    forward1: StrandDirection::Rev,
                    chr2: String::from("chr7"),
                    fp2: 300,
                    st2: 10000,
                    end2: 10400,
                    forward2: StrandDirection::Fwd,
                    sv_type : SVType::default(),
                    name: None,
                    })                  
            
                    ].into_iter()
        );

        let test_tree_chr7 = IntervalTree::from_iter(
                vec![
                (10050..12500, BNDentry{
                        sample: String::from("S1"), 
                        mid : String::from("ID.1.1"), 
                        id: String::from("ID.1.2"),
                        chr2: String::from("chr1"),
                        fp2: 1200,
                        st2: 1000,
                        end2: 1300,
                        forward2: StrandDirection::Fwd,
                        chr1: String::from("chr7"),
                        fp1: 11275,
                        st1: 10050,
                        end1: 12500,
                        forward1: StrandDirection::Fwd,
                        sv_type : SVType::default(),
                        name: Some(String::from("Z--Y")),
                        }),
                (10000..10100, BNDentry{
                        sample: String::from("S2"), 
                        mid : String::from("ID.2.1"), 
                        id: String::from("ID.2.2"),
                        chr2: String::from("chr1"),
                        fp2: 1050,
                        st2: 900,
                        end2: 1200,
                        forward2: StrandDirection::Fwd,
                        chr1: String::from("chr7"),
                        fp1: 10050,
                        st1: 10000,
                        end1: 10100,
                        forward1: StrandDirection::Fwd,
                        sv_type : SVType::default(),
                        name: Some(String::from("Z--Y")),
                        }),
                (10000..10400, BNDentry{
                        sample: String::from("S3"), 
                        mid : String::from("ID.3.1"), 
                        id: String::from("ID.3.2"),
                        chr2: String::from("chr1"),
                        fp2: 1050,
                        st2: 1000,
                        end2: 1100,
                        forward2: StrandDirection::Fwd,
                        chr1: String::from("chr7"),
                        fp1: 10200,
                        st1: 10000,
                        end1: 10400,
                        forward1: StrandDirection::Fwd,
                        sv_type : SVType::default(),
                        name: Some(String::from("Z--Y")),
                        }),
                (100..300, BNDentry{
                        sample: String::from("S4"), 
                        mid : String::from("ID.5.1"), 
                        id: String::from("ID.5.2"),
                        chr2: String::from("chr1"),
                        fp2: 1500,
                        st2: 1400,
                        end2: 1600,
                        forward2: StrandDirection::Fwd,
                        chr1: String::from("chr7"),
                        fp1: 200,
                        st1: 100,
                        end1: 300,
                        forward1: StrandDirection::Fwd,
                        sv_type : SVType::default(),
                        name: Some(String::from("A--B")),
                        }),
                (100..300, BNDentry{
                        sample: String::from("S5"), 
                        mid : String::from("ID.6.1"), 
                        id: String::from("ID.6.2"),
                        chr2: String::from("chr1"),
                        fp2: 1550,
                        st2: 1500,
                        end2: 1600,
                        forward2: StrandDirection::Fwd,
                        chr1: String::from("chr7"),
                        fp1: 200,
                        st1: 100,
                        end1: 300,
                        forward1: StrandDirection::Fwd,
                        sv_type : SVType::default(),
                        name: Some(String::from("A--B")),
                        }),
                (100..300, BNDentry{
                        sample: String::from("S6"), 
                        mid : String::from("ID.7.1"), 
                        id: String::from("ID.7.2"),
                        chr2: String::from("chr1"),
                        fp2: 1600,
                        st2: 1500,
                        end2: 1700,
                        forward2: StrandDirection::Fwd,
                        chr1: String::from("chr7"),
                        fp1: 200,
                        st1: 100,
                        end1: 300,
                        forward1: StrandDirection::Fwd,
                        sv_type : SVType::default(),
                        name: Some(String::from("A--B")),
                        }),
                (100..300, BNDentry{
                        sample: String::from("S7"), 
                        mid : String::from("ID.8.1"), 
                        id: String::from("ID.8.2"),
                        chr2: String::from("chr1"),
                        fp2: 5250,
                        st2: 5000,
                        end2: 5500,
                        forward2: StrandDirection::Fwd,
                        chr1: String::from("chr7"),
                        fp1: 200,
                        st1: 100,
                        end1: 300,
                        forward1: StrandDirection::Fwd,
                        sv_type : SVType::default(),
                        name: None,
                        }),
                (120..500, BNDentry{
                        sample: String::from("S2"), 
                        mid : String::from("ID.9.1"), 
                        id: String::from("ID.9.2"),
                        chr2: String::from("chr1"),
                        fp2: 1100,
                        st2: 1000,
                        end2: 1250,
                        forward2: StrandDirection::Fwd,
                        chr1: String::from("chr7"),
                        fp1: 300,
                        st1: 120,
                        end1: 500,
                        forward1: StrandDirection::Fwd,
                        sv_type : SVType::default(),
                        name: None,
                        }),
                (4800..10400, BNDentry{
                        sample: String::from("S3"), 
                        id : String::from("ID.10.1"), 
                        mid: String::from("ID.10.2"),
                        chr2: String::from("chr1"),
                        fp2: 5200,
                        st2: 4800,
                        end2: 6000,
                        forward2: StrandDirection::Rev,
                        chr1: String::from("chr7"),
                        fp1: 300,
                        st1: 10000,
                        end1: 10400,
                        forward1: StrandDirection::Fwd,
                        sv_type : SVType::default(),
                        name: None,
                        }),     
                ].into_iter()
        );
            
        let schematic = r###"
                CHR1
            
            
            900--------S2:2.1-----------1200->
            1000-----S3:3.1------1100->
            1000-----S1:1.1------------------1300->
            1000------S2:9.1----------1250->
                                                            1400---S4:5.1-----1600->
                                                            1500--S5:6.1----1600->
                                                            1500----S6:7.1----------1700->
                                                                                                    5000-S7:-8.1--5500->
                                                                                                    <-4800---S3:10.1--------6000
            
            ----------------------------------------------------------------------------------------------------------------------------------
            
                CHR7
            
            
                                                                                10050-------------S1:1.2-------------------------------12500->
                                                                                10000----S2:2.2-----10100->
                                                                                10000------S3:3.2-------------10400->
                                                                                10000------S3:10.2------------10400->
            100----S4:5.2-----300->
            100----S5:6.2-----300->                   
            100----S6:7.2-----300->
            100----S7:8.2-----300->
                    120-----------S2:9.2-------500->                 
            "###;                                         
        debug!("Testing the following clustering schematics: {}",schematic);

        // the event 9.1/9.2 raises a typical problem where the cluster of 100--300 is extended till 500 due
        // to one event which defines initially a length of a group but turns later out to be a different pair
        let mut final_tree : FxHashMap<String,IntervalTree<u64,BNDentry>> = FxHashMap::default();
        final_tree.insert(String::from("chr1"), test_tree_chr1);
        final_tree.insert(String::from("chr7"), test_tree_chr7);

        // now we can start and get based on our allowed
        // distance between clusters the ranges which define 
        // a cluster of events:
        let mut contigs : FxHashMap<String, u64> = FxHashMap::default();
        contigs.insert(String::from("chr1"),100000_u64);
        contigs.insert(String::from("chr7"),100000_u64);

        let break_points :  FxHashMap<String,Vec<std::ops::Range<u64>>> = chr_dist_based_groups(&final_tree,range,&contigs);
        let mut break_tree : FxHashMap<String,IntervalTree<u64,String>> = FxHashMap::default();
        for (key,value) in break_points {
            let local_tree = cluster_to_tree(value);
            break_tree.insert(key,local_tree);
        };

        let cluster_pairs = bnd_pairs_clustered(&final_tree,&break_tree,&contigs);
        eprintln!("Cluster pairs: {:?}",cluster_pairs);

        // a above mentioned the problem is essentially now that the cluster from 100--300 was extended till 500 because
        // the event 9 overlaps. In the later step we identified though that these are not having the same matching pair region
        // therefore we need now in a second step now a way to redefine the cluster for a given pair in the
        // same fashion as initially
        // So we take each pair
        let refined_pairs = bnd_filter_pairs(cluster_pairs,range);
        eprintln!("Refined pairs: {:?}",refined_pairs);
        // at this point we finished operations and we collapse now the events for easier processing
        let collapsed_pairs = bnd_collapse_pairs(refined_pairs);
        let multi = false;
        let purge = false;  
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();

        bedpe_collapsed_pairs(collapsed_pairs,multi,purge,Some(path.to_str().unwrap())).expect("ERROR: could not generate the bedpe output!");
        //bedpe_collapsed_pairs(collapsed_pairs,multi,purge,None).expect("ERROR: could not generate the bedpe output!");
        
        let bedpe_test = String::from("test/test_cluster/full_model.bedpe");
        let comparison = is_same_file(&bedpe_test, path.to_str().unwrap());
        assert_eq!(comparison.unwrap(),true);
        std::fs::remove_file(path).unwrap();
    }



    //////////////////////////
    ///   vcf clustering   ///
    //////////////////////////
    

    // a way to define only once the necessary elements to get
    // the tree
    fn get_tree_stranded () -> FullSvTrees {
        let vcf_in   = String::from("test/data/cluster1.vcf");
        let chrom_in = String::from("test/test_cluster/chroms.txt");
        let contigs  = parse_chrom_file(&chrom_in);
        vcf_parsing_tree(
            &vcf_in,
            1,
            &contigs,
            false,
            false,
        )
    }

    // a way to define only once the necessary elements to get
    // the tree
    fn get_tree_unstranded () -> FullSvTrees {
        let vcf_in   = String::from("test/data/cluster1.vcf");
        let chrom_in = String::from("test/test_cluster/chroms.txt");
        let contigs  = parse_chrom_file(&chrom_in);
        vcf_parsing_tree(
            &vcf_in,
            1,
            &contigs,
            true,
            false,
        )
    }

    #[test]
    fn vcf_parsing_0(){
        // this one is a tricky corner case good for testing
        // the start position is 24 and uncertainty is larger with 32bp
        let vcf_in   = String::from("test/vcf_cluster/corner_case1.vcf");
        let chrom_in = String::from("test/vcf_cluster/corner_case1.length");
        let contigs  = parse_chrom_file(&chrom_in);
        let _tree = vcf_parsing_tree(
            &vcf_in,
            1,
            &contigs,
            true,
            false,
        );
    }

    #[test]
    fn vcf_parsing_1(){
        // coerner case with one entry being first base (1) inv 
        let vcf_in   = String::from("test/vcf_cluster/corner_case0.vcf");
        let chrom_in = String::from("test/vcf_cluster/corner_case2.length");
        let contigs  = parse_chrom_file(&chrom_in);
        let _tree = vcf_parsing_tree(
            &vcf_in,
            1,
            &contigs,
            true,
            false,
        );
    }


    #[test]
    fn vcf_parsing_2(){
        // this one is a weird corner case good for testing
        // dont even understand why it is failing at the moment
        let vcf_in   = String::from("test/vcf_cluster/corner_case2.vcf");
        let chrom_in = String::from("test/vcf_cluster/corner_case2.length");
        let contigs  = parse_chrom_file(&chrom_in);
        let _tree = vcf_parsing_tree(
            &vcf_in,
            1,
            &contigs,
            true,
            false,
        );
    }

    #[test]
    fn vcf_parsing_3(){
        // For this one, one of the positions is at the max length
        // of the contig and will have to be corrected for
        let vcf_in   = String::from("test/vcf_cluster/corner_case3.vcf");
        let chrom_in = String::from("test/vcf_cluster/corner_case3.length");
        let contigs  = parse_chrom_file(&chrom_in);
        let _tree = vcf_parsing_tree(
            &vcf_in,
            1,
            &contigs,
            false,
            false,
        );
    }

    #[test]
    #[should_panic]
    fn vcf_parsing_4(){
        // this one is a faulty vcf with a 0 position, should fail
        let vcf_in   = String::from("test/vcf_cluster/corner_case4.vcf");
        let chrom_in = String::from("test/vcf_cluster/corner_case2.length");
        let contigs  = parse_chrom_file(&chrom_in);
        let _tree = vcf_parsing_tree(
            &vcf_in,
            1,
            &contigs,
            true,
            false,
        );
    }

    #[test]
    fn vcf2tree_test_interval_tree() {
        // this test is not resulting in any clustered
        // events, everything becomes a single cluster
        // Additionally we have only 1 file and not multiple vcfs
        let tree = get_tree_stranded();
       
        let value = tree.bnd_pair.get("13").unwrap();
        println!("{:?}",value);

        let chrom_entries : Vec<Entry<'_, u64, BNDentry>> = value.find(123465..123466).collect();
        println!("{:?}",chrom_entries);
        let test1 = BNDentry{
                sample: String::from("cluster1"),
                id: String::from("bnd2.2"),
                mid: String::from("bnd2.1"),
                chr1: String::from("13"),
                fp1: 123465_u64,
                st1: 123465_u64,
                end1: 123466_u64,
                forward1: StrandDirection::Fwd,
                chr2: String::from("2"),
                fp2: 321691_u64,
                st2: 321691_u64,
                end2: 321692_u64,
                forward2: StrandDirection::Rev,
                sv_type: SVType::BndPair,
                name: None,
        };
        assert_eq!(chrom_entries[0].data(),&test1);
    }


    #[test]
    fn locus_cluster_tree_1() {
        // here we have simply 2 entries with a 1 base-pair shift
        // and by extending by +1 bp we allow a new grouping extension
        let tree = IntervalTree::from_iter(vec![
            (3320474..3320475,
            BNDentry { 
                sample: String::from("test2"), 
                id: String::from("NGSAI_GENEFUSION_PID.1629.2"), 
                mid: String::from("NGSAI_GENEFUSION_PID.1629.1"), 
                chr1: String::from("000232F"), 
                fp1: 3320474, 
                st1: 3320474, 
                end1: 3320475, 
                forward1: StrandDirection::Rev, 
                chr2: String::from("000232F"), 
                fp2: 3321632, 
                st2: 3321632, 
                end2: 3321633, 
                forward2: StrandDirection::Rev, 
                sv_type: SVType::BndPair, 
                name: None 
            }),
            (3320476..3320477,
            BNDentry { 
                sample: String::from("test2"), 
                id: String::from("NGSAI_GENEFUSION_PID.1628.2"), 
                mid: String::from("NGSAI_GENEFUSION_PID.1628.1"), 
                chr1: String::from("000232F"), 
                fp1: 3320476, 
                st1: 3320476, 
                end1: 3320477, 
                forward1: StrandDirection::Rev, 
                chr2: String::from("000232F"), 
                fp2: 3321632, 
                st2: 3321632, 
                end2: 3321633, 
                forward2: StrandDirection::Rev, 
                sv_type: SVType::BndPair, 
                name: None 
            })]);
        let local_max = 3320477 +1;
        let entries : Vec<Entry<'_, u64, BNDentry>> = tree.find(0..local_max).collect(); 
        let test :  Vec<std::ops::Range<u64>> = locus_dist_based_groups(entries,1);
        let results = vec![3320474..3320477];
        assert_eq!(test,results);
    }

    #[test]
    fn locus_cluster_tree_2() {
        // here we have simply 2 entries with a 1 base-pair shift
        // which should not be extended as we allow 0 bp extension
        // We create 2 independent entries instead
        let tree = IntervalTree::from_iter(vec![
            (3320474..3320475,
            BNDentry { 
                sample: String::from("test2"), 
                id: String::from("NGSAI_GENEFUSION_PID.1629.2"), 
                mid: String::from("NGSAI_GENEFUSION_PID.1629.1"), 
                chr1: String::from("000232F"), 
                fp1: 3320474, 
                st1: 3320474, 
                end1: 3320475, 
                forward1: StrandDirection::Rev, 
                chr2: String::from("000232F"), 
                fp2: 3321632, 
                st2: 3321632, 
                end2: 3321633, 
                forward2: StrandDirection::Rev, 
                sv_type: SVType::BndPair, 
                name: None 
            }),
            (3320476..3320477,
            BNDentry { 
                sample: String::from("test2"), 
                id: String::from("NGSAI_GENEFUSION_PID.1628.2"), 
                mid: String::from("NGSAI_GENEFUSION_PID.1628.1"), 
                chr1: String::from("000232F"), 
                fp1: 3320476, 
                st1: 3320476, 
                end1: 3320477, 
                forward1: StrandDirection::Rev, 
                chr2: String::from("000232F"), 
                fp2: 3321632, 
                st2: 3321632, 
                end2: 3321633, 
                forward2: StrandDirection::Rev, 
                sv_type: SVType::BndPair, 
                name: None 
            })]);
        let local_max = 3320477 +1;
        let entries : Vec<Entry<'_, u64, BNDentry>> = tree.find(0..local_max).collect(); 
        let test :  Vec<std::ops::Range<u64>> = locus_dist_based_groups(entries,0);
        let results = vec![3320474..3320475,3320476..3320477];
        assert_eq!(test,results);
    }

    #[test]
    fn vcf2tree_test_breakpoints() {
        // this test is not resulting in any clustered
        // events, everything becomes a single cluster
        // Additionally we have only 1 file and not multiple vcfs
        let tree = get_tree_stranded();
        let contigs  = parse_chrom_file("test/test_cluster/chroms.txt");
        //lets test now the breakpoints accordingly
        let bndpair_break_points : FxHashMap<String,Vec<std::ops::Range<u64>>> = chr_dist_based_groups(&tree.bnd_pair,1,&contigs);
        let mut test2 : FxHashMap<String,Vec<std::ops::Range<u64>>> = FxHashMap::default();
        test2.insert(String::from("13"), vec![123456_u64..123457_u64,123465_u64..123466_u64,123472_u64..123473_u64]);
        test2.insert(String::from("2"),  vec![321680_u64..321681_u64,321691_u64..321692_u64,321698_u64..321699_u64]);
        test2.insert(String::from("17"), vec![198981_u64..198982_u64,198992_u64..198993_u64]);
        assert_eq!(bndpair_break_points,test2);
    }

    #[test]
    fn cluster2tree_1() {
        // tests vector of ranges to interval tree
        let tree     = get_tree_stranded();
        let contigs  = parse_chrom_file("test/test_cluster/chroms.txt");
        let bndpair_break_points : FxHashMap<String,Vec<std::ops::Range<u64>>> = chr_dist_based_groups(&tree.bnd_pair,1,&contigs);
        let mut bnd_break_tree   : FxHashMap<String,IntervalTree<u64,String>>  = FxHashMap::default();
        for (key,value) in bndpair_break_points {
            let local_tree = cluster_to_tree(value);
            bnd_break_tree.insert(key,local_tree);
        };
    }
    #[test]
    fn vcf2tree_test_clusterpair1() {
        // here the distance is too small and no clusters are generated
        let tree = get_tree_stranded();
        let contigs  = parse_chrom_file("test/test_cluster/chroms.txt");
        let bndpair_break_points : FxHashMap<String,Vec<std::ops::Range<u64>>> = chr_dist_based_groups(&tree.bnd_pair,1,&contigs);

        let mut bnd_break_tree   : FxHashMap<String,IntervalTree<u64,String>>  = FxHashMap::default();
        for (key,value) in bndpair_break_points {
            let local_tree = cluster_to_tree(value);
            bnd_break_tree.insert(key,local_tree);
        };
        // the next step now defines pairs of grouped events
        let bnd_cluster_pairs   = bnd_pairs_clustered(&tree.bnd_pair,&bnd_break_tree,&contigs);
        let mut test3: FxHashMap<ClusterPairs,Vec<BNDentry>> = FxHashMap::default();
        let cl_pair1 = ClusterPairs{
            chr1: String::from("2"),
            range1: 321691..321692,
            forward1: StrandDirection::Rev, 
            chr2: String::from("13"),
            range2: 123465..123466, 
            forward2: StrandDirection::Fwd
        };
        let info1 = vec![
            BNDentry{
                sample: String::from("cluster1"), 
                id: String::from("bnd2.1"), 
                mid:  String::from("bnd2.2"), 
                chr1:  String::from("2"), 
                fp1: 321691, 
                st1: 321691, 
                end1: 321692, 
                forward1: StrandDirection::Rev, 
                chr2:  String::from("13"), 
                fp2: 123465, 
                st2: 123465, 
                end2: 123466, 
                forward2: StrandDirection::Fwd, 
                sv_type: SVType::BndPair, 
                name: None 
            }];
        test3.insert(cl_pair1.clone(), info1);
        debug!("{:?}",&bnd_cluster_pairs);
        assert_eq!(bnd_cluster_pairs.get(&cl_pair1),test3.get(&cl_pair1));
    }

    #[test]
    fn vcf2tree_test_clusterpair2() {
        // here we increase distance in which we search
        // and thereby create clusters
        let tree = get_tree_stranded();
        let contigs  = parse_chrom_file("test/test_cluster/chroms.txt");
        let bndpair_break_points : FxHashMap<String,Vec<std::ops::Range<u64>>> = chr_dist_based_groups(&tree.bnd_pair,8,&contigs);

        let mut bnd_break_tree   : FxHashMap<String,IntervalTree<u64,String>>  = FxHashMap::default();
        for (key,value) in bndpair_break_points {
            let local_tree = cluster_to_tree(value);
            bnd_break_tree.insert(key,local_tree);
        };
        eprintln!("Tree: {:?}",bnd_break_tree);
        // the next step now defines pairs of grouped events
        let bnd_cluster_pairs   = bnd_pairs_clustered(
            &tree.bnd_pair,
            &bnd_break_tree,
            &contigs
        );
        let mut test3: FxHashMap<ClusterPairs,Vec<BNDentry>> = FxHashMap::default();
        let cl_pair1 = ClusterPairs{
            chr1: String::from("2"),
            range1: 321691..321699,
            forward1: StrandDirection::Rev, 
            chr2: String::from("13"),
            range2: 123456..123473, 
            forward2: StrandDirection::Fwd
        };
        let info1 = vec![
            BNDentry{
                sample: String::from("cluster1"), 
                id: String::from("bnd2.1"), 
                mid:  String::from("bnd2.2"), 
                chr1:  String::from("2"), 
                fp1: 321691, 
                st1: 321691, 
                end1: 321692, 
                forward1: StrandDirection::Rev, 
                chr2:  String::from("13"), 
                fp2: 123465, 
                st2: 123465, 
                end2: 123466, 
                forward2: StrandDirection::Fwd, 
                sv_type: SVType::BndPair, 
                name: None 
            }];
        test3.insert(cl_pair1.clone(), info1);
        eprintln!("{:?}",&bnd_cluster_pairs);
        assert_eq!(bnd_cluster_pairs.get(&cl_pair1),test3.get(&cl_pair1));
    }

    #[test]
    fn collapse_bndpairs_1() {
        //collapsing bnd pair entries
        let mut test3: FxHashMap<ClusterPairs,Vec<BNDentry>> = FxHashMap::default();
        let mut result3: FxHashMap<ClusterPairs,Vec<BNDentry>> = FxHashMap::default();
        let cl_pair1 = ClusterPairs{
            chr1: String::from("000232F"),
            range1: 3321631..3321633,
            forward1: StrandDirection::Rev, 
            chr2: String::from("000232F"),
            range2: 3320476..3320479, 
            forward2: StrandDirection::Rev
        };
        let info1 = vec![
            BNDentry{
                sample: String::from("test2"), 
                id: String::from("NGSAI_GENEFUSION_PID.1629.1"), 
                mid:  String::from("NGSAI_GENEFUSION_PID.1629.2"), 
                chr1:  String::from("000232F"), 
                fp1: 3321631, 
                st1: 3321631, 
                end1: 3321632, 
                forward1: StrandDirection::Rev, 
                chr2:  String::from("000232F"), 
                fp2: 3320477, 
                st2: 3320477, 
                end2: 3320478, 
                forward2: StrandDirection::Rev, 
                sv_type: SVType::BndPair, 
                name: None 
            }];
        let cl_pair2 = ClusterPairs{
            chr1: String::from("000232F"),
            range1: 3320476..3320479,
            forward1: StrandDirection::Rev, 
            chr2: String::from("000232F"),
            range2: 3321631..3321633, 
            forward2: StrandDirection::Rev
        };
        let info2 = vec![
            BNDentry{
                sample: String::from("test2"), 
                id: String::from("NGSAI_GENEFUSION_PID.1629.2"), 
                mid:  String::from("NGSAI_GENEFUSION_PID.1629.1"), 
                chr2:  String::from("000232F"), 
                fp2: 3321631, 
                st2: 3321631, 
                end2: 3321632, 
                forward2: StrandDirection::Rev, 
                chr1:  String::from("000232F"), 
                fp1: 3320477, 
                st1: 3320477, 
                end1: 3320478, 
                forward1: StrandDirection::Rev, 
                sv_type: SVType::BndPair, 
                name: None 
            }];    

        let res_pair1 = ClusterPairs{
            chr1: String::from("000232F"),
            range1: 3321631..3321633,
            forward1: StrandDirection::Rev, 
            chr2: String::from("000232F"),
            range2: 3320476..3320479, 
            forward2: StrandDirection::Rev
        };
        let res1 = vec![
            BNDentry{
                sample: String::from("test2"), 
                id: String::from("NGSAI_GENEFUSION_PID.1629.1"), 
                mid:  String::from("NGSAI_GENEFUSION_PID.1629.2"), 
                chr1:  String::from("000232F"), 
                fp1: 3321631, 
                st1: 3321631, 
                end1: 3321632, 
                forward1: StrandDirection::Rev, 
                chr2:  String::from("000232F"), 
                fp2: 3320477, 
                st2: 3320477, 
                end2: 3320478, 
                forward2: StrandDirection::Rev, 
                sv_type: SVType::BndPair, 
                name: None 
            }];
        test3.insert(cl_pair1.clone(), info1);
        test3.insert(cl_pair2.clone(), info2);
        result3.insert(res_pair1.clone(), res1);
        let result = bnd_collapse_pairs(test3);
        assert_eq!(result,result3);
    }

    
    #[test]
    fn vcf2tree_test_refinement1() {
        // now we test the refinement function 
        // which if necessary breaks down ranges
        let tree = get_tree_stranded();
        let contigs  = parse_chrom_file("test/test_cluster/chroms.txt");
        let bndpair_break_points : FxHashMap<String,Vec<std::ops::Range<u64>>> = chr_dist_based_groups(&tree.bnd_pair,8,&contigs);

        let mut bnd_break_tree   : FxHashMap<String,IntervalTree<u64,String>>  = FxHashMap::default();
        for (key,value) in bndpair_break_points {
            let local_tree = cluster_to_tree(value);
            bnd_break_tree.insert(key,local_tree);
        };
        // the next step now defines pairs of grouped events
        let bnd_cluster_pairs   = bnd_pairs_clustered(&tree.bnd_pair,&bnd_break_tree,&contigs);
        let bnd_refined_pairs   = bnd_filter_pairs(bnd_cluster_pairs,1);
        let mut test3: FxHashMap<ClusterPairs,Vec<BNDentry>> = FxHashMap::default();
        let cl_pair1 = ClusterPairs{
            chr1: String::from("2"),
            range1: 321691..321693,
            forward1: StrandDirection::Rev, 
            chr2: String::from("13"),
            range2: 123465..123467, 
            forward2: StrandDirection::Fwd
        };
        let info1 = vec![
            BNDentry{
                sample: String::from("cluster1"), 
                id: String::from("bnd2.1"), 
                mid:  String::from("bnd2.2"), 
                chr1:  String::from("2"), 
                fp1: 321691, 
                st1: 321691, 
                end1: 321692, 
                forward1: StrandDirection::Rev, 
                chr2:  String::from("13"), 
                fp2: 123465, 
                st2: 123465, 
                end2: 123466, 
                forward2: StrandDirection::Fwd, 
                sv_type: SVType::BndPair, 
                name: None 
            }];
        let cl_pair2 = ClusterPairs{
            chr2: String::from("2"),
            range2: 321691..321693,
            forward2: StrandDirection::Rev, 
            chr1: String::from("13"),
            range1: 123465..123467, 
            forward1: StrandDirection::Fwd
        };
        let info2 = vec![
            BNDentry{
                sample: String::from("cluster1"), 
                id: String::from("bnd2.2"), 
                mid:  String::from("bnd2.1"), 
                chr2:  String::from("2"), 
                fp2: 321691, 
                st2: 321691, 
                end2: 321692, 
                forward2: StrandDirection::Rev, 
                chr1:  String::from("13"), 
                fp1: 123465, 
                st1: 123465, 
                end1: 123466, 
                forward1: StrandDirection::Fwd, 
                sv_type: SVType::BndPair, 
                name: None 
            }];    

        test3.insert(cl_pair1.clone(), info1);
        test3.insert(cl_pair2.clone(), info2);
        eprintln!("{:?}",&bnd_refined_pairs);
        assert_eq!(bnd_refined_pairs.get(&cl_pair1),test3.get(&cl_pair1));
        assert_eq!(bnd_refined_pairs.get(&cl_pair2),test3.get(&cl_pair2));
    }

    #[test]
    fn vcf2tree_test_refinement2() {
        // here we increase distance and do unstranded
        // and thereby create clusters at the later steps, too
        let tree = get_tree_unstranded();
        let contigs  = parse_chrom_file("test/test_cluster/chroms.txt");
        let bndpair_break_points : FxHashMap<String,Vec<std::ops::Range<u64>>> = chr_dist_based_groups(&tree.bnd_pair,8,&contigs);

        let mut bnd_break_tree   : FxHashMap<String,IntervalTree<u64,String>>  = FxHashMap::default();
        for (key,value) in bndpair_break_points {
            let local_tree = cluster_to_tree(value);
            bnd_break_tree.insert(key,local_tree);
        };
        // the next step now defines pairs of grouped events
        let bnd_cluster_pairs   = bnd_pairs_clustered(&tree.bnd_pair,&bnd_break_tree,&contigs);
        let bnd_refined_pairs   = bnd_filter_pairs(bnd_cluster_pairs,10);
        let mut test3: FxHashMap<ClusterPairs,Vec<BNDentry>> = FxHashMap::default();
        let cl_pair1 = ClusterPairs{
            chr1: String::from("2"),
            range1: 321691..321700,
            forward1: StrandDirection::Unknown, 
            chr2: String::from("13"),
            range2: 123465..123474, 
            forward2: StrandDirection::Unknown
        };
        let info1 = vec![
            BNDentry{
                sample: String::from("cluster1"), 
                id: String::from("bnd2.1"), 
                mid:  String::from("bnd2.2"), 
                chr1:  String::from("2"), 
                fp1: 321691, 
                st1: 321691, 
                end1: 321692, 
                forward1: StrandDirection::Unknown, 
                chr2:  String::from("13"), 
                fp2: 123465, 
                st2: 123465, 
                end2: 123466, 
                forward2: StrandDirection::Unknown, 
                sv_type: SVType::BndPair, 
                name: None 
            },
            BNDentry{ 
                sample: String::from("cluster1"), 
                id: String::from("bnd4.1"), 
                mid: String::from("bnd4.2"), 
                chr1: String::from("2"), 
                fp1: 321698, 
                st1: 321698, 
                end1: 321699, 
                forward1: StrandDirection::Unknown, 
                chr2: String::from("13"), 
                fp2: 123472, 
                st2: 123472, 
                end2: 123473, 
                forward2: StrandDirection::Unknown, 
                sv_type: SVType::BndPair, 
                name: None 
            }
            ];

        let cl_pair2 = ClusterPairs{
            chr2: String::from("2"),
            range2: 321691..321700,
            forward2: StrandDirection::Unknown, 
            chr1: String::from("13"),
            range1: 123465..123474, 
            forward1: StrandDirection::Unknown
        };
        let info2 = vec![
            BNDentry{
                sample: String::from("cluster1"), 
                id: String::from("bnd2.2"), 
                mid:  String::from("bnd2.1"), 
                chr2:  String::from("2"), 
                fp2: 321691, 
                st2: 321691, 
                end2: 321692, 
                forward2: StrandDirection::Unknown, 
                chr1:  String::from("13"), 
                fp1: 123465, 
                st1: 123465, 
                end1: 123466, 
                forward1: StrandDirection::Unknown, 
                sv_type: SVType::BndPair, 
                name: None 
            },
            BNDentry{ 
                sample: String::from("cluster1"), 
                id: String::from("bnd4.2"), 
                mid: String::from("bnd4.1"), 
                chr2: String::from("2"), 
                fp2: 321698, 
                st2: 321698, 
                end2: 321699, 
                forward2: StrandDirection::Unknown, 
                chr1: String::from("13"), 
                fp1: 123472, 
                st1: 123472, 
                end1: 123473, 
                forward1: StrandDirection::Unknown, 
                sv_type: SVType::BndPair, 
                name: None 
            }
            ];
    
        test3.insert(cl_pair1.clone(), info1);
        test3.insert(cl_pair2.clone(), info2);
        debug!("{:?}",&bnd_refined_pairs);
        assert_eq!(bnd_refined_pairs.get(&cl_pair1),test3.get(&cl_pair1));
        assert_eq!(bnd_refined_pairs.get(&cl_pair2),test3.get(&cl_pair2));
    }

}