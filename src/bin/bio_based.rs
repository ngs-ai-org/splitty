
//! This is a simple test program to better understand the tree structure and interval trees in general.
//! It does not serve any purpose beyond that and cant be used to take any input data.
//! It will essentially replicate this structure:
//!              CHR1
//!         
//! ```        
//! 900--------S2:2.1-----------1200->
//! 1000-----S3:3.1------1100->
//! 1000-----S1:1.1------------------1300->
//! 1000------S2:9.1----------1250->
//!                                              1400---S4:5.1-----1600->
//!                                                 1500--S5:6.1----1600->
//!                                                 1500----S6:7.1----------1700->
//!                                                                                         5000-S7:-8.1--5500->
//!                                                                                        <-4800---S3:10.1--------6000
//! 
//! ----------------------------------------------------------------------------------------------------------------------------------
//! 
//!   CHR7
//! 
//! 
//!                                                                     10050-------------S1:1.2-------------------------------12500->
//!                                                                     10000----S2:2.2-----10100->
//!                                                                     10000------S3:3.2-------------10400->
//!                                                                     10000------S3:10.2------------10400->
//! 100----S4:5.2-----300->
//! 100----S5:6.2-----300->                   
//! 100----S6:7.2-----300->
//! 100----S7:8.2-----300->
//!       120-----------S2:9.2-------500->    
//! ```


use rustc_hash::FxHashMap;
use bio::data_structures::interval_tree::{*};
use std::iter::FromIterator;

extern crate genefusion;
use genefusion::lib::common::{*};
//use delta_biomed_bio::hts_lib_based::{*};
 

fn main() {
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
        
                ]
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
            ]
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
        // a above mentioned the problem is essentially now that the cluster from 100--300 was extended till 500 because
        // the event 9 overlaps. In the later step we identified though that these are not having the same matching pair region
        // therefore we need now in a second step now a way to redefine the cluster for a given pair in the
        // same fashion as initially
        // So we take each pair
        let refined_pairs = bnd_filter_pairs(cluster_pairs,range);

        // at this point we finished operations and we collapse now the events for easier processing
        let collapsed_pairs = bnd_collapse_pairs(refined_pairs);
         
        
        let multi = false;
        let purge = false;
        println!("{}\n\n\n\n",schematic);
        println!("# elements from BND pairs clustered into logical groups and refined by proper pairs");
        println!("#CHR1\tSTART1\tEND1\tCHR2\tSTART2\tEND2\tELEMENTS");
        for (key,value) in collapsed_pairs {
            // this will squash entries if there are >1 per sample
            // might be a bit debatable in a biological sense but
            // might be absolutely necessary for a proper output format eventually
            let samples : Vec<String> = match purge {
                false => {
                    value.into_iter().map(|x| format!("{}:{}",x.sample,x.id)).collect()
                }
                true  => {
                    let mut tmp = FxHashMap::default();
                    for sa in value {
                        tmp.insert(sa.sample,sa.id);
                    }
                    tmp.into_iter().map(|(key,value)| format!("{}:{}",key,value)).collect()
                }
            };
            match multi {
                false => {
                    print!("{}\t{}\t{}\t{}\t{}\t{}",key.chr1,key.range1.start,key.range1.end,key.chr2,key.range2.start,key.range2.end);
                    for n in samples {
                        print!("\t{}",n);
                    }
                        println!();
                    },
                    true => {
                        if samples.len() > 1 {
                            print!("{}\t{}\t{}\t{}\t{}\t{}",key.chr1,key.range1.start,key.range1.end,key.chr2,key.range2.start,key.range2.end);
                            for n in samples {
                                print!("\t{}",n);
                            }
                            println!();
                        }
                    }
            }
        } 

}

