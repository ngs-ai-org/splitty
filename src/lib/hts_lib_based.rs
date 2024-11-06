use rustc_hash::FxHashMap;
use std::str;
use rust_htslib::{bam, bam::Read};

use std::error::Error;
use linear_map::LinearMap;
use std::collections::HashMap;
use rust_htslib::bcf::header::Header;
use rust_htslib::bcf::Read as BcfRead;
use rust_htslib::bcf::{Format, Writer};
use bio::io::fasta::IndexedReader;
use std::str::from_utf8;
use bio::data_structures::interval_tree::{*};
use std::ops::Deref;
use itertools::Itertools;  // itertools = "0.10"
use std::path::Path;
use std::convert::TryInto;
use std::convert::TryFrom;
use log::debug;
//use common::{*};
use crate::lib::common::{*};


/// returns a median value for a vector of
/// unsigned 32-bit numbers
/// 
/// Unittest: TRUE
///
/// ```
/// use crate::genefusion::lib::hts_lib_based::median_u32;
/// let mut test = vec![5_u32,6_u32,10_u32,20_u32,3_u32];
/// assert_eq!(median_u32(&mut test),6);
/// ```
pub fn median_u32(
    numbers: &mut std::vec::Vec<u32>
) -> u32 {
    numbers.sort_unstable();
    let mid = numbers.len() / 2;
    numbers[mid]
}

/// just a convenience function to easier deal 
/// with immutable and mutable borrowing
/// Returns the HeaderView
///
///Unittest: FALSE
///
pub fn get_header(
    file:&rust_htslib::bcf::Reader
) -> rust_htslib::bcf::header::HeaderView {
    let header = file.header();
    header.clone()
}

/// this checks if our alignments are indeed what we 
/// expect and removes entries that are secondary alignment 
/// or paired or unmapped 
///
///Unittest: TRUE
///
pub fn read_qc_module(
    entry: &bam::Record,
    black_list: &Option<FxHashMap<String, i32>>,
    accept_multi: &bool
) -> bool {
    let name = &str::from_utf8(entry.qname()).unwrap().to_string();
    let on_black_list = match black_list {
        Some(x) => {
            x.contains_key(name)
        },
        None => false
    };
    if !accept_multi{
        !(entry.is_secondary() || entry.mapq()==0 || entry.is_paired() || entry.is_unmapped() || on_black_list)
    }else{
        !(entry.is_secondary() || entry.is_paired() || entry.is_unmapped() || on_black_list)
    }
}

/// this checks if our alignments are indeed what we 
/// expect and removes return false for entries which
/// are not according to specs
/// 
///Unittest: TRUE
///
pub fn read_qc_module_paired(
    entry: &bam::Record,
    black_list: &Option<FxHashMap<String, i32>>,
    accept_multi: &bool
) -> bool {
    let name = &str::from_utf8(entry.qname()).unwrap().to_string();
    if !entry.is_paired() {
        panic!("ERROR: Alignments are not paired-reads!");
    }
    let on_black_list = match black_list {
        Some(x) => {
            x.contains_key(name)
        },
        None => false
    };
    debug!("Is secondary: {:?}",entry);
    if !accept_multi {
        !(entry.is_secondary() || entry.mapq()==0 || entry.is_unmapped() || on_black_list)
    }else{
        !(entry.is_secondary() || entry.is_unmapped() || on_black_list)
    }
}



/// This uses a BAM indexed reader to analyze the primary and SA alignments 
/// It verifies that qc is good, gets the cigar and evaluates it
/// Calls the function AlignInfos and eval_supp_align to extract the
/// information from the reads and organizes it into a vector of 
/// full fusion evidence - it contains essentially the information organized by read
/// All coordinates are 0-based
///
/// Unittest = TRUE
///
pub fn analyze_alignments_from_bam( 
    mut bam: bam::Reader, // the bam file with the alignments
    black: Option<FxHashMap<String, i32>>, // a potential list of IDs which should be ignored
    similarity: f32, // similarity of read mapped to ref genome
    paired: bool, // if paired or unpaired reads
    multi: bool, // accept multi mapping reads with MAPQ=0
    // maximum allowed ratio of short to longer clipping
    // in case that both sides are clipped
    // If none, the larger one is used to define the breakpoint
    ratio: Option<f32>
) -> (FxHashMap<String, Vec<FullFusionEvidence>> , FxHashMap<String,u8> ){
    
    let mut chromosomes: FxHashMap<String,u8> = FxHashMap::default();
    let header      = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    let mut record  = bam::Record::new();

    // Read --> chrA --> chrB --> Vect
    let mut read_based_organization : FxHashMap<String, Vec<FullFusionEvidence>> = FxHashMap::default();

    while let Some(result) = bam.read(&mut record) {
        match result {
            Ok(_) => {
                
                // here comes now the QC part
                if paired {
                    match read_qc_module_paired(&record, &black, &multi) {
                        true => (),
                        false => continue,
                    }
                }else{
                    match read_qc_module(&record, &black, &multi) {
                        true => (),
                        false => continue,
                    }
                }
                // if entry has no supplementary alignment  "SA" then we
                // can directly continue
                let sa_entry = match record.aux(b"SA") {
                    Ok(rust_htslib::bam::record::Aux::String(sa_s)) => sa_s,
                    _ => continue,
                };
                // now as this is boolean we can test directly
                // this returns now the structure CigarInfo
                let cigar_format = format!("{}", record.cigar());

                
                // one is via the NM tag and the other one is rather new and in the
                // cigar string decoded in "=" and "X" for matches and mismatches
                let mismatch     = match record.aux(b"NM") {
                    Ok(rust_htslib::bam::record::Aux::I32(nm_i)) => nm_i as u32 ,
                    _ => 0
                };

                let match_start: u64 = record.pos().try_into().expect("ERROR: could not convert genomic position to u64!");
                // I really dislike 0/1 and prefer rather 1/-1
                let match_strand = match record.is_reverse() {
                    true  => StrandDirection::Rev,
                    false => StrandDirection::Fwd,
                };
                // get the target ID from the header information
                // as we get here only the entry number but not it's name
                let tid = str::from_utf8(header_view.tid2name(record.tid() as u32)).unwrap();
                // now the primary entry is always more precise and reliable
                // then the supplementary one "SA" . We get accordingly
                // for both the information where they predict a fusion position
                let cigar_infos = eval_cigar(
                    &cigar_format, 
                    mismatch.into(), 
                    &match_start, 
                    &match_strand,
                    ratio
                );
                debug!("Cigar information: {:?}",&cigar_infos);
                let success_infos = match cigar_infos {
                    None    => continue,
                    Some(x) => x
                };
                
                if success_infos.similarity <  similarity {
                    continue;
                }

                let prim_infos = AlignInfos {
                    query:          str::from_utf8(record.qname()).unwrap().to_string(),
                    q_length:       success_infos.seq_length as u64, 
                    // now this is strange I wanted to use instead here
                    // record.seq_len() as u64,
                    // but it seems to report sequence length without clipping ??? 
                    cigar:          cigar_format,
                    chrom:          tid.to_string(),
                    fp_target:      success_infos.fp_target,
                    upstream:       success_infos.upstream,
                    downstream:     success_infos.downstream,
                    similarity:     success_infos.similarity,
                    seq_length:     success_infos.seq_length,
                    match_length_s: success_infos.match_length_s,
                    match_length:   success_infos.match_length,
                    fp_query:       success_infos.fp_query,
                    aln_start:      success_infos.aln_start,
                    aln_end:        success_infos.aln_end,
                    match_strand:   success_infos.match_strand,
                };
                debug!("Primary info: {:?}",&prim_infos);

                // here we get all SA information
                let supp_infos = eval_supp_align(
                    sa_entry,
                    str::from_utf8(record.qname()).unwrap(),
                    Some(&prim_infos.seq_length.try_into().unwrap()),
                    ratio
                );
                debug!("Getting back the supp_infos {:?}",&supp_infos);
                // so previously I had here another check for similarity
                // but currently I believe that it might not be a good idea
                // because SA derived CIGAR are just really not that precise and
                //supp_infos.retain(|x| x.similarity > identity);
                // now we check again if SA entries remaining, otherwise we skip
                if supp_infos.is_empty() {
                    continue;
                }
                chromosomes.insert(prim_infos.chrom.clone(), 1);
                for check in supp_infos.iter(){
                    chromosomes.insert(check.chrom.clone(),1);
                }
                // at this point we have now all information
                // which are necessary:
                // primary alignment information and already the supplementary ones
                // The latter are a vector and we need to consider now each case
                // individually for the next steps
                let result: Vec<FullFusionEvidence> = evaluate_and_predict(prim_infos, supp_infos);
                debug!("Gathered for a read all info {:?}",&result);
                for combi in result.into_iter() {
                    debug!("Iterating over results, evaluating now {:?}",&combi);
                    // here we have to clone as we move otherwise part of combi
                    read_based_organization
                        .entry(combi.query_name.clone())
                        .or_default()
                        .push(combi);
                }
            }
            Err(_) => panic!("BAM parsing failed..."),
        }
    }
    (read_based_organization, chromosomes)
}


/// This function takes a position based
/// fusion result and re-analyzes it.
/// Based on encountered split-pairs (SP) 
/// it increases for a riven entry it's support
/// which allows to judge better the validity of a
/// discovered fusion. d
///
/// Unittests: TRUE
///
pub fn add_sp_to_sr_fusions(
    bam: & mut bam::IndexedReader, 
    pos_based_results: FxHashMap<String,Vec<PosBasedFusEvidence>>, 
    range: u32
) -> FxHashMap<String,Vec<PosBasedFusEvidence>> {

    let range_sp = u64::from(range);
    // ideally I would have liked to modify the existing pos_based_results
    // via static borrowing but if we iterate over it rust does not
    // allow a modification in place easily
    // As the input/output is not anymore as complex 
    // at this stage we therefore populate a new
    // Object instead
    let mut added_sp_result : FxHashMap<String,Vec<PosBasedFusEvidence>>  = FxHashMap::default();
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);

    // we do have an organization where
    // the key is a fusion string of chromosome a and chromosome b 
    for (chr_ab,results) in pos_based_results{
        for fusion in results {
            debug!("Analyzing chr_ab : {:?}  with fusion {:?}",&chr_ab,&fusion);
            let mut fusion_copy = fusion.to_owned();
            // now we need to evaluate for PR and SA
            // information if they have reads which
            // reside in the ROI
            let start_a : u64 ;
            let end_a   : u64 ;
            let start_b : u64 ;
            let end_b   : u64 ;
            let chr_a   : &str = &fusion_copy.chrom_pr ;
            let chr_b   : &str = &fusion_copy.chrom_sa;
            // this depends now on orientation
            if fusion_copy.orient_pr_upstream == true {
                if fusion_copy.strand_pr == StrandDirection::Fwd {
                    start_a = if range_sp > fusion_copy.target_break_pr { 0_u64 }else{ fusion_copy.target_break_pr - range_sp };
                    end_a   = fusion_copy.target_break_pr;
                }else if fusion_copy.strand_pr == StrandDirection::Rev {
                    start_a = fusion_copy.target_break_pr ;
                    end_a   = fusion_copy.target_break_pr + range_sp ;
                }else {
                    panic!("ERROR: could not determine PR strand in SP function!");
                }
            }else if fusion_copy.orient_pr_upstream == false {
                if fusion_copy.strand_pr == StrandDirection::Fwd {
                    start_a = fusion_copy.target_break_pr;
                    end_a   = fusion_copy.target_break_pr + range_sp ;
                }else if fusion_copy.strand_pr == StrandDirection::Rev {
                    start_a = if range_sp > fusion_copy.target_break_pr { 0_u64 }else{ fusion_copy.target_break_pr  - range_sp };
                    end_a   = fusion_copy.target_break_pr;
                }else {
                    panic!("ERROR: could not determine PR strand in SP function!");
                }
            }else{
                panic!("ERROR: could not determine PR orientation in SP function!");
            }
            // this depends now on orientation
            if fusion_copy.orient_sa_upstream == true {
                if fusion_copy.strand_sa == StrandDirection::Fwd {
                    start_b = if range_sp > fusion_copy.target_break_sa { 0_u64 }else{ fusion_copy.target_break_sa - range_sp };
                    end_b   = fusion_copy.target_break_sa;
                }else if fusion_copy.strand_sa == StrandDirection::Rev {
                    start_b = fusion_copy.target_break_sa ;
                    end_b   = fusion_copy.target_break_sa + range_sp ;
                }else {
                    panic!("ERROR: could not determine SA strand in SP function!");
                }
            }else if fusion_copy.orient_sa_upstream == false {
                if fusion_copy.strand_sa == StrandDirection::Fwd {
                    start_b = fusion_copy.target_break_sa;
                    end_b   = fusion_copy.target_break_sa + range_sp ;
                }else if fusion_copy.strand_sa == StrandDirection::Rev {
                    start_b = if range_sp > fusion_copy.target_break_sa { 0_u64 }else{ fusion_copy.target_break_sa  - range_sp };
                    end_b   = fusion_copy.target_break_sa;
                }else {
                    panic!("ERROR: could not determine SA strand in SP function!");
                }
            }else{
                //eprintln!("");
                panic!("ERROR: could not determine SA orientation in SP function!");
            }
            // now we have our ranges and we can get
            // the reads at that given position
            // It should be sufficient to get everything on
            // one side and then check if the mate-position indicates
            // a match for the second side
            debug!("Get a pileup of reads for chr {:?} start: {:?} end: {:?}",&chr_a,&start_a,&end_a);
            bam.fetch((chr_a,start_a,end_a)).expect("ERROR: could not fetch region in SP analysis!");
           
            // now to avoid that they are on the same
            // chromosome but not at all associated with the fusion
            // we need to start with our ranges to limit the region of interest
            // Now we defined the start/end purely by coordinates
            // on the chromosome where start is < than end
            // Currently we define an evidence if the read alignment is overlapping
            // our defined region which is for A already by fetching the case
            // With mpos we get the alignment start of the mate
            // this might be though pretty far away if the read is splice-aligned
            // and we cant take that position unfortunately
            // Instead what we do is get all IDs --> MIDs from position A
            // and as well all MIDs-->IDs from the second one
            let mut a_b_relation : FxHashMap<String,bool> = FxHashMap::default();
            let mut b_a_relation : FxHashMap<String,bool> = FxHashMap::default();
            let mut a_reads : FxHashMap<String,bool> = FxHashMap::default();
            let mut b_reads : FxHashMap<String,bool> = FxHashMap::default();
            // lets get the coverage in the adjacent region
            let mut depth_a = vec![];
            for pile in bam.pileup().map(|x| x.expect("ERROR: could not parse BAM file")){
                if (pile.pos() as u64 >= start_a) & ((pile.pos() as u64 )< end_a ) {
                    depth_a.push(pile.depth());
                }
            };

            debug!("Got depth {:?}",&depth_a);
            let d_a = match depth_a.len(){
                0 => 0,
                _ => median_u32(&mut depth_a),
            };
            debug!("Median depth {:?}",&d_a);

            debug!("Get a pileup of reads for chr {:?} start: {:?} end: {:?}",&chr_a,&start_a,&end_a);
            bam.fetch((chr_a,start_a,end_a)).expect("ERROR: could not fetch region in SP analysis!");
            for element in bam.records(){
                let read = match element{
                    Ok(x) => x,
                    _ => panic!("ERROR: could not access read in SP analysis!")
                };
                let name = str::from_utf8(read.qname()).unwrap().to_string();
                a_reads.insert(name.clone(), true);
                
                // this test is not only to make it faster as is skips directly
                // if not a pair mapped
                // Turns out as well that we can get a segfault if 
                // the mate is unmapped and we try to get the mtid()
                if read.is_mate_unmapped() {
                    continue
                }
                // now we need again the mate-target ID as a string
                // rather than a number
                let mtid = str::from_utf8(header_view.tid2name(read.mtid() as u32)).expect("ERROR: could not get mtarget ID in SP analysis!");
                // if the mate is not matching onto the b-target we skip directly
                if mtid != chr_b {
                    continue
                }else{
                    a_b_relation.insert(name.clone(), true);
                }
            }
            debug!("Found the following potentials reads for SP : {:?} with total of reads: {:?}",a_b_relation,a_b_relation.len());
            debug!("Get a pileup of reads for chr {:?} start: {:?} end: {:?}",&chr_b,&start_b,&end_b);
            // now we do the same again but with the b-region:
            bam.fetch((chr_b,start_b,end_b)).expect("ERROR: could not fetch region in SP analysis!");
            let mut depth_b = vec![];
            for pile in bam.pileup().map(|x| x.expect("ERROR: could not parse BAM file")){
                if (pile.pos() as u64 >= start_b) & ((pile.pos() as u64 )< end_b ) {
                    depth_b.push(pile.depth());
                }
            };
            let d_b = match depth_b.len(){
                0 => 0,
                _ => median_u32(&mut depth_b),
            };
            debug!("Median depth {:?}",&d_b);

            bam.fetch((chr_b,start_b,end_b)).expect("ERROR: could not fetch region in SP analysis!");
            for element in bam.records(){
                let read = match element{
                    Ok(x) => x,
                    _ => panic!("ERROR: could not access read in SP analysis!")
                };
                let name = str::from_utf8(read.qname()).unwrap().to_string();
                b_reads.insert(name.clone(), true);
                if read.is_mate_unmapped() {
                    continue
                }
                // now we need again the mate-target ID as a string
                // rather than a number
                let mtid = str::from_utf8(header_view.tid2name(read.mtid() as u32)).unwrap();
                // if the mate is not matching onto the b-target we skip directly
                if mtid != chr_a {
                    continue
                }else{
                    b_a_relation.insert(name.clone(), true);
                }
            }
            debug!("Found the following potentials reads for SP : {:?} with total of reads: {:?}",b_a_relation,b_a_relation.len());
            // If we have removed already reads which are clearly identified as derived from a 
            // clean transcript we should now find for a given gene-fusion a similar number of
            // reads indicate at GeneA a fusion as well as for GeneB. But sometimes we might encounter
            // regions which are highly problematic and therefore we could multiple suggested fusion 
            // outgoing from the same fusion point. This is predominantly not a good sign.
            // Therefore we take the ratio of reads indicative of our fusion and the ones in that region in
            // total. This might indicate then an imbalance towards one edge. 
            debug!("a: {:?} \n b: {:?}",a_b_relation.keys().sorted().len(),b_a_relation.keys().sorted().len());
            for i in a_b_relation.keys(){
                if b_a_relation.contains_key(i){
                    fusion_copy.support_sp += 1;
                }
            };
            let ratio_a = match a_reads.len(){
                0 => 0 as f32,
                _ => (fusion_copy.support_sp as f32 )/ (d_a as f32) 
            };
            let ratio_b = match b_reads.len(){
                0 => 0 as f32,
                _ => (fusion_copy.support_sp as f32 )/ (d_b as f32)
            };
            debug!("Reads supporting : {}, ratio_a {}, ratio_b {}",&fusion_copy.support_sp,&ratio_a,&ratio_b );
            let mut ratios : Vec<f32> = vec![ratio_a,ratio_b];
            ratios.sort_by(|a,b| a.partial_cmp(b).unwrap());
            let delta = ratios[0];
            fusion_copy.sp_ratio = delta;
            // now we can add our entry to the new table
            // entries were either identical to original
            // or modified in the process if evidence was
            // found 
            if added_sp_result.contains_key(&chr_ab) {
                let current_results = added_sp_result.get_mut(&chr_ab).expect("ERROR: could not catch existing elements in SP analysis!");
                current_results.push(fusion_copy);
            }else{
                added_sp_result.insert(chr_ab.to_owned(), vec![fusion_copy]);
            }
        }
    };
    debug!("Results of added split pairs: {:?}",&added_sp_result);
    added_sp_result
}

/*
#[derive(Debug)]
pub struct AlleleInfo {
    allele1: String,
    allele2: String,
}
*/

/// converts VCF to fasta. It expects a pair of VCF entries
/// together with the information how it should handle it.
/// Then it extracts the DNA sequence accordingly.
///
///Unittest: TRUE
/// 
/// ```rust
/// use genefusion::lib::hts_lib_based::vcfpair2fa_record;
/// use rust_htslib::bcf::{Format, Header, Writer};
/// use bio::io::fasta::IndexedReader;
/// use std::fs::File;
/// use std::io::prelude::*;
/// use std::str::from_utf8;
/// use tempfile::NamedTempFile;
/// pretty_env_logger::init();
/// let tmp = NamedTempFile::new().unwrap();
/// let path = tmp.path();
///
/// 
/// // create dummy files
/// let mut file_seq =  File::create("foo.fa").expect("ERROR:could not create seq file!");
/// file_seq.write_all(b">chr1\nGTAGGCTGAAAA\nCCCGATTGACTA\nchr2\nCCCGATTGACTA\nGTAGGCTGAAAA").expect("ERROR:could not write seq file!");
/// let assembly = String::from("foo.fa");
/// 
/// // index as well
/// let mut file_idx =  File::create("foo.fa.fai").expect("ERROR:could not create idx file!");
/// file_idx.write_all(b"chr1\t24\t6\t12\t13\nchr2\t24\t38\t12\t13").expect("ERROR:could not write idx file!");
/// 
/// 
/// let header = Header::new();
/// let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
/// 
/// // FF record test
/// let mut record1 = vcf.empty_record();
/// let alleles: &[&[u8]] = &[b"C", b"C]chr2:6]"];
/// record1.set_alleles(alleles).expect("Failed to set alleles");
/// record1.set_pos(6);
/// let mut record2 = vcf.empty_record();
/// let alleles: &[&[u8]] = &[b"T", b"T]chr1:6]"];
/// record2.set_alleles(alleles).expect("Failed to set alleles");
/// record2.set_pos(6);
/// let fasta_record = vcfpair2fa_record(&record1,&record2,&assembly,10);
/// let empty_record = bio::io::fasta::Record::new();
/// // delete file
/// std::fs::remove_file("foo.fa").unwrap();
/// // delete file
/// std::fs::remove_file("foo.fa.fai").unwrap();
/// 
/// ```
pub fn vcfpair2fa_record (
    entry1: &rust_htslib::bcf::Record,
    entry2: &rust_htslib::bcf::Record,
    assembly: &str,
    range: u32,
) -> bio::io::fasta::Record {
    // here we get botht the a and b nucleotide
    // for both of the alleles provided
    let alleles1  = entry1.alleles();
    let alleles2  = entry2.alleles();
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

    let allele1_desconstr : Bnd2ndAllele = deconstruct_2nd_allele(allele_2,true);
    let allele2_desconstr : Bnd2ndAllele = deconstruct_2nd_allele(allele_1,true);
    // we have the following pattern now
    //
    // FF = A + B 
    // FR = A + B rev comp
    // RF = rev comp A + B 
    // RR = A comp + B rev comp

    // FR = A start-start+x ; B start-x
    let pos_11 : u64 = allele1_desconstr.pos; 
    let pos_21 : u64 = allele2_desconstr.pos; 

    // pos_11 -= 1;
    // pos_21 -= 1;
    // so now we need the start and end but obviously need to verify as well that 
    // not becoming negative as this will create overflow
    let pos_12 : u64 = match allele1_desconstr.second_direction {
        StrandDirection::Fwd => {
            if pos_11 >= range as u64 {
                pos_11 - range as u64
            }else {
                0_u64
            }
        },
        StrandDirection::Rev =>  pos_11 + range as u64,
        _ => panic!("ERROR: could not determine the sense of vcf entry!"),
    };
    let pos_22 : u64 = match allele2_desconstr.second_direction {
        StrandDirection::Fwd => pos_21 + range as u64,
        StrandDirection::Rev => {
            if pos_21 >= range as u64 {
                pos_21 - range as u64
            }else {
                0_u64
            }
        },
        _ => panic!("ERROR: could not determine the sense of vcf entry!"),
    };
    //dbg!(&allele1_desconstr);
    //dbg!(entry1.pos());
    //dbg!(&allele2_desconstr);
    let id1        = String::from_utf8(entry1.id()).expect("ERROR: could not convert ID to utf8!"); 
    let id2        = String::from_utf8(entry2.id()).expect("ERROR: could not convert ID to utf8!"); 
    let read_id    = format!("{}--{}",id1,id2);
    let read_desc  = format!("{}:{}-{}-{:?};{}:{}-{}-{:?}",
        &allele1_desconstr.chr,
        pos_11,
        pos_12,
        allele1_desconstr.second_direction,
        &allele2_desconstr.chr,
        pos_21,
        pos_22,
        allele1_desconstr.second_direction
    ); 

    // here we define FWD fix as direction as we need afterwards to correct
    // the coordinates depending on the order, rather on single elements
    // secondary direction has to be used!
    
    // second_direction has to be used as it actually uses the opposite allele information
    // above!
    // additionally, depending on the order and direction the coordinates might be off by one
    // as the coordinates start at 0 and are exclusive, so the moment we shift direction we
    // need to correct that accordingly
    let fasta_record = match (
        allele1_desconstr.second_direction,
        allele2_desconstr.second_direction 
    ) {
        (StrandDirection::Fwd,StrandDirection::Fwd) => {
            let result = vec![
                fasta_extract_and_correct(
                    assembly,
                    &allele1_desconstr.chr,
                    pos_11+1,
                    pos_12+1,
                    DnaReturnSense::Fwd,
                ),
                fasta_extract_and_correct(
                    assembly,
                    &allele2_desconstr.chr,
                    pos_21,
                    pos_22,
                    DnaReturnSense::Fwd,
                ).to_lowercase()
            ];
            let mut sequence = String::new();
            sequence        += &result[0];
            sequence        += &result[1];
            bio::io::fasta::Record::with_attrs(
                &read_id,
                Some(&read_desc),
                sequence.as_bytes()
            )
        },
        (StrandDirection::Fwd,StrandDirection::Rev) => {
            let result = vec![
                fasta_extract_and_correct(
                    assembly,
                    &allele1_desconstr.chr,
                    pos_11+1,
                    pos_12+1,
                    DnaReturnSense::Fwd,
                ),
                fasta_extract_and_correct(
                    assembly,
                    &allele2_desconstr.chr,
                    pos_21+1,
                    pos_22+1,
                    DnaReturnSense::RevC,
                ).to_lowercase()
            ];
            let mut sequence = String::new();
            sequence        += &result[0];
            sequence        += &result[1];
            bio::io::fasta::Record::with_attrs(
                &read_id,
                Some(&read_desc),
                sequence.as_bytes()
            )
        },
        (StrandDirection::Rev,StrandDirection::Fwd) => {
            let result = vec![
                fasta_extract_and_correct(
                    assembly,
                    &allele1_desconstr.chr,
                    pos_11,
                    pos_12,
                    DnaReturnSense::RevC,
                ),
                fasta_extract_and_correct(
                    assembly,
                    &allele2_desconstr.chr,
                    pos_21,
                    pos_22,
                    DnaReturnSense::Fwd,
                ).to_lowercase()
            ];
            let mut sequence = String::new();
            sequence        += &result[0];
            sequence        += &result[1];
            bio::io::fasta::Record::with_attrs(
                &read_id,
                Some(&read_desc),
                sequence.as_bytes()
            )
        },
        (StrandDirection::Rev,StrandDirection::Rev) => {
            let result = vec![
                fasta_extract_and_correct(
                    assembly,
                    &allele1_desconstr.chr,
                    pos_11,
                    pos_12,
                    DnaReturnSense::RevC,
                ),
                fasta_extract_and_correct(
                    assembly,
                    &allele2_desconstr.chr,
                    pos_21+1,
                    pos_22+1,
                    DnaReturnSense::RevC,
                ).to_lowercase()
            ];
            let mut sequence = String::new();
            sequence        += &result[0];
            sequence        += &result[1];
            bio::io::fasta::Record::with_attrs(
                &read_id,
                Some(&read_desc),
                sequence.as_bytes()
            )
        },
        _=> panic!("ERROR: combination of directions not valid!"),
    };
    fasta_record
}

/// # Position-based gene-fusion to BND entry
/// this simply takes basic information of a gene-fusion
/// entry and will push write it as a vcf entry.
/// to simplify things, the A entry = true, B entry = false. 
/// Similarly for the direction true= forward and false=reverse.
/// This avoids making a structure as we have only 2 possible 
/// entries for both
/// 
/// Unittest: FALSE
///
pub fn push_pos_vcf (
    is_a: bool,
    forward: StrandDirection,
    forward_2nd: StrandDirection,
    vcf: &mut rust_htslib::bcf::Writer,
    rid: std::option::Option<u32>,
    entry: &mut PrinciplePosBasedFusionOutput,
    counter: u32,
    faidx: & mut bio::io::fasta::IndexedReader<std::fs::File>
){

    let mut ref_seq_a = Vec::new();   
    let mut ref_seq_b = Vec::new();
    let mut cipos     = entry.fp_fuzzy;

    // here we get botht the a and b nucleotide
    // for both of the alleles provided
    faidx.fetch(&entry.a_chrom,entry.a_fp ,entry.a_fp  +1).expect("ERROR: could not fetch interval on reference, either not present or coordinates out of scope ");
    // now it can happen that due to inprecision we have 
    // events which go beyond the coordinates of a chromosome
    // if that happens we take the last matching position and add
    // it to CIPOS.
    if faidx.read(&mut ref_seq_a).is_err() {
        // now if it is beyond the size of chromosome 
        // we need to get size first
        // now lets figure out it's size
        // and then get the extreme
        let mut tmp = Vec::new();   
        faidx.fetch_all(&entry.a_chrom).expect("ERROR: could not fetch entire chromosome");
        faidx.read(&mut tmp).expect("ERROR: could not read sequence from reference");
        // now we get the size
        let c_size = tmp.len();
        // now we get the last position
                    
        faidx.fetch(&entry.a_chrom, (c_size -1) as u64, c_size as u64 ).expect("ERROR: could not fetch interval on reference, either not present or coordinates out of scope ");
        faidx.read(&mut ref_seq_a).expect("ERROR: could not read sequence from reference");
        cipos = match cipos {
            Some(x) => Some(x + (entry.a_fp.abs_diff(c_size as u64) as u32)),
            None         => Some((entry.a_fp.abs_diff(c_size as u64)) as u32),
        };
        eprintln!("INFO: position found outside of boundaries for chr: {} position: {} selecting: {} + cipos {}",
            &entry.a_chrom,
            &entry.a_fp,
            &c_size,
            &cipos.unwrap()
        ); 
        entry.a_fp = c_size as u64; 
        
    }
    let nuc_ref_a: String = String::from(from_utf8(&ref_seq_a).unwrap());

    // same for the b entry
    faidx.fetch(&entry.b_chrom,entry.b_fp ,entry.b_fp +1 ).expect("ERROR: could not fetch interval on reference, either not present or coordinates out of scope ");
    if faidx.read(&mut ref_seq_b).is_err() {
        
       // now lets figure out it's size
       // and then get the extreme
        let mut tmp = Vec::new();   
        faidx.fetch_all(&entry.b_chrom).expect("ERROR: could not fetch entire chromosome");
        faidx.read(&mut tmp).expect("ERROR: could not read sequence from reference");
        // now we get the size
        let c_size = tmp.len();
        // now we get the last position
        faidx.fetch(&entry.b_chrom, (c_size -1) as u64, c_size as u64 ).expect("ERROR: could not fetch interval on reference, either not present or coordinates out of scope ");
        faidx.read(&mut ref_seq_b).expect("ERROR: could not read sequence from reference");
        cipos = match cipos {
            Some(x) => Some(x + ((entry.b_fp.abs_diff(c_size as u64) ) as u32)),
            None         => Some((entry.b_fp.abs_diff(c_size as u64 )) as u32),
        };
        eprintln!("INFO: position found outside of boundaries for chr: {} position: {} selecting: {} + cipos {}",
            &entry.a_chrom,
            &entry.a_fp,
            &c_size,
            &cipos.unwrap()
        ); 
        entry.b_fp = c_size as u64;
    
    }
    let nuc_ref_b: String = String::from(from_utf8(&ref_seq_b).unwrap());

    // here we define the names that we use accordingly
    let name1 = format!("NGSAI_GENEFUSION_PID.{}.1",counter);
    let name2 = format!("NGSAI_GENEFUSION_PID.{}.2",counter);

    let allele_2 : String;
    let alleles  : [&[u8]; 2];
    let mut record = vcf.empty_record();
    let fusion_genes = entry.fusion_genes.iter().map(|x| x.as_bytes()).collect::<Vec<_>>();

    // next we set attributes, such as the type, uncertainty
    // support and flags
    record.set_rid(rid);
    record.push_info_string(b"SVTYPE", &["BND".as_bytes()]).expect("ERROR: could not set SVTYPE to BND");
    // now we defined above the CIPOS and can populate it if 
    // necessary
    // Note: all integer in VCF accepted for INFO field must be i32
    if let Some(value) = cipos{
        record.push_info_integer(b"CIPOS", &[- (value as i32),value as i32 ]).expect("ERROR: could not set CIPOS to BND");
    };
    record.push_info_integer(b"SDP", &[(entry.sr_support  as i32 + entry.sp_support  as i32)]).expect("ERROR: could not set DP to BND");
    record.push_info_integer(b"SP",  &[entry.sp_support   as i32 ]).expect("ERROR: could not set SP to BND");
    record.push_info_integer(b"SR",  &[entry.sr_support   as i32]).expect("ERROR: could not set SR to BND");
    record.push_info_float(b"FR"  , &[entry.sp_ratio]).expect("ERROR: could not set FR to BND");
    record.push_info_string(b"GENE_FUSION",&fusion_genes).expect("ERROR: could not set GENE_FUSION ID!");
    if entry.fp_fuzzy.is_none() || entry.fp_fuzzy.unwrap() > 10 || (entry.sr_support + entry.sp_support) == 1 || entry.sp_ratio < 0.1 {
        record.push_filter("SPLITTY_FAILED".as_bytes()).expect("ERROR: SPLITTY_PASSED not defined in VCF header!");
    }
    // now we construct in BND style the 2nd allele
    // depending whether we analyze the a or b entry
    if is_a {
        record.set_pos((entry.a_fp) as i64);
        record.set_id(name1.as_bytes()).expect("ERROR: could not set ID for vcf entry correctly");
        record.push_info_string(b"MATEID", &[name2.as_bytes()]).expect("ERROR: could not set MATEID");
        let allele2_info = Bnd2ndAllele {
                nucleotide: nuc_ref_a.clone(),
                chr: entry.b_chrom.clone(),
                pos: entry.b_fp,
                prim_direction: forward,
                second_direction: forward_2nd,
        };
        allele_2 = construct_2nd_allele(&allele2_info,true);
        alleles  = [nuc_ref_a.as_bytes(), allele_2.as_bytes()];
    }else {
        record.set_pos((entry.b_fp) as i64);
        record.set_id(name2.as_bytes()).expect("ERROR: could not set ID for vcf entry correctly");
        record.push_info_string(b"MATEID", &[name1.as_bytes()]).expect("ERROR: could not set MATEID");
        let allele2_info = Bnd2ndAllele {
            nucleotide: nuc_ref_b.clone(),
            chr: entry.a_chrom.clone(),
            pos: entry.a_fp,
            prim_direction: forward,
            second_direction: forward_2nd,
        };
        allele_2 = construct_2nd_allele(&allele2_info,true);
        alleles  = [nuc_ref_b.as_bytes(), allele_2.as_bytes()];
    }
    record.set_alleles(&alleles).expect("Failed to set alleles");
    vcf.write(&record).unwrap();
}

/// this function takes our vector of final output
/// and writes it in a vcf compatible format to stdout
/// it adds additionally a header with meta information
/// to identify more easily later the context
/// Needs a reference FASTA file to obtain nucleotide information.
/// The information is returned in a 0-based position and 
/// the fusion point indicated is always the last nucleotide 
/// which is still "normal" with the new sequence adjacent to it.
///
///Unittest: FALSE
///
pub fn write_pos_based_vcf_file(
    results: &mut[PrinciplePosBasedFusionOutput], 
    // version: &str, 
    // author: &str, 
    // command: &str,
    infos: &VersionInfo,
    bam_header: HashMap<String, Vec<LinearMap<String, String>>>,
    ref_file: &str,
    vcf_file: &str,
) -> Result<(), Box<dyn Error>>{

    //////////////////////////////
    //// preparing vcf header ////
    //////////////////////////////
    let mut vcf_header    = Header::new();
    let header_sv_line    = r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">"#;
    let header_ci_line    = r#"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">"#;
    let header_dp_line    = r#"##INFO=<ID=SDP,Number=1,Type=Integer,Description="Support Read Depth (SP+SR) of segment containing breakend">"#;
    let header_sp_line    = r#"##INFO=<ID=SP,Number=1,Type=Integer,Description="Supporting split-pairs of gene-fusion">"#;
    let header_sr_line    = r#"##INFO=<ID=SR,Number=1,Type=Integer,Description="Supporting split-reads of gene-fusion">"#;
    let header_ratio_line = r#"##INFO=<ID=FR,Number=1,Type=Float,Description="SP to adjacent coverage ratio, lowest of both partners">"#;
    let header_mateid_line= r#"##INFO=<ID=MATEID,Number=1,Type=String,Description="the ID of the 2nd part from the BND events">"#;
    let header_filter_line= r#"##FILTER=<ID=SPLITTY_FAILED, Description="splitty events which are likely FP events : SDP == 1 | FR <=0.1 | CIPOS >10 | CIPOS not determined ">"#;
    let header_gf_line    = r#"##INFO=<ID=GENE_FUSION,Number=.,Type=String, Description="Suggested fusion partner genes">"#;
    let header_source_line= format!("{}:{}",r#"##source=splitty"# ,infos.version);
    let header_author_line= format!("{}{}",r#"##source_author="# ,infos.author);
    let header_cmd_line   = format!("{}{}",r#"##command="# ,infos.command);
    vcf_header.push_record(header_source_line.as_bytes());
    vcf_header.push_record(header_author_line.as_bytes());
    vcf_header.push_record(header_cmd_line.as_bytes());
    vcf_header.push_record(header_sv_line.as_bytes());
    vcf_header.push_record(header_ci_line.as_bytes());
    vcf_header.push_record(header_dp_line.as_bytes());
    vcf_header.push_record(header_sp_line.as_bytes());
    vcf_header.push_record(header_sr_line.as_bytes());
    vcf_header.push_record(header_gf_line.as_bytes());
    vcf_header.push_record(header_ratio_line.as_bytes());
    vcf_header.push_record(header_mateid_line.as_bytes());
    vcf_header.push_record(header_filter_line.as_bytes());
    // The hasmap of the header contains all the typical fields
    // and we can get now the chromosomes for the VCF header
    let chroms = match bam_header.get("SQ"){
        Some(x) => x,
        None => panic!("ERROR: could not get chromosomes from BAM header for VCF header!")
    };
    // We have for each entry a key-value pair
    // with the name and length of the chromosome
    for entry in chroms {
        let mut header_chr : Option<&String> = None;
        let mut header_size: Option<&String> = None;
        for (key,val) in entry{
            match key.as_str() {
                "SN" =>  {
                    header_chr = Some(val);
                },
                "LN" => {
                    header_size = Some(val);
                },
                _ => panic!("ERROR: unexpected field in header chromosomes encountered!")
            }
        
        }
        vcf_header.push_record(format!("##contig=<ID={},length={},assembly={}>", header_chr.unwrap(), header_size.unwrap(),ref_file).as_bytes());
    };

    let mut vcf   = Writer::from_path(vcf_file,&vcf_header, true, Format::Vcf).expect("ERROR: could not write vcf to provided file");
    // next we prepare the reference indexed fasta file
    let mut faidx = IndexedReader::from_file(&ref_file).expect("ERROR: could not open reference FASTA file!");
	
    ///////////////////////
    //// populate vcf  ////
    ///////////////////////
    //  Note:
    // This could likely be further optimized. 
    // Initially I thought that I could not simplify
    // or collapse that in some functions, but now
    // towards the end I feel that should actually be
    // feasible and be done.
    let mut counter = 0;
    //let mut alleles: &[&[u8]];
    let mut rid : Option<u32> ; 
    for entry in results.iter_mut() {
        counter += 1;
        match (entry.a_strand,entry.b_strand){
            (StrandDirection::Fwd,StrandDirection::Fwd) => {
                rid = Some(vcf.header().name2rid(entry.a_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_pos_vcf(true,StrandDirection::Fwd,StrandDirection::Fwd,&mut vcf,rid,entry,counter,&mut faidx);
                rid = Some(vcf.header().name2rid(entry.b_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_pos_vcf(false,StrandDirection::Fwd,StrandDirection::Fwd,&mut vcf,rid,entry,counter,&mut faidx);
            },
        (StrandDirection::Rev,StrandDirection::Rev) => {
                rid = Some(vcf.header().name2rid(entry.a_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_pos_vcf(true,StrandDirection::Rev,StrandDirection::Rev,&mut vcf,rid,entry,counter,&mut faidx);
                rid = Some(vcf.header().name2rid(entry.b_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_pos_vcf(false,StrandDirection::Rev,StrandDirection::Rev,&mut vcf,rid,entry,counter,&mut faidx);
            },
        (StrandDirection::Fwd,StrandDirection::Rev) => {
                rid = Some(vcf.header().name2rid(entry.a_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_pos_vcf(true,StrandDirection::Fwd,StrandDirection::Rev,&mut vcf,rid,entry,counter,&mut faidx);
                rid = Some(vcf.header().name2rid(entry.b_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_pos_vcf(false,StrandDirection::Rev,StrandDirection::Fwd,&mut vcf,rid,entry,counter,&mut faidx);
        },
        (StrandDirection::Rev,StrandDirection::Fwd) => {
            rid = Some(vcf.header().name2rid(entry.a_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
            push_pos_vcf(true,StrandDirection::Rev,StrandDirection::Fwd,&mut vcf,rid,entry,counter,&mut faidx);
            rid = Some(vcf.header().name2rid(entry.b_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
            push_pos_vcf(false,StrandDirection::Fwd,StrandDirection::Rev,&mut vcf,rid,entry,counter,&mut faidx);
        }
            _ => eprintln!("ERROR: VCF encountered non-compatible directions!")
        }   
    };
    Ok(())
}

/// #  gene-fusion to BND entry
/// this simply takes basic information of a gene-fusion
/// entry and will push write it as a vcf entry.
/// to simplify things, the A entry = true, B entry = false. 
/// Takes 0 based information and writes VCF in 1-based
/// Similarly for the direction true= forward and false=reverse.
/// This avoids making a structure as we have only 2 possible 
/// entries for both
/// 
/// Unittest: FALSE
///
pub fn push_read_vcf (
    is_a: bool,
    forward: StrandDirection,
    forward_2nd: StrandDirection,
    vcf: &mut rust_htslib::bcf::Writer,
    rid: std::option::Option<u32>,
    entry: &PrincipleReadBasedFusionOutput,
    counter: u32,
    faidx: & mut bio::io::fasta::IndexedReader<std::fs::File>
){
    let mut ref_seq_a = Vec::new();   
    let mut ref_seq_b = Vec::new();

    debug!("Preparing VCF for position a {} and b {}", entry.a_fp, entry.b_fp);
    // get the positions which match nuc a
    faidx.fetch(&entry.a_chrom,entry.a_fp ,entry.a_fp +1).expect("ERROR: could not fetch interval on reference ");
    faidx.read(&mut ref_seq_a).expect("ERROR: could not read sequence from reference");
    let nuc_ref_a: String = String::from(from_utf8(&ref_seq_a).unwrap());

    // get the positions which match nuc a
    faidx.fetch(&entry.b_chrom,entry.b_fp ,entry.b_fp +1).expect("ERROR: could not fetch interval on reference ");
    faidx.read(&mut ref_seq_b).expect("ERROR: could not read sequence from reference");
    let nuc_ref_b: String = String::from(from_utf8(&ref_seq_b).unwrap());

    // generate our identifiers
    let name1 = format!("NGSAI_GENEFUSION_FID.{}.1",counter);
    let name2 = format!("NGSAI_GENEFUSION_FID.{}.2",counter);

    let mut record = vcf.empty_record();
    let allele_2 : String;
    let alleles  : [&[u8]; 2];
    let fusion_genes = entry.fusion_genes.iter().map(|x| x.as_bytes()).collect::<Vec<_>>();
    // this determines the fusion position on the fragment
    let frag_position: i32 = (entry.fusion_point +1 ).try_into().expect("ERROR: could not convert fusion point on fragment into i32");

    // generate and add attributes, concerning type and
    // uncertainty 
    record.set_rid(rid);
    record.push_info_string(b"SVTYPE", &["BND".as_bytes()]).expect("ERROR: could not set SVTYPE to BND");
    if entry.fp_fuzzy.is_some(){
        record.push_info_integer(b"CIPOS", &[-(entry.fp_fuzzy.unwrap() as i32),entry.fp_fuzzy.unwrap()  as i32 ]).expect("ERROR: could not set CIPOS to BND");
    }
    record.push_info_string(b"GENE_FUSION",&fusion_genes).expect("ERROR: could not set GENE_FUSION ID!");
    record.push_info_string(b"FRAG_ID",&[entry.query_name.as_bytes()]).expect("ERROR: could not set FRAG_ID!");
    if entry.fp_fuzzy.is_none() || entry.fp_fuzzy.unwrap() > 10  {
        record.push_filter("SPLITTY_FAILED".as_bytes()).expect("ERROR: SPLITTY_PASSED not defined in VCF header!");
    }
    record.push_info_integer(b"FRAG_POS", &[frag_position]).expect("ERROR: Could not pass the fusion position of the fragment!");

    // now depending whether it is the first or the second one
    // in our pair we need to provide different information and
    // generate our 2nd allele
    if is_a {
        record.set_pos((entry.a_fp) as i64);
        record.set_id(name1.as_bytes()).expect("ERROR: could not set ID for vcf entry correctly");
        record.push_info_string(b"MATEID", &[name2.as_bytes()]).expect("ERROR: could not set MATEID");
        let allele2_info = Bnd2ndAllele {
                nucleotide: nuc_ref_a.clone(),
                chr: entry.b_chrom.clone(),
                pos: entry.b_fp ,
                prim_direction: forward,
                second_direction: forward_2nd,
        };
        allele_2 = construct_2nd_allele(&allele2_info,true);
        alleles  = [nuc_ref_a.as_bytes(), allele_2.as_bytes()];
    }else {
        record.set_pos((entry.b_fp) as i64);
        record.set_id(name2.as_bytes()).expect("ERROR: could not set ID for vcf entry correctly");
        record.push_info_string(b"MATEID", &[name1.as_bytes()]).expect("ERROR: could not set MATEID");
        let allele2_info = Bnd2ndAllele {
            nucleotide: nuc_ref_b.clone(),
            chr: entry.a_chrom.clone(),
            pos: entry.a_fp ,
            prim_direction: forward,
            second_direction: forward_2nd,
        };
        allele_2 = construct_2nd_allele(&allele2_info,true);
        alleles  = [nuc_ref_b.as_bytes(), allele_2.as_bytes()];
    }
    // finally we can write now the 2nd allele and the 1st allele
    // to the vcf writer
    record.set_alleles(&alleles).expect("Failed to set alleles");
    vcf.write(&record).unwrap();

    
}

/// this function takes our vector of final output
/// and writes it in a vcf compatible format to stdout
/// it adds additionally a header with meta information
/// to identify more easily later the context
/// Needs a reference FASTA file to obtain nucleotide information.
/// The information is returned in a 1-based position and 
/// the fusion point indicated is always the last nucleotide 
/// which is still "normal" with the new sequence adjacent to it.
///
/// Unittest: FALSE
///
pub fn write_read_based_vcf_file(
    results: &[PrincipleReadBasedFusionOutput], 
    // version: &str, 
    // author: &str, 
    // command: &str,
    infos: &VersionInfo,
    bam_header: HashMap<String, Vec<LinearMap<String, String>>>,
    ref_file: &str,
    vcf_file: &str,
) -> Result<(), Box<dyn Error>>{

    //////////////////////////////
    //// preparing vcf header ////
    //////////////////////////////
    let mut vcf_header    = Header::new();
    let header_sv_line      = r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">"#;
    let header_ci_line      = r#"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">"#;
    let header_source_line= format!("{}:{}",r#"##source=splitty"# ,infos.version);
    let header_mateid_line  = r#"##INFO=<ID=MATEID,Number=1,Type=String,Description="the ID of the 2nd part from the BND events">"#;
    let header_gf_line      = r#"##INFO=<ID=GENE_FUSION,Number=.,Type=String, Description="Suggested fusion partner genes">"#;    
    let header_fragid_line  = r#"##INFO=<ID=FRAG_ID,Number=1,Type=String, Description="Name of the queried fragment">"#;    
    let header_author_line= format!("{}{}",r#"##source_author="# ,infos.author);
    let header_filter_line  = r#"##FILTER=<ID=SPLITTY_FAILED, Description="splitty events which are likely FP events : CIPOS >10 || CIPOS not determined ">"#;
    let header_fragpos_line = r#"##INFO=<ID=FRAG_POS,Number=1,Type=Integer, Description="Position of the fusion on the fragment">"#;    
    let header_cmd_line   = format!("{}{}",r#"##command="# ,infos.command);
    
    vcf_header.push_record(header_source_line.as_bytes());
    vcf_header.push_record(header_author_line.as_bytes());
    vcf_header.push_record(header_cmd_line.as_bytes());
    vcf_header.push_record(header_mateid_line.as_bytes());
    vcf_header.push_record(header_filter_line.as_bytes());
    vcf_header.push_record(header_gf_line.as_bytes());
    vcf_header.push_record(header_fragid_line.as_bytes());
    vcf_header.push_record(header_sv_line.as_bytes());
    vcf_header.push_record(header_ci_line.as_bytes());
    vcf_header.push_record(header_fragpos_line.as_bytes());

    // The hasmap of the header contains all the typical fields
    // and we can get now the chromosomes for the VCF header
    let chroms = match bam_header.get("SQ"){
        Some(x) => x,
        None => panic!("ERROR: could not get chromosomes from BAM header for VCF header!")
    };
    // We have for each entry a key-value pair
    // with the name and length of the chromosome
    for entry in chroms {
        let mut header_chr : Option<&String> = None;
        let mut header_size: Option<&String> = None;
        for (key,val) in entry{
            match key.as_str() {
                "SN" =>  {
                    header_chr = Some(val);
                },
                "LN" => {
                    header_size = Some(val);
                },
                _ => panic!("ERROR: unexpected field in header chromosomes encountered!")
            }
        
        }
        vcf_header.push_record(format!("##contig=<ID={},length={},assembly={}>", header_chr.unwrap(), header_size.unwrap(),ref_file).as_bytes());
    };

    let mut vcf        = Writer::from_path(vcf_file,&vcf_header, true, Format::Vcf).expect("ERROR: could not write vcf to provided file");
    // next we prepare the reference indexed fasta file
    let mut faidx = IndexedReader::from_file(&ref_file).expect("ERROR: could not open reference FASTA file!");
	

    ///////////////////////
    //// populate vcf  ////
    ///////////////////////
    //  Note:
    // This could likely be further optimized. 
    // Initially I thought that I could not simplify
    // or collapse that in some functions, but now
    // towards the end I feel that should actually be
    // feasible and be done.
    let mut counter : u32 = 0;
    //let mut alleles: &[&[u8]];
    //let mut allele_2 : String;
    let mut rid : Option<u32> ; 
    for entry in results.iter() {
        counter += 1;
        debug!("Analyzing now {:?}", entry);
        match (entry.a_strand,entry.b_strand){
        (StrandDirection::Fwd,StrandDirection::Fwd) => {
                rid = Some(vcf.header().name2rid(entry.a_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_read_vcf(true,StrandDirection::Fwd,StrandDirection::Fwd,&mut vcf,rid,entry,counter,&mut faidx);
                rid = Some(vcf.header().name2rid(entry.b_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_read_vcf(false,StrandDirection::Fwd,StrandDirection::Fwd,&mut vcf,rid,entry,counter,&mut faidx);
            },
        (StrandDirection::Rev,StrandDirection::Rev) => {
                rid = Some(vcf.header().name2rid(entry.a_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_read_vcf(true,StrandDirection::Rev,StrandDirection::Rev,&mut vcf,rid,entry,counter,&mut faidx);
                rid = Some(vcf.header().name2rid(entry.b_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_read_vcf(false,StrandDirection::Rev,StrandDirection::Rev,&mut vcf,rid,entry,counter,&mut faidx);
            },
        (StrandDirection::Fwd,StrandDirection::Rev) => {
                rid = Some(vcf.header().name2rid(entry.a_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_read_vcf(true,StrandDirection::Fwd,StrandDirection::Rev,&mut vcf,rid,entry,counter,&mut faidx);
                rid = Some(vcf.header().name2rid(entry.b_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_read_vcf(false,StrandDirection::Rev,StrandDirection::Fwd,&mut vcf,rid,entry,counter,&mut faidx);
        },
        (StrandDirection::Rev,StrandDirection::Fwd) => {
                rid = Some(vcf.header().name2rid(entry.a_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_read_vcf(true,StrandDirection::Rev,StrandDirection::Fwd,&mut vcf,rid,entry,counter,&mut faidx);
                rid = Some(vcf.header().name2rid(entry.b_chrom.as_bytes()).unwrap());//.expect("ERROR: chromosome does not exist in header!");
                push_read_vcf(false,StrandDirection::Fwd,StrandDirection::Rev,&mut vcf,rid,entry,counter,&mut faidx);
        }
            _ => eprintln!("ERROR: VCF encountered non-compatible directions!")
        }   
    };
    Ok(())
}

/// This functions works on vcf/bcf records and verifies if CIPOS
/// is defined. If that's the case it returns a vector with the entries
/// within an Option, otherwise result is None
/// CIPOS according to VCF spec is negative for first and positive for 2nd entry
/// We return 2 positive values
///
/// Unittest: TRUE
///
pub fn get_cipos(
    record: & rust_htslib::bcf::Record
) -> Option<Vec<u32>> {
    let info = match record.info(b"CIPOS").integer() {
        Ok(x) => {
            match x {
               Some(y) =>  Some(y),
               None => None
            }
            //x.map(|y| y)
        },
        Err(_x)     => None
    };
    match info {
    None => None,
    Some(x) => {
        let cipos = x.deref();
        let tmp  : Vec<i32> = cipos.iter().copied().collect();
        let tmp2 : Vec<u32> = tmp.iter().map(|x| u32::try_from(x.abs()).unwrap() ).collect();
        Some(tmp2)
        }
    }
}

/// this function expects a fofn of vcf, one vcf path per line.
/// This is then opened, parsed and a interval tree is generated containing the most
/// key information of all samples
/// It is necessary to provide a hashmap with chromosomes as key and their length as value.
/// It allows to make more efficient structures and is needed to then build clusters, too.
/// Vcf is 1 based and the tree is 0 based, as well as all linked functions.
/// There is though a bug and we must for a range of e.g. 1..1 always use 1..2 and similarly for a 1..2 use 1..3.
/// The documentation says it should work the other way but it just does not.
/// Can be set to be direction agnostic with "ignore_dir=true"
/// It returns a structure of multiple trees
///
/// Unittest: TRUE
///
pub fn vcf_parsing_tree(
    vcf_fofn:   &str ,
    threads:    usize,
    contigs:    &FxHashMap<String, u64>,
    ignore_dir: bool,
    real_fofn:  bool
) -> FullSvTrees {
    
    // first we establish an empty tree for each contig
    let mut bndpair_trees   : FxHashMap<String,IntervalTree<u64,BNDentry>>    = FxHashMap::default();
    let mut bndsingle_trees : FxHashMap<String,IntervalTree<u64,BNDentry>>    = FxHashMap::default();
    let mut indel_trees     : FxHashMap<String,IntervalTree<u64,INDELentry>>  = FxHashMap::default();

    for key in contigs.keys() {
         // make a new IntervalTree
        let tree1 = IntervalTree::new();
        let tree2 = IntervalTree::new();
        let tree3 = IntervalTree::new();
        
        bndpair_trees.insert(key.clone()  , tree1.clone());
        bndsingle_trees.insert(key.clone(), tree2.clone());
        indel_trees.insert(key.clone()    , tree3.clone());
    };
    let mut unknown = 0_u64;
    let fofn : Vec<String> = parse_fofn(vcf_fofn,real_fofn);
    for file in fofn { 
        debug!("Analyzing file: {}",file);       
        // for the following ones, rusts insists that we should remove
        // the mutable part, but this will fail, therefore the exception
        #[allow(unused_mut)]
        let mut vcf_a      = rust_htslib::bcf::Reader::from_path(&file).expect("ERROR: could not open vcf file in fofn !");
        // this is currently extremly stupid as I have to open twice as otherwise I can get the sample ID ????
        #[allow(unused_mut)]
        let mut vcf_b      = rust_htslib::bcf::Reader::from_path(&file).expect("ERROR: could not open vcf file in fofn !");
        vcf_a.set_threads(threads).expect("ERROR: could not set correctly read threads");
        #[allow(unused_mut)]
        let mut vcf_header: Vec<&[u8]> = vcf_b.header().samples();
        let sample : &str;
        if vcf_header.len() > 1 {
            panic!("ERROR: Currently only single sample vcf accepted!");
        }else if vcf_header.is_empty() {
            //eprintln!("WARNING: No sample defined in the vcf file - will use file name instead!");
            sample = match Path::new(&file).file_stem() {
                Some(x) => x.to_str().unwrap(),
                None => panic!("ERROR: could not extract basename from path!"),
            }
        }else{
            sample = str::from_utf8(vcf_header[0]).unwrap();
        }
        // Now we parse all the records for that file
        for entry in vcf_a.records() {
            let record  = entry.expect("ERROR: could not read record of vcf A file!");
            let sample_count = usize::try_from(record.sample_count()).unwrap();   
            if sample_count >1 { 
                panic!("ERROR: Currently only single sample vcf accepted!")
            };
            // depending on it's type we need to add the entry to a different tree and treat as well
            // differently
            let sv_full       = record.info(b"SVTYPE").string().unwrap().unwrap();
            let sv_rec        = str::from_utf8(sv_full[0]).expect("ERROR: could not extract SVTYPE string!");

            let mut sv_type = match sv_rec{
                "INS" => SVType::INS,
                "DEL" => SVType::DEL,
                "BND" => SVType::BndPair,
                _     => SVType::Unknown,
            };
            // let's gather first the BND independent values
            debug!("SV type of entry is {:?}",sv_type);
            let id1 = &record.id() ;
            let id      = str::from_utf8(id1).expect("ERROR: cant extract record id!");
            // we do not accept any non named entry 
            if id == "." {
                panic!("ERROR: all entries must have an annotated entry ID!");
            }
            let rid             = record.rid().expect("ERROR: could not get reference id!");
            let chromosome     = str::from_utf8(record.header().rid2name(rid).expect("ERROR: could not convert rid to chromosome name!")).unwrap();
            let chr1         = chromosome.to_string().clone();
            let chr_length     = contigs.get(chromosome).expect("ERROR: chromosome not existing or no length associated!!");
            let fp1             = u64::try_from(record.pos()).expect("ERROR: got negative chromosomal coordinates!");
            let alleles  = record.alleles();

            // here we currently default these but will
            // add the posibility to have as well event types correctly 
            // classified
            if alleles.len() != 2 {
                panic!("ERROR: encountered more than 2 alleles, which is not supported!")
            };
            let allele_2 = str::from_utf8(alleles[1]).expect("ERROR: could not convert allele!");
            let name  : Option<String> = None;
            // if an entry has an annotated CIPOS, this is currently supported then
            // NOTE: I did not find a proper way yet to test if e.g. multipe CIPOS exist and invalidate our test!
            let (mut left, mut right) : (u64,u64) = (0,0);
            let cipos = get_cipos(&record);
            if cipos.is_some() {
                let positions =  cipos.unwrap();
                if positions.len() !=2 {
                    panic!("ERROR: CIPOS encountered with !=2 entries!");
                }else{
                    left  = u64::try_from(positions[0]).unwrap() ;
                    right = u64::try_from(positions[1]).unwrap() ;
                }
            }

            if sv_type == SVType::BndPair {
                // avoid dropping below 0
                let mut st1   = match fp1.cmp(&left){
                    std::cmp::Ordering::Equal => 0,
                    std::cmp::Ordering::Greater => fp1 - left,
                    std::cmp::Ordering::Less => 0
                } ;
                // avoid exceeding length of chrom (-1 as 0 based)
                let end1  = match chr_length.cmp(&(fp1 + right +1)){
                    std::cmp::Ordering::Equal   => chr_length-1,
                    std::cmp::Ordering::Greater => fp1 + right +1,
                    std::cmp::Ordering::Less    => chr_length-1,
                };
                // nowe in case of extreme towards the end we can experience that start is
                // last nucleotide and the end is last nucleotide 
                debug!("start1: {}, end1: {}", st1, end1);
                if end1 ==  chr_length-1 && st1==end1 {
                    st1 = end1-1;
                };
                debug!("start1: {}, end1: {}", st1, end1);
                let add_fields : BNDentry;
                // now if we have no proper pairs then we add it in single breakend events.
                // actually as well something valid if the second part of the event is unknown
                if record.info(b"MATEID").string().is_err() || record.info(b"MATEID").string().unwrap().is_none() {
                    debug!("No MATEID encountered");
                    // this one will be more tricky to get as not one of standard fields
                    sv_type = SVType::BndSingle;
                    let allele_deconstr = deconstruct_2nd_allele(allele_2,false);
                    if ignore_dir {
                        add_fields = BNDentry{
                            sample: sample.to_string(),
                            id: id.to_string(),
                            mid: String::new(),
                            chr1,
                            fp1,
                            st1,
                            end1,
                            forward1: StrandDirection::Unknown,
                            chr2 : String::new(),
                            fp2  : 0,
                            st2  : 0,
                            end2 : 0,            
                            forward2: StrandDirection::Unknown,
                            sv_type,
                            name
                        };
                    }else {
                        add_fields = BNDentry{
                            sample: sample.to_string(),
                            id: id.to_string(),
                            mid: String::new(),
                            chr1,
                            fp1,
                            st1,
                            end1,
                            forward1: allele_deconstr.prim_direction,
                            chr2 : String::new(),
                            fp2  : 0,
                            st2  : 0,
                            end2 : 0,            
                            forward2: StrandDirection::Unknown,
                            sv_type,
                            name
                            };
                    };
                    debug!("BND entry contained the following fields: {:?}",add_fields);
                    let local_tree = bndsingle_trees.get_mut(chromosome).expect("ERRROR: could not find chromosome in tree!");
                    local_tree.insert(st1..end1,add_fields);                 
                }else{
                    debug!("MATEID encountered");
                    let mate_id1       = record.info(b"MATEID").string().unwrap().unwrap();
                    let mate_id        = str::from_utf8(mate_id1[0]).expect("ERROR: could not extract MateID string!");
                    let allele_deconstr = deconstruct_2nd_allele(allele_2,true);
                    // the 1st entry needs not to be corrected
                    // as the API of vcf is returning 0-based
                    // and mine as well
                    let chr2         = allele_deconstr.chr.clone();
                    let chr2_length     = contigs.get(&chr2).expect("ERROR: chromosome not existing or no length associated!!");
                    
                    debug!("What is second allele: {:?}",&allele_deconstr);
                    // here we avoid going below 0 (which we cant with the type)
                    let mut st2   = match allele_deconstr.pos.cmp(&left){
                        std::cmp::Ordering::Equal => 0,
                        std::cmp::Ordering::Greater => allele_deconstr.pos - left,
                        std::cmp::Ordering::Less => 0
                    } ;
                    // avoid exceeding length of chrom (-1 as 0 based)
                    // here we need the length of the chromosome2 and not chromosome1
                    let end2  = match chr2_length.cmp(&(allele_deconstr.pos + right +1)){
                        std::cmp::Ordering::Equal   => chr2_length-1,
                        std::cmp::Ordering::Greater => allele_deconstr.pos + right +1,
                        std::cmp::Ordering::Less    => chr2_length-1,
                    };
                    debug!("start2: {}, end2: {}", st2, end2);
                    // nowe in case of extreme towards the end we can experience that start is
                    // last nucleotide and the end is last nucleotide 
                    if end2 ==  chr2_length-1 && st2==end2 {
                        st2 = end2-1;
                    };
                    debug!("start2: {}, end2: {}", st2, end2);

                    // populate simple comparison structure
                    if ignore_dir {
                        add_fields = BNDentry{
                            sample: sample.to_string(),
                            id: id.to_string(),
                            mid: mate_id.to_string(),
                            chr1,
                            fp1,
                            st1,
                            end1,
                            forward1: StrandDirection::Unknown,
                            chr2,
                            fp2  : allele_deconstr.pos,
                            st2,
                            end2,            
                            forward2: StrandDirection::Unknown,
                            sv_type,
                            name
                        };
                    }else {
                        add_fields = BNDentry{
                            sample: sample.to_string(),
                            id: id.to_string(),
                            mid: mate_id.to_string(),
                            chr1,
                            fp1,
                            st1,
                            end1,
                            forward1 : allele_deconstr.prim_direction,
                            chr2,
                            fp2      : allele_deconstr.pos,
                            st2,
                            end2,            
                            forward2 : allele_deconstr.second_direction,
                            sv_type,
                            name
                        };
                    }
                    debug!("BND entry contained the following fields: {:?}",add_fields);
                    let local_tree = bndpair_trees.get_mut(chromosome).expect("ERRROR: could not find chromosome in tree!");
                    debug!("Inserting start {} and end {} for chr {} with fields {:?} into tree {:?}",st1,st2,chromosome,add_fields,local_tree);
                    local_tree.insert(st1..end1,add_fields);
                    debug!("Inserted now entry {:?}" ,bndpair_trees.get_mut(chromosome));
                    
                } 
            }else if sv_type == SVType::INS || sv_type == SVType::DEL {
                // avoid dropping below 0
                let st   = match fp1.cmp(&left){
                    std::cmp::Ordering::Equal => 0,
                    std::cmp::Ordering::Greater => fp1 - left,
                    std::cmp::Ordering::Less => 0
                } ;
                // avoid exceeding length of chrom (-1 as 0 based)
                let end  = match chr_length.cmp(&(fp1 + right +1)){
                    std::cmp::Ordering::Equal   => chr_length-1,
                    std::cmp::Ordering::Greater => fp1 + right +1,
                    std::cmp::Ordering::Less    => chr_length-1,
                };
                // populate simple comparison structure
                let add_fields = INDELentry{
                    sample: sample.to_string(),
                    id: id.to_string(),
                    chr:chr1,
                    st,
                    end,
                    sv_type,
                    name
                };
                debug!("SV entry contained the following fields: {:?}",add_fields);
                let local_tree = indel_trees.get_mut(chromosome).expect("ERRROR: could not find chromosome in tree!");
                local_tree.insert(st..end,add_fields);
            }else{
                unknown += 1;
            };
        }
    };
    if unknown != 0 {
        eprintln!("WARNING: {} ignored features as not part of recognized features: BND(paired),INS,DEL",unknown);
    }

    FullSvTrees {
        bnd_pair  : bndpair_trees,
        bnd_single: bndsingle_trees,
        indel     : indel_trees,
    }
}

/// This function expects a NGSAI fusion analysis derived file which contains a
/// fragment-based result with a "FRAG_ID". 
/// It then generates a hashmap with a unique identifier and a vcf entry and a 
/// structure with simplified information linking to that record. 
/// 
/// Unittest: FALSE
///
pub fn vcf_a_parsing (
    vcf_file: &str,
    threads: usize,
    simple: &bool
) -> (FxHashMap<String,rust_htslib::bcf::Record>, Vec<FileAbasedComp> ) {
    let mut full_a_vcf : FxHashMap<String,rust_htslib::bcf::Record> = FxHashMap::default();
    let mut simple_a   = Vec::new();
    let mut vcf_a      = rust_htslib::bcf::Reader::from_path(vcf_file).expect("ERROR: could not open vcf file A !");
    //let header         = vcf_a.header();

    vcf_a.set_threads(threads).expect("ERROR: could not set correctly read threads");
    for entry in vcf_a.records() {
        let record   = entry.expect("ERROR: could not read record of vcf A file!");
        if  !simple && record.info(b"FRAG_ID").string().unwrap().is_none(){
            panic!("ERROR: file A needs to contain FRAG_ID fields!");
        }
        let sample_count = usize::try_from(record.sample_count()).unwrap();   
        if sample_count != 0 { panic!("ERROR: Currently only single sample vcf accepted!")};
        let id1 = &record.id() ;
        let id  = str::from_utf8(id1).expect("ERROR: cant extract record id!");
        if id == "." {
            panic!("ERROR: all entries must have an annotated entry ID!");
        }
        // this one will be more tricky to get as not one of standard fields
        let mate_id1  = record.info(b"MATEID").string().unwrap().unwrap();
        let mate_id   = str::from_utf8(mate_id1[0]).expect("ERROR: could not extract MateID string!");
        let a_id    = format!("{}|{}",id,mate_id);
        let a_rid      = record.rid().expect("ERROR: could not get reference id!");
        let a_chr1  = str::from_utf8(record.header().rid2name(a_rid).expect("ERROR: could not convert rid to chromosome name!")).unwrap().to_string();
        // converting should be safe here as
        // bcf from htslib already does not allow negative 
        // values
        let a_fp1: u64     = record.pos().try_into().unwrap() ;

        let alleles = record.alleles();
        if alleles.len() != 2 { panic!("ERROR: encountered more than 2 alleles, which is not supported!")};
        let allele_2 = str::from_utf8(alleles[1]).expect("ERROR: could not convert allele!");
        // here we use again our function which will deconstruct the 2nd allele
        let allele2_desconstr : Bnd2ndAllele = deconstruct_2nd_allele(allele_2,true);
        // populate hashmap
        full_a_vcf.insert(a_id.clone(), record);
        // populate simple comparison structure
        let add_fields = FileAbasedComp{
            a_id,
            a_chr1,
            a_fp1,
            a_chr2 : allele2_desconstr.chr ,
            a_fp2  : allele2_desconstr.pos,
            b_matches: Vec::new(),
        };
        simple_a.push(add_fields);

    }
    (full_a_vcf, simple_a)
}

/// This function expects a NGSAI fusion analysis derived file which contains a
/// position based result with a. 
/// It then generates a hashmap with a unique identifier and a vcf entry and a 
/// structure with simplified information linking to that record. 
///
/// Unittest: FALSE
///
pub fn vcf_b_parsing (
    vcf_file: &str, 
    range: u32, 
    mut simple_a: Vec<FileAbasedComp>,
    threads: usize,
    simple: &bool
) -> (FxHashMap<String,rust_htslib::bcf::Record>, Vec<FileAbasedComp> ){
    let mut full_b_vcf : FxHashMap<String,rust_htslib::bcf::Record> = FxHashMap::default();
    let mut vcf_b      = rust_htslib::bcf::Reader::from_path(vcf_file).expect("ERROR: could not open vcf file A !");
    //let header         = vcf_a.header();

    vcf_b.set_threads(threads).expect("ERROR: could not set correctly read threads");
    for entry in vcf_b.records() {
        let record   = entry.expect("ERROR: could not read record of vcf B file!");
        if !simple && record.info(b"SP").integer().unwrap().is_none() {
            panic!("ERROR: file B needs to contain SP fields!");
        }
        let sample_count = usize::try_from(record.sample_count()).unwrap();   
        if sample_count != 0 { panic!("ERROR: Currently only single sample vcf accepted!")};
        
        let id1 = &record.id() ;
        let id  = str::from_utf8(id1).expect("ERROR: cant extract record id!");
        if id == "." {
            panic!("ERROR: all entries must have an annotated entry ID!");
        }     
        // this one will be more tricky to get as not one of standard fields
        let mate_id1  = record.info(b"MATEID").string().unwrap().unwrap();
        let mate_id   = str::from_utf8(mate_id1[0]).expect("ERROR: could not extract MateID string!");
        let b_id      = format!("{}|{}",id,mate_id);
        let b_rid     = record.rid().expect("ERROR: could not get reference id!");
        let b_chr1    = str::from_utf8(record.header().rid2name(b_rid).expect("ERROR: could not convert rid to chromosome name!")).unwrap().to_string();
        // converting should be safe here as
        // bcf from htslib already does not allow negative 
        // values
        let b_fp1: u64     = record.pos().try_into().unwrap();
        
        let alleles = record.alleles();
        if alleles.len() != 2 { panic!("ERROR: encountered more than 2 alleles, which is not supported!")};
        let allele_2 = str::from_utf8(alleles[1]).expect("ERROR: could not convert allele!");
        let allele2_desconstr : Bnd2ndAllele = deconstruct_2nd_allele(allele_2,true);
        let b_chr2 = allele2_desconstr.chr;
        let b_fp2     = allele2_desconstr.pos;
        // populate hashmap
        full_b_vcf.insert(b_id.clone(), record);
        for entry in & mut simple_a {
            if b_chr1 != entry.a_chr1 || b_chr2 != entry.a_chr2 {
                continue
            }else{
                let a_range1 = (entry.a_fp1 - range as u64)..(entry.a_fp1 + range as u64);
                let a_range2 = (entry.a_fp2 - range as u64)..(entry.a_fp2 + range as u64);
                if a_range1.contains(&b_fp1) && a_range2.contains(&b_fp2){
                    entry.b_matches.push(b_id.clone());
                }

            }

        }
    }
    (full_b_vcf,simple_a)
}


/// This function now takes again the hashmap of both files
/// and the A entry with matching B entries. It the 
/// goes over our simplify structured and we do 3 things: 
/// - for entries with a matching B combine both and remove B entry in B hashmap
/// - for entries without a matching B add A entry solely
/// - for left B entries write them accordingly as well
/// If both A and B have an entry with "SPLITTY_FAILED" the combined one
/// receives only a "MED" evidence score abd flag is always removed. Otherwise we always judge this
/// as "HIGH" for all other combinations. 
/// For single entries with "SPLITTY_FAILED" this FLAG remains and evidence becomes "LOW"
/// For single entries without "SPLITTY_FAILED" evidence becomes "MED"
///
///Unittest: FALSE
///
pub fn vcf_combine_a_and_b(
    simpl_vcf_a_modified: Vec<FileAbasedComp>, 
    full_vcf_a: & mut FxHashMap<String,rust_htslib::bcf::Record>, 
    full_vcf_b: & mut FxHashMap<String,rust_htslib::bcf::Record>,
    vcf: & mut rust_htslib::bcf::Writer,
    simple: &bool,
    simple_flag: &Option<&str>,
    info_fields: &Vec<&str>,
    //info_frmt: &Vec<&str>
){
    let mut seen_b : FxHashMap<String,bool> = FxHashMap::default();
    for entry in simpl_vcf_a_modified {
        let record_a = full_vcf_a.get_mut(&entry.a_id).expect("ERROR: could not recover vcf entry from full_vcf_a");
        vcf.translate(record_a);
        // if no matching B entry, directly push into new vcf
        if entry.b_matches.is_empty() {
            if *simple {
                record_a.push_info_string(simple_flag.unwrap().as_bytes(),&["0".as_bytes()]).expect("ERROR: could not set CONF field!");
            }else{
                // if single entries have SPLITTY_FAILED then we keep it
                // and add additionally evidence to low
                if record_a.has_filter("SPLITTY_FAILED".as_bytes()) {
                    record_a.push_info_string(b"EVI",&["LOW".as_bytes()]).expect("ERROR: could not set CONF field!");    
                }else{
                    record_a.push_info_string(b"EVI",&["MED".as_bytes()]).expect("ERROR: could not set CONF field!");
                    record_a.remove_filter("SPLITTY_FAILED".as_bytes(),false).expect("ERROR: could not remove SPLITTY_FAILED filter!");

                };
            }
        // if there is a matching b entry    
        } else  if *simple {
            record_a.push_info_string(simple_flag.unwrap().as_bytes(),&["1".as_bytes()]).expect("ERROR: could not set CONF field!");
        }else{
            let conf_a = match record_a.has_filter("SPLITTY_FAILED".as_bytes()) {
                false  => "MED",
                true => "LOW",
            };
            // if we have evidence beyond n=1 then
            // sum them up.
            // Important: if one changes the order it
            // will become 0 for all of them
            let mut sdp = Vec::new();
            let mut sp  = Vec::new();
            let mut sr  = Vec::new();
            let mut fr  = Vec::new();
            let mut ev: Vec<Vec<u8>> = Vec::new();
            // this is for our new pass-through fields
            let mut i_fields: Vec<Vec<i32>> = Vec::new();

            let mut conf_b = "LOW" ;
            if record_a.info(b"SDP").integer().unwrap().is_some() {
                panic!("ERROR: please change File1 and File2");
            }
            // now lets have a look at the entries from b which match
            for b_entry in entry.b_matches {
                let record_b = full_vcf_b.get_mut(&b_entry).expect("ERROR: could not recover vcf entry from full_vcf_a");
                vcf.translate(record_b);
                // here it could happend that we have once low and once med
                // and would override old value
                conf_b = match record_b.has_filter("SPLITTY_FAILED".as_bytes()) {
                    false  => "MED",
                    true => "LOW",
                };
              
                ev.push(record_b.id());
                let mut tmp_fields : Vec<i32> = Vec::new();
                for i_field in info_fields.iter(){
                    tmp_fields.push(record_b.info(i_field.as_bytes()).integer().unwrap().unwrap()[0]);
                };
                i_fields.push(tmp_fields);
                if record_b.info(b"SR").integer().unwrap().is_some() {
                    sdp.push(record_b.info(b"SDP").integer().unwrap().unwrap()[0]);
                    sp.push(record_b.info(b"SP").integer().unwrap().unwrap()[0]);
                    sr.push(record_b.info(b"SR").integer().unwrap().unwrap()[0]);
                    fr.push(record_b.info(b"FR").float().unwrap().unwrap()[0]);
                }else {
                    panic!("ERROR: please change File1 and File2");
                }
                seen_b.insert(b_entry,true);             
            }
            // now we sort and keep the highest one for split reads as it is the 
            // most informative for this purpose
            let my_max = sr.iter().position_max().unwrap();
              // now we collect the interesting values
            let id1 = &ev[my_max] ;
            let id      = str::from_utf8(&id1).expect("ERROR: cant extract record id!");
            
            // now we should not just add them up as they are predominantly the same reads
            record_a.push_info_string(b"TOPREAD",&[id.as_bytes()]).expect("ERROR: could not push TOPREAD info!");
            record_a.push_info_integer(b"SDP",&[sdp[my_max]]).expect("ERROR: could not push SDO info!");
            record_a.push_info_integer(b"SP",&[sp[my_max]]).expect("ERROR: could not push SP info!");
            record_a.push_info_integer(b"SR",&[sr[my_max]]).expect("ERROR: could not push SR info!");
            record_a.push_info_float(b"FR",&[fr[my_max]]).expect("ERROR: could not push FR info!");
            for (n,element) in info_fields.iter().enumerate() {
                record_a.push_info_integer(element.as_bytes(),&[i_fields[my_max][n]]).expect("ERROR: could not push SR info!");
            }
            let final_conf = match(conf_a,conf_b){
                ("LOW","LOW") => "MED",
                _ => "HIGH",
            };
            // remove now the obsolete filter field as it is actually stronger evidence now
            record_a.remove_filter("SPLITTY_FAILED".as_bytes(),true).expect("ERROR: could not remove SPLITTY_FAILED filter!");
            record_a.push_info_string(b"EVI",&[final_conf.as_bytes()]).expect("ERROR: could not set CONF field!");
        
        }
        vcf.write(record_a).expect("ERROR: could not write vcf entry!");

    }
    // now iterate over b and check if already seen 
    // otherwise write, too
    // I added here sorted as they were otherwise un-ordered
    // which was a bit of a nuisance
    if !simple {
        for (key,record) in full_vcf_b.iter_mut().sorted_by_key(|x| x.0) {
            // now if the entry was actually not in
            // a we need to write b as well
            if !seen_b.contains_key(key) {
                vcf.translate(record);
                // if single entries have SPLITTY_FAILED then we keep it
                // and add additionally evidence to low
                if record.has_filter("SPLITTY_FAILED".as_bytes()) {
                    record.push_info_string(b"EVI",&["LOW".as_bytes()]).expect("ERROR: could not set CONF field!");    
                }else{
                    record.push_info_string(b"EVI",&["MED".as_bytes()]).expect("ERROR: could not set CONF field!");
                    record.remove_filter("SPLITTY_FAILED".as_bytes(),false).expect("ERROR: could not remove SPLITTY_FAILED filter!");

                }
                vcf.write(record).expect("ERROR: could not write vcf entry!");
            }
        }
    }
}

/// this function prepares a valid vcf writer
/// specifically for single BND entries.
/// It expects the reference FASTA sequence
/// to properly prepare the header with the needed
/// information, as well as the command, author and 
/// version of the used program.
/// It provides the possibility to add a list of Strings
/// which can be pushed additionally into the header
/// depending on the usecase, e.g.
/// ##INFO=<ID=Contig_ID,Number=1,Type=String, Description="Name of the queried fragment">
///
///Unittest: FALSE
///
pub fn prep_bnd_vcf (
    infos: &VersionInfo,
    add_lines : Option<Vec<&str>>,
    ref_file: &str,
    vcf_out: &str,
) -> rust_htslib::bcf::Writer {

    // first, construct the header with all chromosomes and length.
    // I add as well the number of the used assembly file as this is often
    // very handy later
    let mut vcf_header    = rust_htslib::bcf::Header::new();

    let header_sv_line      = r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">"#;
    let header_ci_line      = r#"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">"#;
    let header_source_line= format!("{}:{}",r#"##source=splitty"# ,infos.version);    
    let header_author_line= format!("{}{}",r#"##source_author="# ,infos.author);
    let header_cmd_line   = format!("{}{}",r#"##command="# ,infos.command);
    let header_mateid_line  = r#"##INFO=<ID=MATEID,Number=1,Type=String,Description="the ID of the 2nd part from paired BND events">"#;
   

    vcf_header.push_record(header_source_line.as_bytes());
    vcf_header.push_record(header_author_line.as_bytes());
    vcf_header.push_record(header_cmd_line.as_bytes());
    vcf_header.push_record(header_sv_line.as_bytes());
    vcf_header.push_record(header_ci_line.as_bytes());
    vcf_header.push_record(header_mateid_line.as_bytes());
    if add_lines.is_some() {
        for element in add_lines.unwrap() {
            vcf_header.push_record(element.as_bytes());
        }
    }
    let ref_fasta =  IndexedReader::from_file(&ref_file).expect("ERROR: could not open reference FASTA file!");
    let idx_fasta = ref_fasta.index;

    // We have for each entry a key-value pair
    // with the name and length of the chromosome
    for entry in idx_fasta.sequences() {
        //dbg!(&entry);
        let key     = entry.name;
        let value   = entry.len.to_string();
        //eprintln!("Name: {:?} Length: {:?}",&key,&value);
        vcf_header.push_record(format!("##contig=<ID={},length={},assembly={}>", key, value, ref_file).as_bytes());
    };


    // important, writer initiation always has to come after the entire header
    // was constructed. Not a well documented and logic
    let vcf_writer        = rust_htslib::bcf::Writer::from_path(vcf_out,&vcf_header,true,Format::Vcf).expect("ERROR: could not generate VCF file!");
    let header_view = vcf_writer.header();    
    eprintln!("Reference FASTA file contained {} chromosomes",header_view.contig_count());

    vcf_writer
}

/// takes a VCF single BND writer, verifies that all is good
/// and then pushes information from AnnotatedRead e.g. (muERV)
/// into said writer. Does not decide if should be sorted or not
/// and takes it as it comes
///
///Unittest: FALSE
///
pub fn write_single_bnd_vcf(
    ref_file    : &str,
    vcf_out     : &mut rust_htslib::bcf::Writer,
    annot_reads : std::vec::Vec<AnnotatedRead>,
    event_type  : &str,
){
    let mut faidx =  IndexedReader::from_file(&ref_file).expect("ERROR: could not open reference FASTA file!");
    // just incremental counter for the IDs within a given file
    let mut counter : u32 = 0;
    // will be reference ID
    let mut rid : Option<u32> ;


    // this is the vessel for the later fetched sequence
    let mut ref_seq_a = Vec::new();   

    for event in annot_reads{
        //dbg!(&event);
        counter += 1;
        let mut record = vcf_out.empty_record();
        // double check that chromosome actually exists in our VCF header
        rid = Some(vcf_out.header().name2rid(event.ref_chr.as_bytes()).expect("ERROR: chromosome does not exist in header!"));
        record.set_rid(rid);
        // generate our identifiers
        let name = format!("NGSAI_{}_ID.{}",event_type,counter);
        record.set_pos(event.integration_site.try_into().unwrap());
        record.set_id(name.as_bytes()).expect("ERROR: could not set ID for vcf entry correctly");
        // generate and add attributes, concerning type and
        // uncertainty 
        record.push_info_string(b"SVTYPE", &["BND".as_bytes()]).expect("ERROR: could not set SVTYPE to BND");
        // next let's push some useful additional information
        record.push_info_string(b"Contig_ID",&[event.read_name.as_bytes()]).expect("ERROR: could not set FRAG_ID!");
        // source of the event
        record.push_info_integer(b"ReadFragRef",&[event.ref_start as i32,event.ref_end as i32]).expect("ERROR: could not set read coordinates for REF!");
        record.push_info_integer(b"ReadFragAlien",&[event.alien_start as i32,event.alien_end as i32]).expect("ERROR: could not set read coordinates for ALIEN!");
        // get the contig and positions which matches the integration point
        faidx.fetch(
            &event.ref_chr, 
            event.integration_site -1 ,
            event.integration_site  )
            .expect("ERROR: could not fetch interval on reference ");
        faidx.read(&mut ref_seq_a).expect("ERROR: could not read sequence from reference");
        // now we have the first allele sequence
        let allele_1 = String::from(from_utf8(&ref_seq_a).unwrap());

        let allele2_info = Bnd2ndAllele {
                nucleotide: allele_1.clone(),
                chr: String::from("NA"),
                pos: 0,
                prim_direction: event.integration_dir,
                second_direction: StrandDirection::Unknown,
        };
        let allele_2 = construct_2nd_allele(&allele2_info,false);
        //dbg!(&allele_2);
        let alleles  : [&[u8]; 2] = [allele_1.as_bytes(), allele_2.as_bytes()];

        // finally we can write now the 2nd allele and the 1st allele
        // to the vcf writer
        record.set_alleles(&alleles).expect("Failed to set alleles");
        //dbg!(&record);
        vcf_out.write(&record).expect("ERROR: could not write VCF record, mal-formatted!");
    
    }
}

/// This function is preparing the output vcf based on the
/// input ones. 
/// It does:
/// - generate a 1st header based on file A
/// - add all the FLAGS which we want finally as well as header and commands
/// - compares the fileA and fileB header to assure that the chromosomes are identical
/// Note: for the latter the name and the length is compared but not the assembly
/// Added a possibility to do "simple" comparison with less header entries and possibility
/// to add a FLAG which is next then populated instead.
///
///Unittest: FALSE
///
pub fn prepare_out_vcf(
    file_o: &str,
    file_a: &str,
    file_b: &str,
    infos: &VersionInfo,
    // version: &str, 
    // author: &str, 
    // command: &str,
    simple: &bool,
    simple_flag: &Option<&str>,
    info_fields: &Vec<&str>,
    //info_frmt: &Vec<&str>
) -> rust_htslib::bcf::Writer {

    //////////////////////////////
    //// preparing vcf header ////
    //////////////////////////////
    let vcf_a      = rust_htslib::bcf::Reader::from_path(file_a).expect("ERROR: could not open vcf file A !");
    let vcf_b      = rust_htslib::bcf::Reader::from_path(file_b).expect("ERROR: could not open vcf file B !");
    let mut vcf_header    = Header::from_template(vcf_a.header());
    
    // this are already the fields which we expect, if empty we will generate
    // essentially lots of whitespace in vcf but no harm done
    let header_sv_line      = r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">"#;
    let header_ci_line      = r#"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">"#;
    let header_source_line= format!("{}:{}",r#"##source=splitty"# ,infos.version);
    let header_mateid_line  = r#"##INFO=<ID=MATEID,Number=1,Type=String,Description="the ID of the 2nd part from the BND events">"#;
    let header_gf_line      = r#"##INFO=<ID=GENE_FUSION,Number=.,Type=String, Description="Suggested fusion partner genes">"#;    
    let header_topread_line = r#"##INFO=<ID=TOPREAD,Number=.,Type=String, Description="PID with strongest support">"#;    
    let header_author_line= format!("{}{}",r#"##source_author="# ,infos.author);
    let header_cmd_line   = format!("{}{}",r#"##command="# ,infos.command);
    vcf_header.push_record(header_source_line.as_bytes());
    vcf_header.push_record(header_author_line.as_bytes());
    vcf_header.push_record(header_cmd_line.as_bytes());
    vcf_header.push_record(header_topread_line.as_bytes());
    vcf_header.push_record(header_mateid_line.as_bytes());
    vcf_header.push_record(header_gf_line.as_bytes());
    vcf_header.push_record(header_sv_line.as_bytes());
    vcf_header.push_record(header_ci_line.as_bytes());

    // adding the provided INFO fields
    for (_n,field) in info_fields.iter().enumerate(){
        //let frmt = info_frmt[n];
        let header_inf_line = format!("##INFO=<ID={field},Number=1,Type=Integer,Description=\"none\">");
        vcf_header.push_record(header_inf_line.as_bytes());
    } 

    // the next lines are specific for splitty input
    if !simple {
        let header_fragid_line= r#"##INFO=<ID=FRAG_ID,Number=1,Type=String, Description="Name of the queried fragment">"#;    
        let header_dp_line    = r#"##INFO=<ID=SDP,Number=1,Type=Integer,Description="Support Read Depth (SP+SR) of segment containing breakend">"#;
        let header_sp_line    = r#"##INFO=<ID=SP,Number=1,Type=Integer,Description="Supporting split-pairs of gene-fusion">"#;
        let header_sr_line    = r#"##INFO=<ID=SR,Number=1,Type=Integer,Description="Supporting split-reads of gene-fusion">"#;
        let header_ratio_line = r#"##INFO=<ID=FR,Number=1,Type=Float,Description="SP to adjacent coverage ratio, lowest of both partners">"#;
        let header_conf_line  = r#"##INFO=<ID=EVI, Number=1, Type=String,Description="EVIDENCE: HIGH, MIDDLE, LOW">"#;
        let header_filter_line= r#"##FILTER=<ID=SPLITTY_FAILED, Description="events seen only in one of 2 strategies and there failing, too">"#;
        vcf_header.push_record(header_fragid_line.as_bytes());
        vcf_header.push_record(header_dp_line.as_bytes());
        vcf_header.push_record(header_sp_line.as_bytes());
        vcf_header.push_record(header_sr_line.as_bytes());
        vcf_header.push_record(header_ratio_line.as_bytes());
        vcf_header.push_record(header_conf_line.as_bytes());
        vcf_header.push_record(header_filter_line.as_bytes());
    }else{
        // just adding the flag for the support accordingly
        if simple_flag.is_some(){
            let tag = simple_flag.unwrap();
            let header_support = format!("##INFO=<ID={},Number=1,Type=Integer,Description=\"support from orthogonal data\">",tag);
            vcf_header.push_record(header_support.as_bytes());
        } else{
            panic!("ERROR: no flag was found, cant continue with vcf header contruction!");
        }
    }
    // now we need from both input vcf the chromosomes 
    // in case they differ, which they must not!
    // for reference in case : https://github.com/rust-bio/rust-htslib/issues/343
    let mut vcf_contigs : FxHashMap<String,linear_map::LinearMap<String,String>> = FxHashMap::default();

    for entries in vcf_a.header().header_records() {
        // here and afterwards key is not used but I cant add "_" because
        // it is defined inside the structure
        #[allow(unused_variables)]
        if let rust_htslib::bcf::header::HeaderRecord::Contig{key,mut values} = entries {
                values.remove("IDX").expect("ERROR:could not remove ID entry!");
                vcf_contigs.insert(values.get("ID").expect("ERROR: could not get ID of contig!").to_string(), values);
        };
    };

    for entries in vcf_b.header().header_records() {
        #[allow(unused_variables)]
        if let rust_htslib::bcf::header::HeaderRecord::Contig{key,mut values} = entries {
            values.remove("IDX").expect("ERROR:could not remove ID entry!");
            let id     = values.get("ID").expect("ERROR: could not get ID of contig!").to_string();
            let length = values.get("length").expect("ERROR: could not get length of contig!").to_string(); 
            if let Some(existant) = vcf_contigs.get(&id){
                if *existant.get("length")
                .expect("ERROR: could not get length of contig!")
                .to_string() != length 
                {
                    panic!("ERROR: length of similar contig IDs differed!");
                }else{
                    continue
                }
            }else{
                vcf_contigs
                    .insert(
                        values.get("ID")
                        .expect("ERROR: could not get ID of contig!")
                        .to_string(),
                        values
                    );
            }
        };
    };

    for (_,value) in vcf_contigs.iter() {
        let id       = value.get("ID").expect("ERROR: could not get ID of contig!").to_string();
        let length   = value.get("length").expect("ERROR: could not get ID of contig!").to_string();
        let assembly = value.get("assembly").expect("ERROR: could not get ID of contig!").to_string();
        vcf_header.push_record(format!("##contig=<ID={},length={},assembly={}>", id, length, assembly).as_bytes());
    }
    Writer::from_path(file_o,&vcf_header, true, Format::Vcf).expect("ERROR: could not write vcf to provided file")
}

#[cfg(test)]
mod tests {
    

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    //use bio::io::fasta::IndexedReader;  
    //use std::str::from_utf8;

    /////////////////////////////////////////
    ///    helper functions      ////////////
    /////////////////////////////////////////
    /// 
    #[test]
    fn get_cipos_some1(){
        let mut vcf = rust_htslib::bcf::Reader::from_path(&"test/vcf2fasta/sim_adipose.hg19_subset_p.ff.vcf").expect("Error opening file.");
        for res in vcf.records() {
            let record = res.unwrap();
            if record.pos() == 54074844 {
                let result : Option<Vec<u32>> = get_cipos(&record);
                assert_eq!(vec![0,0],result.unwrap());
            }
        }
    }
    #[test]
    fn get_cipos_some2(){
        let mut vcf = rust_htslib::bcf::Reader::from_path(&"test/vcf_cluster/single.vcf").expect("Error opening file.");
        for res in vcf.records() {
            let record = res.unwrap();
            if record.pos() == 23127258 {
                let result : Option<Vec<u32>> = get_cipos(&record);
                assert_eq!(vec![7432,7432],result.unwrap());
            }
        }
    }
    #[test]
    fn get_cipos_none(){
        let mut vcf = rust_htslib::bcf::Reader::from_path(&"test/data/cluster1.vcf").expect("Error opening file.");
        for res in vcf.records() {
            let record = res.unwrap();
            if record.pos() == 321680 {
                let result : Option<Vec<u32>> = get_cipos(&record);
                assert_eq!(None,result);
            }
        }
    }
    #[test]
    fn paired_qc_module_1(){
        // checking if record is suitable
        // with no black-list and not accepting multi-mapping reads
        let mut paired_bam = bam::Reader::from_path("test/test_paired/single_paired_1.bam").expect("ERROR: Could not open BAM file ");
        //let record     = bam::Record::new();
        let mut test  : Vec<bool> = Vec::new();
        // now we have only passing reads and final result must 
        // be true for passing 
        for record in paired_bam.records() {
            let result = record.unwrap();
            test.push(read_qc_module_paired(&result, &None,&false));

        };
        assert_eq!(vec![true,true,true], test);
    }
    #[test]
    fn paired_qc_module_2(){
        // checking if record is suitable
        // with a black-list and not accepting multi-mapping reads
        // expecting to return false as the blacklist contains element
        let mut b_list: FxHashMap<String,i32> = FxHashMap::default();
        b_list.insert(String::from("EU432099.1_56_624_1:0:0_3:0:0_4"),1);
        let black_list = Some(b_list);
        let mut paired_bam = bam::Reader::from_path("test/test_paired/single_paired_1.bam").expect("ERROR: Could not open BAM file ");
        //let mut record     = bam::Record::new();
        let mut test  : Vec<bool> = Vec::new();
        // now all reads will be on the black-list
        // and test must switch to false
        for record in paired_bam.records() {
            let result = record.unwrap();
            test.push(read_qc_module_paired(&result, &black_list,&false));

        };
        assert_eq!(vec![false,false,false], test);
    }
    #[test]
    fn paired_qc_module_3(){
        // checking if record is suitable
        // with a black-list and not accepting multi-mapping reads
        // expecting to return true as the blacklist contains elements other than the ones in BAM
        let mut b_list: FxHashMap<String,i32> = FxHashMap::default();
        b_list.insert(String::from("EU432099.1_56_624_1:0:0_3:0:0_5"),1);
        let black_list = Some(b_list);
        let mut paired_bam = bam::Reader::from_path("test/test_paired/single_paired_1.bam").expect("ERROR: Could not open BAM file ");
        //let mut record     = bam::Record::new();
        let mut test  : Vec<bool> = Vec::new();
        // now all reads will be on the black-list
        // and test must switch to false
        for record in paired_bam.records() {
            let result = record.unwrap();
            test.push(read_qc_module_paired(&result, &black_list,&false));

        };
        assert_eq!(vec![true,true,true], test);
    }
    #[test]
    fn paired_qc_module_4(){
        // checking if record is suitable
        // this one are all multi-mappers and should return false
        let mut paired_bam = bam::Reader::from_path("test/test_paired/single_paired_mm.bam").expect("ERROR: Could not open BAM file ");
        //let mut record     = bam::Record::new();
        let mut test  : Vec<bool> = Vec::new();
        // now all reads will be on the black-list
        // and test must switch to false
        for record in paired_bam.records() {
            let result = record.unwrap();
            test.push(read_qc_module_paired(&result, &None,&false));

        };
        assert_eq!(vec![false,false,false], test);
    }
    #[test]
    fn paired_qc_module_5(){
        // checking if record is suitable
        // this one are all multi-mappers by accepting multi-mappers should succeed
        let mut paired_bam = bam::Reader::from_path("test/test_paired/single_paired_mm.bam").expect("ERROR: Could not open BAM file ");
        //let mut record     = bam::Record::new();
        let mut test  : Vec<bool> = Vec::new();
        // now all reads will be on the black-list
        // and test must switch to false
        for record in paired_bam.records() {
            let result = record.unwrap();
            test.push(read_qc_module_paired(&result, &None,&true));

        };
        assert_eq!(vec![true,true,true], test);
    }
    #[test]
    fn paired_qc_module_6(){
        // checking if record is suitable
        // with no black-list and not accepting multi-mapping reads
        // having a secondary entry
        let mut paired_bam = bam::Reader::from_path("test/test_paired/single_paired_sec.bam").expect("ERROR: Could not open BAM file ");
        //let mut record     = bam::Record::new();
        let mut test  : Vec<bool> = Vec::new();
        // now we have only passing reads and final result must 
        // be true for passing 
        for record in paired_bam.records() {
            let result = record.unwrap();
            test.push(read_qc_module_paired(&result, &None,&false));

        };
        assert_eq!(vec![false,true,true], test);
    }
    #[test]
    fn read_qc_module_1(){
        // checking if record is suitable
        // with no black-list and not accepting multi-mapping reads
        let mut fragment_bam = bam::Reader::from_path("test/test_fragment/sim_adipose.hg19_subset_p.ff.bam").expect("ERROR: Could not open BAM file ");
        //let mut record     = bam::Record::new();
        let mut test  : Vec<bool> = Vec::new();
        // now we have only passing reads and final result must 
        // be true for passing 
        for record in fragment_bam.records() {
            let result = record.unwrap();
            test.push(read_qc_module(&result, &None,&false));

        };
        assert_eq!(vec![true,true], test);
    }
    #[test]
    fn read_qc_module_2(){
        // checking if record is suitable
        // with a black-list and not accepting multi-mapping reads
        // expecting to return false as the blacklist contains element
        let mut b_list: FxHashMap<String,i32> = FxHashMap::default();
        b_list.insert(String::from("DKK1|ENSG00000107984.5--C1orf220|ENSG00000213057.4"),1);
        let black_list = Some(b_list);
        let mut fragment_bam = bam::Reader::from_path("test/test_fragment/sim_adipose.hg19_subset_p.ff.bam").expect("ERROR: Could not open BAM file ");
        //let mut record     = bam::Record::new();
        let mut test  : Vec<bool> = Vec::new();
        // now all reads will be on the black-list
        // and test must switch to false
        for record in fragment_bam.records() {
            let result = record.unwrap();
            test.push(read_qc_module(&result, &black_list,&false));

        };
        assert_eq!(vec![false,false], test);
    }
    #[test]
    fn read_qc_module_3(){
        // checking if record is suitable
        // with a black-list and not accepting multi-mapping reads
        // expecting to return true as the blacklist contains elements other than the ones in BAM
        let mut b_list: FxHashMap<String,i32> = FxHashMap::default();
        b_list.insert(String::from("DKK1|ENSG00000107984.5--C1orf220|ENSG00000213057.5"),1);
        let black_list = Some(b_list);
        let mut fragment_bam = bam::Reader::from_path("test/test_fragment/sim_adipose.hg19_subset_p.ff.bam").expect("ERROR: Could not open BAM file ");
        //let mut record     = bam::Record::new();
        let mut test  : Vec<bool> = Vec::new();
        // now all reads will be on the black-list
        // and test must switch to false
        for record in fragment_bam.records() {
            let result = record.unwrap();
            test.push(read_qc_module(&result, &black_list,&false));

        };
        assert_eq!(vec![true,true], test);
    }

    #[test]
    fn splitp_2_splitr_addition(){
        // note : I counted the reads on a which are SP and there are indeed 63 and 68 on the other side
        let mut paired_bam = bam::IndexedReader::from_path("test/test_paired/reads_aligned.bam").expect("ERROR: Could not open BAM file ");
        let mut input : FxHashMap<String,Vec<PosBasedFusEvidence>> = FxHashMap::default();
        let info = PosBasedFusEvidence {
            // name of chr primary match
            chrom_pr: String::from("EU432099.1:0-530"),
            // name of chr SA
            chrom_sa: String::from("EU432099.1:531"),
            // position of break on target from primary match
            // .0-based
            target_break_pr: 530_u64,
            // position of break on target from SA
            // .0-based
            target_break_sa: 1_u64,
            // strand of primary match
            strand_pr: StrandDirection::Rev,
            // strand of SA
            strand_sa: StrandDirection::Rev,
            // placement of primary match in fusion can be 5 or 3
            orient_pr_upstream: false,
            // placement of SA in fusion can be 5 or 3
            orient_sa_upstream: true,
            // collection of target_break_pr values
            precision_pr: vec![530,530],
            // collection of target_break_sa values
            precision_sa: vec![1,1],
            // fusion position of primary alignment
            // .0-based
            aln_start_pr: vec![475,475],
            // fusion position of SA
            // .0-based
            aln_start_sa: vec![93,93],
            // Split-read support is an increasing number of seen
            // evidence supporting that position
            support_sr: 33_u32,
            // Split-pair support is an increasing number of seen
            // evidence supporting that position
            support_sp: 0_u32,
            // The ratio indicates the number of reads spanning a fusion
            // compared to the median coverage in our inquired region. 
            // This ratio is reported for both sites and the lowest is given in
            // this field.
            // Especially for positions where we have ambiguity for 
            // mapping or lots of noise that ratio will be lower and help to identify problems.
            sp_ratio: 0_f32,
            // provides the number of uncertain positions
            // around the fp points
            cipos : Some(0_u32)
        };
        input.insert(String::from("EU432099.1:0-530EU432099.1:531"),vec![info]);
        let range : u32 = 500;
        let result: FxHashMap<String,Vec<PosBasedFusEvidence>>;
        result = add_sp_to_sr_fusions(&mut paired_bam,input,range);
        let test_info = PosBasedFusEvidence {
            // name of chr primary match
            chrom_pr: String::from("EU432099.1:0-530"),
            // name of chr SA
            chrom_sa: String::from("EU432099.1:531"),
            // position of break on target from primary match
            // .0-based
            target_break_pr: 530_u64,
            // position of break on target from SA
            // .0-based
            target_break_sa: 1_u64,
            // strand of primary match
            strand_pr: StrandDirection::Rev,
            // strand of SA
            strand_sa: StrandDirection::Rev,
            // placement of primary match in fusion can be 5 or 3
            orient_pr_upstream: false,
            // placement of SA in fusion can be 5 or 3
            orient_sa_upstream: true,
            // collection of target_break_pr values
            precision_pr: vec![530,530],
            // collection of target_break_sa values
            precision_sa: vec![1,1],
            // fusion position of primary alignment
            // .0-based
            aln_start_pr: vec![475,475],
            // fusion position of SA
            // .0-based
            aln_start_sa: vec![93,93],
            // Split-read support is an increasing number of seen
            // evidence supporting that position
            support_sr: 33_u32,
            // Split-pair support is an increasing number of seen
            // evidence supporting that position
            support_sp: 60_u32,
            // The ratio indicates the number of reads spanning a fusion
            // compared to the median coverage in our inquired region. 
            // This ratio is reported for both sites and the lowest is given in
            // this field.
            // Especially for positions where we have ambiguity for 
            // mapping or lots of noise that ratio will be lower and help to identify problems.
            sp_ratio: 1.7647059_f32,
            // provides the number of uncertain positions
            // around the fp points
            cipos : Some(0_u32)
        };

        let mut test : FxHashMap<String,Vec<PosBasedFusEvidence>>  = FxHashMap::default();
        test.insert(String::from("EU432099.1:0-530EU432099.1:531"),vec![test_info]);
        assert_eq!(result,test);
    }

    #[test]
    fn splitp_2_splitr_addition_tooshort(){
        // so this one is literally so short that we have no
        // read from SP that will support it
        let mut paired_bam = bam::IndexedReader::from_path("test/test_paired/reads_aligned.bam").expect("ERROR: Could not open BAM file ");
        let mut input : FxHashMap<String,Vec<PosBasedFusEvidence>> = FxHashMap::default();
        let info = PosBasedFusEvidence {
            // name of chr primary match
            chrom_pr: String::from("EU432099.1:0-530"),
            // name of chr SA
            chrom_sa: String::from("EU432099.1:531"),
            // position of break on target from primary match
            // .0-based
            target_break_pr: 530_u64,
            // position of break on target from SA
            // .0-based
            target_break_sa: 1_u64,
            // strand of primary match
            strand_pr: StrandDirection::Rev,
            // strand of SA
            strand_sa: StrandDirection::Rev,
            // placement of primary match in fusion can be 5 or 3
            orient_pr_upstream: false,
            // placement of SA in fusion can be 5 or 3
            orient_sa_upstream:true,
            // collection of target_break_pr values
            precision_pr: vec![530,530],
            // collection of target_break_sa values
            precision_sa: vec![1,1],
            // fusion position of primary alignment
            // .0-based
            aln_start_pr: vec![475,475],
            // fusion position of SA
            // .0-based
            aln_start_sa: vec![93,93],
            // Split-read support is an increasing number of seen
            // evidence supporting that position
            support_sr: 33_u32,
            // Split-pair support is an increasing number of seen
            // evidence supporting that position
            support_sp: 0_u32,
            // The ratio indicates the number of reads spanning a fusion
            // compared to the median coverage in our inquired region. 
            // This ratio is reported for both sites and the lowest is given in
            // this field.
            // Especially for positions where we have ambiguity for 
            // mapping or lots of noise that ratio will be lower and help to identify problems.
            sp_ratio: 0_f32,
            // provides the number of uncertain positions
            // around the fp points
            cipos : Some(0_u32)
        };
        input.insert(String::from("EU432099.1:0-530EU432099.1:531"),vec![info]);
        let range : u32 = 50;
        let result: FxHashMap<String,Vec<PosBasedFusEvidence>> ;
        result = add_sp_to_sr_fusions(&mut paired_bam,input,range);
        let test_info = PosBasedFusEvidence {
            // name of chr primary match
            chrom_pr: String::from("EU432099.1:0-530"),
            // name of chr SA
            chrom_sa: String::from("EU432099.1:531"),
            // position of break on target from primary match
            // .0-based
            target_break_pr: 530_u64,
            // position of break on target from SA
            // .0-based
            target_break_sa: 1_u64,
            // strand of primary match
            strand_pr: StrandDirection::Rev,
            // strand of SA
            strand_sa: StrandDirection::Rev,
            // placement of primary match in fusion can be 5 or 3
            orient_pr_upstream: false,
            // placement of SA in fusion can be 5 or 3
            orient_sa_upstream: true,
            // collection of target_break_pr values
            precision_pr: vec![530,530],
            // collection of target_break_sa values
            precision_sa: vec![1,1],
            // fusion position of primary alignment
            // .0-based
            aln_start_pr: vec![475,475],
            // fusion position of SA
            // .0-based
            aln_start_sa: vec![93,93],
            // Split-read support is an increasing number of seen
            // evidence supporting that position
            support_sr: 33_u32,
            // Split-pair support is an increasing number of seen
            // evidence supporting that position
            support_sp: 0_u32,
            // The ratio indicates the number of reads spanning a fusion
            // compared to the median coverage in our inquired region. 
            // This ratio is reported for both sites and the lowest is given in
            // this field.
            // Especially for positions where we have ambiguity for 
            // mapping or lots of noise that ratio will be lower and help to identify problems.
            sp_ratio: 0_f32,
            // provides the number of uncertain positions
            // around the fp points
            cipos : Some(0_u32)
        };

        let mut test : FxHashMap<String,Vec<PosBasedFusEvidence>>  = FxHashMap::default();
        test.insert(String::from("EU432099.1:0-530EU432099.1:531"),vec![test_info]);
        assert_eq!(result,test);
    }

    #[test]
    fn alignments_from_bam_paired1(){
        // this is a single paired read pair which is both, a SR and SP
        let paired_bam = bam::Reader::from_path("test/test_paired/single_paired_1.bam").expect("ERROR: Could not open BAM file ");
        let (results,_ignore) = analyze_alignments_from_bam(paired_bam,None,0.95,true,false,None);
        let mut truth : FxHashMap<String, Vec<FullFusionEvidence>>  = FxHashMap::default();
        //: 
        let info = vec![
            FullFusionEvidence { 
                cigar_pr: String::from("56M94H"), 
                cigar_sa: String::from("60S90M"), 
                chrom_pr: String::from("EU432099.1:0-530"), 
                chrom_sa: String::from("EU432099.1:531"), 
                query_name: String::from("EU432099.1_56_624_1:0:0_3:0:0_4"), 
                query_length: 150, 
                target_break_pr: 529, 
                target_break_sa: 4, 
                query_break_pr: 94, 
                query_break_sa: 89, 
                strand_pr: StrandDirection::Rev, 
                strand_sa: StrandDirection::Rev, 
                orient_pr_upstream: false, 
                orient_sa_upstream: true, 
                aln_start_pr: 474, 
                aln_start_sa: 4, 
                precision: None, 
                support: 0 
            }, 
        FullFusionEvidence { 
            cigar_pr: String::from("60S90M"), 
            cigar_sa: String::from("56M94S"), 
            chrom_pr: String::from("EU432099.1:531"), 
            chrom_sa: String::from("EU432099.1:0-530"), 
            query_name: String::from("EU432099.1_56_624_1:0:0_3:0:0_4"), 
            query_length: 150, 
            target_break_pr: 4, 
            target_break_sa: 529, 
            query_break_pr: 89, 
            query_break_sa: 94, 
            strand_pr: StrandDirection::Rev, 
            strand_sa: StrandDirection::Rev, 
            orient_pr_upstream: true, 
            orient_sa_upstream: false, 
            aln_start_pr: 4, 
            aln_start_sa: 474, 
            precision: None, 
            support: 0 
        }];
        truth.insert(String::from("EU432099.1_56_624_1:0:0_3:0:0_4"),info);
        assert_eq!(truth,results);
    }
    // fn alignments_from_bam_paired2(){
    //     // this is a single paired read pair which is both, a SR and SP
    //     // it is at position 0 which makes it more tricky
    //     let mut paired_bam = bam::Reader::from_path("test/test_paired/single_paired_2.bam").expect("ERROR: Could not open BAM file ");
    //     let (results,ignore) = analyze_alignments_from_bam(paired_bam,None,0.95,true,false);
    //     let mut test : FxHashMap<String, Vec<FullFusionEvidence>>  = FxHashMap::default();
    //     //: 
    //     let info = vec![
    //         FullFusionEvidence { 
    //             cigar_pr: String::from("56M94H"), 
    //             cigar_sa: String::from("60S90M"), 
    //             chrom_pr: String::from("EU432099.1:0-530"), 
    //             chrom_sa: String::from("EU432099.1:531"), 
    //             query_name: String::from("EU432099.1_56_624_1:0:0_3:0:0_4"), 
    //             query_length: 150, 
    //             target_break_pr: 529, 
    //             target_break_sa: 4, 
    //             query_break_pr: 94, 
    //             query_break_sa: 89, 
    //             strand_pr: StrandDirection::Rev, 
    //             strand_sa: StrandDirection::Rev, 
    //             orient_pr: 3, 
    //             orient_sa: 5, 
    //             aln_start_pr: 474, 
    //             aln_start_sa: 92, 
    //             precision: None, 
    //             support: 0 
    //         }, 
    //     FullFusionEvidence { 
    //         cigar_pr: String::from("60S90M"), 
    //         cigar_sa: String::from("56M94S"), 
    //         chrom_pr: String::from("EU432099.1:531"), 
    //         chrom_sa: String::from("EU432099.1:0-530"), 
    //         query_name: String::from("EU432099.1_56_624_1:0:0_3:0:0_4"), 
    //         query_length: 150, 
    //         target_break_pr: 4, 
    //         target_break_sa: 528, 
    //         query_break_pr: 89, 
    //         query_break_sa: 94, 
    //         strand_pr: StrandDirection::Rev, 
    //         strand_sa: StrandDirection::Rev, 
    //         orient_pr: 5, 
    //         orient_sa: 3, 
    //         aln_start_pr: 93, 
    //         aln_start_sa: 474, 
    //         precision: None, 
    //         support: 0 
    //     }];
    //     test.insert(String::from("EU432099.1_56_624_1:0:0_3:0:0_4"),info);
    //     assert_eq!(test,results);
    // }
}

