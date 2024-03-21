use rustc_hash::FxHashMap;
use std::str;
use regex::Regex;
use std::fs::File;
use crossbeam::channel::{unbounded};
use rayon;
use log::{debug};

extern crate noodles;
use noodles::bam as noodle_bam;
use noodles::sam as noodle_sam;
use indexmap::{IndexMap};


use crate::lib::common::{*};


//use noodles_sam::record::Flags;
/// this checks if our alignments are indeed what we 
/// expect and removes return false for entries which
/// are not according to specs.
/// It is not clear if one should panic on non-paired reads
/// at this stage and this might change in the future.
/// Noodles specific
///
/// Unittest: FALSE
///
pub fn read_qc_module_paired_noodles(
    entry: &noodle_bam::record::Record,
    black_list: &Option<FxHashMap<String, i32>>
) -> bool {
    let name = str::from_utf8(
        entry.read_name().expect("ERROR: could not extract read name!")
        .to_bytes()).expect("ERROR: could not convert ref name to string!").to_string();
    let flags = entry.flags();
    if !flags.is_paired() {
        eprintln!("ERROR: Alignments are not paired-reads!");
    }
    let on_black_list = match black_list {
        Some(x) => {
            x.contains_key(&name)
        },
        None => false
    };
    !(!flags.is_paired() || flags.is_secondary() || flags.is_unmapped() || on_black_list)
}


/// very basic function which
/// takes a bam record that contains a
/// SA entry and verifies if the SA elements
/// stem from same target ID or not
/// Needs the source, if CHESS or GENCODE
///
///Unittest: FALSE
///
pub fn analyze_basic_sa_for_fusion_noodles(
    seq_local: &IndexMap<String,noodle_sam::header::ReferenceSequence> ,
    sa_entry: &str,
    record: &noodle_bam::record::Record,
    source: &str, sim: &f32
) -> bool {

    let name   = str::from_utf8(
        record.read_name().expect("ERROR: could not extract read name!")
        .to_bytes()).expect("ERROR: could not convert ref name to string!").to_string();

    // first identify the source accordingly
    let transcriptome_source : GtfSource = match source {
        "CHESS" => GtfSource::CHESS,
        "GENCODE" => GtfSource::GENCODE,
        _ => GtfSource::UNKNOWN,
    };

    // this method might be a bit of over-kill here.
    // It takes a SA entry and splits it in its sub-entries. 
    // Additionally it evaluates the entry and returns information 
    // about the position.
    // We are cheating here and providing a length of 1 to avoid
    // having to pass the real proper length which would need more
    // time-consuming calculations
    // Similarly, I am not even providing the possibility
    // to alter a ratio if both sides are clipped - not the goal
    // here - just put to None instead
    let supp_infos = eval_supp_align(sa_entry,&name,None,None);
    debug!("Noodles supp_infos {:?}",&supp_infos);
    // depending on the source we need to
    // potentially collapse differently transcript to gene
    // level which we need here to be sure it is a different
    // gene and not an isoform
    // so we check each entry and simultaneously remove the ones 
    // which have not enough similarity
    // check if chromosomes of both are different:
    let mut chr1  = record
        .reference_sequence_id()
        .map(|id| {
            seq_local.get_index(i32::from(id) as usize)
            .map(|(_,rs)| rs.name())
            .expect("ERROR: missing reference sequence!")
        }).unwrap_or("*");
    


    let mut chr2 : Vec<String> = Vec::new();
    let mut tmp_chr2 : String ;
    match transcriptome_source {
        GtfSource::CHESS   => {
            // CHESS we simply remove the .[0-9] at the end
            let pattern = Regex::new(r#"^(CHS\.[0-9]+)\.[0-9]+"#).unwrap();
            chr1 = match pattern.captures(chr1){
                Some(x) => x.get(1).unwrap().as_str(),
                _  => panic!("ERROR: Gene name match not captured!"),
            };
            for entry in supp_infos {
                println!("Similarity entry {} , sim-cutoff{}",&entry.similarity,&sim);
                if &entry.similarity < sim {
                    continue
                }
                if pattern.is_match(&entry.chrom) {
                    tmp_chr2 = match pattern.captures(&entry.chrom){
                        Some(x) => x.get(1).unwrap().as_str().to_string(),
                        _  => panic!("ERROR: Gene name match found but not captured!"),
                    };
                    chr2.push(tmp_chr2);
                }
            }
        },
        // Note:
        // For gencode I realized this wont be as easy and we need a translation
        // table :/
        GtfSource::GENCODE => panic!("ERROR: Sorry, Gencode not yet supported !"),
        GtfSource::UNKNOWN => panic!("ERROR: your target is none of the compatible source !"),
    }
    // now we have a curated vector and the question is
    // whether the entry is similar to the primary match
    // or whether it is different
    // The first match is normally the one with the highest score
    // So we go by entry and if anything is identical to 1st chr
    // then we return immediately false
    // otherwise we return true if one is different
    debug!("Noodles identified chr1 {} chr2 {:?}",&chr1,&chr2);
    let mut diff_chrom = false;
    for entry in chr2 {
        if entry == chr1 {
            break
        }else{
            diff_chrom = true;
        }
    }
    diff_chrom

}


/// This function analyzes paired-read alignments and checks if they 
/// provide evidence for a fusion derived event. Importantly it expects
/// it to be derived from a transcriptome based alignment (might works otherwise but not tested).
/// Currently only CHESS is supported as it needs the information of transcript-ID --> gene-ID relationship.
/// This is simple in case of CHESS as it is a logical substring of the ID whereas for 
/// Gencode this requires a real translation table.
/// 
/// Below are 2 examples, one with a SA alignment indicating a fusion read
/// and a second one with a SP from different regions indicating similarly a fusion event
///
///Unittest: TRUE
/// 
/// ```rust
/// use crate::genefusion::lib::noodles_based::analyze_paired_alignments_from_trans_bam_noodles;
/// pretty_env_logger::init();
/// let identity = 98.0_f32;
/// let threads = 2;
/// let source = String::from("CHESS");
/// let bam_file = "test/test_idFusionReads/Sim1-R0P0_1_readsOnTCata_miniSA.bam";
/// assert_eq!(
///     analyze_paired_alignments_from_trans_bam_noodles(
///         bam_file,
///         identity,
///         &source,
///         threads
///     ),
///     vec!["CCDC27|ENSG00000162592.4--CHST11|ENSG00000171310.6:2450525_1_29196_5022_133"]
/// );
/// ```
/// 
/// ```rust
/// use crate::genefusion::lib::noodles_based::analyze_paired_alignments_from_trans_bam_noodles;
/// pretty_env_logger::init();
/// let identity = 98.0_f32;
/// let threads = 2;
/// let source = String::from("CHESS");
/// let bam_file = "test/test_idFusionReads/Sim1-R0P0_1_readsOnTCata_miniSP.bam";
/// assert_eq!(
///     analyze_paired_alignments_from_trans_bam_noodles(
///         bam_file,
///         identity,
///         &source,
///         threads
///     ),
///     vec!["ABCG1|ENSG00000160179.14--GPSM3|ENSG00000213654.5:5202722_0_754_566_262"]
/// );
/// ```
pub fn analyze_paired_alignments_from_trans_bam_noodles(
    bam_file: &str,
    similarity: f32,
    source: &str,
    threads: usize,
) -> Vec<String> {

    let transcriptome_source : GtfSource = match source {
        "CHESS" => GtfSource::CHESS,
        "GENCODE" => GtfSource::GENCODE,
        _ => GtfSource::UNKNOWN,
    };

    let mut reader = File::open(bam_file)
        .map(noodle_bam::Reader::new)
        .expect("ERROR: could not open and read bam file!");
    
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    eprintln!("Current thread-pool number: {:?}",rayon::current_num_threads());
    
    // for noodles we always need first to initialize the header
    // end then in the next step the sequences otherwise we get
    // strange errors - this is though intentional
    match reader.read_header(){
        Ok(_) => (),
        Err(_) => panic!("ERROR: could not read header!"),
    };

    let sequences : IndexMap<String, noodle_sam::header::ReferenceSequence> = match reader.read_reference_sequences() {
        Ok(n) =>  n,
        Err(_) => panic!("ERROR: could not get ref sequence!"),
    };
    
    let (snd, rxv) = unbounded();
    rayon::scope( |child| {
        for record in reader.records(){
            let snd_local = snd.clone();
            let seq_local = &sequences;
            let transcriptome_source_local = &transcriptome_source;
            child.spawn( move |_| {
                let read = record.expect("ERROR: could not extract BAM record!");
                let blacklist : Option<FxHashMap<String, i32>> = None;
                let tmp_result   = str::from_utf8(
                    read.read_name().expect("ERROR: could not extract read name!")
                    .to_bytes()).expect("ERROR: could not convert ref name to string!").to_string();
                let flags = read.flags();

                // here comes now the QC part.
                // As we have currently no blacklist, we provide an empty one
                let qc_result   : bool   = read_qc_module_paired_noodles(&read, &blacklist);
                let mut sa_fusion = false;

                // if the qc fails we skip this part
                // and as well the next then
                if qc_result {
                    // if we have a SA entry it is already becoming
                    // interesting otherwise we need to check read1 vs read2
                    // this is though tricky in noodles
                    // find here some hints: https://github.com/zaeleus/noodles/issues/15
                    let data : FxHashMap<_,_>    = read
                        .data()
                        .fields()
                        .map(|result| 
                            result.map(
                                |field| (field.tag().clone(), field )
                            )
                        )
                        .collect::<Result<_,_>>()
                        .expect("ERROR: could not access data fields of BAM entry!");
                    if let Some(sa_entry) = data.get(&noodle_sam::record::data::field::tag::Tag::OtherAlignments){
                        let value = sa_entry.value().as_str().expect("ERROR: could not convert SA field correctly!");
                        sa_fusion = analyze_basic_sa_for_fusion_noodles(seq_local,value,&read, source ,&similarity);
                        debug!("Noodles sa values {} ",&sa_fusion);
                    };
                    //  now we want to check if the mate is potentially mapped 
                    // onto another chromosome
                    // this takes some processing and we only do it if the
                    // above stage already "failed" and if the second one is 
                    // actually mapped - otherwise does not make sense
                    let mut diff_chr = false;
                    if qc_result && (!flags.is_mate_unmapped() & !sa_fusion ) {
                        // check if chromosomes of both are different:
                        // this is pretty crazy mapping of mapping which I was
                        // not familiar with but got help : https://github.com/zaeleus/noodles/issues/55
                        let mut chr1  = read
                            .reference_sequence_id()
                            .map(|id| {
                                seq_local.get_index(i32::from(id) as usize)
                                .map(|(_,rs)| rs.name())
                                .expect("ERROR: missing reference sequence!")
                            }).unwrap_or("*");
                        let mut chr2  = read
                            .mate_reference_sequence_id()
                            .map(|id| {
                                seq_local.get_index(i32::from(id) as usize)
                                .map(|(_,rs)| rs.name())
                                .expect("ERROR: missing reference sequence!")
                            }).unwrap_or("*");
                        if chr1 != "*" && chr2 != "*"{                       
                            match transcriptome_source_local {
                            GtfSource::CHESS   => {
                                    // CHESS we simply remove the .[0-9] at the end
                                    let pattern = Regex::new(r#"^(CHS\.[0-9]+)\.[0-9]+"#).unwrap();
                                    chr2 = match pattern.captures(chr2){
                                        Some(x) => x.get(1).unwrap().as_str(),
                                        _  => {
                                            eprintln!("ERROR: chr1 {} or chr2 {} were not according to expected format!",&chr1,&chr2);
                                            panic!("ERROR: Gene name match not captured!")
                                        },
                                    };
                                },
                            _ => panic!("ERROR no CHESS trasnscripts provided!")
                            }         
                            
                            match transcriptome_source_local {
                                GtfSource::CHESS   => {
                                    // CHESS we simply remove the .[0-9] at the end
                                    let pattern = Regex::new(r#"^(CHS\.[0-9]+)\.[0-9]+"#).unwrap();
                                    chr1 = match pattern.captures(chr1){
                                    Some(x) => x.get(1).unwrap().as_str(),
                                    
                                        _  => {
                                            eprintln!("ERROR: chr1 {} or chr2 {} were not according to expected format!",&chr1,&chr2);
                                            panic!("ERROR: Gene name match not captured!")
                                        },
                                    };
                                },
                                _ => panic!("ERROR no CHESS trasnscripts provided!")
                            }               
                            if chr1 != chr2 {
                                diff_chr = true;
                            }
                        }
                    }
                    if qc_result && ( sa_fusion | diff_chr ) {
                        snd_local.send(tmp_result).expect("ERROR: thread could not communicate result!");
                    }
                }
            });
        }
    });
    drop(snd);
	let mut thread_results : Vec<std::string::String> = rxv.iter().collect();
    thread_results.sort_unstable();
    thread_results.dedup();
    thread_results
}

