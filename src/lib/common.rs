use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader};
use std::str;
use std::path::Path;
use regex::Regex;
use std::fs::File;
use std::cmp::Ordering;
use std::convert::TryFrom;
use integer_sqrt::IntegerSquareRoot;
use statistical::standard_deviation;
use std::io;
use std::error::Error;
use chrono::{DateTime, Local};
use bio::data_structures::interval_tree::{*};
use std::str::from_utf8;
use bio::alphabets::dna::revcomp;
use bio::alphabets::dna::complement;
use std::io::Read as IoRead;
use log::debug;
use std::io::Write;
use std::io::stdout;
extern crate noodles;


/// # The currently accepted GTF formats.
/// Since gtf/gff are really badly defined as a format,
/// it is almost impossible to accommodate all possibilities.
/// Therefore we have here the currently accepted standard format 
/// which we are able to parse. 
#[derive(Debug, Clone, Copy)]
pub enum GtfSource {
    /// Which is similar to this:
    /// 
    /// ```text
    /// chr1    BestRefSeq      transcript      69091   70008   .       +       .       ID=CHS.6.1;geneID=CHS.6;gene_name=OR4F5
    /// chr1    Gnomon  transcript      91169   120798  .       -       .       ID=CHS.7.1;geneID=CHS.7;gene_name=LOC100996442
    /// ```
    /// 
    /// For the attribute section we expect as separator the "=" and no 
    /// quotes
    CHESS,
    /// Which is similar to this:
    /// 
    /// ```text
    /// chr1    HAVANA  transcript      129081  133566  .       -       .       gene_id "ENSG00000238009.2"; transcript_id "ENST00000453576.2"; gene_type "lincRNA"; gene_status "NOVEL"; gene_name "RP11-34P13.7"; transcript_type "lincRNA"; transcript_status "KNOWN"; transcript_name "RP11-34P13.7-004"; level 2; havana_gene "OTTHUMG00000001096.2"; havana_transcript "OTTHUMT00000003689.1";
    /// chr1    HAVANA  transcript      131025  134836  .       +       .       gene_id "ENSG00000233750.3"; transcript_id "ENST00000442987.3"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "CICP27"; transcript_type "processed_pseudogene"; transcript_status "KNOWN"; transcript_name "CICP27-001"; level 1; ont "PGO:0000004"; tag "pseudo_consens"; havana_gene "OTTHUMG00000001257.3"; havana_transcript "OTTHUMT00000003691.3";
    /// ```
    /// 
    /// Here we have as separator an empty space (which is pretty dangerous) and the value is in brackets
    GENCODE,
    /// This is accepted if we parse the first time all entries
    /// and have no prior knowledge what system is underlying.
    UNKNOWN,
}

#[derive(Debug)]
/// The typcial GTF structure in a tab-separated fashion.
/// Except that we do not encode the strand into +/- in a string but 
/// use instead -1 and 1
pub struct Gtf {
    /// chromosome
    pub chromosome: String,
    /// source of annotation, e.g. ENCODE
    pub source: String,
    /// type of feature, e.g. gene
    pub feature: String,
    /// start of entry
    pub start: u64,
    /// end of entry
    pub end: u64,
    /// score which is most of time "."
    pub score: String,
    /// strand is normally "+/-" but we encode in -1/1
    pub strand: StrandDirection,
    /// frame 1/2/3 mosty "."
    pub frame: String,
    /// the wild-card which is different in almost every gtf
    pub attribute: FxHashMap<String,String>,
}

#[derive(Debug,Default,PartialEq,Eq)]
/// A structure to organize information for the break-point identification. 
/// All coordinates are 0 based
pub struct BreakPoints {
    /// break-point on the read/cDNA
    /// 0-based
    pub break_query: u32,
    /// break-point on the target/chromosome
    /// 0-based
    pub break_target: u64,
    /// alignment start
    /// 0-based
    pub aln_start: u64,
    /// alignment end
    /// 0-based
    pub aln_end: u64,
    /// strand
    pub strand: StrandDirection,
}

#[derive(Debug,Default,PartialEq)]
/// a bunch of information which we
/// calculate based on the information from a cigar
/// this is the base for all SV analysis then.
/// 0-based
pub struct CigarInfos {
    /// position of match on genome/contig
    /// 0-based
    pub fp_target: u64,
    /// number of clipped position on upstream
    /// 0-based
    pub upstream: u32,
    /// number of clipped position on downstream
    /// .0-based
    pub downstream: u32,
    /// gap compressed similarity score
    pub similarity: f32,
    /// length of the sequence (includes INs but no DELs)
    pub seq_length: u32,
    /// match length only M,X and = (no splice,clip,INDELS)
    pub match_length: u32,
    /// length with splicing of match (includes INDELS and N, and CLIP)
    pub match_length_s: u32,
    /// position of match on query
    /// .0-based
    pub fp_query: u32,
    /// start of alignment
    /// .0-based
    pub aln_start: u64,
    /// end of alignment
    /// .0-based
    pub aln_end: u64,
    /// strand of match
    pub match_strand: StrandDirection,
}


#[derive(Debug)]
/// mostly identical to CigarInfos
/// but keeps as well the contig/chromosome name
/// which is an additional information which needs
/// to be correctly linked
/// .0-based
pub struct AlignInfos {
    /// read/CDNA id
    pub query: String,
    /// bp length of the query
    pub q_length: u64,
    /// cigar (might be corrected)
    pub cigar: String,
    /// chromosome/target of match
    pub chrom: String,
    /// position of match on genome/contig
    /// .0-based
    pub fp_target: u64,
    /// number of clipped position upstream
    pub upstream: u32,
    /// number of clipped position downstream
    pub downstream: u32,
    /// gap compressed similarity score
    pub similarity: f32,
    /// length of the sequence (includes INs but no DELs and CLIP)
    pub seq_length: u32,
    /// length with splicing of match (includes INDELS and N, and CLIP)
    pub match_length_s: u32,
    /// length with splicing of match (only matches and mismatchs)
    pub match_length: u32,
    /// position of match on query
    /// .0-based
    pub fp_query: u32,
    /// start of alignment
    /// .0-based
    pub aln_start: u64,
    /// end of alignment
    /// .0-based
    pub aln_end: u64,
    /// strand of match
    pub match_strand: StrandDirection,
}

/*
#[derive(Debug)]
/// matching features of fusion elements
struct ScoredFeatures{
    /// the name , e.g. Gene name
    feature_name  : String ,
    /// the score which describes if lying
    /// within or overlapping
    feature_score : u8 ,
    /// the number of matching basepairs
    /// which is not splice-aware!
    feature_match : u64 ,
}
*/

#[derive(Debug,Clone,Copy)]
pub struct VersionInfo <'a>{
    /// the used program/sub-program
    pub program  : &'a str,
    /// the version of the program
    pub version  : &'a str,
    /// the author
    pub author : &'a str,
    /// the executed command
    pub command : &'a str,
}


#[derive(Debug,Clone,PartialEq,Eq)]
/// this is the richest information
/// which skippy obtains after evaluating
/// all information related to a given reads
/// .0-based
pub struct FullFusionEvidence {
    /// cigar of primary match
    pub cigar_pr: String,
    /// cigar of SA
    pub cigar_sa: String,
    /// name of chr primary match
    pub chrom_pr: String,
    /// name of chr SA
    pub chrom_sa: String,
    /// name of the query read/cDNA
    pub query_name: String,
    /// length in bp of query above
    pub query_length: u64,
    /// position of break on target from primary match
    /// .0-based
    pub target_break_pr: u64,
    /// position of break on target from SA
    /// .0-based
    pub target_break_sa: u64,
    /// position of break on query from primary match
    /// .0-based
    pub query_break_pr: u64,
    /// position of break on query from SA
    /// .0-based
    pub query_break_sa: u64,
    /// strand of primary match
    pub strand_pr: StrandDirection,
    /// strand of SA
    pub strand_sa: StrandDirection,
    /// placement of primary match in fusion
    pub orient_pr_upstream: bool,
    /// placement of SA in fusion
    pub orient_sa_upstream: bool,
    /// fusion position of primary alignment
    /// .0-based
    pub aln_start_pr: u64,
    /// fusion position of SA
    /// .0-based
    pub aln_start_sa: u64,
    /// fuzzyness of precision as bps in both direction
    /// of the middle fusion point indicated
    pub precision: Option<u32>,
    /// support, 0 in initial single entry derived
    /// support, 1 if primary alignment twice observed
    pub support: u32
}


#[derive(Debug,Clone,Default,PartialEq)]
/// this structure organizes 
/// fusion evidence based on coordinates
/// which derive from multiple read-based evidence
/// potentially.
/// For the positions we do provide an exact number
/// and in the precision_pr/sa we accumulate the exact number of
/// each read evidence to have later a value to calculate
/// precision based on
/// .0-based
pub struct PosBasedFusEvidence {
    /// name of chr primary match
    pub chrom_pr: String,
    /// name of chr SA
    pub chrom_sa: String,
    /// position of break on target from primary match
    /// .0-based
    pub target_break_pr: u64,
    /// position of break on target from SA
    /// .0-based
    pub target_break_sa: u64,
    /// strand of primary match
    pub strand_pr: StrandDirection,
    /// strand of SA
    pub strand_sa: StrandDirection,
    /// placement of primary match in fusion 
    pub orient_pr_upstream: bool,
    /// placement of SA in fusion 
    pub orient_sa_upstream: bool,
    /// collection of target_break_pr values
    pub precision_pr: Vec<u64>,
    /// collection of target_break_sa values
    pub precision_sa: Vec<u64>,
    /// fusion position of primary alignment
    /// .0-based
    pub aln_start_pr: Vec<u64>,
    /// fusion position of SA
    /// .0-based
    pub aln_start_sa: Vec<u64>,
    /// Split-read support is an increasing number of seen
    /// evidence supporting that position
    pub support_sr: u32,
    /// Split-pair support is an increasing number of seen
    /// evidence supporting that position
    pub support_sp: u32,
    /// The ratio indicates the number of reads spanning a fusion
    /// compared to the median coverage in our inquired region. 
    /// This ratio is reported for both sites and the lowest is given in
    /// this field.
    /// Especially for positions where we have ambiguity for 
    /// mapping or lots of noise that ratio will be lower and help to identify problems.
    pub sp_ratio: f32,
    /// provides the number of uncertain positions
    /// around the fp points
    pub cipos : Option<u32>
}

/// this structure is specifically for the usage of annotated reads
/// It presumes that a long read/cDNA is composed of 2 elements and one
/// of them is mapping onto a human genome.
/// This helps then for e.g. virus integration to understand better
/// integration sites but might be used, e.g. for transgene integration too.
/// It organizes information which are necessary
/// to understand all basic information for propagation
/// of such information, e.g. into a VCF file subsequently
#[derive(Debug,PartialEq,Eq, Clone)]
pub struct AnnotatedRead {
    /// the reference start coordinate on the read
    pub ref_start  : u64,
    /// the reference end coordinate on the read
    pub ref_end    : u64,
    /// the strand on the reference chromosome 
    pub ref_strand : StrandDirection,
    /// the reference chromosome on which it maps
    pub ref_chr    : String,
    /// the reference genomic mapping start coordinates
    pub ref_gstart : u64,
    /// the reference genomic mapping end coordinates
    pub ref_gend   : u64,
    /// the name of the read/contig
    pub read_name    : String,
    /// the alien start coordinate on the read
    pub alien_start  : u64,
    /// the alien end coordinate on the read
    pub alien_end    : u64,
    /// integration site in the genome
    pub integration_site : u64,
    /// integration of the integration event
    pub integration_dir  : StrandDirection,
}

#[derive(Default,Debug,PartialEq)]
pub struct OverlapResult {
    pub match_length : u64,
    /// match score is a scoring how well 
    /// the annotated feature in the GTF corresponds to the 
    /// presented alignment. 
    /// If an element is stranded, the direction of the annotation
    /// and the alignment is very important and provides  a bonus of 10.
    /// If a feature start is identical to an alignment start 
    /// score increases by 10 - similar if end sections are identical.
    /// If the alignment start/end is located within an annotated element 
    /// but start/end are not identical the score increases by 5 instead of
    /// 10. 
    /// If is extends beyond the score increases only by 1.
    /// The maximum score therefore can be 40 which would be both sides annotated
    /// in proper matching edges (almost impossible from a conceptional point).
    /// The lowest score being 2. Or 0 if nothing is provided
    pub match_score  : u8,
    pub feature_name : String,
}


#[derive(Debug)]
enum Keeper  {
    Yes,
    No,
    Donnu,
    Doubled,
}



#[derive(Debug,PartialEq)]
#[derive(Clone)]
/// this has many parallels to the richer structure "FullFusionEvidence"
/// but drops some of the information as not anymore needed 
/// Important: the coordinates are always sorted, meaning start is supposed
/// to be smaller (upstream) than end (downstream)
/// Therefore as well the need of the strand annotation
/// Positions are 0-based
pub struct PrincipleReadBasedFusionOutput {
    /// given cDNA
    pub query_name: String,
    /// length of cDNA
    pub length: u64,
    /// fusion point on cDNA
    /// .0-based
    pub fusion_point: u64,
    /// GeneA--GeneB or chrA--chrB without annotation provided
    pub fusion_genes: Vec<String>,
    /// fuzzyness of correct fusion point resolution in bp
    pub fp_fuzzy: Option<u32>,
    /// distance between intra-chromosomal events
    pub fusion_distance: Option<u64>,
    /// chrom of geneA
    pub a_chrom: String,
    /// begin of alignment of geneA
    /// .0-based
    pub a_start: u64,
    /// end of alignment of geneA
    /// .0-based
    pub a_end: u64,
    /// fusion point of geneA
    /// .0-based
    pub a_fp: u64,
    /// strand of alignment of geneA 1/-1
    pub a_strand: StrandDirection,
    /// begin of alignment of geneB
    pub b_chrom: String,
    /// end of alignment of geneB
    /// .0-based
    pub b_start: u64,
    /// end of alignment of geneB
    /// .0-based
    pub b_end: u64,
    /// strand of alignment of geneB 1/-1
    pub b_strand: StrandDirection,
    /// fusion point of geneA
    /// .0-based
    pub b_fp: u64
}


#[derive(Debug,Clone)]
/// this has many parallels to the richer structure "PosBasedFusEvidence"
/// but drops some of the information as not anymore needed in final output
/// Important: the coordinates are always sorted, meaning start is supposed
/// to be smaller (upstream) than end (downstream)
/// Therefore as well the need of the strand annotation
pub struct PrinciplePosBasedFusionOutput {
    /// GeneA--GeneB or chrA--chrB without annotation provided
    pub fusion_genes: Vec<String>,
    /// fuzzyness of correct fusion point resolution in bp
    pub fp_fuzzy: Option<u32>,
    /// distance between intra-chromosomal events
    pub fusion_distance: Option<u64>,
    /// chrom of geneA
    pub a_chrom: String,
    /// begin of alignment of geneA
    /// .0-based
    pub a_start: u64,
    /// end of alignment of geneA
    /// .0-based
    pub a_end: u64,
     /// fusion point of geneA
     /// .0-based
     pub a_fp: u64,
    /// strand of alignment of geneA 1/-1
    pub a_strand: StrandDirection,
    /// begin of alignment of geneB
    /// .0-based
    pub b_chrom: String,
    /// end of alignment of geneB
    /// .0-based
    pub b_start: u64,
    /// end of alignment of geneB
    /// .0-based
    pub b_end: u64,
    /// strand of alignment of geneB 1/-1
    pub b_strand: StrandDirection,
    /// fusion point of geneA
    /// .0-basedvariance
    pub b_fp: u64,
    /// number of SR support
    pub sr_support: u32,
    /// number of SP support
    pub sp_support: u32,
    /// The ratio indicates the number of reads spanning a fusion
    /// compared to the median coverage in our inquired region. 
    /// This ratio is reported for both sites and the lowest is given in
    /// this field.
    /// Especially for positions where we have ambiguity for 
    /// mapping or lots of noise that ratio will be lower and help to identify problems.
    pub sp_ratio: f32
}


/// A simple structure to compare 2 vcf files of the type e
/// We have a concatenated ID+MATEID together with the coordinates
/// of the events and the chromosomes.
/// In the matches vector we can then collect entries which match 
/// to later get the corresponding vcf records
#[derive(Debug,PartialEq,Eq,Clone,Default)]
pub struct FileAbasedComp {
    /// the concatenated ID and Mate-ID
    pub a_id: String,
    /// the chromosome of the primary event
    pub a_chr1: String,
    /// the fusion point of the primary event
    pub a_fp1: u64,
    /// the chromosome of the secondary event
    pub a_chr2: String,
    /// the fusion point of the secondary event
    pub a_fp2: u64,
    /// matches within the defined range, vector of concatenated ID and Mate-ID
    pub b_matches: Vec<String>
}


/// this describes possible directions
/// of the BND event. It can obviously be only
/// forward or reverse but the issue is often
/// that we cant determine the direction 
/// from a single BND entry as we need the 
/// 2nd one to complement the direction of the 
/// mate event
#[derive(Debug,Clone,Copy,Hash,Eq, PartialEq,Default)]
pub enum StrandDirection {
    Fwd,
    Rev,
    #[default]
    Unknown,
}


/// A simple structure to hold essential  BND information
/// We have a  SAMPLE and ID as well as MATEID together with the coordinates
/// of the events and the chromosomes.
/// In the matches vector we can then collect entries which match 
/// to later get the corresponding vcf records
/// Positions are 0 based
#[derive(Debug,Clone,Default,PartialEq,Eq)]
pub struct BNDentry {
    // a sample information
    pub sample: String,
    // the BND ID
    pub id: String,
    // the matching mate BND ID
    pub mid: String,
    // the chromosome of the primary event
    pub chr1: String,
    // the fusion point of the primary event
    pub fp1: u64,
    // the start point of primary event
    pub st1: u64,
    // the end point of the primary event
    pub end1: u64,
    // orientation of 1st event
    pub forward1: StrandDirection,
    // the chromosome of the secondary event
    pub chr2: String,
    // the fusion point of the secondary event
    pub fp2: u64,
    // the fusion start point of the 2nd event
    pub st2: u64,
    // the fusion end point of the 2nd event
    pub end2: u64,
    // the orientation of the second event
    pub forward2: StrandDirection,
    // the type
    pub sv_type: SVType,
    // a name one can provide, e.g. 
    // in case of gene fusions the name of the genes
    pub name: Option<String>,

}

/// This is the equivalent of the BNDentry structure
/// which holds information that is INDEL related
#[derive(Debug,Clone,Default)]
pub struct INDELentry {
    // a sample information
    pub sample: String,
    // the BND ID
    pub id: String,
    // the chromosome of the event
    pub chr: String,
    // the start point of the event
    pub st: u64,
    // the end point of the event
    pub end: u64,
    // the type
    pub sv_type: SVType,
    // a name one can provide, e.g. 
    // in case of gene fusions the name of the genes
    pub name: Option<String>,

}

#[derive(Debug,Clone,Eq,Hash,PartialEq,Copy,Default)]
pub enum SVType  {
    BndPair,
    BndSingle,
    INS,
    DEL,
    TRA,
    INV,
    #[default]
    Unknown,
}


/// this enum allows functions to take all type
/// of SV trees potentially, when they are compatible
pub enum SvTypeEntry {
    Bnd(BNDentry),
    InDel(INDELentry),
}


/// this structure bundles all possible SVs in the proper structured
/// trees
#[derive(Debug,Clone)]
pub struct FullSvTrees {
    // a tree for which we expect paired entries which are
    // linked through the ID-MATEID pairs
    pub bnd_pair   : FxHashMap<String,IntervalTree<u64,BNDentry>>,
    // a tree for which we have single BND events and dont know
    // the 2nd part of the event
    pub bnd_single : FxHashMap<String,IntervalTree<u64,BNDentry>> ,
    // classical deletion or insertion events
    pub indel      : FxHashMap<String,IntervalTree<u64,INDELentry>>,
}

/// this structure bundles all possible SVs in the proper structured
/// trees
#[derive(Debug,Clone)]
pub struct BedpeTrees {
    // a tree for which we expect paired entries which are
    // linked through the ID-MATEID pairs
    pub paired   : FxHashMap<String,IntervalTree<u64,ClusterPairs>>,
    // This is everything which is not paired
    pub single : FxHashMap<String,IntervalTree<u64,SvClusters>> ,
}

/// cluster pairs describes a structure which is used as keys
/// in a hash table of entries which have a similar BND cluster ID 
/// and BND mate cluster id
/// Holds as well the direction, which is though currently not used.
#[derive(Debug,Clone,Default,Hash,Eq, PartialEq)]
pub struct ClusterPairs {
    // the first chromosome/contig of the cluster of the BND event
    pub chr1: String,
    // the range from the 1st cluster of the BND event
    pub range1: std::ops::Range<u64>,
    // the orientation of the 1st event
    pub forward1: StrandDirection,
    // the 2nd chromosome/contig of the cluster of the BND event
    pub chr2: String,
    // the range from the 2nd cluster of the BND event
    pub range2:std::ops::Range<u64>,
    // the orientation of the matching mate event
    pub forward2: StrandDirection,
}

///  describes a structure which is used as keys
/// in a hash table of entries which have a similar SV cluster IDs
/// Holds as well the direction, which is though currently not used.
#[derive(Debug,Clone,Default,Hash,Eq, PartialEq,)]
pub struct SvClusters {
    // the first chromosome/contig of the cluster of the SV event
    pub chr: String,
    // the range from the 1st cluster of the SV event
    pub range: std::ops::Range<u64>,
    // the orientation of the of the SV event
    pub forward1: StrandDirection,
    // type of the SVcluster
    pub svtype: SVType,
}


/// Structure which describes all information encoded
/// in a second allele of a BND entry which has a mate
#[derive(Debug,Clone,Eq, PartialEq)]
pub struct Bnd2ndAllele {
    // the reference nucelotide which is altered
    pub nucleotide: String,
    // the mate chromosome where the 2nd entry is located
    pub chr: String,
    // the position on the mate chromosome 0-based
    pub pos: u64,
    // the direction of the primary event
    pub prim_direction: StrandDirection,
    // the direction of the associated event
    pub second_direction: StrandDirection
}

#[derive(Debug,Clone,Hash,Eq, PartialEq,)]
// just the difference sense which we allow for DNA
// which is the normal sense, complementary, reverse and reverse complementary
pub enum DnaReturnSense {
    Fwd,
    FwdC,
    Rev,
    RevC
}


/// this function parses a BEDPE files for the organization
/// into pairs of ranges. These pairs are predominantly thought
/// to be used in the comparison with gene-fusion or BND events
/// in general but could potentially used in different scenarios as well.
/// According to the BEDPE specs one can use for the scoring   to 1000, inclusive.
/// But many use as well floats. Here we will accept float but they will actually be
/// converted into u32 with a warning being displayed.
/// The function returns a hashmap containing all chromosomes and for each
/// an interval tree containing the range and the ClusterPair feature observed
/// 
/// One issue though is that it is pretty arbiturary which element is first and
/// which is second and order might stem from input program which we dont control.
/// Therefore for any given event we need to adapt this accordingly.
/// 
/// Another issue is that there is a need to treat the last nucleotide of a contig slightly
/// different as the API to find a range will always expect for the range a +1 extension
/// Therefore we check in the BEDPE if we have the max of the chromosome and extend in that case the
/// range in the tree by +1
/// 
/// Unittest: TRUE
///
/// ```rust
/// use genefusion::lib::common::{*};
/// use genefusion::lib::hts_lib_based::{*};
/// use bio::data_structures::interval_tree::{*};
/// use rust_htslib::bcf::header::Header;
/// use rust_htslib::bcf::{Read as BcfRead};
/// use rust_htslib::bcf::{Format};
/// use std::fs::File;
/// use std::io::prelude::*;
/// use rustc_hash::FxHashMap;
/// pretty_env_logger::init();
/// 
/// let mut chromosomes : FxHashMap<String,u64> = FxHashMap::default();
/// chromosomes.insert(String::from("chr1"),100000);
/// chromosomes.insert(String::from("chr2"),100000);
/// chromosomes.insert(String::from("chr4"),190214555);
/// chromosomes.insert(String::from("KI270466.1"),1233);
/// // create dummy files
/// let mut file_bedpe_ff =  File::create("test_ff.bedpe").expect("ERROR:could not create BEDPE file!");
/// file_bedpe_ff.write_all(b"chr1\t5\t10\tchr2\t15\t20\tPA1214:bnd_SMAP1265_1;PA2566-R2P8:bnd_SMAP1247_1;_2\t2\t+\t+\tBND_PAIR\n").expect("ERROR:could not write BEDPE file!");
/// /// second case with last position
/// file_bedpe_ff.write_all(b"KI270466.1\t1231\t1232\tchr4\t49711479\t49711480\ttest_file:NGSAI_GENEFUSION_FID.4240.2;\t1\t.\t.\tBND_PAIR\n").expect("ERROR:could not write BEDPE file!");
/// 
/// let tree      = svtree_from_bedpe(&String::from("test_ff.bedpe"),None,&chromosomes);
/// let result1   = tree.paired.get("chr1").unwrap();
/// let result2   = tree.paired.get("chr2").unwrap();
/// let result3   = tree.paired.get("chr4").unwrap();
/// let result4   = tree.paired.get("KI270466.1").unwrap();
/// for r in result1.find(6..7){
///     assert_eq!(r.interval().start,5);
///     assert_eq!(r.interval().end,10);
/// };
/// for r in result2.find(17..18){
///     assert_eq!(r.interval().start,15);
///     assert_eq!(r.interval().end,20);
/// };
/// for r in result3.find(49711479..49711480){
///     assert_eq!(r.interval().start,49711479);
///     assert_eq!(r.interval().end,49711480);
/// };
/// 
/// for r in result4.find(1232..1233){
///     assert_eq!(r.interval().start,1231);
///     assert_eq!(r.interval().end,1233);
/// };
/// 
/// ```
pub fn svtree_from_bedpe(
    my_path: &str, 
    cutoff : std::option::Option<&str>,
    contigs: &FxHashMap<String, u64>,
) ->  BedpeTrees  {
    // get filtering score
    let min_score : Option<u8> = cutoff
        .map(|_| cutoff.unwrap().parse::<u8>().expect("ERROR: could not parse the cut-off correctly!"));

    assert!(
        Path::new(my_path).exists(),
        "ERROR: BEDPE file {:?} does not exists!",
        my_path
    );

    let mut pair_trees  : FxHashMap<String,IntervalTree<u64,ClusterPairs>>  = FxHashMap::default();
    let mut single_trees: FxHashMap<String,IntervalTree<u64,SvClusters>>    = FxHashMap::default();
    

    for key in contigs.keys() {
        // make a new IntervalTree
        let tree1 = IntervalTree::new();
        let tree2   = IntervalTree::new();
        pair_trees.insert(key.clone(), tree1);
        single_trees.insert(key.clone(), tree2);
    };
    eprintln!("INFO: Chromosomes in interval tree from BEDPE file : {}", pair_trees.len());
    let input  = File::open(my_path).expect("Unable to open BEDPE file");
    let reader = BufReader::new(input);
    let mut n_ranges: usize = 0;
    for (_line_number,line) in reader.lines().enumerate() {
        let l = line.expect("Unable to read line");
        // now we split our BEDPE by tab and fill a vector with it
        let elements: Vec<&str> = l.split('\t').collect();
        // if we have < 10 elements in vector then BEDPE is not valid
        if elements.len() < 11  {
            panic!("ERROR: BEDPE does not contain minimal fields + last with type annotation! ");
        };
        // now the 11th needs to determine the type, if it is a BND with pair
        // then we expect both being accordingly filled 
        let sv_type = match elements[10] {
            "BND_PAIR" => SVType::BndPair,
            "INS"      => SVType::INS,
            "DEL"      => SVType::DEL,
            _          => panic!("ERROR: for 11th column only these types are supported: BND_PAIR, INS, DEL !"),
        };
        // the strand is normally encoded in "+" and "-" but can be
        // as well not defined
        let chrom1  = elements[0].to_string();
        let start1  = elements[1].parse::<u64>().unwrap();
        let c_length1 = contigs.get(&chrom1).expect("ERROR: could not get chromosome length from chromosome file!");
        let mut end1    = elements[2].parse::<u64>().unwrap();
        // verify now that we have not last bp of chrom and if case correct
        if &(end1 +1)  == c_length1 {
            end1 +=1;
        };

        let strand1 = match elements[8]{
            "+" => StrandDirection::Fwd,
            "-" => StrandDirection::Rev,
            "." => StrandDirection::Unknown,
            _ => panic!("ERROR: your strand information is not valid!"),
        };

        let _fname  = elements[6].to_string();
        let score   = elements[7].parse::<u8>().unwrap();
        if min_score.is_some() && score <  min_score.unwrap()  {
            continue
        }

        if sv_type == SVType::BndPair {
            let chrom2  = elements[3].to_string();
            let start2  = elements[4].parse::<u64>().unwrap();
            let mut end2    = elements[5].parse::<u64>().unwrap();
            let c_length2 = contigs.get(&chrom2).expect("ERROR: could not get chromosome length from chromosome file!");
            // verify now that we have not last bp of chrom and if case correct
            if &(end2 +1)  == c_length2 {
                end2 +=1;
            };

            let strand2 = match elements[9]{
                "+" => StrandDirection::Fwd,
                "-" => StrandDirection::Rev,
                "." => StrandDirection::Unknown,
                _ => panic!("ERROR: your strand information is not valid!"),
            };
        
            let record = ClusterPairs {
                chr1: chrom1.clone(),
                range1: std::ops::Range {
                    start: start1,
                    end: end1
                },
                forward1: strand1,
                chr2: chrom2.clone(),
                range2: std::ops::Range {
                    start: start2,
                    end: end2,
                },
                forward2: strand2,
            };

            // now depending on sorting or other aspects
            // a +/- event can be perceived as a -/+ event as
            // well. E.g. whenever we have not stranded sequencing
            // which is **predominantly** the case.
            // Therefore we need to report both direction on both s
            // e.g. a chr2:FWD and chr6:REV
            // should be encoded in   chr2:chr6 FWD:REV
            // and similarly as well  chr6:chr2 REV:FWD

            let record2 = ClusterPairs {
                chr1: chrom2.clone(),
                range1: std::ops::Range {
                    start: start2,
                    end: end2
                },
                forward1: strand2,
                chr2: chrom1.clone(),
                range2: std::ops::Range {
                    start: start1,
                    end: end1,
                },
                forward2: strand1,
            };


            let mut local_tree1 = pair_trees.get_mut(&record.chr1).expect("ERRROR: could not find chromosome in tree!");
            local_tree1.insert(&record.range1,record.clone());
            local_tree1.insert(&record.range1,record2.clone());
            local_tree1 = pair_trees.get_mut(&record.chr2).expect("ERRROR: could not find chromosome in tree!");
            local_tree1.insert(&record.range2,record.clone());
            local_tree1.insert(&record.range2,record2.clone());
            n_ranges +=1;
        
        }else if sv_type == SVType::INS || sv_type == SVType::DEL {
            // not necessary to write out again, but more bulletproof 
            // if I add later other types rather then just else
            let record = SvClusters {
                chr: chrom1,
                range: std::ops::Range {
                    start: start1,
                    end: end1
                },
                forward1: strand1,
                svtype: sv_type,
            };
            let local_tree1 = single_trees.get_mut(&record.chr).expect("ERRROR: could not find chromosome in tree!");
            local_tree1.insert(&record.range,record.clone());
            n_ranges +=1;
        }
    }
    eprintln!("INFO: Successfully parsed BEDPE entries: {}",n_ranges);
    //eprintln!("INFO: Entries in intervall tree: paired BND {}; single {};",pair_trees.len(),single_trees.len());
    BedpeTrees {
        paired: pair_trees,
        single: single_trees,
    }
}


/// this function takes a fofn file 
/// with one entry per line and returns simply an array
/// IF fofn is false it assumes it is only a single file
/// instead. This is a work-around as the tree building
/// generates otherwise issues if done per file base
///
/// Unittest: TRUE
///
pub fn parse_fofn(
    my_file: &str,
    fofn   : bool
) -> Vec<String> {
    let mut files : Vec<String> = Vec::new();
    if fofn {
        eprintln!("INFO: FOFN {} provided, reading entries...", my_file);
        assert!(
            Path::new(my_file).exists(),
            "ERROR: ID list file {:?} does not exist!",
            my_file
        );

        let input  = File::open(my_file).expect("Unable to open filter file");
        let reader = BufReader::new(input);
            
        for line in reader.lines() {
            let l = line.expect("ERROR: could not read line!");
            // now we split by empty space
            let e: Vec<&str> = l.split(' ').collect();
            if e.len() != 1 {
                panic!("ERROR: your FOFN contains delimited entries!");
            }
            files.push(e[0].to_string());
        }
    }else{
        files.push(my_file.to_string());
    }
    files
}


/// this function takes a txt file 
/// with one entry per line and returns simply a Hashmap
/// to have the collection of chromosomes for a white/black-list
/// 
/// Unittest: TRUE
///
pub fn parse_chroms_txt(
    my_file: &str,
) -> FxHashMap<String,u8> {
    let mut chroms : FxHashMap<String,u8> = FxHashMap::default();
    assert!(
        Path::new(my_file).exists(),
        "ERROR: List with chromosomes {:?} does not exist!",
        my_file
    );
    let input  = File::open(my_file).expect("Unable to open filter file");
    let reader = BufReader::new(input);
    for line in reader.lines() {
        let l = line.expect("ERROR: could not read line!");
        // now we split by empty space
        let e: Vec<&str> = l.split(' ').collect();
        if e.len() != 1 {
            panic!("ERROR: your chrom file contains whitespace!");
        }
        chroms.insert(e[0].to_string(), 1);
        
    }
    chroms
}


/// this function takes a fofn file 
/// with one entry per line and returns simply a hashmap
/// of the chromosome name and size
/// 
///Unittest: TRUE
///
pub fn parse_chrom_file(
    my_file: &str
) -> FxHashMap<String,u64> {

    let mut chromosomes : FxHashMap<String,u64> = FxHashMap::default();
    assert!(
        Path::new(my_file).exists(),
        "ERROR: Chromosome file {:?} does not exist!",
        my_file
    );

    let input  = File::open(my_file).expect("Unable to open filter file");
    let reader = BufReader::new(input);
        
    for l in reader.lines() {
        let line = l.expect("ERROR: could not read line!");
        let elements: Vec<&str> = line.split('\t')
            .collect();
        if elements.len() != 2 {
            panic!("ERROR: Chromosomes must be tab-separated \"Chrom\"\t \"size in bp\"!");
        }
        if chromosomes.contains_key(elements[0]) {
            panic!("ERROR: chromosome twice encountered!");
        }else{
            chromosomes.insert(elements[0].to_string(), elements[1].parse::<u64>().unwrap());
        }
    }
    chromosomes
}

/// adapted from here https://users.rust-lang.org/t/efficient-way-of-checking-if-two-files-have-the-same-content/74735
/// very useful for tests with external files and to verify that the results is identical
/// to a previously manually generated result file
pub fn is_same_file(
    file1: &Path, 
    file2: &Path
) -> Result<bool, std::io::Error> {
    println!("INFO: comparing file1 {:?} and file2 with each other {:?}", file1.to_str(), file2.to_str());
    let f1 = File::open(file1).expect("ERROR: could not open file");
    let f2 = File::open(file2).expect("ERROR: could not open file");
    

    // Use buf readers since they are much faster
    let f1r = BufReader::new(f1);
    let f2r = BufReader::new(f2);

    // Do a byte to byte comparison of the two files
    for (b1, b2) in f1r.bytes().zip(f2r.bytes()) {
        if b1.unwrap() != b2.unwrap() {
            return Ok(false);
        }
    }
   Ok(true)
}

/// this takes an indexed FASTA reader and return a sequence as a string
/// It does account as well for direction and will reverse the sequence
/// if sense if Rev. Additionally, it will return a shortened sequence
/// if it realized that it is out of boundaries.
/// Careful! This might be not intended by the user! And will only print a warning in that case.
/// Furthermore, if start and stop are reversed it will sort them first.
/// Returns the sequence as a String
///
/// Unittest: TRUE
/// 
/// Example:
/// ```rust
/// use genefusion::lib::common::fasta_extract_and_correct;
/// use genefusion::lib::common::DnaReturnSense;
/// use bio::io::fasta::IndexedReader;
/// use std::fs::File;
/// use std::io::prelude::*;
/// pretty_env_logger::init();
/// 
/// // create dummy files
/// let mut file_seq =  File::create("foo.fa").expect("ERROR:could not create seq file!");
/// file_seq.write_all(b">chr1\nGTAGGCTGAAAA\nCCCC").expect("ERROR:could not write seq file!");
/// 
/// // index as well
/// let mut file_idx =  File::create("foo.fa.fai").expect("ERROR:could not create idx file!");
/// file_idx.write_all(b"chr1\t16\t6\t12\t13").expect("ERROR:could not write idx file!");
/// 
/// let contig:  &str = &String::from("chr1");
/// let start: u64 = 0;  // start is 0-based, inclusive
/// let stop : u64 = 10; // stop is 0-based, exclusive
/// let far  : u64 = 20; // this is too far
/// let maxn : u64 = 16; // maximum
/// 
/// // now  reversing complementing the sequence
/// let result = fasta_extract_and_correct(&String::from("foo.fa"),contig,start,stop,DnaReturnSense::RevC);
/// assert_eq!(result, String::from("TTCAGCCTAC"));
/// 
/// 
/// // cleaning up the fasta and index file
/// std::fs::remove_file("foo.fa").expect("ERROR: could not remove the temporary fasta file!");
/// std::fs::remove_file("foo.fa.fai").expect("ERROR: could not remove the temporary fasta index file!");
/// ```
pub fn fasta_extract_and_correct (
    faidx : &str,
    contig: &str,
    pos1: u64,
    pos2: u64,
    sense: DnaReturnSense,
)-> String {
    // the vector into which we read the sequence
    let mut seq = Vec::new();   
    let mut fasta_reader = bio::io::fasta::IndexedReader::from_file(&faidx).expect("ERROR: could not open reference FASTA file!");
    // first we fetch, then we read
    // we order the positions as the function goes otherwise wild
    let (start,end) : (u64,u64) ;
    match pos1.cmp(&pos2) {
        Ordering::Less => {
            start = pos1;
            end   = pos2;
        },
        Ordering::Greater => {
            start = pos2;
            end   = pos1;
        },
        _ => {
            panic!("ERROR: extracting the positions from FASTA resulted in non compliant coordinates!");
        },
    };
    //println!("Extracting from chr {} position {} to {} in sense {:?}",&contig,&start,&end,&sense);
    fasta_reader.fetch(
        contig,
        start,
        end 
    ).expect(
        "ERROR: could not fetch interval on reference, either not present or coordinates out of scope "
    );
    if fasta_reader.read(&mut seq).is_err() {
       // now lets figure out it's size
       // and then get the extreme
        let mut tmp = Vec::new();   
        fasta_reader.fetch_all(contig).expect("ERROR: could not fetch entire chromosome");
        fasta_reader.read(&mut tmp).expect("ERROR: could not read sequence from reference");
        // now we get the size
        let c_size = tmp.len();
        // now we get the last position
                    
        fasta_reader.fetch(contig,start, c_size as u64 ).expect("ERROR: could not fetch interval on reference, either not present or coordinates out of scope ");
        fasta_reader.read(&mut seq).expect("ERROR: could not read sequence from reference");
        eprintln!("INFO: position found outside of boundaries for chr: {} position: {} selecting: {}",
            contig,
            start,
            &c_size,
        ); 
    }
    let nuc = String::from(from_utf8(&seq).unwrap());
    // now comes the different operations
    // depending on the sense one wants
    match sense{
        DnaReturnSense::Fwd  => nuc,
        DnaReturnSense::Rev  => nuc.chars().rev().collect::<String>(),
        DnaReturnSense::RevC => {
            let tmp = revcomp(nuc.as_bytes());
            String::from_utf8(tmp).expect("ERROR: Found invalid UTF-8")
        },
        DnaReturnSense::FwdC => {
            let tmp = nuc.as_bytes().iter().map(|x| complement(*x)).collect();
            String::from_utf8(tmp).expect("ERROR: Found invalid UTF-8")
        },
    }
}

/// this function will read a provided BED file expecting
/// 2 lines/read and being sorted by the read-name.
/// It verifies the order on the read, the direction and
/// genomic mappings. Finally it returns an annotated read structure
/// which contains the genomic integration position and the direction.
/// The overlap options defines the maximum overlap between the alien and
/// normal feature in base-pairs. This is often necessary as mapping
/// is not a perfect process and we can potentially generate overlaps.
/// If there is unresolved ambiguity concerning the overlap then it will print a 
/// warning and continue with the next entry.
/// If there is no regex match for the alien and normal annotation it will exit hard.
///
///Unittest: TRUE
///
/// ```rust
/// use genefusion::lib::common::{*};
/// use std::fs::File;
/// use std::io::prelude::*;
/// use std::str::from_utf8;
/// use tempfile::NamedTempFile;
/// pretty_env_logger::init();
/// 
/// let tmp = NamedTempFile::new().unwrap();
/// let path = tmp.path();
/// // This contains all possible combinations:
/// let mut file_A =  File::create("A.bed").expect("ERROR:could not create seq file!");
/// file_A.write_all(b"m64143_211228_145337/1508655/ccs\t0\t1250\tmouse_merged,muERV_merged\t0\t.\n\
///     m64143_211228_145337/1508655/ccs\t1275\t8961\thuman_merged\t-\tchr9@132857445-132865123\n").expect("ERROR:could not write BED file!");
/// let input    = String::from("A.bed");
/// let result1  = read_2lines_bed(&input,"human","muERV",&1_u64);
/// let result_A = AnnotatedRead {
///     ref_start: 1275_u64,
///     ref_end: 8961_u64,
///     ref_strand: StrandDirection::Rev,
///     ref_chr: String::from("chr9"),
///     ref_gstart: 132857445_u64,
///     ref_gend: 132865123_u64,
///     read_name: String::from("m64143_211228_145337/1508655/ccs"),
///     alien_start: 0_u64,
///     alien_end: 1250_u64,
///     integration_site: 132865123_u64,
///     integration_dir: StrandDirection::Fwd,
/// };
/// assert_eq!(result1[0],result_A);
/// // delete file
/// std::fs::remove_file(input).unwrap();
/// ```
pub fn read_2lines_bed (
    bed: &str,
    regex_human: &str,
    regex_mouse: &str,
    overlap: &u64
) -> Vec<AnnotatedRead> {
    let h_regex                       = Regex::new(regex_human).unwrap();
    let m_regex                       = Regex::new(regex_mouse).unwrap();
    let input                          = File::open(bed).expect("Unable to open BED file");
    let reader              = BufReader::new(input);
    let mut collection : Vec<AnnotatedRead>  = Vec::new();

    // now we read line by line and work then on pairs of lines
    let mut n        : i64 = 0 ;
    let mut previous : Vec<String> = Vec::new();
    for l in reader.lines(){
        n += 1;
        let line = l.expect("ERROR: could not read line!");
        let current: Vec<String> = line.split('\t').map(|x| x.to_string()).collect::<Vec<String>>();
        if current.len() != 6 {
            panic!("ERROR: Did not encounter a valid 6 column tab-separated BED file !");
        }
        if n % 2 != 0 {
            previous = current;
        }else{
            // now we got a pair of entries and we can continue
            if previous[0] != current[0] {
                panic!("ERROR: encountered pairs which had different read names. Verify that file is sorted by read-name!")
            }
            //compiled_info.read_name =  previous[0].clone();

            let current_a     : bool;
            let current_human : bool;
            // let's now decide which ( current or previous ) is the beginning and end
            // of the read annotation 
            //if previous[1] <= current[2] {
            if (current[1].parse::<u64>().unwrap() + overlap) >= previous[2].parse::<u64>().unwrap() {    
                // previous entry is the A one 
                //eprintln!("INFO: entry A is ");
                //dbg!(&previous);
                current_a = false;
            }else if (previous[1].parse::<u64>().unwrap() + overlap)  >= current[2].parse::<u64>().unwrap() {
                // current entry is the A one 
                //eprintln!("INFO: entry A is ");
                current_a = true;
            }else{
                eprintln!("WARNING: could not determine for {} which one is A and which is B ",&previous[0]);
                continue;
            }

            // next we need to define unambigiously which one is human 
            // and which entry is mouse-virus
            if h_regex.is_match(&current[3]) && m_regex.is_match(&previous[3]) {
                // current entry is the human 
                //eprintln!("INFO: human is ");
                //dbg!(&current);
                current_human = true;
            }else if h_regex.is_match(&previous[3]) && m_regex.is_match(&current[3]) {
                // current entry is the A one 
                //eprintln!("INFO: mouse is ");
                current_human = false;
            }else{
                eprintln!("ERROR: could not find conclusive matches for human regex {} and mouse regex {} ",h_regex,m_regex);
                panic!();
            }

            // now since we know who is who and where, we can analyze these
            // and put it in our provided structure
            let integration_position : u64 ;
            let integration_sense : StrandDirection ;
            // now it can happen that we have 2 regions annotated where it maps
            // this is then obviously where it becomes complicated as we need to propagate
            // the information and can only later decide which one is key
            // FOR THE MOMENT THESE ARE IGNORED!!
            let pos_regex = Regex::new(r"@|-").unwrap();
            let coordinate : Vec<String> ;
            let mut compiled_info = AnnotatedRead {
                ref_start  : 0,
                ref_end    : 0,
                ref_strand : StrandDirection::Unknown,
                ref_chr    : String::from("Unknown"),
                ref_gstart : 0,
                ref_gend   : 0,
                read_name    : current[0].clone(),
                alien_start  : 0,
                alien_end    : 0,
                integration_site : 0,
                integration_dir : StrandDirection::Unknown,
            };
            // analyzed current one is actually human
            if current_human {
                let pos_infos : Vec<String> = current[5].split(',').map(|x| x.to_string()).collect();  
                //dbg!(&pos_infos);
                if pos_infos.len() == 1usize  {
                    coordinate = pos_regex.split(&pos_infos[0]).map(|x| x.to_string()).collect();
                    //dbg!(&coordinate);
                }else{
                    eprintln!("WARNING: Position information unclear!");
                    continue;
                }
                if coordinate.len() != 3usize {
                    panic!("ERROR: got more than 3 coordinates entries, expecting chromosome, start and end!");
                }
                let pos_chr   : &str = &coordinate[0];
                let pos_start : u64  = coordinate[1].parse::<u64>().unwrap();
                let pos_end   : u64  = coordinate[2].parse::<u64>().unwrap();
                //dbg!(&pos_infos);
                compiled_info.ref_chr      = pos_chr.to_string();
                compiled_info.ref_gstart   = pos_start;
                compiled_info.ref_gend     = pos_end;
                compiled_info.ref_start    = current[1].parse::<u64>().unwrap();
                compiled_info.ref_end      = current[2].parse::<u64>().unwrap();
                compiled_info.alien_start  = previous[1].parse::<u64>().unwrap();
                compiled_info.alien_end    = previous[2].parse::<u64>().unwrap();
                // the human one is located in the first part
                if current_a {
                    // now we know orientations and need
                    // to understand how where the integration position is
                    if current[4] == "+" {
                        // integration 3' of alignment end
                        compiled_info.ref_strand = StrandDirection::Fwd;
                        integration_sense        = StrandDirection::Fwd;        
                        integration_position     = compiled_info.ref_gend;                
                    } else if current[4] == "-" {
                        // integration 5' of alignment start
                        compiled_info.ref_strand = StrandDirection::Rev;
                        integration_sense        = StrandDirection::Rev;
                        integration_position     = compiled_info.ref_gstart;         
                    }else{
                        eprintln!("WARNING: Encountered elements which could not be differentiated via +/-");
                        continue
                    }
                }else{
                    // now we know orientations and need
                    // to understand how where the integration position is
                    if current[4] == "+" {
                        // integration 5' of alignment start
                        compiled_info.ref_strand = StrandDirection::Fwd;
                        integration_sense        = StrandDirection::Rev;
                        integration_position     = compiled_info.ref_gstart;         
                    } else if current[4] == "-" {
                        // integration 3' of alignment end
                        compiled_info.ref_strand = StrandDirection::Rev;
                        integration_sense        = StrandDirection::Fwd;
                        integration_position     = compiled_info.ref_gend;                
                    }else{
                        eprintln!("WARNING: Encountered elements which could not be differentiated via +/-");
                        continue
                    }
                }
                compiled_info.integration_dir  = integration_sense;
                compiled_info.integration_site = integration_position;

            }else if !current_human {
                let pos_infos : Vec<String> = previous[5].split(',').map(|x| x.to_string()).collect();  
                //dbg!(&pos_infos);
                if pos_infos.len() == 1usize  {
                    coordinate = pos_regex.split(&pos_infos[0]).map(|x| x.to_string()).collect();
                    //dbg!(&coordinate);
                }else{
                    eprintln!("WARNING: Position information unclear!");
                    continue;
                }
                if coordinate.len() != 3usize {
                    panic!("ERROR: got more than 3 coordinates entries, expecting chromosome, start and end!");
                }
                let pos_chr   : &str = &coordinate[0];
                let pos_start : u64  = coordinate[1].parse::<u64>().unwrap();
                let pos_end   : u64  = coordinate[2].parse::<u64>().unwrap();
                //dbg!(&pos_infos);

                compiled_info.ref_chr      = pos_chr.to_string();
                compiled_info.ref_gstart   = pos_start;
                compiled_info.ref_gend     = pos_end;
                compiled_info.ref_start    = previous[1].parse::<u64>().unwrap();
                compiled_info.ref_end      = previous[2].parse::<u64>().unwrap();
                compiled_info.alien_start  = current[1].parse::<u64>().unwrap();
                compiled_info.alien_end    = current[2].parse::<u64>().unwrap();

                if !current_a {
                    // this means the "previous" which is human will be a

                    // now we know orientations and need
                    // to understand how where the integration position is
                    if previous[4] == "+" {
                        // integration 3' of alignment end
                        compiled_info.ref_strand = StrandDirection::Fwd;
                        integration_sense        = StrandDirection::Fwd;
                        integration_position     = compiled_info.ref_gend;    
                    } else if previous[4] == "-" {
                        // integration 5' of alignment start
                        compiled_info.ref_strand = StrandDirection::Rev;
                        integration_sense        = StrandDirection::Rev;
                        integration_position     = compiled_info.ref_gstart;                
                    }else{
                        eprintln!("WARNING: Encountered elements which could not be differentiated via +/-");
                        continue
                    }
                }else{
                    // now we know orientations and need
                    // to understand how where the integration position is
                    if previous[4] == "+" {
                        // integration 5' of alignment start
                        compiled_info.ref_strand = StrandDirection::Fwd;
                        integration_sense        = StrandDirection::Rev;
                        integration_position     = compiled_info.ref_gstart;                
                    } else if previous[4] == "-" {
                        // integration 3' of alignment end
                        compiled_info.ref_strand = StrandDirection::Rev;
                        integration_sense        = StrandDirection::Fwd;
                        integration_position     = compiled_info.ref_gend;                

                    }else{
                        eprintln!("WARNING: Encountered elements which could not be differentiated via +/-");
                        continue
                    }
                }
                compiled_info.integration_dir  = integration_sense;
                compiled_info.integration_site = integration_position;          
            }
            collection.push(compiled_info);
        } 
    }
    collection
} 

/// this function is key for BND entries with an associated mate.
/// It will take the 2nd allele (multiple ones not supported so far)
/// and will deconstruct the information and return a structure with
/// the proper information organized.
/// Important: this will return though in a 0 based fashion the
/// coordinates to be easier compliant with the VCF rust function
/// which returns 1 based VCF positions in a 0 based fashion.
/// Will only work for entries with one additional BND allele
///
///Unittest: TRUE
///
/// ```rust
/// use genefusion::lib::common::{*};
/// pretty_env_logger::init();
/// 
/// let ff = String::from("G]chr16:18804692]");
/// let ff_check = Bnd2ndAllele {
///     nucleotide: String::from("G"),
///     chr: String::from("chr16"),
///     pos: 18804691_u64,
///     prim_direction: StrandDirection::Fwd,
///     second_direction: StrandDirection::Fwd,
/// };
/// let ff_result = deconstruct_2nd_allele(&ff,true);
/// assert_eq!(ff_result,ff_check);
/// ```
pub fn deconstruct_2nd_allele(
    allele_2:&str,
    paired: bool,
) -> Bnd2ndAllele {

    let new_allele = match paired {
        true => {
            let bracket1    : Vec<&str> = allele_2.split(']').collect();
            let bracket2    : Vec<&str> = allele_2.split('[').collect();
            let allele_info : &str;
            let allele_nuc  : &str;
            // now only one of the above should result into 
            // something with more than 0 entries
            let (forward1,forward2) = match bracket1.len(){
                3 => {
                    allele_info = bracket1[1];
                    // this can be now either FF or RF
                    if !bracket1[0].is_empty(){
                        // FF
                        allele_nuc  = bracket1[0];
                        (StrandDirection::Fwd,StrandDirection::Fwd)
                    }else{
                        // RF
                        allele_nuc  = bracket1[2];
                        (StrandDirection::Rev,StrandDirection::Fwd)
                    }
                },
                1 => {
                    allele_info = bracket2[1];

                    // this can be now either RR or FR
                    if !bracket2[0].is_empty(){      
                        // FR      
                        allele_nuc  = bracket2[0];
                        (StrandDirection::Fwd,StrandDirection::Rev)
                    }else{
                        // RR
                        allele_nuc  = bracket2[2];
                        (StrandDirection::Rev,StrandDirection::Rev)
                    }
                },

                _ => panic!("ERROR: could not determine direction of 2nd allele")
            };
            let location2 : Vec<&str>= allele_info.split(':').collect();
            if location2.len() != 2 { panic!("ERROR: allele2 splitting gave not chromosome and position!");}
            let tmp_chrom = location2[0];
            let tmp_pos    = (location2[1]).parse::<u64>().expect("ERROR: could not extract position correctly from 2nd allele!")-1;
            Bnd2ndAllele {
                nucleotide: allele_nuc.to_string(),
                chr: tmp_chrom.to_string(),
                pos: tmp_pos,
                prim_direction: forward1,
                second_direction: forward2,
            }
        },
        false => {
            let side    : Vec<&str> = allele_2.split('.').collect();
            let allele_nuc  : &str;
            let forward1 : StrandDirection;
            if side.len() != 2 {
                panic!("ERROR: single BND entry not correctly formatted!");
            }else if side[0].is_empty() {
                // Reverse
                allele_nuc  = side[1];
                forward1 = StrandDirection::Rev;
            }else if side[1].is_empty() {
                // Forwards
                allele_nuc  = side[0];
                forward1 = StrandDirection::Fwd;
            }else{
                panic!("ERROR: could not deconstruct single BND event!");
            };
            Bnd2ndAllele {
                nucleotide: allele_nuc.to_string(),
                chr: String::new(),
                pos: 0,
                prim_direction: forward1,
                second_direction: StrandDirection::Unknown,
            }
        }
    };
    new_allele

}

/// this function constructs the 2nd allele field for the VCF 4.2 specification
/// it is for single 2nd alleles only and one needs to verify beforehand that 
/// no multiallelic sites are being generated/investigated.
/// The input is a structure which holds the information necessary already in a
/// proper organized way.
/// It is though 0 based and the final VCF needs to be 1 based again
///
/// Unittest: TRUE
///
/// Example:
/// ```rust
/// use genefusion::lib::common::{*};
/// pretty_env_logger::init();
/// 
/// let ff = String::from("G]chr16:18804692]");
/// let ff_check = Bnd2ndAllele {
///     nucleotide: String::from("G"),
///     chr: String::from("chr16"),
///     pos: 18804691_u64,
///     prim_direction: StrandDirection::Fwd,
///     second_direction: StrandDirection::Fwd,
/// };
/// let ff_result = construct_2nd_allele(&ff_check,true);
/// assert_eq!(ff_result,ff);
/// ```
pub fn construct_2nd_allele(
    information: &Bnd2ndAllele,
    paired: bool,
) ->  String {
    match paired {
        true => {
            match (information.prim_direction,information.second_direction){
                (StrandDirection::Fwd,StrandDirection::Fwd) => format!("{}]{}:{}]",information.nucleotide,information.chr,information.pos+1),
                (StrandDirection::Fwd,StrandDirection::Rev) => format!("{}[{}:{}[",information.nucleotide,information.chr,information.pos+1),
                (StrandDirection::Rev,StrandDirection::Fwd) => format!("]{}:{}]{}",information.chr,information.pos+1,information.nucleotide),
                (StrandDirection::Rev,StrandDirection::Rev) => format!("[{}:{}[{}",information.chr,information.pos+1,information.nucleotide),
                _ => panic!("ERROR: unexpected direction pairs in paired BND events encountered!"),
            }
        },
        false => {
            match (information.prim_direction,information.second_direction){
                (StrandDirection::Fwd,StrandDirection::Unknown) => format!("{}.",information.nucleotide),
                (StrandDirection::Rev,StrandDirection::Unknown) => format!(".{}",information.nucleotide),
                _ => panic!("ERROR: unexpected direction pairs for single BND events encountered"),
            }
        }
    }
}


/// this function expects pairs of BND events
/// and will collapse them. It assures that the directions
/// match and will that way simplify the further analysis, sorting and
/// reporting of these events with a lower memory foot-pring
///
///Unittest: TRUE
///
pub fn bnd_collapse_pairs (
    cluster_pairs: FxHashMap<ClusterPairs,Vec<BNDentry>>
)->  FxHashMap<ClusterPairs,Vec<BNDentry>>{
    // collection of what we keep
    let mut collapsed_pairs : FxHashMap<ClusterPairs, Vec<BNDentry>> = FxHashMap::default();
    // so essentially we go through all and if one is 
    // found which is in its opposite direction present as well
    // then we add it to collapsed pairs.
    debug!("Collapsing following pairs: {:?}",&cluster_pairs);
    for (key, value) in &cluster_pairs {
        // here we generate the complementary element
        let comp_key = ClusterPairs{
            chr1    : key.chr2.clone(),
            range1  : key.range2.clone(),
            forward1: key.forward2,
            chr2    : key.chr1.clone(),
            range2  : key.range1.clone(),
            forward2: key.forward1,
        };
        // now we check if that complementary pair actually exits
        if cluster_pairs.contains_key(&comp_key){
            // now we want to add it, but only if the complementary 
            // is not already there, otherwise we would not collapse
            if !collapsed_pairs.contains_key(&comp_key){
                collapsed_pairs.insert(key.clone(),value.clone());
            }
        }else {
            eprintln!("for key: {:?} no matching comp_key {:?} ", &key,&comp_key);
        }
    }  
    debug!("Collappsed pairs {:?}",&collapsed_pairs);  
    collapsed_pairs
}


/// This function takes collapsed BND pairs
/// and writes them to STDOUT as BEDPE format (https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format).
/// This is a 0 based position system.
/// Since ranges somehow need always to extend beyond +1 we need here to correct the end
/// coordinates before printing
/// Name represents a ";"-separated collection of samples in which
/// said BND event appeared and the score represents the number of samples
/// that were observed
///
///Unittest: FALSE
///
pub fn bedpe_collapsed_pairs(
    collapsed_pairs:FxHashMap<ClusterPairs, Vec<BNDentry>>,
    multi:bool,
    purge:bool,
    out_vcf: Option<&str>
) -> Result<(), Box<dyn Error>> {

    //shamelessly taken 1:1 from https://stackoverflow.com/a/42216134/11255396
    let mut out_writer = match out_vcf {
        Some(x) => {
            let path = Path::new(x);
            Box::new(File::create(path).unwrap()) as Box<dyn Write>
        }
        None => Box::new(stdout()) as Box<dyn Write>,
    };

    for (key,value) in collapsed_pairs {
        debug!("Key {:?}, value {:?}",key,value);
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
        let (dir1,dir2) = match (key.forward1,key.forward2){
            (StrandDirection::Fwd,StrandDirection::Fwd) => ("+","+"),
            (StrandDirection::Rev,StrandDirection::Rev) => ("-","-"),
            (StrandDirection::Fwd,StrandDirection::Rev) => ("+","-"),
            (StrandDirection::Rev,StrandDirection::Fwd) => ("-","+"),
            (StrandDirection::Unknown,StrandDirection::Unknown) => (".","."),
            _=> panic!("ERROR: no mixed orientation models allowed!"),
        };
        match multi {
            false => {
                write!(out_writer,"{}\t{}\t{}\t{}\t{}\t{}\t",key.chr1,key.range1.start,key.range1.end-1,key.chr2,key.range2.start,key.range2.end-1).expect("ERROR: could not write bedpe ");
                let mut score = 0_i16;      
                for n in samples {
                    score+=1;
                    write!(out_writer,"{};",n).expect("ERROR: could not write bedpe ");
                }
                write!(out_writer,"\t{}",score).expect("ERROR: could not write bedpe ");
                write!(out_writer,"\t{}\t{}",dir1,dir2).expect("ERROR: could not write bedpe ");
                writeln!(out_writer,"\tBND_PAIR").expect("ERROR: could not write bedpe ");
            },
            true => {
                if samples.len() > 1 {
                    write!(out_writer,"{}\t{}\t{}\t{}\t{}\t{}\t",key.chr1,key.range1.start,key.range1.end-1,key.chr2,key.range2.start,key.range2.end-1).expect("ERROR: could not write bedpe ");
                    let mut score = 0_i16;
                    for n in samples {
                        score+=1;
                        write!(out_writer,"{};",n).expect("ERROR: could not write bedpe ");
                    }
                    write!(out_writer,"\t{}",score).expect("ERROR: could not write bedpe ");
                    write!(out_writer,"\t{}\t{}",dir1,dir2).expect("ERROR: could not write bedpe ");
                    writeln!(out_writer,"\tBND_PAIR").expect("ERROR: could not write bedpe ");
                }
            }
        }
    }      
    Ok(())
}


/// This function takes collapsed BND pairs
/// and writes them to STDOUT as BEDPE format (https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format).
/// This is a 0 based position system.
/// Name represents a ";"-separated collection of samples in which
/// said BND event appeared and the score represents the number of samples
/// that were observed
///
///Unittest: FALSE
///
pub fn bedpe_sv_clusters(
    sv_clusters: FxHashMap<SvClusters, Vec<INDELentry>>,
    multi:bool,
    purge:bool,
    type_sv: &str,
    out_vcf: Option<&str>
) {
    // shamelessly taken 1:1 from https://stackoverflow.com/a/42216134/11255396
    let mut out_writer = match out_vcf {
        Some(x) => {
            let path = Path::new(x);
            Box::new(File::create(path).unwrap()) as Box<dyn Write>
        }
        None => Box::new(io::stdout()) as Box<dyn Write>,
    };
    for (key,value) in sv_clusters {
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
        let (dir1,dir2) = match key.forward1{
            StrandDirection::Fwd => ("+","."),
            StrandDirection::Rev => ("-","."),
            StrandDirection::Unknown => (".","."),
        };
        match multi {
            false => {
                write!(out_writer,"{}\t{}\t{}\t.\t-1\t-1\t",key.chr,key.range.start,key.range.end,).expect("ERROR: could not write bedpe ");
                let mut score = 0_i16;      
                for n in samples {
                    score+=1;
                    write!(out_writer,"{};",n).expect("ERROR: could not write bedpe ");
                }
                write!(out_writer,"\t{}",score).expect("ERROR: could not write bedpe ");
                write!(out_writer,"\t{}\t{}",dir1,dir2).expect("ERROR: could not write bedpe ");
                writeln!(out_writer,"\t{}",type_sv).expect("ERROR: could not write bedpe ");
            },
            true => {
                if samples.len() > 1 {
                    write!(out_writer,"{}\t{}\t{}\t.\t-1\t-1\t",key.chr,key.range.start,key.range.end).expect("ERROR: could not write bedpe ");
                    let mut score = 0_i16;
                    for n in samples {
                        score+=1;
                        write!(out_writer,"{};",n).expect("ERROR: could not write bedpe ");
                    }
                    write!(out_writer,"\t{}",score).expect("ERROR: could not write bedpe ");
                    write!(out_writer,"\t{}\t{}",dir1,dir2).expect("ERROR: could not write bedpe ");
                    writeln!(out_writer,"\t{}",type_sv).expect("ERROR: could not write bedpe ");
                }
            }
        }
    }      
}


/// This functions takes cluster pairs and will refine them.
/// Essentially it will take of a pair the elements in a range and
/// re-define the range since in a first iteration it might have 
/// contained elements which were overlapping in the region but
/// did not share a common partner for the 2nd event location.
///
///Unittest: TRUE
///
pub fn bnd_filter_pairs(
    cluster_pairs: FxHashMap<ClusterPairs,Vec<BNDentry>>,
    range:  u64,
) ->  FxHashMap<ClusterPairs,Vec<BNDentry>> {
    let mut refined_pairs : FxHashMap<ClusterPairs, Vec<BNDentry>> = FxHashMap::default();
    debug!("Analyzing cluster pairs, length {:?}",cluster_pairs.len());
    for (mut key,value) in cluster_pairs {
            debug!("Key: {:?} value {:?}",key, value);
            // first we make again now a tree which contains only elements
            // within that range
            debug!("Length of value: {:?}",value.len());
            let tmp_tree : IntervalTree<u64,BNDentry> = value.iter()
                .map(|x| (
                    // needs to span before always
                    x.st1
                    ..
                    x.end1,
                    x.clone()
                )).collect::<IntervalTree<_, _>>();
            let local_max = key.range1.end;// + 1 ;
            debug!("Local tree {:?}",tmp_tree);
            // here we have now a collection of entries falling within our interval
            let tree_entries : Vec<Entry<'_, u64, BNDentry>> = tmp_tree.find(0..local_max).collect(); 
            debug!("tree entry within filter: {:?}, in total elements {:?}",tree_entries,tree_entries.len());
            // this tree we can use now to try building again a more
            // refined range
            let sub_points :  Vec<std::ops::Range<u64>> = locus_dist_based_groups(tree_entries,range);
            debug!("sub_points entry within filter: {:?}",sub_points);
            // And again it might now happen that we have within a large stretch
            // after the break-down >1 element (saw that in real data) and we 
            // need to treat each individually now.
            for sub_elements in sub_points {
                key.range1.start = sub_elements.start;
                key.range1.end   = sub_elements.end +1;
                // now we do the same for the element on the 2nd locus
                //debug!("Length of value: {:?}",value.len());
                let tmp_tree2 : IntervalTree<u64,BNDentry> = value.iter()
                    .map(|x| (
                        x.st2
                        ..
                        x.end2,
                        x.clone()
                    )).collect::<IntervalTree<_, _>>();
                // here we need +1 otherwise it always returns an empty one
                let local_max2 = key.range2.end + 1 ;
                debug!("Local tree2 {:?} with start : {:?} and end : {:?}",tmp_tree2,key.range2.start,key.range2.end);
                
                let tree_entries2 : Vec<Entry<'_, u64, BNDentry>> = tmp_tree2.find(0..local_max2).collect();
                debug!("tree2 entry within filter: {:?}, in total elements {:?}",tree_entries2,tree_entries2.len());

                // this tree we can use now to try building again a more
                // refined range
                let sub_points2 :  Vec<std::ops::Range<u64>> = locus_dist_based_groups(tree_entries2,range);
                debug!("sub_points2 entry within filter: {:?}",sub_points2);
                for sub_elements2 in sub_points2 {
                    key.range2.start =  sub_elements2.start;
                    key.range2.end   =  sub_elements2.end+1;
                    // okay, so now we have actually the problem that
                    // the values have to be again evaluated, too - not only the keys.
                    let mut modified_value : Vec<BNDentry> = Vec::new();
                    for bnd in &value {
                        if bnd.chr1 == key.chr1 && 
                            key.range1.contains(&(bnd.st1)) && 
                            key.range1.contains(&bnd.end1) && 
                            bnd.chr2 == key.chr2 && 
                            key.range2.contains(&(bnd.st2)) && 
                            key.range2.contains(&bnd.end2)
                        {
                            modified_value.push(bnd.clone());
                        }
                    }
                    if !modified_value.is_empty() {
                        debug!("Pushing new entry key: {:?} value {:?}",key,modified_value);
                        refined_pairs.insert(key.clone(), modified_value);
                    }else{
                        debug!("No new entry to push: {:?}",modified_value);
                    }
                }
                
            }
    }
    refined_pairs
}

/// A function which will split intervals into groups 
/// based on non-overlapping distances between ranges.
/// To optimize the usage, it assumes that a tree is organized
/// by chromosomes and verifies that the chromosomes exists + size matches.
/// In the example below we have a few ranges which should group
//  into 3 groups, 2 with 3 elements and one
/// with a single element
/// 
/// It returns a vector with ranges for each groups.
///
///Unittest: TRUE
/// 
/// Example:
/// 
/// ```bash
///  900--------------------1200
///    1000-----------1100
///    1000-----------------------1300
///                                                 1400---------1600
///                                                    1500------1600
///                                                    1500--------------1700
///                                    
/// 
///                                                         5000----5500
/// ```
pub fn locus_dist_based_groups (
    mut chr_tree :  Vec<Entry<'_, u64, BNDentry>>,
    clust_dist: u64
) -> Vec<std::ops::Range<u64>> {
    let mut chrom_breaks : Vec<std::ops::Range<u64>> = Vec::new();
    // initialize empty the start and end for sliding along
    // the chromosome
    let mut prev_end   : Option<u64>  = None;
    let mut prev_start : Option<u64 > = None;
    // now if there is nothing obtained for a chromosome 
    // we still keep them. Makes it easier later to assume
    // all chromosomes need to exist in all trees
    // sort all entries to allow "walking" along chromosome 
    chr_tree.sort_by(|x1, x2| x1.interval().start.cmp(&x2.interval().start));
    for element in chr_tree{
        // GET FIRST ELEMENT
        if prev_start.is_none(){
            prev_start = Some(element.interval().start );
        }
        if prev_end.is_none(){
            prev_end = Some(element.interval().end );
        }
        // now if the distance between 2 elements grow larger than specified
        // then we want to generate a groups entry and start again defining a new one
        if (
            element.interval().start > prev_end.unwrap()
        ) && (
            element.interval().start - prev_end.unwrap()
        ) > clust_dist {
            let new_range = std::ops::Range {
                start: prev_start.unwrap(),
                end: prev_end.unwrap(),
            };
            chrom_breaks.push(new_range);
            prev_start = Some(element.interval().start);
            prev_end   = Some(element.interval().end);
        // it can happen that the range is still within groups
        // but end is earlier than existing one due to fact that
        // it iterates based on start, not end
        }else if element.interval().end > prev_end.unwrap() {
            prev_end = Some(element.interval().end );
        }
        //println!("prev_Start {}, prev_End {}", prev_start.unwrap(),prev_end.unwrap());
    }
    // wrap up last element
    if prev_start.is_some() && prev_end.is_some() {
        let new_range = std::ops::Range {
            start: prev_start.unwrap(),
            end: prev_end.unwrap(),
        };
        //eprintln!("cluster start {} cluster end {}", new_range.start,new_range.end);
        chrom_breaks.push(new_range);
    }
    debug!("New chrom_breaks: {:?}",chrom_breaks);
    chrom_breaks
}





/// A function which will split intervals into groups 
/// based on non-overlapping distances between ranges.
/// To optimize the usage, it assumes that a tree is organized
/// by chromosomes and verifies that the chromosomes exists + size matches.
/// In the example below we have a few ranges which should group
//  into 3 groups, 2 with 3 elements and one
/// with a single element
/// Added a possibility to select based on the SVtype as well
///  
/// It returns a vector with ranges for each groups.
///
///Unittest: FALSE
/// 
/// Example:
/// 
/// ```bash
///  900--------------------1200
///    1000-----------1100
///    1000-----------------------1300
///                                                 1400---------1600
///                                                    1500------1600
///                                                    1500--------------1700
///                                    
/// 
///                                                         5000----5500
/// ```
pub fn locus_dist_based_groups_indel (
    mut chr_tree :  Vec<Entry<'_, u64, INDELentry>>,
    clust_dist: u64,
    selection: &Option<SVType>,
) -> Vec<std::ops::Range<u64>> {
    let mut chrom_breaks : Vec<std::ops::Range<u64>> = Vec::new();
    // initialize empty the start and end for sliding along
    // the chromosome
    let mut prev_end   : Option<u64>  = None;
    let mut prev_start : Option<u64 > = None;
    // now if there is nothing obtained for a chromosome 
    // we still keep them. Makes it easier later to assume
    // all chromosomes need to exist in all trees
    // sort all entries to allow "walking" along chromosome 
    chr_tree.sort_by(|x1, x2| x1.interval().start.cmp(&x2.interval().start));
    for element in chr_tree{
        // now we can actually skip if e.g. the demanded type is INS and we have an DEL instead
        if selection.is_some() && element.data().sv_type != selection.unwrap(){
            continue
        }
        // GET FIRST ELEMENT
        if prev_start.is_none(){
            prev_start = Some(element.interval().start );
        }
        if prev_end.is_none(){
            prev_end = Some(element.interval().end );
        }
        // now if the distance between 2 elements grow larger than specified
        // then we want to generate a groups entry and start again defining a new one
        // We add to the end position a +1 as the entry giving rise to that event is
        // otherwise later not caught with a range function
        if (
            element.interval().start > prev_end.unwrap()
        ) && (
            element.interval().start - prev_end.unwrap()
        ) > clust_dist {
            let new_range = std::ops::Range {
                start: prev_start.unwrap(),
                end: prev_end.unwrap(),
            };
            chrom_breaks.push(new_range);
            prev_start = Some(element.interval().start);
            prev_end   = Some(element.interval().end);
        // it can happen that the range is still within groups
        // but end is earlier than existing one due to fact that
        // it iterates based on start, not end
        }else if element.interval().end > prev_end.unwrap() {
            prev_end = Some(element.interval().end );
        }
        //println!("prev_Start {}, prev_End {}", prev_start.unwrap(),prev_end.unwrap());
    }
    // wrap up last element
    if prev_start.is_some() && prev_end.is_some() {
        let new_range = std::ops::Range {
            start: prev_start.unwrap(),
            end: prev_end.unwrap(),
        };
        //eprintln!("cluster start {} cluster end {}", new_range.start,new_range.end);
        chrom_breaks.push(new_range);
    }
    chrom_breaks
}




/// function which wraps the locus_dist_based_groups function
/// and executes it on a per chromosome base.
/// 
///Unittest: TRUE
///
pub fn chr_dist_based_groups(
    tree: &FxHashMap<String,IntervalTree<u64,BNDentry>>,
    clust_dist: u64,
    contigs: &FxHashMap<String, u64>
) -> FxHashMap<String,Vec<std::ops::Range<u64>>> {

    // now we generate a tree within a hashmap with the key being the chromosomes
    let mut break_points : FxHashMap<String,Vec<std::ops::Range<u64>>> = FxHashMap::default();
    for (chrom,value) in tree {
        // iterate in sorted (start) order 
        // and define groups
        // Therefore we need to know the max length of the chromosome
        let chrom_length = contigs.get(chrom).expect("ERROR: could not retrive chromosome in provided chromsome lengths!");
        // get all entries for that chromosome
        let chrom_entries : Vec<Entry<'_, u64, BNDentry>> = value.find(0..*chrom_length).collect();
        if chrom_entries.is_empty() {
            continue
        }
        debug!("Got the following chr entries: {:?}",chrom_entries);
        let chrom_breaks = locus_dist_based_groups(chrom_entries,clust_dist);
        //debug!("Got the chrom_breaks: {:?}",chrom_breaks);
        break_points.insert(chrom.to_string(), chrom_breaks);
    };
    break_points
}

/// function which wraps the locus_dist_based_groups function
/// and executes it on a per chromosome base.
/// NOTE: I tried to make a generic type with the above BNDentry similarly
/// but it broke in so many places downstream that they have currently
/// kind of duplicated functions with tiny differences.
/// Needs to be improved and consolidated in the future
/// The Selection field allows to return ranges for events with a 
/// certain type. This allows to generate one big tree but only
/// then filter for certain event types
///
///Unittest: FALSE
///
pub fn chr_dist_based_groups_indel(
    tree      : &FxHashMap<String,IntervalTree<u64,INDELentry>>,
    clust_dist: u64,
    contigs   : &FxHashMap<String, u64>,
    selection : Option<SVType>
) -> FxHashMap<String,Vec<std::ops::Range<u64>>> {

    // now we generate a tree within a hashmap with the key being the chromosomes
    let mut break_points : FxHashMap<String,Vec<std::ops::Range<u64>>> = FxHashMap::default();
    for (chrom,value) in tree {
        // iterate in sorted (start) order 
        // and define groups
        // Therefore we need to know the max length of the chromosome
        let chrom_length = contigs.get(chrom).expect("ERROR: could not retrive chromosome in provided chromsome lengths!");
        // get all entries for that chromosome
        let chrom_entries : Vec<Entry<'_, u64, INDELentry>> = value.find(0..*chrom_length).collect();
        // here we can now filter for our event type of interest
        let chrom_breaks = locus_dist_based_groups_indel(chrom_entries,clust_dist,&selection);
        break_points.insert(chrom.to_string(), chrom_breaks);
    };
    break_points
}


/// This function takes now our full interval tree with all events
/// together with the tree of our cluster intervals.
/// It then checks for each entry falling into a cluster interval where it's corresponding
/// cluster is located. It then returns a Hashmap containing as key a simple ClusterPair structure
/// which involves the chromosome, range and direction and reports all elements which were found in 
/// that particular pair.
/// 
/// Note: All events will be reported twice, once from each member of the pairs oriented.
/// This is by design to allow later to easily generate a matching BND ID for each in vcf entries
///
///Unittest: TRUE
///
pub fn bnd_pairs_clustered (
    full_interval_tree: &FxHashMap<String,IntervalTree<u64,BNDentry>>,
    cluster_tree      : &FxHashMap<String,IntervalTree<u64,String>>,
    contigs           : &FxHashMap<String, u64>
) -> FxHashMap<ClusterPairs,Vec<BNDentry>> {

    let mut range_pairs : FxHashMap<ClusterPairs,Vec<BNDentry>> = FxHashMap::default();
    // now we can use the generated cluster ranges to get the elements
    // of ranges which match each cluster
    for (chrom,value_cluster) in cluster_tree {
        let chrom_length = contigs.get(chrom).expect("ERROR: could not retrive chromosome in provided chromsome lengths!");
        // get all entries for that chromosome
        // and sort them
        let mut clusters : Vec<Entry<'_, u64, String>> = value_cluster.find(0..*chrom_length).collect();
        // here we can quit anything
        // directly for which no entries 
        // are at all present
        if clusters.is_empty(){
            continue
        }else{
            debug!("Chromosome {:?} value_cluster {:?}",chrom,value_cluster);
        }
        clusters.sort_by(|x1, x2| x1.interval().start.cmp(&x2.interval().start));
        // now we go over the cluster entries and get the matching
        // original entries again
        for breaks in clusters {
            let range  = &breaks.interval();
            // this is now a tree with entries for one cluster
            // again we sort first as by default sorting is not done
            let mut in_range : Vec<Entry<'_, u64, BNDentry>>  = full_interval_tree
                .get(chrom)
                .expect("ERROR: could not fetch the chromosome in full tree!")
                .find(range.start..range.end).collect();
            in_range.sort_by(|x1, x2| x1.interval().start.cmp(&x2.interval().start));   
            debug!("What is in range of break {:?}: {:?}",&range, &in_range);     
            // now thanks to the information in the BNDentry we can jump directly to the 
            // cluster which should contain it's counter-part
            for entry in in_range {
                // now we need to simply lookup in which cluster of break-tree our complementary
                // entry falls 
                //let compl_range    = std::ops::Range{ start: entry.value.st2, end: entry.value.end2 };
                debug!("Trying entry {} with start : {} and end : {}",entry.data().chr2,entry.data().st2,entry.data().end2);
                let compl_in_range = cluster_tree
                    .get(&entry.data().chr2)
                    .expect("ERROR: could not fetch the chromosome in break tree!")
                    .find(
                        std::ops::Range{ 
                            start: entry.data().st2, 
                            end: entry.data().end2 
                        }
                        .clone()
                    );
                // so this is now an iterator which must contain only 1 entry
                let mut z = 0;
                for cluster_range in compl_in_range {
                    debug!("Found for entry {:?} a complementary range {:?}", entry ,cluster_range);
                    z+=1;
                    if z > 1 {
                        panic!("ERROR: there should be only 1 matching cluster found!");
                    }
                    // here we generate now our key which is actually
                    // a structure which identifies the chromosome and the
                    // ranges of the **cluster** as well as it's orientation.
                    // This should help to separate as well events if the direction
                    // is not identical.
                    let pair_entry = ClusterPairs {
                        chr1      : entry.data().chr1.clone(),
                        range1    : range.start..range.end,
                        forward1  : entry.data().forward1,
                        chr2      : entry.data().chr2.clone(),
                        range2    : cluster_range.interval().start..cluster_range.interval().end,
                        forward2  : entry.data().forward2,
                    };
                    if range_pairs.contains_key(&pair_entry){
                        range_pairs.get_mut(&pair_entry).unwrap().push(entry.data().clone());
                    }else{
                        range_pairs.insert(pair_entry,vec![entry.data().clone()]);
                    }
                }
            }
        }
    }
    range_pairs
}



/// This function takes now our full interval tree with all events
/// together with the tree of our cluster intervals.
/// It then checks for each entry falling into a cluster interval where it's corresponding
/// cluster is located. It then returns a Hashmap containing as key a simple ClusterPair structure
/// which involves the chromosome, range and direction and reports all elements which were found in 
/// that particular pair.
/// 
/// Note: All events will be reported twice, once from each member of the pairs oriented.
/// This is by design to allow later to easily generate a matching BND ID for each in vcf entries
///
///Unittest: FALSE
///
pub fn indels_clustered (
    full_interval_tree: &FxHashMap<String,IntervalTree<u64,INDELentry>>,
    cluster_tree      : &FxHashMap<String,IntervalTree<u64,String>>,
    contigs           : &FxHashMap<String, u64>,
    selection         : &Option<SVType>,
) -> FxHashMap<SvClusters,Vec<INDELentry>> {

    let mut compiled_clusters : FxHashMap<SvClusters,Vec<INDELentry>> = FxHashMap::default();
    // now we can use the generated cluster ranges to get the elements
    // of ranges which match each cluster
    for (chrom,value_cluster) in cluster_tree {
        let chrom_length = contigs.get(chrom).expect("ERROR: could not retrive chromosome in provided chromsome lengths!");
        // get all entries for that chromosome
        // and sort them
        let mut clusters : Vec<Entry<'_, u64, String>> = value_cluster.find(0..*chrom_length).collect();
        clusters.sort_by(|x1, x2| x1.interval().start.cmp(&x2.interval().start));
        // now we go over the cluster entries and get the matching
        // original entries again
        for breaks in clusters {
            //dbg!(&breaks);
            let range  = &breaks.interval();
            // this is now a tree with entries for one cluster
            // again we sort first as by default sorting is not done
            let mut in_range : Vec<Entry<'_, u64, INDELentry>>  = full_interval_tree
                .get(chrom)
                .expect("ERROR: could not fetch the chromosome in full tree!")
                .find(range.start..range.end).collect();
            in_range.sort_by(|x1, x2| x1.interval().start.cmp(&x2.interval().start));   
            for entry in in_range {
                if selection.is_some() && entry.data().sv_type != selection.unwrap() {
                    continue
                }
                // here we generate now our key which is actually
                // a structure which identifies the chromosome and the
                // ranges of the **cluster** as well as it's orientation.
                // The latter is unknown for INDELS but we might want to 
                // use the SvClusters structure for single BND events, thats
                // why we keep the field
                let sv_entry = SvClusters {
                    chr: entry.data().chr.clone(),
                    range: 
                        range.start
                        ..
                        range.end,
                    forward1: StrandDirection::Unknown,
                    svtype: entry.data().sv_type,
                };
                if compiled_clusters.contains_key(&sv_entry){
                    compiled_clusters.get_mut(&sv_entry).unwrap().push(entry.data().clone());
                }else{
                    compiled_clusters.insert(sv_entry,vec![entry.data().clone()]);
                }
            }
        }
    }
    compiled_clusters
}


/// This simply takes a vector of ranges and translates it into an interval Tree.
/// It adds as a value "cluster" -name 
/// 
///Unittest: TRUE
///
pub fn cluster_to_tree (
    ranges:Vec<std::ops::Range<u64>>
) -> IntervalTree<u64, String> {
    let mut out : IntervalTree<u64, String> = IntervalTree::new();
    let mut n = 0_u32;
    for element in &ranges {
        n+=1;
        out.insert(element.clone(),  format!("Cluster{}",n));
    };
    out
}

/// Calculates the median of u64 values
///
/// Unittest: TRUE
///
/// ```
/// use crate::genefusion::lib::common::median_u64;
/// let mut test = vec![5_u64,6_u64,10_u64,20_u64,3_u64];
/// assert_eq!(median_u64(&mut test),6);
/// ```
pub fn median_u64(
    numbers: &mut std::vec::Vec<u64>
) -> u64 {
    numbers.sort_unstable();
    let mid = numbers.len() / 2;
    numbers[mid]
}


/// The supplementary alignment (SA) field of the long read alignments
/// are really the core to determine any re-arrangements. This function
/// takes the SA field from the BAM alignment as well as the query name 
/// and returns a vector with the alignment information, similarly to 
/// primary alignment information.
/// Again we need to be prudent and aware that SA derived Cigars have a limited
/// information content compared to primary alignment derived Cigars.
/// Note: for speed reasons there are occasions where we dont want
/// to have the CIGAR fully evaluated. Therefore the sequence length is an
/// option 
/// The coordinates of the SA entry are as well 1 based, not 0 based.
/// We correct this in order to be coherent everywhere
/// 
///Unittest: FALSE
///
pub fn eval_supp_align(
    supp: &str, 
    query_name: &str,
    seq_length: Option<&u64>,
    // maximum allowed ratio of short to longer clipping
    // in case that both sides are clipped
    // If none, the larger one is used to define the breakpoint
    ratio: Option<f32>
) -> Vec<AlignInfos> {
    let mut all_sa = Vec::new();
    let mut subresult: Vec<&str> = supp.split(';').collect();
    // as it generates a last empty element we pop it
    subresult.pop();
    for entries in subresult.iter() {
        debug!("Potential supp alignment entries: {:?}",&entries);
        let overview_match: Vec<&str> = entries.split(',').collect();
        //overview_match.pop();
        let target = overview_match[0];
        let start = match overview_match[1].parse::<u64>() {
            Ok(v) => v-1,
            _ 	  => panic!(
                "ERROR: could not parse SA field of {} --> {} --> {}",
                &supp, &entries, overview_match[1]
            ),
        };
        let strand: StrandDirection = match overview_match[2] {
            "+" => StrandDirection::Fwd,
            "-" => StrandDirection::Rev,
            _   => panic!(
                "ERROR: could not parse SA field of {} --> {} --> {}",
                &supp, &entries, overview_match[2]
            ),
        };
        let cigar = overview_match[3];
        debug!("Supp alignment, cigar: {}", &cigar);
        let mismatch = match overview_match[5].parse::<u64>() {
            Ok(v) => v,
            _     => panic!(
                "ERROR: could not parse SA field of {} --> {} --> {}",
                &supp, &entries, overview_match[5]
            ),
        };
        ///////////////////////////
        //// correct SA cigar /////
        ///////////////////////////
        //
        // This is really stupid because the position
        // is in SA correct but the representation often not
        // In that case we would loose the match because similarity
        // would be really low (e.g. large D at end before H)
        // So we do both, a correction for the similarity and
        // positions for getting the positions are from original cigar
        let new_cigar = cigar_fix(cigar,seq_length);
        debug!("Fixed cigar from originally {} to now {}",&cigar,&new_cigar);

        let cig_info = eval_cigar(
            &new_cigar,
            mismatch,
            &start,
            &strand,
            ratio
        );
        debug!("Parsed cigar and got: {:?}", &cig_info);
        match cig_info {
            None => debug!("Got an none CIGAR back"),
            Some(x) => {
                let sa_info  = AlignInfos {
                    query:          query_name.to_string(),
                    q_length:       *seq_length.unwrap_or(&0_u64),
                    cigar:          cigar.to_string(),
                    chrom:          target.to_string(),
                    fp_target:      x.fp_target,
                    upstream:       x.upstream,
                    downstream:     x.downstream,
                    similarity:     x.similarity,
                    seq_length    : x.seq_length,
                    match_length_s: x.match_length_s,
                    match_length:   x.match_length,
                    fp_query:       x.fp_query,
                    aln_start:      x.aln_start,
                    aln_end:        x.aln_end,
                    match_strand:   x.match_strand,
                };
                debug!("Gathered information: {:?}", &sa_info);
                // gathering all SA info we can
                // get situation where we could not determine a potential
                // fusion point, remove these
                all_sa.push(sa_info);
            },
        }
    }
    all_sa
}



/// This function takes a file an returns a hashmap
/// with chromosomes as key and value being vectors
/// with structure as entries. It requests a filter on
/// a certain type of annotation e.g. "gene", "transcript" or "exon"
/// and will otherwise default to "transcript".
/// This was chosen as it becomes otherwise really extensive and that might change
/// in the future to accomodate as well combinations of features.
/// 
/// It expects a valid GTF file with header and 9 columns and will otherwise panic
/// I added the possibility to add chromosomes of interest to speed up reading. 
/// In that case it will skip any chromosomes which is absent in that list
/// and allows to read e.g. first the BAM file and only generate annotation for
/// chromosomes of interest
///
/// *Important*: GTF is normally 1 based and we are working 0 based
/// Therefore I convert here everything to 0 based as well
///
///Unittest: FALSE
///
pub fn parse_gtf_annotation(
    my_path: &str, 
    my_filter: Option<&str>, 
    my_chromosome: Option<Vec<String>>
) -> FxHashMap<String, Vec<Gtf>> {

    eprintln!("INFO: annotation {} provided, parsing entries...", my_path);
    assert!(
        Path::new(my_path).exists(),
        "ERROR: annotation file {:?} does not exists!",
        my_path
    );
    let mut gtf: FxHashMap<String, Vec<Gtf>> = FxHashMap::default();
    let input  = File::open(my_path).expect("Unable to open annotation file");
    let reader = BufReader::new(input);
    let reg_header = Regex::new(r"^#.*").unwrap();
    // now we read line by line, skip header and push
    // for each chromosome a vector with our structure
    let mut gff_system  = GtfSource::UNKNOWN;
    let mut ignored_lines : i32 = 0;
    for (line_number,line) in reader.lines().enumerate() {
        let l = line.expect("Unable to read line");
        // if we have a header or comments skip
        if reg_header.is_match(&l) {
            ignored_lines+= 1;
            continue;
        }
        // now we split our GTF by tab and fill a vector with it
        let elements: Vec<&str> = l.split('\t').collect();
        // if we have not 9 elements in vector then GTF is not valid
        if elements.len() != 9 {
            panic!("ERROR: GTF does not contain 9 fields ");
        }
        // the strand is normally encoded in "+" and "-"
        // which I prefer in 1 and -1 and easier to work with
        let direction: StrandDirection = match elements[6] {
            "+" => StrandDirection::Fwd,
            "-" => StrandDirection::Rev,
            _ => panic!("ERROR: your stand information is not valid!"),
        };
        let mut annot_struct = Gtf {
            chromosome: elements[0].to_string(),
            source:     elements[1].to_string(),
            feature:    elements[2].to_string(),
            start:      elements[3].parse::<u64>().unwrap() -1 ,
            end:        elements[4].parse::<u64>().unwrap(),
            score:      elements[5].to_string(),
            strand:     direction,
            frame:      elements[7].to_string(),
            attribute:  FxHashMap::default(),
        };
         // As we provide the possibility to filter as well for
        // feature to reduce size we need to check if the option is
        // given or not, if match we can skip directly to next one
        let feature = my_filter.unwrap_or("transcript");
        if feature != elements[2] {
            ignored_lines+=1;
            continue;
        };
        // so this we want later to treat correctly
        // but it is a bit intensive so we first check if 
        // we actually keep that entry
        let tmp_attribute = elements[8].to_string();
 
       
        // now everything is controled and we can
        // get our subset which we keep in this particular program
        // we want to capture either the "gene_name or gene_id"

        // the only GFF/GTF rule for that last element
        // is that it has to be ; separated, regex is slow
        // so lets put it in a hash and see if key exists
        let mut tmp_hash : FxHashMap<String,String> = FxHashMap::default();
        let set1    = Regex::new(r#"gene_name "{1}.*?"{1}"#).unwrap();
        let set2    = Regex::new(r#"gene_name="#).unwrap();
        // in line 1 we still need to figure out what is the
        // underlying structure - in line 2 it needs to be determined or
        // we panic
        //eprintln!("Line number {:?}", &ignored_lines);
        if line_number == ignored_lines as usize {
            for s in tmp_attribute.split(';'){
                //let id = cleaner.replace_all(&id, "").to_string();
                // next we want a key : value relationship
                // there are in principle 2 possibilities I want to support
                // gencode or CHESS
                // gencode is: gene_name "myGene"
                // CHESS is  : gene_name=myGene
                // First time encounter check what can be found
                // and set divider accordingly
                match &gff_system {
                    GtfSource::UNKNOWN => {               
                        if set1.is_match(s){
                            gff_system  = GtfSource::GENCODE;
                            eprintln!("INFO: GTF identified as GENCODE-like");
                        }else if set2.is_match(s) {
                            gff_system = GtfSource::CHESS;
                            eprintln!("INFO: GTF identified as CHESS-like");
                        }
                    },
                    GtfSource::GENCODE => (),
                    GtfSource::CHESS   => (),
                };

            }
        }else if let GtfSource::UNKNOWN = gff_system {panic!("ERROR: could not determine your GTF style!")};

        // Here we guarantee now that the chromosome
        // is actually of interest
        if my_chromosome.is_some() {
            let mut encountered = false;
            for chr in my_chromosome.as_ref().unwrap().iter() {
                if chr == &annot_struct.chromosome {
                    encountered = true;
                    continue
                }
            }
            if !encountered {
                //eprintln!("Chr {} not in list of chromosomes, skipping", &annot_struct.chromosome);
                continue
            }
        }
        if annot_struct.end < annot_struct.start {
            panic!("ERROR: start position higher then end position!");
        }
        // now we can start splitting our attribute fields
        match gff_system {
            GtfSource::UNKNOWN => panic!("ERROR: could not determine your GTF style!"),
            GtfSource::GENCODE => {
                for s in tmp_attribute.split(';'){
                    let splited : Vec<&str> = s.split(r#" ""#).collect();
                    if splited.len() != 2 {
                        continue
                    }
                    let new_key   = splited[0].to_string().trim_end_matches(' ').trim_start_matches(' ').to_string();
                    let new_value = splited[1].to_string().trim_end_matches('"').trim_end_matches(' ').trim_start_matches(' ').to_string();
                    tmp_hash.insert(new_key,new_value);
                }
            },
            GtfSource::CHESS   => {
                for s in tmp_attribute.split(';'){
                    let splited : Vec<&str> = s.split('=').collect();
                    if splited.len() != 2 {
                        continue
                    }
                    let new_key   = splited[0].to_string().trim_end_matches(' ').trim_start_matches(' ').to_string();
                    let new_value = splited[1].to_string().trim_end_matches(' ').trim_start_matches(' ').to_string();
                    tmp_hash.insert(new_key,new_value);
                }
            },
        }
        annot_struct.attribute = tmp_hash;
        //annot_struct.name = id;
        // now we add entries to a hashmap which
        // has the chromosome as key and values are
        // a vector where each entry of the vector is the above
        // gtf structure
        gtf.entry(annot_struct.chromosome.to_string())
            .or_default()
            .push(annot_struct);
    }
    gtf
}


/// here we do compressed gap similarity
/// We count the number and length of INDELS
/// and use it together with information of mis-matches.
/// The latter can come originally from both
/// the CIGAR or mismatch tag.
/// Matches have to be collection of matches + mismatches
/// 
///Unittest: TRUE
///
/// ```rust
/// use genefusion::lib::common::{*};
/// pretty_env_logger::init();
/// let cigar = String::from("18M3D2M2D2M1I22M");
/// let sim   = calc_gapcompressed_sim(&cigar,44,&1);
/// let truth = 91.489365_f32;
/// assert_eq!(truth,sim); 
/// 
/// let cigar2 = String::from("544S692M1190N1796M");
/// let sim2   = calc_gapcompressed_sim(&cigar2,2488,&0);
/// let truth2 = 100_f32;
/// assert_eq!(truth2,sim2); 
/// ```
pub fn calc_gapcompressed_sim(
    ciggi: &str, 
    matches: u32, 
    mismatch: &u64
) -> f32 {

    let mut open: u32 = 0;
    let mut crap: u32 = 0;
    let mut full: u32 = 0;
    let re_crap = Regex::new(r"([0-9]+)[D|I]").unwrap();
    for m in re_crap.captures_iter(ciggi) {
        crap += &m[1].parse().unwrap();
        open += 1;
    }
    let re_length = Regex::new(r"([0-9]+)[M|D|I|X|=]").unwrap();
        for m in re_length.captures_iter(ciggi) {
            full += &m[1].parse().unwrap();
        }
    // IdentityGapComp = 100.0 * match / (match + mismatch + delEvents + insEvents);
    debug!("We have full: {}, crap {} , opening {}, matches {} , mismathes {}",&full,&crap,&open,&matches,&mismatch);
    (100f64 * (matches as f64 - *mismatch as f64)
        / (full as f64 - open as f64)) as f32 
}


/// it can happen in SA alignment that for the cigar
/// the alignment is unprecise and we observe e.g. 500M10D50H
/// Obviously a D prior a clip makes no match at all and this needs
/// to be fixed in that case (as far as I can tell only trailing)
/// I discovered that the following occuring unsetteling scenarios:
/// - CIGAR which ends with a deletion + soft/hard-clipe (see issue)[https://github.com/lh3/minimap2/issues/524]
/// - CIGAR which has a clip, larget than the entire sequence "chr2,68346460,-,2172S830M37403D,60,73" for a 3002bp long seq
/// 
/// Unittest: TRUE
///
/// ```rust
/// use genefusion::lib::common::{*};
/// pretty_env_logger::init();
/// 
/// let test1 = String::from("544S692M1190N1796M");
/// let res1  = cigar_fix(&test1,Some(&3032_u64));
/// assert_eq!(test1,res1);
/// 
/// let test2  = String::from("545M245D2487S");
/// let res2   = cigar_fix(&test2,Some(&3277_u64));
/// let truth2 = String::from("545M2732S");
/// assert_eq!(truth2,res2);
///
/// let test3  = String::from("2172S830M37403D");
/// let res3   = cigar_fix(&test3,Some(&3002_u64));
/// let truth3 = String::from("2172S830M");
/// assert_eq!(truth3,res3);
/// ```
pub fn cigar_fix(
    ciggi: &str,
    seq_length: Option<&u64>,
) -> String {

    let false_end = Regex::new(r"(?P<c>.*?)(?P<d>[0-9]+)D(?P<h>[0-9]+)(?P<k>[H|S])$").unwrap();
    let mut tmp: String = String::from(ciggi);
    if false_end.is_match(ciggi) {
        // since we check already if match exist this can be done directly
        let tmp_match = false_end.captures(ciggi).unwrap();
        let nw_cif = tmp_match.name("c").unwrap().as_str();
        let d_n = tmp_match
            .name("d")
            .unwrap()
            .as_str()
            .parse::<u32>()
            .unwrap();
        let h_n = tmp_match
            .name("h")
            .unwrap()
            .as_str()
            .parse::<u32>()
            .unwrap();
        let k_t = tmp_match.name("k").unwrap().as_str();
        let n_n = (d_n + h_n).to_string();
        // okay now we modified our strings and we want to pass now the new cigar back
        let x_t = format!("{}{}{}", nw_cif, n_n, k_t);
        tmp = x_t.as_str().to_owned();
    };
    // Found as well situation with splicing where we have a 206S1118M1488D
    // This is similarly stupid as we have immediately a match similarity of <50%
    let false_end2 = Regex::new(r"(?P<c>.*?)(?P<d>[0-9]+)D$").unwrap();
    if false_end2.is_match(ciggi) {
        // since we check already if match exist this can be done directly
        let tmp_match = false_end2.captures(ciggi).unwrap();
        let nw_cif = tmp_match.name("c").unwrap().as_str();
        let d_n = tmp_match
            .name("d")
            .unwrap()
            .as_str()
            .parse::<u32>()
            .unwrap();
        // okay now we modified our strings and we want to pass now the new cigar back
        let x_t = format!("{}{}H", nw_cif, d_n);
        tmp = x_t.as_str().to_owned();
    };
    // finally we can have situations where the CIGAR is way longer than the real
    // sequence length and it seems then it is the end again which is messed up
    
    ///////////////////////////////
    ////     length of seq   //////
    ///////////////////////////////
    let mut sa_length : u64 = 0;
    let re_length = Regex::new(r"([0-9]+)[H|S|M|I|X|=]").unwrap();
    for m in re_length.captures_iter(&tmp) {
        sa_length += &m[1].parse().unwrap();
    }
    // so in this case the safest (not a good) option is to remove clipping or deletions at
    // the end 
    if seq_length.is_some() && &sa_length  > seq_length.unwrap() {
        debug!("Wrong SA cigar with length {} for original {} and cigar {}",&sa_length,&seq_length.unwrap(),&tmp);
        let rem_end = Regex::new(r"(?P<c>.*?)(?P<d>[0-9]+)H$").unwrap();
        if rem_end.is_match(&tmp) {
            let tmp_match = rem_end.captures(&tmp).unwrap();
            tmp = tmp_match.name("c").unwrap().as_str().to_owned();
        }
    }
    tmp
}


/// this functions turns information about
/// the clipping and the length + start/end 
/// of matches into valid information about 
/// the breakpoints on a genome
/// Coordinates are 0 based
/// 
/// Just a few unit-tests below to test some sanity scenarios
/// In the example the length of clipping is larger than
/// the entire sequence and willr result in None
///
/// Unittest: TRUE
///
/// ```rust
/// use genefusion::lib::common::{*};
/// pretty_env_logger::init();
/// 
/// let test1  = calc_align_properties(7554_u32,2850_u32,3000_u32,5050_u32,7554,StrandDirection::Unknown,None);
/// let truth1 =  None;
/// assert_eq!(test1,truth1);
/// ```
pub fn calc_align_properties(
    // hard- or soft-clipping on 5' end
    clip5: u32,
    // hard- or soft-clipping on 3' end
    clip3: u32,
    // length of the query sequence
    slength: u32,
    // length of the match including splicing
    mlength_spliced: u32,
    // alignment start on target, 0 based
    start: u64,
    // alignment direction on target
    strand: StrandDirection,
    // maximum allowed ratio of short to longer clipping
    // in case that both sides are clipped
    // If none, the larger one is used to define the breakpoint
    ratio: Option<f32>,
) -> Option<BreakPoints> {

    // since our positions are 0 based and not 1 based, we need to account for that
    // e.g. if alignment start is 0 and end is 49, length is still 50
    // so everything which uses either slength or mlength_spliced  and clip for coordinates
    // needs this fix

    // We make them immutable to avoid
    // redefinition which must not happen
    // here normall
    // We further define them as option as
    // it can happen that we dont satisfy our
    // conditions and cant determine the
    // breakpoints

    let mut break_target: Option<u64> = None;
    let mut break_query:  Option<u32> = None;
    let mut align_start:  Option<u64> = None;
    let mut align_end:    Option<u64> = None;

    ///////////////////////
    //// FP on targets ////
    ///////////////////////
    // so here I had in a rare situation following
    debug!("start {}, slength {}, mlength_spliced {}, clip5 {}, clip3 {}, ratio {:?}",&start,&slength,&mlength_spliced,&clip5,&clip3,&ratio);

    // so we need to assure that the clips are not more than the length!
    if (start + slength as u64 > clip5 as u64 + clip3 as u64) &&  
        (
            // simple case 
            (clip5 != 0 && clip3 == 0) | 
            // no ratio cut-off, just take the larger one
            (ratio.is_none() && clip5 > clip3) | 
            // ratio cut-off, verify if matches
            (ratio.is_some() && clip5 > clip3 && (clip3 as f32/clip5 as f32)  < ratio.unwrap())
        ) {
        break_target = Some(start);
        align_start  = Some(start);
        align_end    = Some(start + (mlength_spliced as u64 ) -1);
    } else if (start + slength  as u64 > clip5 as u64 + clip3 as u64) &&  
        (
            (clip5 == 0 && clip3 != 0) | 
            (ratio.is_none() && clip5 < clip3) | 
            (ratio.is_some() && clip5 < clip3 && (clip5 as f32/clip3 as f32) < ratio.unwrap() )
        ) {
        break_target = Some(start + (mlength_spliced as u64 ) -1);
        align_start  = Some(start);
        align_end    = Some(start + (mlength_spliced as u64 ) -1);
    } else {
        debug!("WARNING: No break_target, align_start and end defined, skipped read!")
    };
    debug!("Now we have  align_start {:?}, align_end {:?}, break_target {:?}",&align_start,&align_end,&break_target);
    /////////////////////
    //// FP on query ////
    /////////////////////

    if ((clip5 != 0)    & (clip3 == 0) & (strand == StrandDirection::Fwd)) | 
	((clip5 > clip3) & (strand == StrandDirection::Fwd)) {
        break_query = Some(clip5 );
        debug!("Encountered case {}",1);
    } else if ((clip5 == 0) & (clip3 != 0) & (strand == StrandDirection::Fwd)) | 
	((clip5 < clip3) & (strand == StrandDirection::Fwd)) 
	{
        break_query = Some(slength - clip3 -1);
        debug!("Encountered case {}",2);
    } else if ((clip5 != 0) & (clip3 == 0) & (strand == StrandDirection::Rev))
        | ((clip5 > clip3)  & (strand == StrandDirection::Rev))
    {
        break_query = Some(slength - clip5 -1);
        debug!("Encountered case {}",3);
    } else if ((clip5 == 0) & (clip3 != 0) & (strand == StrandDirection::Rev))
        | ((clip5 < clip3 ) & (strand == StrandDirection::Rev))
    {
        break_query = Some(clip3);
        debug!("Encountered case {}",4);
    } 
    else 	{
        debug!("WARNING: No break query found")
    };
    
    match (break_target,break_query,align_start,align_end) {
        (Some(x),Some(y),Some(z),Some(n)) => Some(BreakPoints {
                        // this is special, if the field name 
                        // equals variable name we can simply 
                        // write it like that
                        break_query:  y,
                        break_target: x,
                        aln_start:    z,
                        aln_end:      n,
                        strand
                }       
            ),
        _ => None
    }
}

/// This function evaluates the cigar of an alignment and calculates
/// a couple of information based on its content. It can evaluate directly 
/// mismatches if they are indicated in a "=" and "X" fashion in the cigar ,
/// otherwise it needs the information of number of mismatches from the 
/// according tag in the BAM alignment provided
/// Positions are 0-based
/// 
/// So if we take e.g. this excelent resource on match similarity 
/// from [Heng li blog](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity)
/// 18M3D2M2D2M1I22M
/// 
/// Note: as the cigar in the blog does not use the X for mismatches
/// we provide a mismatch for the function itself 
///
/// Unittest: TRUE
/// 
/// Example:
///
/// ```rust
/// use genefusion::lib::common::{*};
/// pretty_env_logger::init();
///
/// let cigar4            = String::from("50M25S");
/// let result4           = eval_cigar(&cigar4,0_u64,&0_u64,&StrandDirection::Fwd,None);
/// let mut truth4        = CigarInfos::default();
/// truth4.fp_target      = 49_u64;
/// truth4.aln_start      = 0_u64;
/// truth4.aln_end        = 49;
/// truth4.downstream     = 25_u32;
/// truth4.seq_length     = 75_u32;
/// truth4.match_length   = 50_u32;
/// truth4.match_length_s = 50_u32;
/// truth4.similarity     = 100_f32;
/// truth4.fp_query       = 49_u32;
/// truth4.match_strand   = StrandDirection::Fwd;
/// assert_eq!(result4.unwrap(),truth4);
///
/// let cigar5            = String::from("50S25M");
/// let result5           = eval_cigar(&cigar5,0_u64,&49_u64,&StrandDirection::Rev,None);
/// let mut truth5        = CigarInfos::default();
/// truth5.fp_target      = 49_u64;
/// truth5.aln_start      = 49_u64;
/// truth5.aln_end        = 73;
/// truth5.upstream       = 50_u32;
/// truth5.seq_length     = 75_u32;
/// truth5.match_length   = 25_u32;
/// truth5.match_length_s = 25_u32;
/// truth5.similarity     = 100_f32;
/// truth5.fp_query       = 24_u32;
/// truth5.match_strand   = StrandDirection::Rev;
/// assert_eq!(result5.unwrap(),truth5);
///
/// ```
pub fn  eval_cigar(
    // the cigar which can include as well X for mismatches
    ciggi: &str,
    // number of mismatches (might be needed if cigar contains not info)
    mut seq_mmatch: u64,
    // 0 based alignment position 
    match_start: &u64,
    // direction of alignment
    match_strand: &StrandDirection,
    // maximum allowed ratio of short to longer clipping
    // in case that both sides are clipped
    // If none, the larger one is used to define the breakpoint
    ratio: Option<f32>,
) -> Option<CigarInfos> {

    let mut infos = CigarInfos::default();
    
    // just double-check again that we have not an
    // unmapped read
    if ciggi != "*" {
        let mut seq_match: u32 = 0;

        ///////////////////////////
        ////  clipping infos //////
        ///////////////////////////
        let up_match   = Regex::new(r"^([0-9]+)[H|S].*").unwrap();
        let down_match = Regex::new(r".*?([0-9]+)[H|S]$").unwrap();

        for m in up_match.captures_iter(ciggi) {
            infos.upstream = m[1].parse().unwrap();
        }
        for m in down_match.captures_iter(ciggi) {
            infos.downstream = m[1].parse().unwrap();
        }

        ///////////////////////////
        ////   # of matches  //////
        ///////////////////////////
        let re_match  = Regex::new(r"([0-9]+)M").unwrap();
        let re_match2 = Regex::new(r"([0-9]+)=").unwrap();
        // now if no "M" flag we do rather the match on "="
        if re_match.is_match(ciggi) {
            for m in re_match.captures_iter(ciggi) {
                seq_match += &m[1].parse().unwrap();
            }
        } else if re_match2.is_match(ciggi) {
            for m in re_match2.captures_iter(ciggi) {
                seq_match += &m[1].parse().unwrap();
            }
        } else {
            panic!("ERROR: no match in Cigar detected !");
        }

        ///////////////////////////////
        ////   # of mis-matches  //////
        ///////////////////////////////
        let re_mmatch = Regex::new(r"([0-9]+)X").unwrap();
        if re_mmatch.is_match(ciggi) {
            for m in re_mmatch.captures_iter(ciggi) {
                seq_mmatch += &m[1].parse().unwrap();
                // we add here as well the mismatches to matches
                // to facilitate later the calculations
                seq_match += &m[1].parse().unwrap();
            }
        }

        ///////////////////////////////
        ////     length of seq   //////
        ///////////////////////////////
        let re_length = Regex::new(r"([0-9]+)[H|S|M|I|X|=]").unwrap();
        for m in re_length.captures_iter(ciggi) {
            infos.seq_length += &m[1].parse().unwrap();
        }

        ///////////////////////////////
        ////     length of match //////
        ///////////////////////////////
        let re_length_m = Regex::new(r"([0-9]+)[M|X|=]").unwrap();
        for m in re_length_m.captures_iter(ciggi) {
            infos.match_length += &m[1].parse().unwrap();
        }

        ////////////////////////////////////
        ////   spliced match length   //////
        ////////////////////////////////////
        //let re_length_s = Regex::new(r"([0-9]+)[H|S|M|N|D|I|X|=]").unwrap();
        let re_length_s = Regex::new(r"([0-9]+)[M|N|D|X|=]").unwrap();
        for m in re_length_s.captures_iter(ciggi) {
            infos.match_length_s += &m[1].parse().unwrap();
        }
        
        /////////////////////////
        ////   similarity   /////
        /////////////////////////
        infos.similarity = calc_gapcompressed_sim(
            ciggi,
            seq_match,
            &seq_mmatch
        );

        debug!("Calculated similarity {}",&infos.similarity);
        ////////////////////////
        //// get breakpoint ////
        ////////////////////////
        let breakpoints = calc_align_properties(
            infos.upstream,
            infos.downstream,
            infos.seq_length,
            infos.match_length_s,
            *match_start,
            *match_strand,
            ratio
        );

        debug!("Calculated breakpoint proterties {:?}",&breakpoints);
        match breakpoints{ 
            None => None ,
            Some(x) => {
                infos.fp_target    = x.break_target;
                infos.fp_query     = x.break_query;
                infos.aln_start    = x.aln_start;
                infos.aln_end      = x.aln_end;
                infos.match_strand = x.strand;
                Some(infos)
            }
        }
    }else{
        None
    }
    
}


/// this takes now our accumulated information
/// and tries to provide a proper return of
/// the direction order and position of the
/// gene fusion position
/// removes entries which are too close
/// All coordinates are 0 based
/// 
/// Unittest: TRUE
///
pub fn evaluate_and_predict(
    primary_infos: AlignInfos,
    sa_infos: Vec<AlignInfos>,
) -> Vec<FullFusionEvidence> {

    let mut gathered: Vec<FullFusionEvidence> = Vec::new();
    // we can have >1 entry of SA mappings
    for align in sa_infos.iter() {
        debug!("SA Alignment information to be evaluated: {:?}",&align);
        debug!("PR Alignment information to be evaluated: {:?}",&primary_infos);
        // now we want to determine which one is the 5'part of the
        // aligned piece and which is the 3'part and adjust as well
        // coordinates accordingly

        let mut orient1_up : Option<bool> = None;
        let mut f_begin1   : Option<u64> = Some(primary_infos.aln_start );
        let mut fp_target1 : Option<u64> = Some(primary_infos.fp_target );
        let mut orient2_up : Option<bool> = None;
        let mut f_begin2   : Option<u64> = Some(align.aln_start);
        let mut fp_target2 : Option<u64> = Some(align.fp_target );
        
        // Lets mix up things and evaluate the similarity between clipping information and match length
        // of the corresponding one. It should deal better with ambigious cases
        
        // PRIMARY ALIGNMENT 
        let up_diff = ((primary_infos.upstream as i64) - ( align.match_length as i64)).abs() as u32;
        let down_diff = ((primary_infos.downstream as i64) - ( align.match_length as i64)).abs() as u32;
        debug!("PR: up_diff difference absolute {} down_diff difference absolute {} ", up_diff,down_diff);
        // the assumption is that the one which is smaller indicates which clipping direction
        // is the more likely one

        // PRIMARY ALIGNMENT 
        //
        // upstream positioned
        if (( up_diff>down_diff) & (primary_infos.match_strand ==  StrandDirection::Fwd )) |
        ((up_diff<down_diff) & (primary_infos.match_strand ==  StrandDirection::Rev )){
            
            orient1_up    = Some(true);    
            if primary_infos.match_strand ==  StrandDirection::Fwd {
                debug!("Primary is upstream and FWD",);
                // verified and correct
                f_begin1   = Some(primary_infos.aln_start);
                fp_target1 = Some(fp_target1.unwrap() );
            }else{
                debug!("Primary is upstream and REF",);
                // verified and correct
                f_begin1   = Some(primary_infos.aln_start);
                fp_target1 = Some(fp_target1.unwrap());
            }
        // downstream positioned
        } else if ((up_diff<down_diff) & (primary_infos.match_strand == StrandDirection::Fwd)) |
			      ((up_diff>down_diff ) & (primary_infos.match_strand == StrandDirection::Rev)) {
            orient1_up = Some(false);
            if primary_infos.match_strand == StrandDirection::Fwd {
                debug!("Primary is downstream and FWD",);
                // verified and correct
                f_begin1 = Some(primary_infos.aln_start);
                fp_target1 = Some(primary_infos.aln_start );

            }else{
                debug!("Primary is downstream and REF",);
                // verified and correct
                fp_target1 = Some(primary_infos.aln_end );
                f_begin1   = Some(primary_infos.aln_start);
            }
        }

        let up_diff_sa   = (( align.upstream as i64) - (primary_infos.match_length as i64)  ).abs() as u32;
        let down_diff_sa = (( align.downstream as i64) - (primary_infos.match_length as i64)  ).abs() as u32;
        debug!("SA: up_diff_sa difference absolute {} down_diff_sa difference absolute {} ", up_diff_sa,down_diff_sa);
        // the assumption is that the one which is smaller indicates which clipping direction
        // is the more likely one

        // SA ALIGNMENT 
        //
        // upstream
        if (( up_diff_sa>down_diff_sa) & (align.match_strand == StrandDirection::Fwd )) | 
			((up_diff_sa<down_diff_sa) & (align.match_strand == StrandDirection::Rev)) {
            orient2_up    = Some(true);    
            if align.match_strand == StrandDirection::Fwd {
                debug!("SA is upstream and FWD",);

                // verified and correct
                f_begin2   = Some(align.aln_start);
                fp_target2 = Some(align.aln_end);

            }else{
                debug!("SA is upstream and REF",);

                // verified and correct
                f_begin2   = Some(align.aln_start);
                fp_target2 = Some(align.aln_start);

            }
        // downstream
        } else if ((up_diff_sa<down_diff_sa) & (align.match_strand == StrandDirection::Fwd)) | 
				((  up_diff_sa>down_diff_sa) & (align.match_strand == StrandDirection::Rev) ) {
            orient2_up    = Some(false);
            if align.match_strand == StrandDirection::Fwd {
                debug!("SA is downstream and FWD",);

                // verified and correct
                fp_target2 = Some(fp_target2.unwrap());
                f_begin2   = Some(align.aln_start);
            }else{
                debug!("SA is downstream and REF",);

                // verified and correct
                fp_target2 = Some(fp_target2.unwrap());
                f_begin2   = Some(align.aln_start);
            }

        } 
        // Now if both break information are contradictive trash it
        // as well if it was not resolved
        if orient1_up.is_none() | orient2_up.is_none() |  fp_target1.is_none()  | fp_target2.is_none() | f_begin1.is_none() | f_begin2.is_none()
        {
            
			debug!("INFO: Removed entry primary {:?} with SA {:?} as fusion positions ambiguous", primary_infos,align);
            debug!("INFO: orient1 {:?} orient2 {:?} fp_target1 {:?} fp_target2 {:?} f_begin1 {:?} f_begin2 {:?}",orient1_up,orient2_up,fp_target1,fp_target2,f_begin1,f_begin2);
            continue;
        }else if  orient1_up.unwrap() == orient2_up.unwrap() {
            eprintln!("WARNING: opposing information discovered in evaluting sense of fusion - this should not happen!");
            debug!("INFO: Removed entry primary {:?} with SA {:?} as fusion positions ambiguous", primary_infos,align);
            debug!("INFO: orient1 {:?} orient2 {:?} fp_target1 {:?} fp_target2 {:?} f_begin1 {:?} f_begin2 {:?}",orient1_up,orient2_up,fp_target1,fp_target2,f_begin1,f_begin2);
            continue;
        }
        
        // as this is reported from the other side around we need to correct
        // the coordinate
        let populator = FullFusionEvidence {
            cigar_pr:        primary_infos.cigar.to_owned(),
            cigar_sa:        align.cigar.to_owned(),
            chrom_pr:        primary_infos.chrom.to_owned(),
            chrom_sa:        align.chrom.to_owned(),
            query_name:      primary_infos.query.to_owned(),
            query_length:    primary_infos.q_length.to_owned(),
            target_break_pr: fp_target1.unwrap(),
            target_break_sa: fp_target2.unwrap(),
            query_break_pr:  primary_infos.fp_query as u64,
            query_break_sa:  align.fp_query as u64,
            strand_pr:       primary_infos.match_strand,
            strand_sa:       align.match_strand,
            orient_pr_upstream:       orient1_up.unwrap(),
            orient_sa_upstream:       orient2_up.unwrap(),
            aln_start_pr:    f_begin1.unwrap(),
            aln_start_sa:    f_begin2.unwrap(),
            // for both the below values we have
            // at this point no the information
            // available and need to evaluate this later
            support: 0,
            precision: None,
        };
        gathered.push(populator);
    }
    gathered
}


/// this function takes a read-based organization
/// which must contain the FullFusionEvidence structure
/// it then evaluates all entries and removes if necessary
/// entry which have insufficient evidence.
/// Furthermore we take combine the information to have
/// from both the information in principle coming from the
/// primary alignment rather than using SA derived information
/// This provides as well a score to evaluate precision of determined fusion point:
/// 0 = no complementary entry (e.g. sim. too low for one of them)
/// 1 = both determined
/// Additionally offering now the option to keep all entries 
/// for suggested SA entries, not just the best. This might be favorable in some
/// situation (e.g. integration analysis)
/// 
/// ## Improvements:
/// I have the feeling that this section could actually 
/// benefit from some improvement. Functional wise it does what it 
/// is supposed to do but code could def. be optimized and I 
/// have some values being cloned way too often but could not 
/// get is better optimized for the moment.
///
/// Unittest: FALSE
///
pub fn prune_read_based_organization(
    my_hash: FxHashMap<String, Vec<FullFusionEvidence>>,
    keep_all: bool,
) -> FxHashMap<String, Vec<FullFusionEvidence>> {

    let mut combined_entries: FxHashMap<String, Vec<FullFusionEvidence>> = FxHashMap::default();
    
    for (key, read) in my_hash.iter() {
        debug!("Analyzing key {} read {:?}",&key,&read);
        // we want to assure that each entry
        // has its counterpart, otherwise we put the score at 0
        let read_copy = read;
        let mut gathered: Vec<FullFusionEvidence> = Vec::new();
        let mut already_done : Option<usize> = None;
        let mut precision:    Option<u32> = None;
        let mut best_precision : Vec<u32> =Vec::new();
        //let v: Vec<FullFusionEvidence> = read.into_iter().unique().collect();
        for (index, check) in read.iter().enumerate() {
            let chr_a = &check.chrom_pr;
            let chr_b = &check.chrom_sa;
            let mut get_match: bool = false;
            // if this position was already treated and
            // found a match ignore it to avoid
            // double entries
            if already_done == Some(index) {
                continue;
            }
            for (index2, double_check) in read_copy.iter().enumerate() {
                // if we look at identical entries skip to next
                if index == index2 {
                    continue;
                }

                // now occasionally we can actually get multiple primary and multiple SA ones
                // to avoid that multiple primary ones are paired with each other we check that they 
                // are not actually identical
                if check.cigar_pr == double_check.cigar_pr {
                    //eprintln!("check {:?} vs check2 {:?} skipped as same primary cigar {} vs {}", &check,&double_check,&pos,&pos2);
                    continue;
                }

                // this still leaves potentially open a scenario where
                // we could mix entries if e.g. 2 SA alignments both
                // point to a diff. location on same chromosome 
                // This cant be though prevented and should then be determined by the precision
                debug!("query break pr : {} 2nd query_break_pr {} : diff {}", &check.query_break_pr,&double_check.query_break_pr,check.query_break_pr.abs_diff(double_check.query_break_pr));
                if (&double_check.chrom_pr == chr_b) & (&double_check.chrom_sa == chr_a) {

                    precision = Some(
                        u32::try_from(check.query_break_pr.abs_diff(double_check.query_break_pr )/2)
                            .expect("ERROR: could not determine precision as conversion into u32 failed"));
                    // now we concatenate the evidence and take from
                    // each entry the primary alignment information which
                    // contains the precise alignment information
                    let tmp = FullFusionEvidence {
                        cigar_pr:        check.cigar_pr.to_owned(),
                        cigar_sa:        double_check.cigar_pr.to_owned(),
                        chrom_pr:        check.chrom_pr.to_owned(),
                        chrom_sa:        double_check.chrom_pr.to_owned(),
                        query_name:      check.query_name.to_owned(),
                        query_length:    check.query_length,
                        target_break_pr: check.target_break_pr,
                        target_break_sa: double_check.target_break_pr,
                        query_break_pr:  check.query_break_pr,
                        query_break_sa:  double_check.query_break_pr,
                        strand_pr:       check.strand_pr,
                        strand_sa:       double_check.strand_pr,
                        orient_pr_upstream:       check.orient_pr_upstream,
                        orient_sa_upstream:       double_check.orient_pr_upstream,
                        aln_start_pr:    check.aln_start_pr,
                        aln_start_sa:    double_check.aln_start_pr,
                        precision,
                        // here now support is 1 as we got both evidence
                        support: 1,
                    };
                    //dbg!(&tmp);
                    get_match = true;
                    gathered.push(tmp);
                    best_precision.push(precision.unwrap());
                    already_done = Some(index2);
                }
            }
            // now if we did not find a partner for that match
            // we keep it but now precision is really lower due to
            // limited information from SA cigar
            if !get_match {
                let tmp2 = FullFusionEvidence {
                    cigar_pr:        check.cigar_pr.to_owned(),
                    cigar_sa:        check.cigar_sa.to_owned(),
                    chrom_pr:        check.chrom_pr.to_owned(),
                    chrom_sa:        check.chrom_sa.to_owned(),
                    query_name:      check.query_name.to_owned(),
                    query_length:    check.query_length.to_owned(),
                    target_break_pr: check.target_break_pr,
                    target_break_sa: check.target_break_sa,
                    query_break_pr:  check.query_break_pr,
                    query_break_sa:  check.query_break_sa,
                    strand_pr:       check.strand_pr,
                    strand_sa:       check.strand_sa,
                    orient_pr_upstream:       check.orient_pr_upstream,
                    orient_sa_upstream:       check.orient_sa_upstream,
                    aln_start_pr:    check.aln_start_pr,
                    aln_start_sa:    check.aln_start_sa,
                    precision,
                    support: 0,
                };
                gathered.push(tmp2);
            }
        }
        // lets only keep best entries and purge the rest 
        best_precision.sort_unstable();
        if !best_precision.is_empty() && !keep_all {
            let precision_to_keep = best_precision[0];
            gathered.retain(|x| x.precision <= Some(precision_to_keep));
        }
        // we still might have duplicates unfortunately which we need 
        // to remove again
        // construct unique entries
        let mut final_purge : FxHashMap<String,FullFusionEvidence> = FxHashMap::default();
        for final_candidates in gathered.into_iter(){
            // due to many:many entries we can have the exact 
            // same entry but pa/sa switched
            // therefore we construct here a sorted identifier based on cigars and
            // collapse essentially identical entries
            let mut uniq_key_v : Vec<String> = vec![final_candidates.cigar_pr.clone(),final_candidates.cigar_sa.clone()];
            uniq_key_v.sort();
            let uniq_key   = uniq_key_v.join(":");
            final_purge.insert(uniq_key,final_candidates);
        };
        let mut final_candidates : Vec<FullFusionEvidence> = Vec::new();
        for entry in final_purge.values().cloned() {
            final_candidates.push(entry);
        }
        debug!("Final candidate {:?}", &final_candidates);
        combined_entries.insert(key.to_string(), final_candidates);
    }
    combined_entries
}



/// This function needs the attributes parsed from a GTF and checks if
/// it can obtain a gene name or a gene ID. If that fails it defaults 
/// to NA
/// 
/// Unittest: TRUE
///
pub fn get_gene_name(
    attributes: &FxHashMap<String,String>
) -> String {
    if attributes.contains_key("gene_name"){
        attributes.get("gene_name").unwrap().to_owned()
    }else if attributes.contains_key("geneID"){
        attributes.get("geneID").unwrap().to_owned()
    }else if attributes.contains_key("gene_id"){
        attributes.get("gene_id").unwrap().to_owned()
    }else{
        //eprintln!("I got no name from {:?}", attributes);
        
        String::from("NA")
    }
}

/// this takes the start and end of each annotated fusion feature
/// and verifies if we have in the annotation a corresponding feature
/// if that is the case it mutates the annotated fusion gene definition
/// It calculates for each potential candidate a score 
/// which indicates the positioning and a similarity 
/// in percent which describes the percent covered 
/// Currently, only the score is propagated in the resulting output.
/// If only one of the features is annotated then it reflects that ones score,
/// if both are annotated it reflects the average of both.
/// 
/// If we cant distinguish which gene is underlying as >1 have identical scores,
/// then we generate for each an independent output. 
/// One might want to indicate in the future how many alternatives for 1 elements 
/// were potentially obtained.
///
///Unittest: FALSE
///
pub fn add_fusion_annotation_read_based(
    out: PrincipleReadBasedFusionOutput, 
    annot: &Option<FxHashMap<String, Vec<Gtf>>>, 
    stranded: bool
) -> Vec<PrincipleReadBasedFusionOutput>{

    let mut my_return :  Vec<PrincipleReadBasedFusionOutput> = Vec::new();
    let annot2    = annot.as_ref().unwrap();
    // first lets verify if the chromosome is actually in annotation
    // if not we panic here
    if !annot2.contains_key(&out.a_chrom) | !annot2.contains_key(&out.b_chrom){
        debug!("WARNING: No annotation available for one of chromosomes from entry {:?}",out);
    }
    // there are now essentially 3 scenarios to test
	// - candidate region is outside gene annotation region --> drop it
	// - candidate is within gene annotation region --> keep it
	// - candidate is on any side overlapping --> keep it 

    let mut sug_feature_a : Vec<OverlapResult> = Vec::new();
    let mut sug_feature_b : Vec<OverlapResult> = Vec::new();
    if annot2.contains_key(&out.a_chrom) {
        
    // GeneA
        for feature in annot2[&out.a_chrom].iter() {     
            let mut result_a : OverlapResult =  calc_match_overlap(
                out.a_start,
                out.a_end, 
                out.a_strand,
                feature.start,
                feature.end,
                feature.strand,
                stranded
            );
            if result_a.match_score == 0 {
                continue;
            }
            // now we get a matching gene name field
            result_a.feature_name = get_gene_name(&feature.attribute);
            debug!("feature name: {:?}",&result_a.feature_name );

            if sug_feature_a.is_empty()  {
                sug_feature_a.push(result_a);    
            }else{
                
                // remove all entries with lesser score
                let rm_index  = sug_feature_a.iter().position(|x| x.match_score < result_a.match_score);
                sug_feature_a =  match rm_index {
                    Some(value) => {
                        sug_feature_a.remove(value);
                        sug_feature_a
                    },
                    None => sug_feature_a,
                };
                
                if sug_feature_a.is_empty(){
                    sug_feature_a.push(result_a);    
                    continue;
                }
                let mut to_push = Keeper::Donnu ;

                for (_pos,candidates) in sug_feature_a.iter().enumerate() {
                    // if it is worse than the prev. found skip
                    if result_a.match_score < candidates.match_score {
                        //eprintln!("A: {} Inferior to feature {:?} ",&result_a.feature_name,&candidates.feature_name);
                        to_push = Keeper::No;
                        break;
                    }
                    // if score is identical but same feature name skip
                    if (result_a.match_score == candidates.match_score) & ( result_a.feature_name == candidates.feature_name) {
                        //eprintln!("A: {} identical to feature {:?} ",&result_a.feature_name,&candidates.feature_name);
                        to_push = Keeper::Doubled;
                    break;
                    // if different features keep them because we need them if we cant distinguish 
                    }else if (result_a.match_score == candidates.match_score) & ( result_a.feature_name != candidates.feature_name) {
                        //eprintln!("A: {} similar to feature {:?} but diff name ",&result_a.feature_name,&candidates.feature_name);
                        match &mut to_push {
                            Keeper::No      => to_push = Keeper::No,
                            Keeper::Doubled => to_push = Keeper::No,
                            _ => to_push = Keeper::Yes,
                        };
                    // now if it is better we remove that old candidate in favor for new
                    }
                }
                if let Keeper::Yes = to_push { sug_feature_a.push(result_a)};
            }
        };
    }
    
    if annot2.contains_key(&out.b_chrom) {
    // GeneB
        for feature in annot2[&out.b_chrom].iter() {
            // first if no overlap at all skip directly
            let mut result_b : OverlapResult =  calc_match_overlap(
                out.b_start,
                out.b_end,
                out.b_strand,
                feature.start,
                feature.end,
                feature.strand,
                stranded
            );

            if result_b.match_score == 0 {
                continue;
            }
            // now we get a matching gene name field
            result_b.feature_name = get_gene_name(&feature.attribute);
            //dbg!(&result_b);
            if sug_feature_b.is_empty() {
                sug_feature_b.push(result_b);    
            }else{
                // remove all entries with lesser score
                let rm_index = sug_feature_b.iter().position(|x| x.match_score < result_b.match_score);
                sug_feature_b =  match rm_index {
                    Some(value) => {
                        sug_feature_b.remove(value);
                        sug_feature_b
                    },
                    None => sug_feature_b,
                };          // now if the list is again empty we can just add
                if sug_feature_b.is_empty(){
                    sug_feature_b.push(result_b);    
                    continue;
                }
                let mut to_push = Keeper::Donnu ;

                for (_pos,candidates) in sug_feature_b.iter().enumerate() {
                    // if it is worse than the prev. found skip
                    if result_b.match_score < candidates.match_score {
                        //eprintln!("A: {} Inferior to feature {:?} ",&result_b.feature_name,&candidates.feature_name);
                        to_push = Keeper::No;
                        break;
                    }
                    // if score is identical but same feature name skip
                    if (result_b.match_score == candidates.match_score) & ( result_b.feature_name == candidates.feature_name) {
                        //eprintln!("A: {} identical to feature {:?} ",&result_b.feature_name,&candidates.feature_name);
                        to_push = Keeper::Doubled;
                        break;
                    // if different features keep them because we need them if we cant distinguish 
                    }else if (result_b.match_score == candidates.match_score) & ( result_b.feature_name != candidates.feature_name) {
                        //eprintln!("A: {} similar to feature {:?} but diff name ",&result_b.feature_name,&candidates.feature_name);
                        match &mut to_push {
                            Keeper::No      => to_push = Keeper::No,
                            Keeper::Doubled => to_push = Keeper::No,
                            _ => to_push = Keeper::Yes,
                        };
                    // now if it is better we remove that old candidate in favor for new
                    }
                }
                
                if let Keeper::Yes = to_push {sug_feature_b.push(result_b)}
            }
        };
    }
    // at this step we have an array with only 1 element
    // therefore we can pick without issues the first element
    // Later it will get populated with annotations and cand
    // contain more than one element
    if out.fusion_genes.len() > 1_usize {
        panic!("ERROR: Input event contained more than 1 annotated fusion element!");
    }
    // here we need to generate now all combinations
    // it is not necessary to specify each line all conditions but 
    // makes it easier to read and understand logic
    if (sug_feature_a.is_empty() ) & (sug_feature_b.is_empty() ){
        my_return.push(out);
    }else if (sug_feature_a.is_empty()) & (!sug_feature_b.is_empty() ) {
        let mut new_name;
        let reg_a  = Regex::new(r"^.*--").unwrap();
        let get_a  = reg_a.captures(&out.fusion_genes[0]).unwrap();
        let name_a = get_a.get(1).map_or("NA", |m| m.as_str());
        let mut tmp_out  = out.clone();
        tmp_out.fusion_genes.pop();
        for b_solutions in sug_feature_b.iter(){
            new_name = format!("{}--{}",&name_a,b_solutions.feature_name);
            tmp_out.fusion_genes.push(new_name);
            //tmp_out.annot_score  = b_solutions.match_score;
        }
        my_return.push(tmp_out);

    }else if (!sug_feature_a.is_empty()) & (sug_feature_b.is_empty() ){
        let mut new_name;
        let reg_b  = Regex::new(r"--.*$").unwrap();
        let get_b  = reg_b.captures(&out.fusion_genes[0]).unwrap();
        let name_b = get_b.get(1).map_or("NA", |m| m.as_str());
        let mut tmp_out  = out.clone();
        tmp_out.fusion_genes.pop();
        for a_solutions in sug_feature_a.iter(){
            new_name = format!("{}--{}",a_solutions.feature_name,&name_b);
            tmp_out.fusion_genes.push(new_name);
            //tmp_out.annot_score  = a_solutions.match_score;
        }
        my_return.push(tmp_out);
    }else if (!sug_feature_a.is_empty()) & (!sug_feature_b.is_empty() ){
        let mut tmp_out      = out;
        tmp_out.fusion_genes.pop();
        for a_solutions in sug_feature_a.iter(){
            for b_solutions in sug_feature_b.iter(){
                tmp_out.fusion_genes.push(format!("{}--{}",a_solutions.feature_name,b_solutions.feature_name));
                //tmp_out.annot_score  = (a_solutions.match_score + b_solutions.match_score)/2;
            }
        }
        my_return.push(tmp_out);
    }else{
        eprintln!("ERROR: feature evaluation not successful for {:?}",out);
    }
    
    my_return
}



/// this takes the start and end of each annotated fusion feature
/// and verifies if we have in the annotation a corresponding feature
/// if that is the case it mutates the annotated fusion gene definition
/// It calculates for each potential candidate a score 
/// which indicates the positioning and a similarity 
/// in percent which describes the percent covered 
/// Currently, only the score is propagated in the resulting output.
/// If only one of the features is annotated then it reflects that ones score,
/// if both are annotated it reflects the average of both.
/// 
/// If we cant distinguish which gene is underlying as >1 have identical scores,
/// then we generate for each an independent output. 
/// One might want to indicate in the future how many alternatives for 1 elements 
/// were potentially obtained.
///
///Unitest: FALSE
///
pub fn add_fusion_annotation_pos_based(
    out: PrinciplePosBasedFusionOutput, 
    annot: &Option<FxHashMap<String, Vec<Gtf>>>, 
    stranded: bool
) -> Vec<PrinciplePosBasedFusionOutput>{

    let mut my_return :  Vec<PrinciplePosBasedFusionOutput> = Vec::new();
    let annot2    = annot.as_ref().unwrap();
    // first lets verify if the chromosome is actually in annotation
    // if not we panic here
    if !annot2.contains_key(&out.a_chrom) | !annot2.contains_key(&out.b_chrom){
        debug!("WARNING: No annotation available for one of chromosomes from entry {:?}",out);
    }
    // there are now essentially 3 scenarios to test
	// - candidate region is outside gene annotation region --> drop it
	// - candidate is within gene annotation region --> keep it
	// - candidate is on any side overlapping --> keep it 


    let mut sug_feature_a : Vec<OverlapResult> = Vec::new();
    let mut sug_feature_b : Vec<OverlapResult> = Vec::new();

    if annot2.contains_key(&out.a_chrom) {  
    // GeneA
        for feature in annot2[&out.a_chrom].iter() {     
            let mut result_a : OverlapResult =  calc_match_overlap(
                out.a_start,
                out.a_end, 
                out.a_strand,
                feature.start,
                feature.end,
                feature.strand,
                stranded
            );
            if result_a.match_score == 0 {
                continue;
            }
            // now we get a matching gene name field
            result_a.feature_name = get_gene_name(&feature.attribute);
            //dbg!(&result_a);

            if sug_feature_a.is_empty()  {
                sug_feature_a.push(result_a);    
            }else{
                
                // remove all entries with lesser score
                let rm_index  = sug_feature_a.iter().position(|x| x.match_score < result_a.match_score);
                sug_feature_a =  match rm_index {
                    Some(value) => {
                        sug_feature_a.remove(value);
                        sug_feature_a
                    },
                    None => sug_feature_a,
                };
                
                if sug_feature_a.is_empty(){
                    sug_feature_a.push(result_a);    
                    continue;
                }
                let mut to_push = Keeper::Donnu ;

                for (_pos,candidates) in sug_feature_a.iter().enumerate() {
                    // if it is worse than the prev. found skip
                    if result_a.match_score < candidates.match_score {
                        //eprintln!("A: {} Inferior to feature {:?} ",&result_a.feature_name,&candidates.feature_name);
                        to_push = Keeper::No;
                        break;
                    }
                    // if score is identical but same feature name skip
                    if (result_a.match_score == candidates.match_score) & ( result_a.feature_name == candidates.feature_name) {
                        //eprintln!("A: {} identical to feature {:?} ",&result_a.feature_name,&candidates.feature_name);
                        to_push = Keeper::Doubled;
                    break;
                    // if different features keep them because we need them if we cant distinguish 
                    }else if (result_a.match_score == candidates.match_score) & ( result_a.feature_name != candidates.feature_name) {
                        //eprintln!("A: {} similar to feature {:?} but diff name ",&result_a.feature_name,&candidates.feature_name);
                        match &mut to_push {
                            Keeper::No      => to_push = Keeper::No,
                            Keeper::Doubled => to_push = Keeper::No,
                            _ => to_push = Keeper::Yes,
                        };
                    // now if it is better we remove that old candidate in favor for new
                    }
                }
                if let Keeper::Yes = to_push { sug_feature_a.push(result_a)};
            }
        };
    }
    
    if annot2.contains_key(&out.b_chrom) {
    // GeneB
        for feature in annot2[&out.b_chrom].iter() {
            // first if no overlap at all skip directly
            let mut result_b : OverlapResult =  calc_match_overlap(
                out.b_start,
                out.b_end,
                out.b_strand,
                feature.start,
                feature.end,
                feature.strand,
                stranded
            );

            if result_b.match_score == 0 {
                continue;
            }
            // now we get a matching gene name field
            result_b.feature_name = get_gene_name(&feature.attribute);
            //dbg!(&result_b);
            if sug_feature_b.is_empty() {
                sug_feature_b.push(result_b);    
            }else{
                // remove all entries with lesser score
                let rm_index = sug_feature_b.iter().position(|x| x.match_score < result_b.match_score);
                sug_feature_b =  match rm_index {
                    Some(value) => {
                        sug_feature_b.remove(value);
                        sug_feature_b
                    },
                    None => sug_feature_b,
                };          // now if the list is again empty we can just add
                if sug_feature_b.is_empty(){
                    sug_feature_b.push(result_b);    
                    continue;
                }
                let mut to_push = Keeper::Donnu ;

                for (_pos,candidates) in sug_feature_b.iter().enumerate() {
                    // if it is worse than the prev. found skip
                    if result_b.match_score < candidates.match_score {
                        debug!("A: {} Inferior to feature {:?} ",&result_b.feature_name,&candidates.feature_name);
                        to_push = Keeper::No;
                        break;
                    }
                    // if score is identical but same feature name skip
                    if (result_b.match_score == candidates.match_score) & ( result_b.feature_name == candidates.feature_name) {
                        debug!("A: {} identical to feature {:?} ",&result_b.feature_name,&candidates.feature_name);
                        to_push = Keeper::Doubled;
                        break;
                    // if different features keep them because we need them if we cant distinguish 
                    }else if (result_b.match_score == candidates.match_score) & ( result_b.feature_name != candidates.feature_name) {
                        debug!("A: {} similar to feature {:?} but diff name ",&result_b.feature_name,&candidates.feature_name);
                        match &mut to_push {
                            Keeper::No      => to_push = Keeper::No,
                            Keeper::Doubled => to_push = Keeper::No,
                            _ => to_push = Keeper::Yes,
                        };
                    // now if it is better we remove that old candidate in favor for new
                    }
                }
                
                if let Keeper::Yes = to_push {sug_feature_b.push(result_b)}
            }
        };
    }
    // at this step we have an array with only 1 element
    // therefore we can pick without issues the first element
    // Later it will get populated with annotations and cand
    // contain more than one element
    if out.fusion_genes.len() > 1_usize {
        panic!("ERROR: Input event contained more than 1 annotated fusion element!");
    }
    // here we need to generate now all combinations
    // it is not necessary to specify each line all conditions but 
    // makes it easier to read and understand logic
    if (sug_feature_a.is_empty() ) & (sug_feature_b.is_empty() ){
        my_return.push(out);
    }else if (sug_feature_a.is_empty()) & (!sug_feature_b.is_empty() ) {
        let mut new_name;
        let reg_a  = Regex::new(r"^.*--").unwrap();
        let get_a  = reg_a.captures(&out.fusion_genes[0]).unwrap();
        let name_a = get_a.get(1).map_or("NA", |m| m.as_str());
        let mut tmp_out  = out.clone();
        tmp_out.fusion_genes.pop();
        for b_solutions in sug_feature_b.iter(){
            new_name = format!("{}--{}",&name_a,b_solutions.feature_name);
            tmp_out.fusion_genes.push(new_name);
        };
        my_return.push(tmp_out);

    }else if (!sug_feature_a.is_empty()) & (sug_feature_b.is_empty() ){
        let mut new_name;
        let reg_b  = Regex::new(r"--.*$").unwrap();
        let get_b  = reg_b.captures(&out.fusion_genes[0]).unwrap();
        let name_b = get_b.get(1).map_or("NA", |m| m.as_str());
        let mut tmp_out  = out.clone();
        tmp_out.fusion_genes.pop();
        for a_solutions in sug_feature_a.iter(){  
            new_name = format!("{}--{}",a_solutions.feature_name,&name_b);
            tmp_out.fusion_genes.push(new_name);
        }
        my_return.push(tmp_out);
    }else if (!sug_feature_a.is_empty()) & (!sug_feature_b.is_empty() ){
        let mut tmp_out      = out;
        tmp_out.fusion_genes.pop();
        for a_solutions in sug_feature_a.iter(){
            for b_solutions in sug_feature_b.iter(){
                tmp_out.fusion_genes.push(format!("{}--{}",a_solutions.feature_name,b_solutions.feature_name));
            }
        }
        my_return.push(tmp_out);
    }else{
        eprintln!("ERROR: feature evaluation not successful for {:?}",out);
    }
    //dbg!(&my_return);
    my_return
}

/// function which takes the final results
/// and annotation if available and formats the  output accordingly. 
/// All coordinates are 0 based.
/// 
/// Notion: if no annotation is provided, gene fusion names will default to 
/// the chromosome e.g. chr3--chr5 , if there is annnotation provided
/// it will become GeneA--GeneB or in case of missing annotation for a region
/// GeneA--NA or NA--NA
/// 
/// Unittest: TRUE
///
pub fn format_and_annotate_read_based_fusions(
    read_collection: FxHashMap<String, Vec<FullFusionEvidence>>,
    gtf: Option<FxHashMap<String, Vec<Gtf>>>,
    stranded: bool,
    max_dist: u64
) -> Vec<PrincipleReadBasedFusionOutput> {

    let mut output : Vec<PrincipleReadBasedFusionOutput> = Vec::new();
    for (cdna, entries) in read_collection.iter() {
        debug!("Analyzing read {:?}",entries);
        // after annotation it is possible that we generate more 
        // than one entry again per item as we cant always resolve 
        // underlying genes perfectly
        // therefore we build yet another vector here
        for attributes in entries.iter() {
            //dbg!(attributes);
            // we first initiate it emptry and then populate it
            let mut result = PrincipleReadBasedFusionOutput {
                query_name      : cdna.to_string(),
                length          : attributes.query_length,
                fusion_point    : (attributes.query_break_pr + attributes.query_break_sa) / 2 ,
                fusion_genes    : vec![String::from("NA")],
                fp_fuzzy        : attributes.precision,
                fusion_distance : None, 
                a_chrom         : String::from("NA"),
                a_start:      0,
                a_end:        0,
                a_strand:     StrandDirection::Unknown,
                a_fp:         0,
                b_chrom:      String::from("NA"),
                b_start:      0,
                b_end:        0,
                b_strand:     StrandDirection::Unknown,
                b_fp:         0,
            };
            
            // now we define the order which one is
            // 5' and 3' for the output formatting
            result.fusion_genes[0] = match gtf.is_some() {
                true => format!("{}--{}","NA","NA"),
                false => format!("{}--{}",attributes.chrom_sa.to_owned(), attributes.chrom_pr.to_owned())
            };
            match attributes.orient_pr_upstream {
                true => {
                    result.a_chrom     = attributes.chrom_pr.to_owned();
                    result.a_strand    = attributes.strand_pr;
                    if attributes.aln_start_pr < attributes.target_break_pr {
                        result.a_start = attributes.aln_start_pr;
                        result.a_end   = attributes.target_break_pr;
                    } else {
                        result.a_start = attributes.target_break_pr;
                        result.a_end   = attributes.aln_start_pr;
                    }
                    result.b_chrom     = attributes.chrom_sa.to_owned();
                    result.b_strand    = attributes.strand_sa;
                    if attributes.aln_start_sa < attributes.target_break_sa {
                        result.b_start = attributes.aln_start_sa;
                        result.b_end   = attributes.target_break_sa;
                    } else {
                        result.b_start = attributes.target_break_sa;
                        result.b_end   = attributes.aln_start_sa;
                    }
                }
                false => {
                    result.a_chrom     = attributes.chrom_sa.to_owned();
                    result.a_strand    = attributes.strand_sa;
                    if attributes.aln_start_sa < attributes.target_break_sa {
                        result.a_start = attributes.aln_start_sa;
                        result.a_end   = attributes.target_break_sa;
                    } else {
                        result.a_start = attributes.target_break_sa;
                        result.a_end   = attributes.aln_start_sa;
                    }
                    result.b_chrom     = attributes.chrom_pr.to_owned();
                    result.b_strand    = attributes.strand_pr;
                    if attributes.aln_start_pr < attributes.target_break_pr {
                        result.b_start = attributes.aln_start_pr;
                        result.b_end   = attributes.target_break_pr;
                    } else {
                        result.b_start = attributes.target_break_pr;
                        result.b_end   = attributes.aln_start_pr;
                    }
                }
            };
            if result.a_chrom == result.b_chrom {
                // A located 5' of B
                if result.a_end < result.b_start {
                    result.fusion_distance = Some(result.b_start - result.a_end );
                // A located 3' of B
                }else if result.b_end < result.a_start {
                    result.fusion_distance = Some(result.a_start - result.b_end );
                // overlapping
                }else if (( result.a_end > result.b_start ) & ( result.a_end < result.b_end )) 
                    | (( result.b_end > result.a_start ) & ( result.b_end < result.a_end )) {
                        result.fusion_distance = Some(0);
                }else{
                    debug!("INFO: could not correctly assign distance!");
                }
            }
            //dbg!(&result);
            // now we do want to check as well if they are not located
            // too close to each other and ignore the event otherwise
            // this is e.g. for read-through events very important
            if result.a_chrom == result.b_chrom {
                // ---A---B--->
                if result.a_end < result.b_start && result.b_start - result.a_end <= max_dist {
                    continue
                }
                // --B--A--->
                if result.b_end < result.a_start && result.a_start - result.b_end <= max_dist {
                    continue
                }
                // if overlapping -- not sure that can even happen here
                if (result.a_end >= result.b_start ) & (result.a_start <= result.b_start){
                    continue
                }
                if (result.b_end >= result.a_start ) & (result.b_start <= result.a_start){
                    continue
                }

            }

            // now we are just filling out the fusion point 
            // information
            if result.a_strand == StrandDirection::Rev {
                result.a_fp = result.a_start
            }else if result.a_strand == StrandDirection::Fwd {
                result.a_fp = result.a_end
            }else{
                panic!("ERROR: GeneA strand not correctly annotated!");
            }
            if result.b_strand == StrandDirection::Rev {
                result.b_fp = result.b_end
            }else if result.b_strand == StrandDirection::Fwd {
                result.b_fp = result.b_start
            }else{
                panic!("ERROR: GeneB strand not correctly annotated!");
            }
            //println!("result: {:?}", result);
            // this is now a complete output already
            // lets now see if we have a potential annotation
            // and accordingly annotate if possible
            if gtf.is_some() {
                let mut annotated = add_fusion_annotation_read_based(result,&gtf,stranded); 
                output.append(&mut annotated);
            }else{
                output.push(result);
            }
        }
    }
    output
}


/// function which takes the final results
/// and annotation if available and formats the  output accordingly. 
/// We need at least 2 events supporting it, otherwise precision
/// is NA
///
///Unittest: FALSE
///
pub fn format_and_annotate_pos_based_fusions(
    read_collection: FxHashMap<String, Vec<PosBasedFusEvidence>>,
    gtf: Option<FxHashMap<String, Vec<Gtf>>>,
    stranded: bool,
    max_dist: u64
) -> Vec<PrinciplePosBasedFusionOutput> {

    let mut output : Vec<PrinciplePosBasedFusionOutput> = Vec::new();
    for (_cdna, entries) in read_collection.iter() {
        // after annotation it is possible that we generate more 
        // than one entry again per item as we cant always resolve 
        // underlying genes perfectly
        // therefore we build yet another vector here
        for attributes in entries.iter() {
            // currently we have for the alignment an array of information
            // at this point we take the median as a conservative value
            let mut aln_start_array_pr = attributes.aln_start_pr.clone();
            let mut aln_start_array_sa = attributes.aln_start_sa.clone();
            let aln_start_pr_med : u64 = median_u64(&mut aln_start_array_pr);
            let aln_start_sa_med : u64 = median_u64(&mut aln_start_array_sa);

            // here we combine precision as difference for both, primary and secondary events
            let precision_array_i1 = attributes.precision_pr.clone();
            let precision_array_i2 = attributes.precision_sa.clone();
            let mut precision_array_f : Vec<f64>= precision_array_i1.into_iter().flat_map(|v| vec![v.abs_diff(attributes.target_break_pr) as f64]).collect();
            let mut precision_array_f2 : Vec<f64>= precision_array_i2.into_iter().flat_map(|v| vec![v.abs_diff(attributes.target_break_sa) as f64]).collect();
            precision_array_f.append(&mut precision_array_f2);
            debug!("Precision from {:?}",precision_array_f);
            // here we want though again an int as float really does not 
            // make sense for bps eventually
            // Caution though: if we have <2 data-points standard_deviation will fail 
            
            let aln_dev : Option<u32> = if precision_array_f.len() > 2 {
                let stdev = standard_deviation(&precision_array_f,None).round() as u32;
                debug!("Precision {:?} from {:?}",stdev,precision_array_f);
                Some(stdev)
            }else{
                debug!("Precision absent from {:?}",precision_array_f);
                None 
            };
            
            // we first initiate it emptry and then populate it
            let mut result = PrinciplePosBasedFusionOutput {
                fusion_genes    : vec![String::from("NA")],
                fp_fuzzy        : aln_dev,
                fusion_distance : None, 
                a_chrom         : String::from("NA"),
                a_start         : 0,
                a_end           : 0,
                a_strand        : StrandDirection::Unknown,
                a_fp            : 0,
                b_chrom         : String::from("NA"),
                b_start         : 0,
                b_end           : 0,
                b_strand        : StrandDirection::Unknown,
                b_fp            : 0,
                sp_support      : attributes.support_sp,
                sr_support      : attributes.support_sr,
                sp_ratio        : attributes.sp_ratio,
            };
            // now we only define the order which one is
            // 5' and 3' for the output formatting
            match attributes.orient_pr_upstream {
                true => {
                    result.fusion_genes[0]= format!(
                        "{}--{}",
                        "NA",
                        "NA"
                    );
                    result.a_chrom     = attributes.chrom_pr.to_owned();
                    result.a_strand    = attributes.strand_pr;

                    if aln_start_pr_med < attributes.target_break_pr {
                        result.a_start = aln_start_pr_med;
                        result.a_end   = attributes.target_break_pr;
                    } else {
                        result.a_start = attributes.target_break_pr;
                        result.a_end   = aln_start_pr_med;
                    }
                    result.b_chrom     = attributes.chrom_sa.to_owned();
                    result.b_strand    = attributes.strand_sa;
                    if aln_start_sa_med < attributes.target_break_sa {
                        result.b_start = aln_start_sa_med;
                        result.b_end   = attributes.target_break_sa;
                    } else {
                        result.b_start = attributes.target_break_sa;
                        result.b_end   = aln_start_sa_med;
                    }
                }
                false => {
                    result.fusion_genes[0] = format!(
                        "{}--{}",
                        attributes.chrom_sa.to_owned(),
                        attributes.chrom_pr.to_owned()
                    );
                    result.a_chrom     = attributes.chrom_sa.to_owned();
                    result.a_strand    = attributes.strand_sa;
                    if aln_start_sa_med < attributes.target_break_sa {
                        result.a_start = aln_start_sa_med;
                        result.a_end   = attributes.target_break_sa;
                    } else {
                        result.a_start = attributes.target_break_sa;
                        result.a_end   = aln_start_sa_med;
                    }
                    result.b_chrom     = attributes.chrom_pr.to_owned();
                    result.b_strand    = attributes.strand_pr;
                    if aln_start_pr_med < attributes.target_break_pr {
                        result.b_start = aln_start_pr_med;
                        result.b_end   = attributes.target_break_pr;
                    } else {
                        result.b_start = attributes.target_break_pr;
                        result.b_end   = aln_start_pr_med;
                    }
                }
                
            };
            if result.a_chrom == result.b_chrom {
                // A located 5' of B
                if result.a_end < result.b_start {
                    result.fusion_distance = Some(result.b_start - result.a_end );
                // A located 3' of B
                }else if result.b_end < result.a_start {
                    result.fusion_distance = Some(result.a_start - result.b_end );
                // overlapping
                }else if (( result.a_end > result.b_start ) & ( result.a_end < result.b_end )) 
                    | (( result.b_end > result.a_start ) & ( result.b_end < result.a_end )) {
                        result.fusion_distance = Some(0);
                }else{
                    eprintln!("INFO: could not correctly assign distance!");
                }
            }
            //dbg!(&result);
            // now we do want to check as well if they are not located
            // too close to each other and ignore the event otherwise
            // this is e.g. for read-through events very important
            if result.a_chrom == result.b_chrom {
                // ---A---B--->
                if result.a_end < result.b_start && result.b_start - result.a_end <= max_dist {
                    continue
                }
                // --B--A--->
                if result.b_end < result.a_start && result.a_start - result.b_end <= max_dist {
                    continue
                }
                // if overlapping -- not sure that can even happen here
                if (result.a_end >= result.b_start ) & (result.a_start <= result.b_start){
                    continue
                }
                if (result.b_end >= result.a_start ) & (result.b_start <= result.a_start){
                    continue
                }

            }

            // now we are just filling out the fusion point 
            // information
            if result.a_strand == StrandDirection::Rev {
                result.a_fp = result.a_start
            }else if result.a_strand == StrandDirection::Fwd {
                result.a_fp = result.a_end
            }else{
                panic!("ERROR: GeneA strand not correctly annotated!");
            }
            if result.b_strand == StrandDirection::Rev {
                result.b_fp = result.b_end
            }else if result.b_strand == StrandDirection::Fwd {
                result.b_fp = result.b_start
            }else{
                panic!("ERROR: GeneB strand not correctly annotated!");
            }

            // this is now a complete output already
            // lets now see if we have a potential annotation
            // and accordingly annotate if possible
           
            if gtf.is_some() {
                let mut annotated =add_fusion_annotation_pos_based(result,&gtf,stranded); 
                output.append(&mut annotated);
            }else{
                output.push(result);
            }
        }
    }
    output
}




/// A function which evaluates how well a feature matches a given annotation.
/// It simply takes start and end of both as well as the strand. 
/// 
/// Unittest: FALSE
///
/// Example:
///
/// ```rust
/// use genefusion::lib::common::{*};
/// pretty_env_logger::init();
///
/// let q_start           = 500_u64;
/// let q_end             = 1000_u64;
/// let q_strand          = StrandDirection::Fwd;
/// let t_start           = 750_u64;
/// let t_end             = 1250_u64;
/// let t_strand          = StrandDirection::Fwd; 
/// let stranded          = true;
/// let result = calc_match_overlap(q_start,q_end,q_strand,t_start,t_end,t_strand,stranded);
/// let test   = OverlapResult {
///        match_length : 250,
///        match_score  : 16,
///        feature_name : String::from(""),
///    };
/// assert_eq!(result,test);
/// ```
pub fn calc_match_overlap (
    q_start : u64,
    q_end   : u64, 
    q_strand: StrandDirection,
    t_start : u64,
    t_end   : u64,
    t_strand: StrandDirection,
    stranded: bool
) -> OverlapResult {
    debug!("Evaluation overlap for feature q_start {:?} q_end {:?} q_strand {:?} t_start {:?} t_end {:?} t_stramd {:?} stranded {:?}",q_start,q_end,q_strand,t_start,t_end,t_strand,stranded);
        // so now we get check the amount of overlap and
        // distribute scores
        let mut result : OverlapResult = Default::default();
        // first if no overlap at all skip directly
        if (q_end <= t_start) | (q_start >= t_end ){
            return result
        }
        let mut overlap_start : Option<u64> = None;
        let mut overlap_end   : Option<u64> = None;
        // now lets check feature and add score
        // is it overlapping or within a feature (latter preferred)
        // located 3' of start and potentially overlapping 3' 
        if stranded & (q_strand == t_strand ){
            result.match_score += 10;
        }
        // now first the case that both are in same 
        // direction stranded

        if (q_start == t_start) | (q_start > t_start) {
            debug!("Starts are either identical or q>t");
            match q_start.cmp(&t_start){
                Ordering::Equal => {
                    debug!("Starts are identical");
                    result.match_score  += 10; 
                    overlap_start = Some(q_start.to_owned());
                },
                Ordering::Greater => {
                    debug!("Qstart > Tstart");
                    result.match_score  += 1;
                    overlap_start = Some(q_start.to_owned());
                },
                Ordering::Less => {
                    debug!("Qstart < Tstart");
                    panic!("ERROR: smaller encountered where not possible for overlap comparison")
                }
            }

            match q_end.cmp(&t_end) {
                Ordering::Less => {
                    debug!("Qend < Tend");
                    result.match_score += 5;
                    overlap_end = Some(q_end.to_owned());
                },
                Ordering::Greater => {
                    debug!("Qend > Tend");
                    result.match_score += 1;
                    overlap_end = Some(t_end.to_owned());
                },
                Ordering::Equal => {
                    debug!("Qend == Tend");
                    result.match_score +=10;
                    overlap_end = Some(t_end.to_owned());
                }
            }

        // overlapping from 5'
        } else if q_start < t_start {
            result.match_score += 1;
            overlap_start = Some(t_start.to_owned());
            // finishing with annotated feature

            match q_end.cmp(&t_end) {
                Ordering::Greater => {
                    result.match_score += 1;
                    overlap_end = Some(t_end.to_owned());
                },
                Ordering::Equal => {
                    result.match_score += 10; 
                    overlap_end = Some(t_end.to_owned());
                },
                Ordering::Less => {
                    result.match_score += 5;
                    overlap_end = Some(q_end.to_owned());
                }
            }
        };
        if overlap_end.is_none() | overlap_start.is_none(){ 
           return result
        }
        result.match_length = overlap_end.unwrap() - overlap_start.unwrap();
        // now we double check that it resulted in something sensible concerning the length
        // there is an exception which was recently added, that the start position of the annotation is at 0 (contig start) which invalidated the 
        // incomplete test before - therefore we have now the testing if start is 0
        if (result.match_length == overlap_end.unwrap() && overlap_start.unwrap()!=0) | (result.match_length == overlap_start.unwrap()&& overlap_start.unwrap()!=0){
            panic!("ERROR: A:  overlap with annotation wrong, overlap is {:?}",result.match_length);
        }else{
            result
        }
}




/// this function takes the read-based organization of fusion detection
/// and re-organizes it into a position based structure. 
/// The output is a hashmap of hashmap which contains primary and secondary chromosome as keys
/// to optimize speed in looking it up in the position based comparison steps. 
/// To be considered the same fusion event, a second event needs to be in 
/// similar direction for both, target and query.
/// Additionally it must be within a given range of the position which can be supplied as argument.
/// In this case, the median of the visited positions is taken and the new value needs to be located 
/// within that median +/- the provided range which is considered acceptable.
/// The new entry is then added in the array of target and query position for its obtained position
/// and the support is increased by a count of 1.
/// If the input is supposedly stranded, then it will be strict in the way
/// that only events with same orientation and direction will be collapsed.
/// If we have though unstranded information then we collapse events which are 
/// otherwise compatible as then a 5-3 orentation with mapping in 1,-1 is identical to
/// a 3-5 orientation with mapping -1,1 for example
///
/// unittest: FALSE
///
pub fn read_based_2_pos_based_fusions(
    read_based: FxHashMap<String, Vec<FullFusionEvidence>>,
    range: u32,
    stranded:bool
) -> FxHashMap<String,Vec<PosBasedFusEvidence>> {
    let mut pos_based_result : FxHashMap<String,Vec<PosBasedFusEvidence>>  = FxHashMap::default();
    for (_key,value) in read_based {
        for evidence in value {
            let a  =  &evidence.chrom_pr.to_owned();
            let b  =  &evidence.chrom_sa.to_owned();
            let chr_ab = format!("{}{}",a ,b) ;
            // if chrom1 does not exist yet
            if ! pos_based_result.contains_key(&chr_ab) {
                // now we populate our new format and then add this
                // in the new structure
                let new_evidence = PosBasedFusEvidence{
                    chrom_pr     : evidence.chrom_pr,
                    chrom_sa     : evidence.chrom_sa,
                    target_break_pr: evidence.target_break_pr,
                    target_break_sa: evidence.target_break_sa,
                    strand_pr    : evidence.strand_pr,
                    strand_sa    : evidence.strand_sa,
                    orient_pr_upstream    : evidence.orient_pr_upstream,
                    orient_sa_upstream    : evidence.orient_sa_upstream,
                    precision_pr : vec![evidence.target_break_pr],
                    precision_sa : vec![evidence.target_break_sa],
                    aln_start_pr : vec![evidence.aln_start_pr],
                    aln_start_sa : vec![evidence.aln_start_sa],
                    support_sr   : 1,
                    support_sp   : 0,
                    sp_ratio     : 0 as f32,
                    cipos        : evidence.precision
                };
                pos_based_result.insert(chr_ab, vec![new_evidence]);

            
            }else{ 
                // now it becomes much more complicated
                // essentially we first check for any already existing elements
                // whether we have a similar one.
                // If not we will then add a new one 
                let mut found_elements = pos_based_result.get(&chr_ab).expect("ERROR: could not get the prev found position!").to_vec();
                let mut found_match = false;
                for elements in found_elements.iter_mut() {
                    if stranded {
                        if elements.chrom_pr == evidence.chrom_pr &&
                            elements.chrom_sa == evidence.chrom_sa &&
                            elements.strand_pr == evidence.strand_pr &&
                            elements.strand_sa == evidence.strand_sa && 
                            elements.orient_pr_upstream == evidence.orient_pr_upstream &&
                            elements.orient_sa_upstream == evidence.orient_sa_upstream && 
                            u32::try_from(
                                (evidence.target_break_pr.abs_diff(elements.target_break_pr)) * (evidence.target_break_pr.abs_diff(elements.target_break_pr))
                                .integer_sqrt()
                            ).expect("ERROR: could not convert range into u32!") < range && 
                            u32::try_from(
                                (evidence.target_break_sa.abs_diff(elements.target_break_sa)) * (evidence.target_break_sa.abs_diff(elements.target_break_sa))
                                .integer_sqrt()
                            ).expect("ERROR: could not convert range into u32!") < range 
                        {
                            elements.precision_pr.push(evidence.target_break_pr);
                            elements.precision_sa.push(evidence.target_break_sa);
                            elements.aln_start_pr.push(evidence.aln_start_pr);
                            elements.aln_start_sa.push(evidence.aln_start_sa);
                            elements.support_sr +=1;
                            // now we should as well again
                            // correct potentially our one-value fix primary and supplementary 
                            // target break-point information.
                            // Otherwise we might have by chance the 1st value 
                            // being most unprecise and determinig the output value
                            // Therefore we provide the array, as we now
                            // can calculate a median on observation:
                            elements.target_break_pr = median_u64(&mut elements.precision_pr);
                            elements.target_break_sa = median_u64(&mut elements.precision_sa);
                            found_match = true;
                        }
                    }else if 
                        (
                            elements.chrom_pr == evidence.chrom_pr &&
                            elements.chrom_sa == evidence.chrom_sa 
                        ) & (
                            (
                                elements.strand_pr != evidence.strand_pr &&
                                elements.strand_sa != evidence.strand_sa && 
                                elements.orient_pr_upstream != evidence.orient_pr_upstream &&
                                elements.orient_sa_upstream != evidence.orient_sa_upstream 
                            ) | (
                                elements.strand_pr == evidence.strand_pr &&
                                elements.strand_sa == evidence.strand_sa && 
                                elements.orient_pr_upstream == evidence.orient_pr_upstream &&
                                elements.orient_sa_upstream == evidence.orient_sa_upstream 
                            )
                        ) & (
                            u32::try_from(( (evidence.target_break_pr.abs_diff(elements.target_break_pr)) * (evidence.target_break_pr.abs_diff(elements.target_break_pr)) ).integer_sqrt()).expect("ERROR: could not cast range into u32!") < range && 
                            u32::try_from(( (evidence.target_break_sa.abs_diff(elements.target_break_sa)) * (evidence.target_break_sa.abs_diff(elements.target_break_sa)) ).integer_sqrt()).expect("ERROR: could not cast range into u32!") < range 
                        )
                    {
                        elements.precision_pr.push(evidence.target_break_pr);
                        elements.precision_sa.push(evidence.target_break_sa);
                        elements.aln_start_pr.push(evidence.aln_start_pr);
                        elements.aln_start_sa.push(evidence.aln_start_sa);
                        elements.support_sr +=1;
                        // now we should as well again
                        // correct potentially our one-value fix primary and supplementary 
                        // target break-point information.
                        // Otherwise we might have by chance the 1st value 
                        // being most unprecise and determinig the output value
                        // Therefore we provide the array, as we now
                        // can calculate a median on observation:
                        elements.target_break_pr = median_u64(&mut elements.precision_pr);
                        elements.target_break_sa = median_u64(&mut elements.precision_sa);
                        found_match = true;
                    }
                }
                if !found_match {
                    // now we populate our new format and then add this
                    // in the new structure
                    let new_evidence = PosBasedFusEvidence{
                        chrom_pr     : evidence.chrom_pr,
                        chrom_sa     : evidence.chrom_sa,
                        target_break_pr: evidence.target_break_pr,
                        target_break_sa: evidence.target_break_sa,
                        strand_pr    : evidence.strand_pr,
                        strand_sa    : evidence.strand_sa,
                        orient_pr_upstream    : evidence.orient_pr_upstream,
                        orient_sa_upstream    : evidence.orient_sa_upstream,
                        precision_pr : vec![evidence.target_break_pr],
                        precision_sa : vec![evidence.target_break_sa],
                        aln_start_pr : vec![evidence.aln_start_pr],
                        aln_start_sa : vec![evidence.aln_start_sa],
                        support_sr   : 1,
                        support_sp   : 0,
                        sp_ratio     : 0 as f32,
                        cipos        : evidence.precision
                    };
                    found_elements.push(new_evidence);
                }
                debug!("Adding new element for chr {:?} : {:?}",&chr_ab,&found_elements);
                pos_based_result.insert(chr_ab, found_elements);
            }
        }
    }
    pos_based_result
}


/// this function takes our vector of final output
/// and writes it in a tab-separated format to stdout
/// it adds additionally a header with meta information
/// to identify more easily later the context.
/// The information is returned in a 1-based position and 
/// the fusion point indicated is always the last nucleotide 
/// which is still "normal" with the new sequence adjacent to it.
///
/// Unittest: FALSE
///s
pub fn write_read_based_tsv_stdout(
    // 0-based results
    results: &[PrincipleReadBasedFusionOutput], 
    version: &str, 
    author: &str, 
    command: &str
) -> Result<(), Box<dyn Error>>{
    
    eprintln!("INFO: writing output in tab-separated format to stdout");
    let now: DateTime<Local> = Local::now();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::stdout());
    // here we add some empty fields as otherwise the writere realizes that this is wrong    
    writer.write_record(["##","splitty:",version,"","","","","","","","","","","",""])?;
    writer.write_record(["##","author:",author,"","","","","","","","","","","",""])?;
    writer.write_record(["##","date:",&now.to_rfc2822(),"","","","","","","","","","","",""])?;
    writer.write_record(["##","command:", command,"","","","","","","","","","","",""])?;    
    writer.write_record([
        "#contig",
        "contig_length",
        "fusion_point",
        "fusion_genes",
        "precision",
        "chrA",
        "startA",
        "endA",
        "strandA",
        "fpA",
        "chrB",
        "startB",
        "endB",
        "strandB",
        "fpB"
        ])?;
    for entry in results.iter() {
        let fusion_fuzzy = match entry.fp_fuzzy{
            Some(x) => x.to_string(),
            None => String::from("NA"),
        };
        let strand1 = match entry.a_strand {
            StrandDirection::Fwd => "1",
            StrandDirection::Rev => "-1",
            StrandDirection::Unknown => ".",
        };
        let strand2 = match entry.b_strand{
            StrandDirection::Fwd => "1",
            StrandDirection::Rev => "-1",
            StrandDirection::Unknown => ".",
        };
        writer.write_record(&[
            entry.query_name.clone(),
            entry.length.to_string(),
            (entry.fusion_point +1).to_string(),
            entry.fusion_genes.join(";"),
            fusion_fuzzy,
            entry.a_chrom.clone(),
            (entry.a_start +1).to_string(),
            (entry.a_end +1).to_string(),
            strand1.to_string(),
            (entry.a_fp +1).to_string(),
            entry.b_chrom.clone(),
            (entry.b_start +1).to_string(),
            (entry.b_end +1).to_string(),
            strand2.to_string(),
            (entry.b_fp +1).to_string()
            ])?;
    };
    writer.flush()?;
    Ok(())
}



/// this function takes our vector of final output
/// and writes it in a tab-separated format to stdout
/// it adds additionally a header with meta information
/// to identify more easily later the context
/// The information is returned in a 0-based position and 
/// the fusion point indicated is always the last nucleotide 
/// which is still "normal" with the new sequence adjacent to it.
///
/// Unittest: FALSE
///
pub fn write_pos_based_tsv_stdout(
    results: &[PrinciplePosBasedFusionOutput], 
    version: &str, 
    author: &str, 
    command: &str
) -> Result<(), Box<dyn Error>>{
    
    eprintln!("INFO: writing output in tab-separated format to stdout");
    let now: DateTime<Local> = Local::now();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(io::stdout());
    // here we add some empty fields as otherwise the writere realizes that this is wrong    
    writer.write_record(["##","splitty:",version,"","","","","","","","","","","","",""])?;
    writer.write_record(["##","author:",author,"","","","","","","","","","","","",""])?;
    writer.write_record(["##","date:",&now.to_rfc2822(),"","","","","","","","","","","","",""])?;
    writer.write_record(["##","command:", command,"","","","","","","","","","","","",""])?;    
    writer.write_record([
        "#fusion_genes",
        "fp_fuzzy",
        "fusion_distance",
        "chrA",
        "startA",
        "endA",
        "strandA",
        "fpA",
        "chrB",
        "startB",
        "endB",
        "strandB",
        "fpB",
        "sr_support",
        "sp_support",
        "sp_ratio",
        ])?;
    for entry in results.iter() {
        let fusion_distance = match entry.fusion_distance{
            Some(x) => x.to_string(),
            None => String::from("NA"),
        };
        let fusion_fuzzy = match entry.fp_fuzzy{
            Some(x) => x.to_string(),
            None => String::from("NA"),
        };
        let strand1 = match entry.a_strand {
            StrandDirection::Fwd => "1",
            StrandDirection::Rev => "-1",
            StrandDirection::Unknown => ".",
        };
        let strand2 = match entry.b_strand{
            StrandDirection::Fwd => "1",
            StrandDirection::Rev => "-1",
            StrandDirection::Unknown => ".",
        };
        writer.write_record(&[
            entry.fusion_genes.join(";"),
            fusion_fuzzy,
            fusion_distance,
            entry.a_chrom.clone(),
            (entry.a_start +1).to_string(),
            (entry.a_end +1).to_string(),
            strand1.to_string(),
            (entry.a_fp +1).to_string(),
            entry.b_chrom.clone(),
            (entry.b_start +1).to_string(),
            (entry.b_end +1).to_string(),
            strand2.to_string(),
            (entry.b_fp +1).to_string(),
            entry.sr_support.to_string(),
            entry.sp_support.to_string(),
            format!{"{:.3}", entry.sp_ratio},
            ])?;
    };
    writer.flush()?;
    Ok(())
}



#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    //use bio::io::fasta::IndexedReader;  
    use std::fs::File;
    //use std::str::from_utf8;
    use tempfile::NamedTempFile;

    /////////////////////////////////////////
    ///       CIGAR EVALUATION   ////////////
    /////////////////////////////////////////
    /// 
    #[test]
    fn cigar_eval_test1(){
        // this one will fail because we have no clipping at all
        let cigar            = String::from("18M3D2M2D2M1I22M");
        let result           = eval_cigar(&cigar,1_u64,&1_u64,&StrandDirection::Fwd,None);
        assert_eq!(result,None);
        // let mut truth        = CigarInfos::default();
        // truth.seq_length     = 45_u32;
        // truth.match_length   = 44_u32;
        // truth.match_length_s = 49_u32;
        // truth.similarity     = 91.489365_f32;
        // truth.match_strand   = StrandDirection::Fwd;
        // assert_eq!(result.unwrap(),truth);
    }
    #[test]
    fn cigar_eval_test2(){
        let cigar2             = String::from("2172H7M73I19407N116M2836N102M2789N131M2751N101M9693N300M");
        let result2            = eval_cigar(&cigar2,0_u64,&68346459_u64,&StrandDirection::Rev,None);
        let mut truth2         = CigarInfos::default();
        truth2.fp_target       = 68346459_u64;
        truth2.aln_start       = 68346459_u64;
        truth2.aln_end         = 68384691;
        truth2.upstream        = 2172_u32;
        truth2.seq_length      = 3002_u32;
        truth2.match_length    = 757_u32;
        truth2.match_length_s  = 38233_u32;
        truth2.similarity      = 91.314835_f32;
        truth2.fp_query        = 829_u32;
        truth2.match_strand    = StrandDirection::Rev;
        assert_eq!(result2.unwrap(),truth2);
    }
    #[test]
    fn cigar_eval_test6(){
        // NAA38|ENSG00000128534.3--POLR2J|ENSG00000005075.11
        let cigar6            = String::from("1539M2583N129M831S");
        let result6           = eval_cigar(&cigar6,0_u64,&117824209_u64,&StrandDirection::Fwd,None);
        let mut truth6        = CigarInfos::default();
        truth6.fp_target      = 117828459_u64;
        truth6.aln_start      = 117824209_u64;
        truth6.aln_end        = 117828459;
        truth6.upstream       = 0_u32;
        truth6.downstream     = 831_u32;
        truth6.seq_length     = 2499_u32;
        truth6.match_length   = 1668_u32;
        truth6.match_length_s = 4251_u32;
        truth6.similarity     = 100_f32;
        truth6.fp_query       = 1667_u32;
        truth6.match_strand   = StrandDirection::Fwd;
        assert_eq!(result6.unwrap(),truth6);
    }
    #[test]
    fn cigar_eval_test7(){
        // AC011530.4|ENSG00000268434.1--IL4|ENSG00000113520.6
        let cigar7            = String::from("102S120M");
        let result7           = eval_cigar(&cigar7,0_u64,&46288851_u64,&StrandDirection::Rev,None);
        let mut truth7        = CigarInfos::default();
        truth7.fp_target      = 46288851_u64;
        truth7.aln_start      = 46288851_u64;
        truth7.aln_end        = 46288970;
        truth7.upstream       = 102_u32;
        truth7.downstream     = 0_u32;
        truth7.seq_length     = 222_u32;
        truth7.match_length   = 120_u32;
        truth7.match_length_s = 120_u32;
        truth7.similarity     = 100_f32;
        truth7.fp_query       = 119_u32;
        truth7.match_strand   = StrandDirection::Rev;
        assert_eq!(result7.unwrap(),truth7);
    }
    #[test]
    fn cigar_eval_test8(){
        // AC011530.4|ENSG00000268434.1--IL4|ENSG00000113520.6
        let cigar8            = String::from("119H103M");
        let result8           = eval_cigar(&cigar8,0_u64,&132018176_u64,&StrandDirection::Fwd,None);
        let mut truth8        = CigarInfos::default();
        truth8.fp_target      = 132018176_u64;
        truth8.aln_start      = 132018176_u64;
        truth8.aln_end        = 132018278;
        truth8.upstream       = 119_u32;
        truth8.downstream     = 0_u32;
        truth8.seq_length     = 222_u32;
        truth8.match_length   = 103_u32;
        truth8.match_length_s = 103_u32;
        truth8.similarity     = 100_f32;
        truth8.fp_query       = 119_u32;
        truth8.match_strand   = StrandDirection::Fwd;
        assert_eq!(result8.unwrap(),truth8);
    }
    #[test]
    fn cigar_eval_test9(){
        // ATP13A4|ENSG00000127249.10--CAMLG|ENSG00000164615.3
        let cigar9            = String::from("619S125M5770N112M2787N154M225N87M1615N151M3530N62M2116N189M896N158M1133N33M");
        let result9           = eval_cigar(&cigar9,0_u64,&193165994_u64,&StrandDirection::Rev,None);
        let mut truth9        = CigarInfos::default();
        truth9.fp_target      = 193165994_u64;
        truth9.aln_start      = 193165994_u64;
        truth9.aln_end        = 193185136;
        truth9.upstream       = 619_u32;
        truth9.downstream     = 0_u32;
        truth9.seq_length     = 1690_u32;
        truth9.match_length   = 1071_u32;
        truth9.match_length_s = 19143_u32;
        truth9.similarity     = 100_f32;
        truth9.fp_query       = 1070_u32;
        truth9.match_strand   = StrandDirection::Rev;
        assert_eq!(result9.unwrap(),truth9);
    }
    #[test]
    fn cigar_eval_test10(){
        // ATP13A4|ENSG00000127249.10--CAMLG|ENSG00000164615.3
        let cigar10            = String::from("1066H624M");
        let result10           = eval_cigar(&cigar10,0_u64,&134086443_u64,&StrandDirection::Fwd,None);
        let mut truth10        = CigarInfos::default();
        truth10.fp_target      = 134086443_u64;
        truth10.aln_start      = 134086443_u64;
        truth10.aln_end        = 134087066;
        truth10.upstream       = 1066_u32;
        truth10.downstream     = 0_u32;
        truth10.seq_length     = 1690_u32;
        truth10.match_length   = 624_u32;
        truth10.match_length_s = 624_u32;
        truth10.similarity     = 100_f32;
        truth10.fp_query       = 1066_u32;
        truth10.match_strand   = StrandDirection::Fwd;
        assert_eq!(result10.unwrap(),truth10);
    }


    /////////////////////////////////////////////
    ///       Alignment properties   ////////////
    /////////////////////////////////////////////
    #[test]
    fn calc_align_test1(){
        // so here
        let test1  = calc_align_properties(2172_u32,0_u32,3002_u32,38233_u32,68346459_u64,StrandDirection::Rev,None);
        let truth1 = BreakPoints {
                        break_query: 829_u32,
                        break_target: 68346459_u64,
                        aln_start: 68346459_u64,
                        aln_end: 68384691_u64,
                        strand: StrandDirection::Rev 
                    };
        assert_eq!(test1.unwrap(),truth1);
    }
    #[test]
    fn calc_align_test2(){
        // start 0, slength 150, mlength_spliced 58, clip5 92, clip3 0
        let test2  = calc_align_properties(92_u32,0_u32,150_u32,58_u32,0_u64,StrandDirection::Rev,None);
        let truth2 = BreakPoints {
                        break_query: 57_u32,
                        break_target: 0_u64,
                        aln_start: 0_u64,
                        aln_end: 57_u64,
                        strand: StrandDirection::Rev 
                    };
        assert_eq!(test2.unwrap(),truth2);
    }
    #[test]
    fn calc_align_test3(){
        let test3  = calc_align_properties(0_u32,25_u32,75_u32,50_u32,0_u64,StrandDirection::Fwd,None);
        let truth3 = BreakPoints {
                        break_query: 49_u32,
                        break_target: 49_u64,
                        aln_start: 0_u64,
                        aln_end: 49_u64,
                        strand: StrandDirection::Fwd 
                    };
        assert_eq!(test3.unwrap(),truth3);
    }
    #[test]
    fn calc_align_test4(){
        let test4  = calc_align_properties(50_u32,0_u32,75_u32,25_u32,49_u64,StrandDirection::Rev,None);
        let truth4 = BreakPoints {
                        break_query: 24_u32,
                        break_target: 49_u64,
                        aln_start: 49_u64,
                        aln_end: 73_u64,
                        strand: StrandDirection::Rev 
                    };
        assert_eq!(test4.unwrap(),truth4);
    }
    #[test]
    fn calc_align_test5(){
        // same as above but checking if both 5' and 3' clipping
        let test5  = calc_align_properties(50_u32,5_u32,75_u32,25_u32,49_u64,StrandDirection::Rev,None);
        let truth5 = BreakPoints {
                        break_query: 24_u32,
                        break_target: 49_u64,
                        aln_start: 49_u64,
                        aln_end: 73_u64,
                        strand: StrandDirection::Rev 
                    };
        assert_eq!(test5.unwrap(),truth5);
    }
    #[test]
    fn calc_align_test6(){
        // same as above but checking if both 5' and 3' clipping with a 10% inbalance ratio
        let test6  = calc_align_properties(50_u32,5_u32,75_u32,25_u32,49_u64,StrandDirection::Rev,Some(0.9_f32));
        let truth6 = BreakPoints {
                        break_query: 24_u32,
                        break_target: 49_u64,
                        aln_start: 49_u64,
                        aln_end: 73_u64,
                        strand: StrandDirection::Rev 
                    };
        assert_eq!(test6.unwrap(),truth6);
    }
    #[test]
    fn calc_align_test7(){
        // same as above but checking if both 5' and 3' clipping with a <10% inbalance ratio which should result in None
        let test7  = calc_align_properties(50_u32,5_u32,75_u32,25_u32,49_u64,StrandDirection::Rev,Some(0.05_f32));
        let truth7 : Option<BreakPoints> = None;
        assert_eq!(test7,truth7);
    }
    

    /////////////////////////////////////////////
    ///       construct 2nd allele   ////////////
    /////////////////////////////////////////////
    #[test]
    fn construct_2ndallele_test1(){
        let ff = String::from("G]chr16:18804692]");
        let ff_check = Bnd2ndAllele {
            nucleotide: String::from("G"),
            chr: String::from("chr16"),
            pos: 18804691_u64,
            prim_direction: StrandDirection::Fwd,
            second_direction: StrandDirection::Fwd,
        };
        let ff_result = construct_2nd_allele(&ff_check,true);
        assert_eq!(ff_result,ff);
    }
    #[test]
    fn construct_2ndallele_test2(){
        let rr = String::from("[chr19:50193168[A");
        let rr_check = Bnd2ndAllele {
            nucleotide: String::from("A"),
            chr: String::from("chr19"),
            pos: 50193167_u64,
            prim_direction: StrandDirection::Rev,
            second_direction: StrandDirection::Rev,
        };
        let rr_result = construct_2nd_allele(&rr_check,true);
        assert_eq!(rr_result,rr);
    }
    #[test]
    fn construct_2ndallele_test3(){
        let fr = String::from("G[chr16:18804692[");
        let fr_check = Bnd2ndAllele {
            nucleotide: String::from("G"),
            chr: String::from("chr16"),
            pos: 18804691_u64,
            prim_direction: StrandDirection::Fwd,
            second_direction: StrandDirection::Rev,
        };
        let fr_result = construct_2nd_allele(&fr_check,true);
        assert_eq!(fr_result,fr);
    }
    #[test]
    fn construct_2ndallele_test4(){
        let rf = String::from("]chr19:50193168]A");
        let rf_check = Bnd2ndAllele {
            nucleotide: String::from("A"),
            chr: String::from("chr19"),
            pos: 50193167_u64,
            prim_direction: StrandDirection::Rev,
            second_direction: StrandDirection::Fwd,
        };
        let rf_result = construct_2nd_allele(&rf_check,true);
        assert_eq!(rf_result,rf);
    }

    /////////////////////////////////////////////
    ///     deconstruct 2nd allele   ////////////
    /////////////////////////////////////////////
    #[test]
    fn deconstruct_2ndallele_test1(){
        let ff = String::from("G]chr16:18804692]");
        let ff_check = Bnd2ndAllele {
            nucleotide: String::from("G"),
            chr: String::from("chr16"),
            pos: 18804691_u64,
            prim_direction: StrandDirection::Fwd,
            second_direction: StrandDirection::Fwd,
        };
        let ff_result = deconstruct_2nd_allele(&ff,true);
        assert_eq!(ff_result,ff_check);
    }

    #[test]
    fn deconstruct_2ndallele_test2(){
        let rr = String::from("[chr19:50193168[A");
        let rr_check = Bnd2ndAllele {
            nucleotide: String::from("A"),
            chr: String::from("chr19"),
            pos: 50193167_u64,
            prim_direction: StrandDirection::Rev,
            second_direction: StrandDirection::Rev,
        };
        let rr_result = deconstruct_2nd_allele(&rr,true);
        assert_eq!(rr_result,rr_check);
    }

    #[test]
    fn deconstruct_2ndallele_test3(){
        let fr = String::from("G[chr16:18804692[");
        let fr_check = Bnd2ndAllele {
            nucleotide: String::from("G"),
            chr: String::from("chr16"),
            pos: 18804691_u64,
            prim_direction: StrandDirection::Fwd,
            second_direction: StrandDirection::Rev,
        };
        let fr_result = deconstruct_2nd_allele(&fr,true);
        assert_eq!(fr_result,fr_check);
    }

    #[test]
    fn deconstruct_2ndallele_test4(){
        let rf = String::from("]chr19:50193168]A");
        let rf_check = Bnd2ndAllele {
            nucleotide: String::from("A"),
            chr: String::from("chr19"),
            pos: 50193167_u64,
            prim_direction: StrandDirection::Rev,
            second_direction: StrandDirection::Fwd,
        }; 
        let rf_result = deconstruct_2nd_allele(&rf,true);
        assert_eq!(rf_result,rf_check);
    }

    #[test]
    fn deconstruct_2ndallele_test5(){
        let f  = String::from("A.");
        let f_check = Bnd2ndAllele {
            nucleotide: String::from("A"),
            chr: String::new(),
            pos: 0_u64,
            prim_direction: StrandDirection::Fwd,
            second_direction: StrandDirection::Unknown,
        };
        let f_result  = deconstruct_2nd_allele(&f,false);
        assert_eq!(f_result ,f_check);
    }

    #[test]
    fn deconstruct_2ndallele_test6(){
        let r  = String::from(".C");
        let r_check = Bnd2ndAllele {
            nucleotide: String::from("C"),
            chr: String::new(),
            pos: 0_u64,
            prim_direction: StrandDirection::Rev,
            second_direction: StrandDirection::Unknown,
        };
        let r_result  = deconstruct_2nd_allele(&r,false);
        assert_eq!(r_result ,r_check);
    }

    /////////////////////////////////////////////
    ///     extract and manip Seq    ////////////
    /////////////////////////////////////////////
    #[test]
    fn extract_test1(){ 
        let mut file_seq =  File::create("foo1.fa").expect("ERROR:could not create seq file!");
        file_seq.write_all(b">chr1\nGTAGGCTGAAAA\nCCCC").expect("ERROR:could not write seq file!");
        // index as well
        let mut file_idx =  File::create("foo1.fa.fai").expect("ERROR:could not create idx file!");
        file_idx.write_all(b"chr1\t16\t6\t12\t13").expect("ERROR:could not write idx file!");

        let contig:  &str = &String::from("chr1");
        let start: u64 = 0;  // start is 0-based, inclusive
        let stop : u64 = 10; // stop is 0-based, exclusive
        
        let result = fasta_extract_and_correct(&String::from("foo1.fa"),contig,start,stop,DnaReturnSense::Fwd);
        assert_eq!(result, String::from("GTAGGCTGAA"));

        // cleaning up the fasta and index file
        std::fs::remove_file("foo1.fa").expect("ERROR: could not remove the temporary fasta file!");
        std::fs::remove_file("foo1.fa.fai").expect("ERROR: could not remove the temporary fasta index file!");
    }

    #[test]
    fn extract_test2(){ 
        let mut file_seq =  File::create("foo2.fa").expect("ERROR:could not create seq file!");
        file_seq.write_all(b">chr1\nGTAGGCTGAAAA\nCCCC").expect("ERROR:could not write seq file!");
        // index as well
        let mut file_idx =  File::create("foo2.fa.fai").expect("ERROR:could not create idx file!");
        file_idx.write_all(b"chr1\t16\t6\t12\t13").expect("ERROR:could not write idx file!");

        let contig:  &str = &String::from("chr1");
        let start: u64 = 0;  // start is 0-based, inclusive
        let stop : u64 = 10; // stop is 0-based, exclusive
        
        // lets assume somebody provides reversed coordinates
        let result = fasta_extract_and_correct(&String::from("foo2.fa"),contig,stop,start,DnaReturnSense::Fwd);
        assert_eq!(result, String::from("GTAGGCTGAA"));

        // cleaning up the fasta and index file
        std::fs::remove_file("foo2.fa").expect("ERROR: could not remove the temporary fasta file!");
        std::fs::remove_file("foo2.fa.fai").expect("ERROR: could not remove the temporary fasta index file!");
        
    }

    #[test]
    fn extract_test3(){ 
        let mut file_seq =  File::create("foo3.fa").expect("ERROR:could not create seq file!");
        file_seq.write_all(b">chr1\nGTAGGCTGAAAA\nCCCC").expect("ERROR:could not write seq file!");
        // index as well
        let mut file_idx =  File::create("foo3.fa.fai").expect("ERROR:could not create idx file!");
        file_idx.write_all(b"chr1\t16\t6\t12\t13").expect("ERROR:could not write idx file!");

        let contig:  &str = &String::from("chr1");
        let start: u64 = 0;  // start is 0-based, inclusive
        let maxn : u64 = 16; // maximum
        
        // lets try last base
        let result = fasta_extract_and_correct(&String::from("foo3.fa"),contig,start,maxn,DnaReturnSense::Fwd);
        assert_eq!(result, String::from("GTAGGCTGAAAACCCC"));

        // cleaning up the fasta and index file
        std::fs::remove_file("foo3.fa").expect("ERROR: could not remove the temporary fasta file!");
        std::fs::remove_file("foo3.fa.fai").expect("ERROR: could not remove the temporary fasta index file!");
    }

    #[test]
    fn extract_test4(){ 
        let mut file_seq =  File::create("foo4.fa").expect("ERROR:could not create seq file!");
        file_seq.write_all(b">chr1\nGTAGGCTGAAAA\nCCCC").expect("ERROR:could not write seq file!");
        // index as well
        let mut file_idx =  File::create("foo4.fa.fai").expect("ERROR:could not create idx file!");
        file_idx.write_all(b"chr1\t16\t6\t12\t13").expect("ERROR:could not write idx file!");

        let contig:  &str = &String::from("chr1");
        let start: u64 = 0;  // start is 0-based, inclusive
        let far  : u64 = 20; // this is too far
        
        // lets try a out of range coordinate
        let result = fasta_extract_and_correct(&String::from("foo4.fa"),contig,start,far,DnaReturnSense::Fwd);
        assert_eq!(result, String::from("GTAGGCTGAAAACCCC"));

        // cleaning up the fasta and index file
        std::fs::remove_file("foo4.fa").expect("ERROR: could not remove the temporary fasta file!");
        std::fs::remove_file("foo4.fa.fai").expect("ERROR: could not remove the temporary fasta index file!");
    }

    #[test]
    fn extract_test5(){ 
        let mut file_seq =  File::create("foo5.fa").expect("ERROR:could not create seq file!");
        file_seq.write_all(b">chr1\nGTAGGCTGAAAA\nCCCC").expect("ERROR:could not write seq file!");
        // index as well
        let mut file_idx =  File::create("foo5.fa.fai").expect("ERROR:could not create idx file!");
        file_idx.write_all(b"chr1\t16\t6\t12\t13").expect("ERROR:could not write idx file!");

        let contig:  &str = &String::from("chr1");
        let start: u64 = 0;  // start is 0-based, inclusive
        let far  : u64 = 20; // this is too far
        
        // lets try a out of range coordinate and complementary
        let result = fasta_extract_and_correct(&String::from("foo5.fa"),contig,start,far,DnaReturnSense::FwdC);
        assert_eq!(result, String::from("CATCCGACTTTTGGGG"));    
        // cleaning up the fasta and index file
        std::fs::remove_file("foo5.fa").expect("ERROR: could not remove the temporary fasta file!");
        std::fs::remove_file("foo5.fa.fai").expect("ERROR: could not remove the temporary fasta index file!");    
    }

    #[test]
    fn extract_test6(){ 
        let mut file_seq =  File::create("foo6.fa").expect("ERROR:could not create seq file!");
        file_seq.write_all(b">chr1\nGTAGGCTGAAAA\nCCCC").expect("ERROR:could not write seq file!");
        // index as well
        let mut file_idx =  File::create("foo6.fa.fai").expect("ERROR:could not create idx file!");
        file_idx.write_all(b"chr1\t16\t6\t12\t13").expect("ERROR:could not write idx file!");

        let contig:  &str = &String::from("chr1");
        let start: u64 = 0;  // start is 0-based, inclusive
        let stop : u64 = 10; // stop is 0-based, exclusive
        
        // now  reversing the sequence
        let result = fasta_extract_and_correct(&String::from("foo6.fa"),contig,start,stop,DnaReturnSense::Rev);
        assert_eq!(result, String::from("AAGTCGGATG"));
        // cleaning up the fasta and index file
        std::fs::remove_file("foo6.fa").expect("ERROR: could not remove the temporary fasta file!");
        std::fs::remove_file("foo6.fa.fai").expect("ERROR: could not remove the temporary fasta index file!");
    }

    #[test]
    fn extract_test7(){ 
        let mut file_seq =  File::create("foo7.fa").expect("ERROR:could not create seq file!");
        file_seq.write_all(b">chr1\nGTAGGCTGAAAA\nCCCC").expect("ERROR:could not write seq file!");
        // index as well
        let mut file_idx =  File::create("foo7.fa.fai").expect("ERROR:could not create idx file!");
        file_idx.write_all(b"chr1\t16\t6\t12\t13").expect("ERROR:could not write idx file!");

        let contig:  &str = &String::from("chr1");
        let start: u64 = 0;  // start is 0-based, inclusive
        let stop : u64 = 10; // stop is 0-based, exclusive
        
        // now  reversing complementing the sequence
        let result = fasta_extract_and_correct(&String::from("foo7.fa"),contig,start,stop,DnaReturnSense::RevC);
        assert_eq!(result, String::from("TTCAGCCTAC"));
        
        // cleaning up the fasta and index file
        std::fs::remove_file("foo7.fa").expect("ERROR: could not remove the temporary fasta file!");
        std::fs::remove_file("foo7.fa.fai").expect("ERROR: could not remove the temporary fasta index file!");
    }

    ///////////////////////
    /// parsing input /////
    ///////////////////////
    #[test]
    fn fofn_test0(){
        // testing multiple files
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file =  File::create(path).expect("ERROR:could not create  file!");
        file.write_all(b"file1.txt\nfile2.txt\n").expect("ERROR:could not write fofn file!");
        let result = parse_fofn(&path.to_str().unwrap(),true);
        let test = vec![String::from("file1.txt"),String::from("file2.txt")];
        assert_eq!(result,test);
    }
    #[test]
    #[should_panic]
    fn fofn_test1(){
        // testing multiple files, wrong formatted
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file =  File::create(path).expect("ERROR:could not create  file!");
        file.write_all(b"file1.txt wrong\nfile2.txt\tnot wrong\n").expect("ERROR:could not write fofn file!");
        parse_fofn(&path.to_str().unwrap(),true);
    }
    #[test]
    fn fofn_test2(){
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file =  File::create(path).expect("ERROR:could not create  file!");
        file.write_all(b"file1.txt\nfile2.txt\n").expect("ERROR:could not write fofn file!");
        let result = parse_fofn(&path.to_str().unwrap(),false);
        let test = vec![path.to_str().unwrap().to_string()];
        assert_eq!(result,test);
    }
    #[test]
    fn chrs_test0(){
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file =  File::create(path).expect("ERROR:could not create  file!");
        file.write_all(b"chr1\nchr2\n").expect("ERROR:could not write chrom file!");
        let result = parse_chroms_txt(&path.to_str().unwrap());
        let mut chroms : FxHashMap<String,u8> = FxHashMap::default();
        chroms.insert(String::from("chr1"),1);
        chroms.insert(String::from("chr2"),1);
        assert_eq!(result,chroms);
    }
    #[test]
    #[should_panic]
    fn chrs_test1(){
        //  wrong formatted
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file =  File::create(path).expect("ERROR:could not create  file!");
        file.write_all(b"chr1 wrong\nchr2\tnot wrong2\n").expect("ERROR:could not write chrom file!");
        parse_chroms_txt(&path.to_str().unwrap());
    }
    #[test]
    fn chrs2_test0(){
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file =  File::create(path).expect("ERROR:could not create  file!");
        file.write_all(b"chr1\t10\nchr2\t20\n").expect("ERROR:could not write chrom file!");
        let result = parse_chrom_file(&path.to_str().unwrap());
        let mut chroms : FxHashMap<String,u64> = FxHashMap::default();
        chroms.insert(String::from("chr1"),10);
        chroms.insert(String::from("chr2"),20);
        assert_eq!(result,chroms);
    }
    #[test]
    #[should_panic]
    fn chrs2_test1(){
        // no second column
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file =  File::create(path).expect("ERROR:could not create  file!");
        file.write_all(b"chr1\nchr2\n").expect("ERROR:could not write chrom file!");
        parse_chrom_file(&path.to_str().unwrap());
    }
    #[test]
    #[should_panic]
    fn chrs2_test2(){
        // not u64
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file =  File::create(path).expect("ERROR:could not create  file!");
        file.write_all(b"chr1\t-2\nchr2\t-1\n").expect("ERROR:could not write chrom file!");
        parse_chrom_file(&path.to_str().unwrap());
    }
    /////////////////////////////////////////////
    ///     muERV interleaved from BED    ///////
    /////////////////////////////////////////////
    #[test]
    fn muerv_interleaved_test1(){ 
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file_a =  File::create(path).expect("ERROR:could not create seq file!");
        file_a.write_all(b"m64143_211228_145337/1508655/ccs\t0\t1250\tmouse_merged,muERV_merged\t0\t.\n\
            m64143_211228_145337/1508655/ccs\t1275\t8961\thuman_merged\t-\tchr9@132857445-132865123\n").expect("ERROR:could not write BED file!");
        let result1  = read_2lines_bed(&path.to_str().unwrap(),"human","muERV",&1_u64);
        let result_a = AnnotatedRead {
            ref_start: 1275_u64,
            ref_end: 8961_u64,
            ref_strand: StrandDirection::Rev,
            ref_chr: String::from("chr9"),
            ref_gstart: 132857445_u64,
            ref_gend: 132865123_u64,
            read_name: String::from("m64143_211228_145337/1508655/ccs"),
            alien_start: 0_u64,
            alien_end: 1250_u64,
            integration_site: 132865123_u64,
            integration_dir: StrandDirection::Fwd,
        };
        assert_eq!(result1[0],result_a);
        // delete file
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn muerv_interleaved_test2(){ 
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file_a =  File::create(path).expect("ERROR:could not create seq file!");
        file_a.write_all(b"m64143_211228_145337/2622036/ccs\t0\t3205\thuman_merged\t-\tchr2@236659019-236662225\n\
            m64143_211228_145337/2622036/ccs\t3205\t8087\tmouse_merged,muERV_merged\t0\t.\n").expect("ERROR:could not write BED file!");
        let result1  = read_2lines_bed(&path.to_str().unwrap(),"human","muERV",&1_u64);
        let result_a = AnnotatedRead {
           ref_start: 0_u64,
           ref_end: 3205_u64,
           ref_strand: StrandDirection::Rev,
           ref_chr: String::from("chr2"),
           ref_gstart: 236659019_u64,
           ref_gend: 236662225_u64,
           read_name: String::from("m64143_211228_145337/2622036/ccs"),
           alien_start: 3205_u64,
           alien_end: 8087_u64,
           integration_site: 236659019_u64,
           integration_dir: StrandDirection::Rev,
        };
        assert_eq!(result1[0],result_a);
        // delete file
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn muerv_interleaved_test3(){ 
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file_c =  File::create(path).expect("ERROR:could not create seq file!");
        file_c.write_all(b"m64143_211228_145337/6423823/ccs\t37\t7201\thuman_merged\t+\tchr11@8837266-8844392\n\
            m64143_211228_145337/6423823/ccs\t7202\t7904\tmouse_merged,muERV_merged\t0\t.\n").expect("ERROR:could not write BED file!");        
        let result1  = read_2lines_bed(&path.to_str().unwrap(),"human","muERV",&1_u64);
        let result_b = AnnotatedRead {
            ref_start: 37_u64,
            ref_end: 7201_u64,
            ref_strand: StrandDirection::Fwd,
            ref_chr: String::from("chr11"),
            ref_gstart: 8837266_u64,
            ref_gend: 8844392_u64,
            read_name: String::from("m64143_211228_145337/6423823/ccs"),
            alien_start: 7202_u64,
            alien_end: 7904_u64,
            integration_site: 8844392_u64,
            integration_dir: StrandDirection::Fwd,
        };
        assert_eq!(result1[0],result_b);
        // delete file
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn muerv_interleaved_test4(){ 
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file_d =  File::create(path).expect("ERROR:could not create seq file!");
        file_d.write_all(b"m64143_211228_145337/14680935/ccs\t0\t3644\tmouse_merged,muERV_merged\t0\t.\n\
            m64143_211228_145337/14680935/ccs\t3646\t8956\thuman_merged\t+\tchr15@44444731-44450034\n").expect("ERROR:could not write BED file!");
        let result1  = read_2lines_bed(&path.to_str().unwrap(),"human","muERV",&1_u64);
        let result_b = AnnotatedRead {
            ref_start: 3646_u64,
            ref_end: 8956_u64,
            ref_strand: StrandDirection::Fwd,
            ref_chr: String::from("chr15"),
            ref_gstart: 44444731_u64,
            ref_gend: 44450034_u64,
            read_name: String::from("m64143_211228_145337/14680935/ccs"),
            alien_start: 0_u64,
            alien_end: 3644_u64,
            integration_site: 44444731_u64,
            integration_dir: StrandDirection::Rev,
        };
        assert_eq!(result1[0],result_b);
        // delete file
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn muerv_interleaved_test5(){ 
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file_e =  File::create(path).expect("ERROR:could not create seq file!");
        file_e.write_all(b"m64143_211027_074753/7669299/ccs\t21\t3288\thuman_merged\t-\tchr3@149151633-149154890\n\
            m64143_211027_074753/7669299/ccs\t3299\t10407\tmouse_merged,muERV_merged\t0\t.\n").expect("ERROR:could not write BED file!");
        let result1  = read_2lines_bed(&path.to_str().unwrap(),"human","muERV",&1_u64);
        let result_e = AnnotatedRead {
            ref_start: 21_u64,
            ref_end: 3288_u64,
            ref_strand: StrandDirection::Rev,
            ref_chr: String::from("chr3"),
            ref_gstart: 149151633_u64,
            ref_gend: 149154890_u64,
            read_name: String::from("m64143_211027_074753/7669299/ccs"),
            alien_start: 3299_u64,
            alien_end: 10407_u64,
            integration_site: 149151633_u64,
            integration_dir: StrandDirection::Rev,
        };
        assert_eq!(result1[0],result_e);
        // delete file
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn muerv_interleaved_test6(){ 
        let tmp  = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut file_f =  File::create(path).expect("ERROR:could not create seq file!");
        file_f.write_all(b"m64143_211105_211504/3080490/ccs\t0\t4057\thuman_merged\t+\tchr15@48988687-48992741\n\
            m64143_211105_211504/3080490/ccs\t4054\t12097\tmouse_merged,muERV_merged\t0\t.\n").expect("ERROR:could not write BED file!");
        let result6  = read_2lines_bed(&path.to_str().unwrap(),"human","muERV",&10_u64);
        let result7  = read_2lines_bed(&path.to_str().unwrap(),"human","muERV",&1_u64);
        let result_f = AnnotatedRead {
            ref_start: 0_u64,
            ref_end: 4057_u64,
            ref_strand: StrandDirection::Fwd,
            ref_chr: String::from("chr15"),
            ref_gstart: 48988687_u64,
            ref_gend: 48992741_u64,
            read_name: String::from("m64143_211105_211504/3080490/ccs"),
            alien_start: 4054_u64,
            alien_end: 12097_u64,
            integration_site: 48992741_u64,
            integration_dir: StrandDirection::Fwd,
        };
        let empty : Vec<AnnotatedRead>  = Vec::new();
        assert_eq!(result6[0],result_f);
        assert_eq!(result7,empty);
        // delete file
        std::fs::remove_file(path).unwrap();
    }

    ///////////////////////
    // fusion events //////
    ///////////////////////

    #[test]
    fn fusion_annotation(){
        let _read_collection : FxHashMap<String, Vec<PosBasedFusEvidence>> = FxHashMap::default();
        //"EU432099.1:0-530EU432099.1:531": [PosBasedFusEvidence { chrom_pr: "EU432099.1:0-530", chrom_sa: "EU432099.1:531", target_break_pr: 529, target_break_sa: 0, strand_pr: Rev, strand_sa: Rev, orient_pr: 3, orient_sa: 5, precision_pr: [525, 529, 529, 529, 529, 529, 529, 529, 529, 529, 529, 529, 529, 529, 529, 529, 529, 529, 529], precision_sa: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 4], aln_start_pr: [450, 496, 468, 467, 495, 438, 474, 470, 445, 468, 447, 455, 451, 455, 460, 491, 470, 426, 497], aln_start_sa: [69, 115, 87, 86, 114, 57, 93, 89, 64, 87, 66, 74, 70, 74, 79, 110, 89, 45, 116], support_sr: 19, support_sp: 60, sp_ratio: 1.7647059, cipos: Some(0) }]
        let _gtf : Option<&str> = None;
        let _stranded = false;
        let _max_dist = 500_u64;

    }
    
    #[test]
    // test to check simple basic functionality
    fn score_annot_overlap_1(){
        let q_start           = 500_u64;
        let q_end             = 1000_u64;
        let q_strand          = StrandDirection::Fwd;
        let t_start           = 750_u64;
        let t_end             = 1250_u64;
        let t_strand          = StrandDirection::Fwd; 
        let stranded          = true;
        let result = calc_match_overlap(q_start,q_end,q_strand,t_start,t_end,t_strand,stranded);
        let test   = OverlapResult {
            match_length : 250,
            match_score  : 16,
            feature_name : String::from(""),
        };
        assert_eq!(result,test);
    }

    #[test]
    // test based on reported bug from Nikos #64
    fn score_annot_overlap_2(){
        let q_start           = 0_u64;
        let q_end             = 5958_u64;
        let q_strand          = StrandDirection::Rev;
        let t_start           = 0_u64;
        let t_end             = 12793_u64;
        let t_strand          = StrandDirection::Rev; 
        let stranded          = true;
        let result = calc_match_overlap(q_start,q_end,q_strand,t_start,t_end,t_strand,stranded);
        let test   = OverlapResult {
            match_length : 5958,
            match_score  : 25,
            feature_name : String::from(""),
        };
        assert_eq!(result,test);
    }

    #[test]
    // test based on easy fragment fusion event
    fn test_orientation_fwd_rev(){

        let primary_infos =  AlignInfos { 
            query: String::from("transcript/196704;full_length_coverage=33;length=984;NpolyAremoved=0;noPolyA"),
            q_length: 984,
            cigar: String::from("7S79=891N76=5875N7=1X139=1452N95=580S"),
            chrom: String::from("chr7"),
            fp_target: 55200409,
            upstream: 7,
            downstream: 580,
            similarity: 99.74811,
            seq_length: 984,
            match_length_s: 8615,
            match_length:397,
            fp_query: 403,
            aln_start: 55191795,
            aln_end: 55200409,
            match_strand: StrandDirection::Fwd 
        };

        let mut sa_infos : Vec<AlignInfos> = Vec::new();

        let new_info = AlignInfos { 
            query: String::from("transcript/196704;full_length_coverage=33;length=984;NpolyAremoved=0;noPolyA"),
            q_length: 984,
            cigar: String::from("177S403=404S"),
            chrom: String::from("chr7"),
            fp_target: 55796094,
            upstream: 177,
            downstream: 404,
            similarity: 100.0,
            seq_length: 984,
            match_length_s: 403,
            match_length:403,
            fp_query: 404,
            aln_start: 55795692,
            aln_end: 55796094,
            match_strand: StrandDirection::Rev
         };
         sa_infos.push(new_info);

         let mut truth : Vec<FullFusionEvidence> = Vec::new();
         let truth_event = FullFusionEvidence { 
            cigar_pr: String::from("7S79=891N76=5875N7=1X139=1452N95=580S"),
            cigar_sa: String::from("177S403=404S"),
            chrom_pr: String::from("chr7"),
            chrom_sa: String::from("chr7"),
            query_name: String::from("transcript/196704;full_length_coverage=33;length=984;NpolyAremoved=0;noPolyA"),
            query_length: 984,
            target_break_pr: 55200409,// manually checked
            target_break_sa: 55796094,// manually checked
            query_break_pr: 403,
            query_break_sa: 404,
            strand_pr: StrandDirection::Fwd,
            strand_sa: StrandDirection::Rev,
            orient_pr_upstream: true,
            orient_sa_upstream: false,
            aln_start_pr: 55191795, // manually checked
            aln_start_sa: 55795692, // manually checked
            precision: None,
            support: 0 
        };
        
        truth.push(truth_event);

        let result = evaluate_and_predict(primary_infos,sa_infos);
        assert_eq!(truth,result);
    }
    
    #[test]
    // test based on easy fragment fusion event
    fn test_orientation_fwd_fwd(){

        let primary_infos =  AlignInfos { 
            query: String::from("transcript/192912;full_length_coverage=48;length=984;NpolyAremoved=0;noPolyA"),
            q_length: 984,
            cigar: String::from("409S101=511N147=413N131=189N19=177S"),
            chrom: String::from("chr1"),
            fp_target: 156874907,
            upstream: 409,
            downstream: 177,
            similarity: 100.0,
            seq_length: 984,
            match_length_s: 1511,
            match_length: 398,
            fp_query: 409,
            aln_start: 156874907,
            aln_end: 156876417,
            match_strand: StrandDirection::Fwd
        };

        let mut sa_infos : Vec<AlignInfos> = Vec::new();

        let new_info = AlignInfos { 
            query: String::from("transcript/192912;full_length_coverage=48;length=984;NpolyAremoved=0;noPolyA"),
            q_length: 984,
            cigar: String::from("7S243=15342N159=575S"),
            chrom: String::from("chr1"),
            fp_target: 156130774,
            upstream: 7,
            downstream: 575,
            similarity: 100.0,
            seq_length: 984,
            match_length_s: 15744,
            match_length: 402,
            fp_query: 408,
            aln_start: 156115031,
            aln_end: 156130774,
            match_strand: StrandDirection::Fwd
         };
         sa_infos.push(new_info);

         let mut truth : Vec<FullFusionEvidence> = Vec::new();
         let truth_event = FullFusionEvidence { 
            cigar_pr: String::from("409S101=511N147=413N131=189N19=177S"),
            cigar_sa: String::from("7S243=15342N159=575S"),
            chrom_pr: String::from("chr1"),
            chrom_sa: String::from("chr1"),
            query_name: String::from("transcript/192912;full_length_coverage=48;length=984;NpolyAremoved=0;noPolyA"),
            query_length: 984, 
            target_break_pr: 156874907, // manually checked
            target_break_sa: 156130774, // manually checked
            query_break_pr: 409, 
            query_break_sa: 408, 
            strand_pr: StrandDirection::Fwd,
            strand_sa: StrandDirection::Fwd,
            orient_pr_upstream: false,
            orient_sa_upstream: true,
            aln_start_pr: 156874907, // manually checked
            aln_start_sa: 156115031, // manually checked
            precision: None, 
            support: 0
        };
        truth.push(truth_event);

        let result = evaluate_and_predict(primary_infos,sa_infos);
        assert_eq!(truth,result);
    }
    
    #[test]
    fn test_orientation_ref_fwd(){

        let primary_infos =  AlignInfos { 
            query: String::from("transcript/210380;full_length_coverage=33;length=887;NpolyAremoved=0;noPolyA"),
            q_length: 887,
            cigar: String::from("577S216=1X86=7S"),
            chrom: String::from("chr10"),
            fp_target: 59906121,
            upstream: 577,
            downstream: 7,
            similarity: 99.66997,
            seq_length: 887,
            match_length_s: 303,
            match_length: 303,
            fp_query: 309,
            aln_start: 59906121,
            aln_end: 59906423,
            match_strand: StrandDirection::Rev
        };

        let mut sa_infos : Vec<AlignInfos> = Vec::new();

        let new_info = AlignInfos { 
            query: String::from("transcript/210380;full_length_coverage=33;length=887;NpolyAremoved=0;noPolyA"),
            q_length: 887,
            cigar: String::from("310S148=1641N22=1X85=1050N144=177S"),
            chrom: String::from("chr10"),
            fp_target: 43116583,
            upstream: 310,
            downstream: 177,
            similarity: 99.5,
            seq_length: 887,
            match_length_s: 3091,
            match_length: 400,
            fp_query: 310,
            aln_start: 43116583,
            aln_end: 43119673,
            match_strand: StrandDirection::Fwd
         };
         sa_infos.push(new_info);

         let mut truth : Vec<FullFusionEvidence> = Vec::new();
         let truth_event = FullFusionEvidence { 
            cigar_pr: String::from("577S216=1X86=7S"),
            cigar_sa: String::from("310S148=1641N22=1X85=1050N144=177S"),
            chrom_pr: String::from("chr10"),
            chrom_sa: String::from("chr10"),
            query_name: String::from("transcript/210380;full_length_coverage=33;length=887;NpolyAremoved=0;noPolyA"),
            query_length: 887,
            target_break_pr: 59906121,// manually checked
            target_break_sa: 43116583,// manually checked
            query_break_pr: 309,
            query_break_sa: 310,
            strand_pr: StrandDirection::Rev,
            strand_sa: StrandDirection::Fwd,
            orient_pr_upstream: true,
            orient_sa_upstream: false,
            aln_start_pr: 59906121,// manually checked
            aln_start_sa: 43116583,// manually checked
            precision: None,
            support: 0
        };
        truth.push(truth_event);

        let result = evaluate_and_predict(primary_infos,sa_infos);
        assert_eq!(truth,result);
    }

    #[test]
    fn test_orientation_rev_rev(){

        let primary_infos =  AlignInfos { 
            query: String::from("transcript/240115;full_length_coverage=52;length=662;NpolyAremoved=0;noPolyA"),
            q_length: 662,
            cigar: String::from("177S212=130040N86=9096N105=82S"),
            chrom: String::from("chr21"),
            fp_target: 38445409,
            upstream: 177,
            downstream: 82,
            similarity: 100.0,
            seq_length: 662,
            match_length_s: 139539,
            match_length: 403,
            fp_query: 484,
            aln_start: 38445409,
            aln_end: 38584947,
            match_strand: StrandDirection::Rev
        };

        let mut sa_infos : Vec<AlignInfos> = Vec::new();

        let new_info = AlignInfos { 
            query: String::from("transcript/240115;full_length_coverage=52;length=662;NpolyAremoved=0;noPolyA"),
            q_length: 662,
            cigar: String::from("580S75=7S"),
            chrom: String::from("chr21"),
            fp_target: 41508083,
            upstream: 580,
            downstream: 7,
            similarity: 100.0,
            seq_length: 662,
            match_length_s: 75,
            match_length: 75,
            fp_query: 81,
            aln_start: 41508083,
            aln_end: 41508157,
            match_strand: StrandDirection::Rev
         };
         sa_infos.push(new_info);

         let mut truth : Vec<FullFusionEvidence> = Vec::new();
         let truth_event = FullFusionEvidence { 
            cigar_pr: String::from("177S212=130040N86=9096N105=82S"),
            cigar_sa: String::from("580S75=7S"),
            chrom_pr: String::from("chr21"),
            chrom_sa: String::from("chr21"),
            query_name: String::from("transcript/240115;full_length_coverage=52;length=662;NpolyAremoved=0;noPolyA"),
            query_length: 662,
            target_break_pr: 38584947,// manually checked
            target_break_sa: 41508083,// manually checked
            query_break_pr: 484,
            query_break_sa: 81,
            strand_pr: StrandDirection::Rev,
            strand_sa: StrandDirection::Rev,
            orient_pr_upstream: false,
            orient_sa_upstream: true,
            aln_start_pr: 38445409, // manually checked
            aln_start_sa: 41508083, // manually checked
            precision: None,
            support: 0
        };
        truth.push(truth_event);

        let result = evaluate_and_predict(primary_infos,sa_infos);
        assert_eq!(truth,result);
    }
    
    
}
