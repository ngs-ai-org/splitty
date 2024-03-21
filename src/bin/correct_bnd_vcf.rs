use clap::{app_from_crate,crate_name,crate_description,crate_authors,crate_version,Arg};
use std::env;
//use bio::data_structures::interval_tree;
// our library which is within the same project
extern crate genefusion;
use genefusion::lib::common::{*};

// essentially we will do for the 1st vcf a first run population
// a hashmap were the keys are concatenated iD:mateID + the link to the vcf entry 
// together with a structure containging very simplified intormation 
// In a second step we do the same for the 2nd vcf but instead of the 2nd structure
// we compare it directly to the 1st strucure and add it accordingly 
use rust_htslib::bcf::header::Header;
use rust_htslib::bcf::{Read as BcfRead};
use rust_htslib::bcf::{Format, Writer};
use std::str;

#[derive(Debug)]
pub struct AlleleInfo {
    allele1: String,
    allele2: String,
}


#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf::{Format, Header, Writer};
    use tempfile::NamedTempFile;
    
    // "G]chr16:18804692]"  = FF 
    // "[chr19:50193168[A"  = RR 
    // "G[chr16:18804692["  = FR 
    // "]chr19:50193168]A"  = RF 

    #[test]
    fn test_ff_correct() {
        pretty_env_logger::init();

        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record1 = vcf.empty_record();
        let allele1: &[&[u8]] = &[b"G", b"G]chr7:137513870]"];
        record1.set_alleles(allele1).expect("Failed to set alleles");
        let mut record2 = vcf.empty_record();
        let allele2: &[&[u8]] = &[b"G", b"G]chr21:31668501]"];
        record2.set_alleles(allele2).expect("Failed to set alleles");
        let (new_allele1,new_allele2) = correct_alleles(&record1,&record2);
        assert_eq!(new_allele1.allele1,String::from("G"));
        assert_eq!(new_allele1.allele2,String::from("G]chr7:137513870]"));
        assert_eq!(new_allele2.allele1,String::from("G"));
        assert_eq!(new_allele2.allele2,String::from("G]chr21:31668501]"));
    }

    #[test]
    fn test_rr_correct() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record1 = vcf.empty_record();
        let allele1: &[&[u8]] = &[b"A", b"[chr13:29509670[A"];
        record1.set_alleles(allele1).expect("Failed to set alleles");
        let mut record2 = vcf.empty_record();
        let allele2: &[&[u8]] = &[b"G", b"[chr19:39476566[G"];
        record2.set_alleles(allele2).expect("Failed to set alleles");
        let (new_allele1,new_allele2) = correct_alleles(&record1,&record2);
        assert_eq!(new_allele1.allele1,String::from("A"));
        assert_eq!(new_allele1.allele2,String::from("[chr13:29509670[A"));
        assert_eq!(new_allele2.allele1,String::from("G"));
        assert_eq!(new_allele2.allele2,String::from("[chr19:39476566[G"));
    }
     

    #[test]
    fn test_fr_correct() {

        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record1 = vcf.empty_record();
        let allele1: &[&[u8]] = &[b"G", b"G]chr13:21176536]"];
        record1.set_alleles(allele1).expect("Failed to set alleles");
        let mut record2 = vcf.empty_record();
        let allele2: &[&[u8]] = &[b"C", b"[chr19:3053319[C"];
        record2.set_alleles(allele2).expect("Failed to set alleles");
        let (new_allele1,new_allele2) = correct_alleles(&record1,&record2);
        assert_eq!(new_allele1.allele1,String::from("G"));
        assert_eq!(new_allele1.allele2,String::from("G[chr13:21176536["));
        assert_eq!(new_allele2.allele1,String::from("C"));
        assert_eq!(new_allele2.allele2,String::from("]chr19:3053319]C"));
    }

    
    #[test]
    fn test_rf_correct() {

        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf = Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let mut record1 = vcf.empty_record();
        let allele1: &[&[u8]] = &[b"T", b"[chr12:555021[T"];
        record1.set_alleles(allele1).expect("Failed to set alleles");
        let mut record2 = vcf.empty_record();
        let allele2: &[&[u8]] = &[b"T", b"T]chr14:34800260]"];
        record2.set_alleles(allele2).expect("Failed to set alleles");
        let (new_allele1,new_allele2) = correct_alleles(&record1,&record2);
        assert_eq!(new_allele1.allele1,String::from("T"));
        assert_eq!(new_allele1.allele2,String::from("]chr12:555021]T"));
        assert_eq!(new_allele2.allele1,String::from("T"));
        assert_eq!(new_allele2.allele2,String::from("T[chr14:34800260["));
    }
}

/// function takes 2 associated entries and corrects 
/// the mate id orientation in the 2nd allele. It is based on the assumption 
/// that the primary orientation is correct and that the mate
/// one has to be re-oriented based on the 2nd primary entry
pub fn correct_alleles(
    prev_record: &rust_htslib::bcf::Record ,
    cur_record: &rust_htslib::bcf::Record 
) -> (AlleleInfo,AlleleInfo){
            let alleles1        = prev_record.alleles();
            let alleles2        = cur_record.alleles();
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
            let mut allele1_desconstr : Bnd2ndAllele = deconstruct_2nd_allele(allele_1,true);
            let mut allele2_desconstr : Bnd2ndAllele = deconstruct_2nd_allele(allele_2,true);
            
            // as mentioned initially, the information
            // of the mate-direction is currently wrong and we will
            // correct it based on the matching entry primary info
            allele1_desconstr.second_direction = allele2_desconstr.prim_direction;
            allele2_desconstr.second_direction = allele1_desconstr.prim_direction;
            
            // now we need to construct again our allele field accordingly:
            // carefully the function returns now a 0 based system
            let new_allele1 = construct_2nd_allele(&allele1_desconstr,true);
            let new_allele2 = construct_2nd_allele(&allele2_desconstr,true);

            // now we generate a basic information structure
            // with both alleles and return it
            let info1 = AlleleInfo {
                allele1: allele1_desconstr.nucleotide,
                allele2: new_allele1,
            };
            let info2 = AlleleInfo {
                allele1: allele2_desconstr.nucleotide,
                allele2: new_allele2,
            };
            (info1,info2)
}

fn main() {
    pretty_env_logger::init();

    // now the next is not really for any argument 
    // parsing but simply to get the command which 
    // was used to execute as I cant get this from clap
    //let args: Vec<String> = env::args().collect();
    //let args_string = args.join(" ");
    let matches = app_from_crate!()
    .about("This tool verifies that pairs of BND events are behaving as expected. 
            It assumes that the direction of a given event is correct and that the orientation
            of the mate-event is not correctly encoded in the primary event.
            This became necessary as the genefusion pipeline at NGSAI suffered from that particular problem.
            ")
    .arg(Arg::with_name("FILE")
            .short("i")
            .long("input")
            .value_name("FILE")
            .help("the vcf file to parse")
            .takes_value(true)
            .required(true))
    .get_matches();

  
    let file         = matches.value_of("FILE").unwrap();
    
    // so now we have a BND tree which pretty much consists only of 2 event types:
    // "G]chr16:18804692]"  = FF --> exists
    // "[chr19:50193168[A"  = RR --> exists
    // "G[chr16:18804692["  = FR --> absent
    // "]chr19:50193168]A"  = RF --> absent
    // we want to fix that and the easiest might be to go over all entries
    // of a vcf file and read always 2 lines which must consist of a pair of entries
    // We then the the direction of the each ID and subplement the mate with the 
    // corresponding primary information
    // Then we modify for both the associated event information in the allele field and write both entries again
    let mut vcf_a      = rust_htslib::bcf::Reader::from_path(file).expect("ERROR: could not open vcf file in fofn !");
    let vcf_header = vcf_a.header().samples().to_vec();
    if vcf_header.len() > 1 {
        panic!("ERROR: Currently only single sample vcf accepted!");
    }

    // we need to initialize an empty record first
    let mut prev_record = vcf_a.empty_record();
    // we go now essentially always over the 2nd entry
    // keeping the previous cached
    let mut vcf_writer = Writer::from_stdout(&Header::from_template(vcf_a.header()), true, Format::Vcf).expect("ERROR: could not write vcf to provided file");
    for (n,entry) in vcf_a.records().enumerate() {
        // Turns out first entry is 0 not 1,
        // so we need to count  
        if let 0=n%2 { 
            prev_record = entry.expect("ERROR: could not read 1st record of vcf file!"); 
            continue
        }else{
            let mut cur_record = entry.expect("ERROR: could not read 1st record of vcf file!");
            let id_1 = &prev_record.id() ;
            let id_2 = &cur_record.id() ;

            let id1  = str::from_utf8(id_1).expect("ERROR: cant extract record id!");
            if id1 == "." {
                panic!("ERROR: all entries must have an annotated entry ID!");
            };
            let id2  = str::from_utf8(id_2).expect("ERROR: cant extract record id!");
            if id2 == "." {
                panic!("ERROR: all entries must have an annotated entry ID!");
            };
            // now we need to verify that we are indeed comparing proper pairs
            // we expect and ID .1 and ID .2 following each other
            // Otherwise we have either a different ordered file
            // which we cant treat here or we have a missing entry
            let mut name1 : Vec<&str> = id1.split('.').collect();
            let mut name2 : Vec<&str> = id2.split('.').collect();
            name1.pop();
            name2.pop();
            if (name1 != name2) || (name1.is_empty()) {
                eprintln!("Non-matching entries {:?} and {:?}",name1,name2);
                panic!("ERROR: entries are either mal-formated or not sorted by matching IDs!");
            }
            // this function now returns our corrected alleles
            let (new_allele1,new_allele2) = correct_alleles(&prev_record,&cur_record);            
            
            // now we change the allele record of our elements
            prev_record.set_alleles(&[new_allele1.allele1.as_bytes(),new_allele1.allele2.as_bytes()]).expect("Failed to set alleles");
            cur_record.set_alleles( &[new_allele2.allele1.as_bytes(),new_allele2.allele2.as_bytes()]).expect("Failed to set alleles");
            
            // just to be sure we translate to the header
            vcf_writer.translate(&mut prev_record);
            // now we can write
            vcf_writer.write(&prev_record).expect("Could not write record!");   
            vcf_writer.translate(&mut cur_record);
            vcf_writer.write(&cur_record).expect("Could not write record!");
        };
    }

}