
//! # Gene-fusion and integration libraries
//! 
//! Author: Emanuel Schmid-Siegert
//! 
//! This libraries are a collection of functions and structures which
//! help working with integration and fusion events.
//! They are in principle used in the fusionR suite which thrives to identify
//! events of integration or fusion from either paired-reads or mapped cDNA elements.
//! 
//! Unfortunately, it turned out that passing records around using rust-htslib is buggy
//! at the time this program was written and for the sub-program `id_fusion_reads` this resulted
//! either in 
//! - very slow processing of reads
//! - irregular segmentation faults
//! - blowing up RAM usage into >100 GB
//! 
//! Therefore I started to use noodles in parallel to do multi-threaded operations.
//! Unfortunately it has a very different data-structure to hts-lib and therefore
//! needed many similar functions with different input and intermediate results.
//! Subsequently I split here the libraries into :
//!  - common: functions + structures used for noodles and htslib
//!  - hts_lib_based: functions specific for htslib derived input
//!  - noodles_based: functions specific for noodles derived input
//!

/// functions + structures used for noodles and htslib
pub mod lib {
    pub mod common;
    /// functions specific for htslib derived input
    pub mod hts_lib_based;
    /// functions specific for noodles derived input
    pub mod noodles_based;
}