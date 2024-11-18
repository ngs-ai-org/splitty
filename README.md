# Splitty - tool suite for fusion detection

This respository is targeted around the detection of fusion events.
It contains algorithms for fusion detection as well as tools to cluster and filter the VCF-based events in a BND-pair fashion.

Splitty will derive the position of integration/fusion based on alignment information, in particular the information of supplementary alignments of reads which are split. Therefore it is important to use callers which support such a reporting like BWA-MEM2 for short reads or minimap2 for long reads.

Alignments are not solved and can create ambiguity, especially around break-end positions.

This is amplified and reflected in

- alignments adjacent to edges, especially for spliced elements
- scoring, where multiple targets result in identical scores
- situations where mapping of supplementary element is not resulting in a high enough score 

It is subsequently reflected in the fuzzy-ness of the exact reported position and splitty provides infomration (CIPOS) in bp reporting the expected fuziness.

Occasionally (and in particular for integration analysis) it can happen that a predicated site is due to these ambiguities beyond the boundaries of a chromosome/contig. Meaning it is < 0 or > contigSize. If that happens, the program will assign the last valid position and report the difference to the predicted one as fuzzyness, too.

Both, fragment and paired read fusion detection uses annotation to provide additional information of the involved genes (transcripts).
Importantly, if a annotation file is present the notion for known fusions is `GeneA--GeneB` with `NA` for unknown elements e.g. `GeneA--NA`.
If no annotation is present, this though replaced by the chromosome names, e.g. `ChrA--ChrB` or `ChrA--ChrA` and never defaulting to `NA`.

Splitty uses the annotation of the underlying transcripts to determine which transcript is the most likely to be source of the resulting fragment. If multiple transcripts are overlapping, priority is given to the one in the same sense of the transcript and then the larger proportion of the transcript being represented in the alignment block. If this results in identical "scores", splitty will report both (or multiple) transcripts as possible annotations of the same VCF event.

On the other hand, if one of the targets results in multi-mappers and cant be resolved properly then splitty will report one VCF entry (pair) per target as they represent different genomic coordinates and we emphasize coordinates rather than gene-annotations as fusion reporting.
E.g. we might encounter an entry `SA:Z:chr14,42077618,+,3121M279334D415S,4,1486; chr14,35031866,-,1809M65484D1727S,60,0;`
Which described more than one entry and will result potentially into 2 subsequently proposed gene fusions.
Normally the entry with the highest score is kept but sometimes we do encounter identical scores, in this case both are propagated if specified.


## Splitty fragment 

It's main function is to predict and annotate fusion positions based on a contig alignment against a reference genome. There is not necessarily a limit to the length of an aligned contig but it is unlikely to work well on any given short reads.


```bash
splitty-fragment 
classifying cDNA/DNA mappings into gene-fusion events, not for paired-end reads

USAGE:
    splitty fragment [FLAGS] [OPTIONS] --bam <FILE>

FLAGS:
    -k, --keep          keep all suggested entries if multiple SAs available, otherwise only most precisecareful though
                        this might increase the number of FPs and precision might become wrong for these entriesin
                        particular for entries which are |----A-->|---B-->|---C--> architecture
    -u, --unique-off    accepts reads with non-unique mapping mapq=0
    -s, --stranded      implies that used cDNA is stranded this provides better annotation precision
    -h, --help          Prints help information
    -V, --version       Prints version information

OPTIONS:
    -b, --bam <FILE>              long read/cDNA alignment vs reference
    -l, --blacklist <FILE>        black-list of queries to ignore , one ID per line
    -w, --distance <int>          minimum allowed distance for intra-chromosomal fusions [default=10000]
    -g, --gtf <FILE>              Provide annotation of the genome. Use option m to specify selected feature The field
                                  "gene_id" or, if available "gene_name" will be used for reporting Comply either to
                                  CHESS or Gencode standard
    -m, --match <feature-name>    specifies the feature field in a GTF to be used, [defaults: "transcript" ] 
                                  activates automatically option "gtf" too
    -d, --ratio <float>           ratio of accepted 5' to 3' clipping (and other way around), 0.0 will deactivate check
                                  [default: 0.1]
    -r, --reference <int>         reference in fasta format provided, required for vcf output
    -i, --similarity <float>      minimum seq similarity of matches [default:98.0]
    -t, --threads <int>           number of threads for reading BAM [default: 1]
    -f, --vcf <int>               generate vcf encoded fusion call, requires reference

```

If there is actually a paired-read data detected the program will exit as this is not compatible. It has been tested to work on assembled transcript fragments, entire transcripts, assembled integration sites and IsoSeq reads.

`splitty fragment` scans the BAM file and extracts entries which have a `SA` field of supplementary alignment annotated and evaluates whether it classifies as a fusion indication. 
It checks that the minimum sequence similarity (gap-compressed) and that information are coherent - otherwise it will trash said entry. 
In a second step called "pruning", it expects for each contig >= 2 alignments to be found with an `SA` entry and that these are complementary. 
Based on these complementary primary alignments, it then generates a precise estimation of the integration/fusion site. It can as well estimate it from a single non-complementary entry, but then precision is sub-par as SA entries do not contain alignments but only mapping results.

If both primary alignments disagree of fusion position on the query contig we will report this as a `precision` value in the end.
It takes the difference in bp between both suggestions and divides it by 2 providing the fuzzyness in bp on both sides.
If that precision score is not equal to 0 one should consider the following positions as imprecise:

- fusion point on query position
- fusion point on genome from geneA
- fusion point on genome from geneB

At the end splitty generates a BND style VCF with paired entries.
The header of these files contain supplementary information about the version of the program, the date and the used command.

**Note**: currently single entries are permitted and annotated but the information not returned at the very end of the program.
This might have to be changed in the future.

**Note**: As briefly mentioned above, the SA-alignments do generate as well cigars which are being used to estimate the position of the gene fusion for each alignment individually.
These cigars are though often flawed and approximations and do no contain full alignment information.
This is supposedly a feature not a [bug](https://github.com/lh3/minimap2/issues/524).

In some cases, the alignment mayb contain towards both extremities clipping which makes the determination of the sense of the fusion less intuitive. Splitty has 2 ways to deal with this:

1. use a ratio of maximum accepted clipping between both sides - default 10%
2. deactivate and use comparison with 2nd aligned elements

The first would e.g. discard an element for which one side contained 100bp of clipping and the other side >10bp of clipping.
This can be useful to reduce the number of FP elements, especially for very similar reference-to-sequencing situations.
The second case ignores though such a factor, and rather would compare than with the SA alignment the similarity. It generates an absolute delta between clipping and 2nd alignment, e.g. 2nd alignment would be 110 bp, than for the 100bp clipping the delta would be 10bp and for a e.g. 15bp clipping on the second clipping it would result in 85bp. Based on that the 1st would be chosen as fusion direction.
This option might be generally preferable, but has not been extensively tested and compared.

## Paired

This tool is in principle the extension of splitty, meaning it does as well analyze gene fusions from alignments.
The underlying assumptions and therefore result-organization is a bit different though.

```bash
splitty-paired 
detecting integrations and genefusions from mapped paired-reads

USAGE:
    splitty paired [FLAGS] [OPTIONS] --bam <FILE>

FLAGS:
    -u, --unique-off    accepts reads with non-unique mapping mapq=0
    -s, --stranded      implies that used cDNA is stranded this provides better annotation precision
    -h, --help          Prints help information
    -V, --version       Prints version information

OPTIONS:
    -b, --bam <FILE>              long read/cDNA alignment vs reference
    -w, --distance <int>          minimum allowed distance for intra-chromosomal fusions [default=10000]
    -g, --gtf <FILE>              Provide annotation of the genome. Use option m to specify selected feature The field
                                  "gene_id" or, if available "gene_name" will be used for reporting Comply either to
                                  CHESS or Gencode standard
    -m, --match <feature-name>    specifies the feature field in a GTF to be used, [defaults: "transcript" ] 
                                  activates automatically option "gtf" too
    -p, --precision <int>         the required precision in bp for split-read derived fusion/integration points Note:
                                  increasing might help to pool events in close proximity. This might fuse though
                                  independent events and the lowest possible value might be advisable [default: 5]
    -z, --range <int>             the allowed range in bp within which a split pair is considered support for a fusion
                                  point. Applies for both sides of the fusion [default: 500]
    -r, --reference <int>         reference in fasta format provided, required for vcf output
    -i, --similarity <float>      minimum seq similarity of matches [default:98.0]
    -t, --threads <int>           number of threads for reading BAM [default: 1]
    -f, --vcf <int>               generate vcf encoded fusion call, requires reference

```


It does not derive from splice-aware alignment of full- or fragmented-cDNA but from splice-aware alignment of paired-short reads.
This can be either WGS short-reads for integration or RNAseq reads for gene fusions. The latter must be aligned with a splice-aware aligner such as STAR. Each read is analyzed whether it contains indication of a fusion + each read pair is as well evaluated whether the pairs do represent information. It is though not possible to obtain fusion/integration information without split-reads as these are the base in the 1st step of the program.
The second step then consists into re-analyzing the entire BAM for discovered points to evaluate whether supporting information can be found.

Precision is the standard deviation from multiple supports compared to the median value (in both directions).
If not enough support is gathered (because only 1 read found), then it becomes NA instead.

Opposed to long read sequencing, RNAseq is most of the time not stranded and therefore the option for strandness should be ignored.
Multi-threading is limited to multiple threads for reading and alignment should be reduced before to informative alignment-only.

## Utilities

This repo contains furthermore tools which help to deal properly with BND-pair VCF entries.


### id_fusion_reads

This tool allows to select from an alignment all reads which indicate fusion of transcripts.
It is designed to work on paired-read alignment on a collection of transcripts.
In principle it does for each alignment 2 operations:

- verifies if the 2nd in pair aligns to a different transcript
- verifies if the read itself has a supplementary alignment (SA) tag and evaluates it

For the 2nd it similarly evaluates if an SA field then points to a different transcript or is
within the same transcript present.

As alignments on transcripts suffer from potential discordant alignment between isoforms we need
to particularly take case that different targets are really different genes as well.

Therefore one needs to specify the transcriptome source to faciliate this.
For CHESS this is implemented and rather easy as the CHESS ID and transcript ID are easy to collapse
onto a gene-level annotation. For ENSEMBLE/GENCODE this is more complicated as the identifiers do not easily allow this without more knowledge.
Meaning we will need a GTF describing for each transcript the gene-level ID as well.
Currently, GENCODE/ENSEMBLE is not yet implemented.

```bash
USAGE:
    splitty id_fusion_reads [OPTIONS] --bam <FILE> --source <STRING>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -b, --bam <FILE>            long read/cDNA alignment vs reference
    -r, --source <STRING>       define the source of your transcriptome data. This is key to identify if 2 FASTA entries
                                are different transcript but from same gene or if they are indeed 2 different genes and
                                are indicative of a fusion [default: CHESS]  [possible values: CHESS]
    -i, --similarity <float>    minimum seq similarity of matches [default:98.0]
    -t, --threads <int>         number of threads for multi-threading [default: 10]

```

Multi-threading is fully implemented here and will benefit massively from as many CPUs provided as possibe.
Eventually, BAM reading IO-capabilities become the bottleneck in our experience.


### vcf_bnd_merge

This tools is really designed to merge paired read and CDNA-fragment derived information.
It expects 2 VCF files and will only work IMHO on splitty derived output files.
This tool combines 2 vcf files describing BND. SVTYPE needs to be BND and ID + MATEID has to be set. Important: if you
have Fragment and read derived results, please pass them in this order ignoring the "--simple" and "--flag". To simply
check if entries in FileA were supported in FileB use the "--simple" flag. In this case the "--flag" field will be added
as INFO field and be set to 0 if no support was found and 1 if support was detected. One can provide further INFO fields
which will be propagated accordingly in a comma separated list. Currently this requires though elements which have a
single entrie. E.g. CIPOS which provides 2 numbers typically are not supported.

Scoring: if only entryA/B but not both are present, low confidence events will get the label
"SPLITTY_FAILED" . If both evidence are combined the filter label is set to "PASS".
Otherwise no filter field is applied.
IMPORTANT: currently the strand is not verified of the events!

```bash
USAGE:
    vcf_bnd_merge [FLAGS] [OPTIONS] --FileA <vcf> --FileB <vcf> --out <vcf>

FLAGS:
    -s, --simple     simple true/false comparison if fields in A are supported in B one can set a flag name which is
                     then added as TRUE if existent in B
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -a, --FileA <vcf>             the first file used to compare BND events
    -b, --FileB <vcf>             the 2nd file used to compare BND events
    -c, --cipos <int>             the distance within which events are combined [default: 500]
    -f, --flag <string>           name of the flag value in case of "simple"
    -d, --info-fields <string>    a comma separated list of INFO fields to propagate, fomrat of field must be Integer
                                  [default: SR]
    -o, --out <vcf>               the combined out vcf
    -t, --threads <int>           vcf reading threads (files need to be huge to make an impact) [default: 2]

```

If these 2 files are provided it follows the following weighting schema:

seen in both files:

- use fragment based for output VCF
- extend that entry with information based on read-info
- weight: **HIGH**

seen in fragments only:

- use framgent based for output VCF
- use fragment based weighting
  - SPLITTY_FAILED in fragment based; weight **LOW**
  - otherwise; weight **MED**

seen in reads only:

- use read based for VCF
- use read-based for weighting:
  - SPLITTY_FAILED in read based; weight **LOW**
  - otherwise; weight **MED**

### vcf_sv_cluster

This tool compares and extends BND based events for single or multiple samples. It expects a list of vcf files as fofn
or a single vcf which is parsed and organized in an inverval tree. 
SVTYPE can be BND, INS and DEL. For BND ID + MATEID has to be set. 
Events of same type are clustered together based on a defined maximal distance and potentially thereby extend the
original event range. 
The program produces a 0-based BEDPE file with the ranges of the clusters and members that have been observed within.

By adding the `purge` option, it will only keep events which have more than 1 sample annotated and only 1 event pers sample.
This helps to reduce the size of the final file if many events cluster together.
The `multi-only` option will only output in the BEDPE file events which had >1 event in a given cluster.
Otherwise each VCF BND entry will potentially generate one cluster entry.


```bash


USAGE:
    vcf_sv_cluster [FLAGS] [OPTIONS] --chroms <FILE> --fofn <FILE> --vcf <FILE>

FLAGS:
    -m, --multi-only    will only report multi event regions, no singletons
    -p, --purge         will only report one event per sample. Can be used with single files to remove redundant events,
                        too
    -u, --uni-dir       uni-directional, ignoring direction of events
    -h, --help          Prints help information
    -V, --version       Prints version information

OPTIONS:
    -c, --chroms <FILE>    tab-separated file with name in 1st column and chrom length in second column
    -f, --fofn <FILE>      a file with files to parse, one entry per line
    -r, --range <int>      the distance within which events are combined [default: 500]
    -t, --threads <int>    vcf reading threads (files need to be huge to make an impact) [default: 2]
    -v, --vcf <FILE>       a single vcf file to collapse and cluster entries

```

Combinations:

 - multiple samples = report unique names of samples in description and all the events for each sample;  score = number of samples with >= 1 observed event
 - multiple samples + multi-only = report unique names of samples in description and all the events for each sample;  score = number of samples with **> 1** observed event
 - multiple samples + purge = report unique names of samples in description **only 1 representative event** per sample ; score = number of samples with >= 1 observed event
 - multiple samples + purge + multi-only =  report unique names of samples in description **only 1 representative event** per sample ; score = number of samples with **> 1** observed event
 - single sample + multi-only = report unique name of sample in description and all the events for the sample ; score = number of events/reads supporting the sample
 - single sample + purge = report unique name of sample in description **only 1 representative event** for the sample ; score = 1
 - single sample + purge + multi-only = always empty result

### vcf_bnd_2fasta

This tool takes BND events and extracts annotated FASTA files with the in silico combined event.It expects pairs of BND
events and needs the matching reference indexed fasta sequence. It cant currently work with events other than BND and
will ignore these. By default the header of a FASTA entry will contain the name of the BND events. The comment field in
the header contains additionally information about the position of the fusion point (FP) and the range chosen for the in
silico edges.The sequence is split where the first half is upper case and the second part of the sequence is lower-case.
This faciliates the identification of the event point Entries need to be sorted by the identifier of the BND events and
pairs must follow each other, e.g.:

**Important limitation**:
This is an *in silico* reconstruction and cant correctly rebuild all genomic features.
It will fail for the following:

- SNVs in the adjacent proximity
- SVs in the adjacent proximity
- splicing in the adjacent proximity

For the above there is **no** possibility to reconstruct this from a VCF file or even to start understanding if reconstruction is correct or wrong.

```bash
USAGE:
    vcf_bnd_2fasta [OPTIONS] --out <file> --reference <file> --vcf <file>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -o, --out <file>          the in silico sequence of the combined BND events
    -p, --prefix <string>     a prefix for the fasta header, can be e.g. sample name - will default to file name without
                              prefix
    -e, --extent <int>        the range for both BND position till which the in silico sequence will extend [default:
                              500]
    -r, --reference <file>    the already indexed FASTA sequence of the VCF matching reference
    -i, --vcf <file>          VCF file with BND events

```

### vcf_sv_subset

This tool allows to subset a VCF file based on ranges in a BEDPE file.
It's primary application is paired BND events for which the SVTYPE field must be properly labeled and an ID + MATEID
must be set.These pairs must have a '.1' and '.2' ending pair BND ID. Single BND events and INDELS are not yet fully
implemented.Resulting VCF files will be sorted by BND PAIR entry IDs, not genomic positions.

```bash
USAGE:
    vcf_sv_subset [FLAGS] [OPTIONS] --bedpe <FILE> --chroms <FILE> --vcf <FILE>

FLAGS:
    -i, --inverse    will do inverse selection, reporting all events except regions in BEDPE file
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -b, --bedpe <FILE>     BEDPE format file indicating 2 regions per line
    -c, --chroms <FILE>    tab-separated file with name in 1st column and chrom length in second column
    -o, --out <FILE>       outfile name for the VCF
    -s, --score <INT>      if set, allows to filter BEDPE regions based on scores, will keep only entries higher than
                           this (e.g. number of observations)
    -a, --vcf <FILE>       the VCF file to subset
```

### vcf_chromremoval 

This tools will analyze a VCF and only keep events provided in a list of chromosomes.
Importantly this will work both for all type of events (e.g. INS,DEL,SNV,DUPL,..) but as well paired BND events.
For the latter it will remove both pairs if one of them is on the black list. The chromosome section of the VCF header
will **not be altered**.
    
```bash
USAGE:
    vcf_chromremoval [FLAGS] [OPTIONS] --chroms <FILE> --vcf <FILE>

FLAGS:
    -i, --inverse    will do inverse selection, removing all entries matching the provided list 
    -s, --simple     will not do a 2-pass but filter directly. Allows missing VCF-ID. Will fail on paired BND events
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -c, --chroms <FILE>    list of chromosomes to keep, one entry per line
    -o, --out <FILE>       outfile name for the VCF
    -a, --vcf <FILE>       the VCF file to subset
```

## Dependencies

Splitty needs rust for the compilation and furthermore `cmake` and a recent `libclang-dev` environment.
