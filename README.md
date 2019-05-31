# README #

# What is targetsqc #
targetsqc is a quality control tool that alerts about regions of low coverage and missing regions in amplicon-based exome sequencing. It can also filter an annotated VCF file (ANNOVAR format) based on the provided gene list.  
targetsqc is currently in **early alpha version**.  


# How to install targetsqc #
1 . Requisites:  

 * Python 3.5+ (I suggest [Anaconda Python](https://www.continuum.io/downloads))  
 * The [pandas](https://pandas.pydata.org/) library (pandas)  
 
2 . Required data files:  
 
 * Annotated genome data, in gff3 format, of the same version used for read alignment and variant calling. E.g. [GRCh37](ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GRCh37.p13_interim_annotation/interim_GRCh37.p13_top_level_2017-01-13.gff3.gz) (unpack the archive for use)  
 * BED file containing the targets list (AmpliseqExome.<...>.designed.bed). For IonTorrent sequencers, this file can be downloaded from the IonTorrent Server in Configurations > References > Target Regions (\<your server>/configure/references/#target-regions).  
  
3 . Download or clone the source code for targetsqc. The Python module can be ran directly from source.  

# Using targetsqc #
targetsqc can be easily ran within a Galaxy instance.
However, it can also be executed from the command line. The basic usage is:  

~~~~
python targetsqc.py  --bamfile <your_sample_bam_file>  
    --target <ampliseq_targets_bed_file>  
    --annotation <genome_gff3_annotation_file>  
    --wannovar <wANNOVAR_annotated_VCF_file>  
    --genelist <list_of_genes_of_interest>  
~~~~

# Structure of the "genes of interest" file #
Genes of interest should be stored in a text file, one gene name per line.  
Lines can be commented out using the hash symbol (`#`).  
`# List of genes relevant for disease X`  
When the gene name *in the gff3 genome annotation file* and gene name *in the Ampliseq BED file* are different, both gene names must be specified in a single line, separated by a tab.  
`GENOME_NAME <tab> AMPLISEQ_EXOME_NAME`  

# Output files #
targetsqc returns the following files:  

1 . "targetsqc_report_<...>.tsv"  
Tab-separated file containing coverage information.  
The purpose of this file is to allow the user to verify a selection of CDS regions with no amplicons, and possibly failed regions, on [IGV](https://software.broadinstitute.org/software/igv/).  

* The first section informs the parameters used for the analysis.  

* The second section lists the genes of interest, as parsed by the script. Two columns, "Annotation check" and "Exome check" are used to notify about gene names that were not found either in the `gff3` file or the Ampliseq Exome `designed.bed` file.  
  We suggest determining correct gene names by loading the `designed.bed`file on [IGV](https://software.broadinstitute.org/software/igv/) and verifying target names at the locus of interest.  
  _Please note that if a gene is indicated as not found in this list, it is considered missing and will not be processed further!_  
  NEW: Columns will also indicate how much percent of the total CDS region is present in the amplicon panel ("% CDS within the panel") and is covered effectively ("% CDS covered >= {n} reads"). _Please remember that this script is offered as-is, with no guarantee of exactedness!_  
  
* The third section, starting with "These CDS regions had no amplicons:", indicates any region for which there is _no amplicon designed_ in the Exome panel.  
  _Please note that if a region is indicated as not found in this list, it is considered missing and will not be processed further!_  
  
* The fourth and last section, starting with "Some amplicons had low coverage:", lists any region where the parameters "min_coverage", "min_cov_each_strand", and "strand_bias" fall below the quality threshold.  


2 . "targetsqc_output_<...>.csv"  
Comma-separated (wANNOVAR-style) selection of variants located within the genes of interest.  
The variants are selected according to wANNOVAR's "Gene.refgene" column.  

3 . "targetsqc_genes_<...>.bed"  
BED file showing all exonic regions considered for the analysis.  

4 . "targetsqc_failed_<...>.bed"  
BED file showing all regions that failed to pass the quality thresholds.

# LICENSE #
**targetsqc is offered under the GNU General Public License.**  
**Please read it here: https://www.gnu.org/copyleft/gpl.html**  

# DISCLAIMER #
_**targetsqc is tested at the best of our knowledge but is offered as-is, with no guarantee of exactness.**_

