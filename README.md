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
 * BED file containing the targets list (AmpliseqExome.<...>.designed.bed). For IonTorrent sequencers, this file can be downloaded from the IonTorrent Server in Configurations > References > Target Regions (\<your server\>/configure/references/#target-regions).  
  
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



# LICENSE #

**targetsqc is offered under the GNU General Public License.**  
**Please read it here: https://www.gnu.org/copyleft/gpl.html**  

