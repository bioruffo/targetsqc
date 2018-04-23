#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 14:39:58 2017

@author: roberto

"""

# TODO REMEMBER CFTR-37!! The Ion Proton / PGM will trim primer sequences
# for variant analysis, but the BAM still contain these ends!! Don't count them!!
# TIP: just employ TVC post-processed BAM


import pandas as pd
import pysam
import csv
import argparse
from bisect import bisect
from collections import defaultdict
from utils import createdir


def parse_cmd():
    '''
    parsing params from command line
    '''
    parser = argparse.ArgumentParser(description='A coverage analysis tool' +
    'data obtained from amplicon-based exome sequencing')
    parser.add_argument('-v', '--version', action='version',
                        version='1.0.0b1')
    parser.add_argument('-b', '--bamfile', type=str, required=True,
                        help='alignment filename (BAM format)')
    parser.add_argument('-t', '--target', type=str, required=True,
                        help='ampliseq targets filename (BED format)')
    parser.add_argument('-w', '--wannovar', type=str, required=True,
                        help='variant annotation filename (CSV format)')
    parser.add_argument('-a', '--annotation', type=str, required=True,
                        help='genome annotation file (GFF3 format)')
    parser.add_argument('-g', '--genelist', type=str, required=True,
                        help='list of genes of interest (text file).' +
                        ' One gene per line')
    parser.add_argument('-p', '--platform', type=str,
                        nargs='?', default='Proton',
                        choices=set(('Proton', 'PGM')),
                        help='sequencing platform')
    parser.add_argument('-s', '--select', type=bool,
                        nargs='?', default=True,
                        choices=set((True, False)),
                        help='select variants that are within genelist?')
    parser.add_argument('-q', '--quality', type=int,
                        nargs='?', default=1,
                        choices=set((0, 1, 2, 3, 4)))
    parser.add_argument('-m', '--min-mapping-qv', type=int,
                        nargs='?', default=4)
    parser.add_argument('-f', '--min-failed-size', type=int,
                        nargs='?', default=2,
                        help='ignore 1-bp dels')
    parser.add_argument('-u', '--over', type=int,
                        nargs='?', default=10,
                        help='add this many nucleotides on each side of exons' +
                        ' and CDS')
    parser.add_argument('-o', '--outdir', type=str,
                        nargs='?', default='targetsqc_output',
                        help='output dir name')
    parser.add_argument('-n', '--sample-name', type=str,
                        nargs='?', default='sample',
                        help='sample identifier')

    opts = parser.parse_args()
    return opts


class Genedata(object):
    def __init__(self, genesfile, bedfile, annot_genes, over, bedout):
        # Load desired gene names from .txt file
        # self.genes_ex = user-declared gene names to search the exome .bed file
        # self.bridgedict = alternate gene names for the genome annotation file
        self.genes_ex, self.bridgedict = self.load_genes(genesfile)
        self.genes_ex_input = self.genes_ex.copy()
        # self.genes_ann = user-declared gene names to search the genome annotation file
        self.genes_ann = [self.bridge(gene) for gene in self.genes_ex_input]
        print("Loaded {} gene names.".format(len(self.genes_ex)))
        # Loading amplicon data from the .bed file
        self.amplicons, self.notfound_bed = self.load_amplicons(bedfile)
        self.genes_ex = sorted(set(self.genes_ex).difference(self.notfound_bed))
        # Loading transcript data from the genome annotation file
        self.transcripts, self.notfound_transc = self.load_transcripts(annot_genes)
        self.genes_ann = sorted(set(self.genes_ann).difference(self.notfound_transc))
        '''
        self.data[gene] = {'chr': 'chr17', 'CDS': [], 'exon': [],
                           "5'-splice": [], "3'-splice": [], 'missingCDS': []}
        '''
        self.data = dict([(gene, dict()) for gene in self.genes_ex])
        # populate self.data
        for gene in self.genes_ex:
            self.populate_data(gene, over)
        self.write_gene_bed(bedout)

    def load_genes(self, genesfile):
        with open(genesfile, "r") as f:
            rawdata = f.read().split('\n')
        linesnum = 0
        geneslist = []
        bridgedict = {}
        for item in rawdata:
            if item and not item.startswith('#'):
                linesnum += 1
                if '\t' in item:
                    gene_ann, gene_ex = item.split('\t')
                    if gene_ex.strip() == '':
                        gene_ex = gene_ann
                    if gene_ex.strip() not in geneslist:
                        geneslist.append(gene_ex.strip())
                        bridgedict[gene_ex.strip()] = gene_ann.strip()
                else:
                    if item.strip() not in geneslist:
                        geneslist.append(item.strip())
        print("{} valid lines read.".format(linesnum))
        return geneslist, bridgedict

    def bridge(self, gene):
        return self.bridgedict.get(gene, gene)

    def load_amplicons(self, bedfile):
        print("Reading amplicons from bed file...")
        amplicons = pd.read_csv(bedfile, sep='\t', header=None, skiprows=1)

        amplicons['gene'] = amplicons[7].str.extract("GENE_ID=(\w*?);",
                 expand=False)
        amplicons = amplicons[[0, 1, 2, 3, 'gene']]
        amplicons = amplicons[amplicons['gene'].isin(self.genes_ex)]

        notfound_bed = sorted(set(self.genes_ex).difference(
                set(amplicons['gene'])))
        if notfound_bed:
            print("I couldn't locate these {} genes within the exome panel" \
                  "(see report):\n{}".format(len(notfound_bed), ', '.join(notfound_bed)))
        return amplicons, notfound_bed

    def load_transcripts(self, annot_genes):
        print("Reading transcript CDSs from genome annotation...")
        transcripts = pd.read_csv(annot_genes, sep='\t',
                comment="#", header=None, usecols=[0, 2, 3, 4, 8]) # TODO see if this works
        # Using only NC_ contigs
        transcripts = transcripts[transcripts[0].str.startswith("NC_")]
        #  no pseudogenes
        transcripts = transcripts[~transcripts[8].str.contains("pseudo=true")]
        # Using transcript rows, not 'gene'
        transcripts = transcripts[transcripts[2].isin(['exon', 'CDS'])]
        transcripts = transcripts.drop_duplicates()
        transcripts['gene'] = transcripts[8].str.extract("\;gene=(\w*?)\;",
                expand=False)

        transcripts = transcripts[transcripts['gene'].isin(self.genes_ann)]

        notfound_transc = set(self.genes_ann).difference(set(transcripts['gene']))
        if notfound_transc:
            print("I couldn't locate these {} genes within the transcripts file" \
                  "(see report):\n{}".format(len(notfound_transc), ', '.join( \
                   sorted(list(notfound_transc)))))
        return transcripts, notfound_transc

    def get_annot_genenames(self, annot_genes):
        with open(annot_genes, 'r') as f:
            genes = pd.read_csv(f, delimiter='\t', comment='#', header=None, usecols=[0, 2, 8])

        genes = genes[genes[2]=='gene']
        genes = genes[genes[0].str.startswith("NC_000")]
        genes = genes[~genes[8].str.contains("pseudo=true")]
        genes['name'] = genes[8].str.extract(';Name=(.*?);', expand=False)
        genes['synonyms'] = genes[8].str.extract(';gene_synonym=(.*?);', expand=False)
        # some genes are still duplicated: tRNA, and a few genes from chromosomes x and y
        # genes[genes.duplicated('name',keep=False)]

        # Take all genes with synonyms (they're very few)
        syns = genes[~genes['synonyms'].isnull()][['name', 'synonyms']].to_dict(orient='split')['data']
        names = [s[0] for s in syns]
        synonyms = defaultdict(list)
        for name, entry in syns:
            for synonym in entry.split(','):
                if name.upper() not in synonyms[synonym.upper()]:
                    synonyms[synonym.upper()].append(name.upper())
        return names, synonyms

    def populate_data(self, gene, over):
        # print("Populating data for gene: {}".format(gene))
        self.data[gene]['chr'] = self.amplicons[self.amplicons['gene'] == gene].iloc[0,0]
        gene_transcripts = self.transcripts[self.transcripts['gene'] == self.bridge(gene)]
        for item in ['CDS', 'exon']:
            data = gene_transcripts[gene_transcripts[2] == item]
            # lowering the start nucleotide by 1 to include the first nucleotide
            sel = nttuple(data, [3,4], minusone=True)
            boundaries = set_boundaries(sel, 0)
            self.data[gene][item] = boundaries
        # `data` and `sel` now contain exon information
        for i, item in enumerate(["5'-splice", "3'-splice"]):
            splice = [sorted([s[i], s[i]+over*([-1, 1][i])]) for s in sel]
            self.data[gene][item] = set_boundaries(splice, 0)
        self.data[gene]['missingCDS'] = compare_tuplelists(self.data[gene]['CDS'],
                 set_boundaries(nttuple(self.amplicons[self.amplicons['gene'] == gene], [1, 2]),
                                over=0))

    def write_gene_bed(self, bedout):
        print("Writing gene data to", bedout)
        with open(bedout, "w") as f:
            for gene in self.data:
                for category in ['CDS', 'exon', "5'-splice", "3'-splice"]:
                    for item in self.data[gene][category]:
                        f.write('{}\t{}\t{}\t{}\t{}\n'.format(self.data[gene]['chr'],
                                                       item[0],
                                                       item[1],
                                                       category,
                                                       gene))


class NGSData(object):
    def __init__(self, wannovarfile, bamfile, genes):
        self.genes = genes
        self.wannovar = wannovarfile
        print("Reading coverage data from:", str(bamfile))
        self.bam = pysam.AlignmentFile(bamfile, "rb")


    def select(self, output = "selector_output.csv"):
        global out
        out = self
        print("Loading data from:", str(self.wannovar))
        self.data = pd.read_csv(self.wannovar)
        # Compatibility fix to accept both old wANNOVAR column label "Gene.refgene"
        # and the new label "Gene.RefGene"
        try:
            genenamecol = [colname for colname in self.data.columns.values \
                           if colname.lower() == 'gene.refgene'][0]
        except IndexError:
            print("### WARNING ### Could not identify gene names column in wANNOVAR.")
            raise
        print("Performing selection of variants in the desired genes...")
        self.selected = self.data[[any(x in self.genes.genes_ex for x in \
                                   genes.split(',')) if type(genes) == str else False \
                                   for genes in self.data[genenamecol]]]
        print("Saving variants from selected genes to:", output)
        self.selected.to_csv(output, index=False)


    def analyze_coverage(self, gene, run_params):
        min_mapping_qv, min_coverage, min_cov_each_strand, strand_bias, \
                min_failed_size, over = run_params
        flattened = dict()
        for item in ['CDS', 'exon', "5'-splice", "3'-splice"]:
            flattened[item] = flatten(self.genes.data[gene][item])
        covlist = []
        failedlist = []
        for __, amplicon in self.genes.amplicons[self.genes.amplicons['gene'] == gene].iterrows():
            failedchrom = None
            failedreasons = set()
            failedwhere = set()
            chrom, start, stop = amplicon[:3]
            samstart = start - 1
            samstop = stop + 1
            pileiter = self.bam.pileup(chrom, samstart, samstop)
            for column in pileiter:
                location = column.pos + 1
                if location >= start and location <= stop:
                    r_directs = [(not read.alignment.is_reverse) for read in column.pileups
                                 if read.query_position is not None
                                 and read.alignment.mapq >= min_mapping_qv]
                    fcount = r_directs.count(True)
                    rcount = r_directs.count(False)
                    rowdata = [chrom, location, fcount, rcount, amplicon['gene']]
                    reason = []
                    if (fcount + rcount) < min_coverage:
                        reason.append('min_coverage')
                    if min(fcount, rcount) < min_cov_each_strand:
                        reason.append('min_cov_each_strand')
                    if ((fcount + rcount) > 0) and ((max(fcount, rcount) / (fcount + rcount)) > strand_bias):
                        reason.append('strand_bias')
                    where = self.get_where(location, flattened)
                    rowdata.extend([(reason == []), reason, where])
                    if reason:
                        failedreasons.update(reason)
                        failedwhere.update(where)
                        if not failedchrom:
                            failedchrom = chrom
                            failedlist.append([amplicon['gene'], failedchrom, location])
                        else:
                            pass # we're inside a failed region
                    else:
                        if failedchrom:
                            failedlist[-1].extend([location-1, (location - failedlist[-1][-1]),
                                                  sorted(failedreasons), sorted(failedwhere)])
                            if failedlist[-1][4] < min_failed_size:
                                failedlist.pop(-1)

                        failedchrom = None
                        failedreasons = set()
                        failedwhere = set()
                    covlist.append(rowdata)
            # closing last pileup region, if open
            if failedchrom:
                failedlist[-1].extend([min(location, stop),
                                       (min(location, stop) - failedlist[-1][-1] + 1),
                                       sorted(failedreasons), sorted(failedwhere)])
                if failedlist[-1][4] < min_failed_size:
                    failedlist.pop(-1)
        return covlist, failedlist

    def get_where(self, location, flattened):
        where = []
        for item in ["CDS", "exon"]:
            if bisect(flattened[item], location) % 2:
                where.append(item)
                break
        for item in ["5'-splice", "3'-splice"]:
            if bisect(flattened[item], location) % 2:
                where.append(item)
        if where == []:
                where.append("intron")
        return where


def report(bam, run_params, genes, data,
           reportfile = "selector_report.tsv",
           failedbed = "selector_failed.bed"):
    reportfile = reportfile
    min_mapping_qv, min_coverage, min_cov_each_strand, \
                  strand_bias, min_failed_size, over = run_params
    missing = []
    failedlist = []
    for gene in genes.genes_ex:
        if gene not in genes.notfound_bed:
            cov, fail = data.analyze_coverage(gene, run_params)
            failedlist.extend(fail)
        if genes.data[gene].get('missingCDS', None):
            for item in genes.data[gene]['missingCDS']:
                missing.append([gene, genes.data[gene]['chr'], str(item[0]), str(item[1]), str(item[1]-item[0]+1)])
            # TODO "coverage" and "failed" can be used in a per-nt basis
            #covlist.extend(cov)

    #coverage = pd.DataFrame(covlist)
    #failed = coverage[coverage[5] == False]
    if failedlist:
        print("Some genes had areas of low coverage (see report).")

    if missing:
        print("Some genes had CDS regions without a designed amplicon (see report).")

    print("Writing report...{}".format(["", " THERE ARE WARNINGS!"][genes.notfound_bed != [] \
                                                                    or failedlist != [] \
                                                                    or missing != []]))
    with open(reportfile, 'w') as r:
        with open(failedbed, 'w') as f:
            r.write("Analyzed file:\t{}\n\n".format(bam))

            r.write("Parameters used:\n\t{}\n\n".format("\n\t".join(["Minimum coverage\t" + str(min_coverage),
                        "Minimum cov. each strand\t" + str(min_cov_each_strand),
                        "Maximum strand bias\t" + str(strand_bias),
                        "Minimum mapping QV\t"+str(min_mapping_qv),
                        "Minimum failed region size\t"+str(min_failed_size)])))

            r.write("Genes of interest:\nExome name\tAnnotation name\tExome check\tAnnotation check\n")
            for gene in genes.genes_ex_input:
                r.write('\t'.join([gene,
                                   genes.bridge(gene),
                                   ["", "Exome: not found"][gene in genes.notfound_bed],
                                   ["", "Annotation: not found", ""][gene in genes.notfound_transc]]) \
                        + '\n')

            if missing:
                r.write("WARNING:\nThese CDS regions had no amplicons:\n")
                r.write("\tGene\tChrom\tStart\tEnd\tExtension\n")
                for miss in missing:
                    r.write('\t' + '\t'.join(miss) + '\n')
                    n, c, s, e, *_ = miss
                    f.write("{}\t{}\t{}\tmissing\t{}\n".format(c, s, e, n))
                r.write('\n')
            else:
                r.write("All CDS were represented by amplicons.\n\n")

            if failedlist:
                r.write("WARNING:\nSome amplicons had low coverage:\n\n")
                r.write("\tGene\tChrom\tStart\tEnd\tExtension\tReason\tLocation\n")
                for item in failedlist:
                    n, c, s, e, ex, re, *_,  = item
                    r.write(((7*'\t{}')+'\n').format(*item[:-2], ', '.join(item[-2]), ', '.join(item[-1])))
                    f.write("{}\t{}\t{}\t{}\t{}\n".format(c, s, e, '-'.join(re), n))

            else:
                r.write("All positions had sufficient coverage.")


def set_boundaries(tpls, over):
    boundaries = dict()
    for start, end in tpls:
        start = max(0, start - over)
        end = end + over
        fhit = start
        rhit = end
        for bstart, bend in list(boundaries.items()):
            if bstart <= start and bend >= start:
                fhit = bstart
                boundaries.pop(bstart)
            if bstart <= end and bend >= end:
                rhit = bend
                if bstart in boundaries:
                    boundaries.pop(bstart)
        boundaries[fhit] = rhit
    return sorted(boundaries.items())


def nttuple(dataframe, columns, minusone = False):
    col1, col2 = columns
    return [(int(tpl[col1] - int(minusone)), int(tpl[col2])) \
            for idx, tpl in dataframe[columns].iterrows()]


def flatten(tuplelist):
    return [item for tpl in tuplelist for item in tpl]


def compare_tuplelists(query, reference):
    '''
    return, as a list of tuples, all regions of `query` that are not in `reference`.
    `reference` should be an amplicon block (not list of amplicons, but processed blocks).
    '''
    flatref = flatten(reference)
    external = []
    for item in query:
        startbis, endbis = map(lambda x: bisect(flatref, x), item)
        if startbis % 2:
            # start is within one block
            if startbis == endbis:
                # start and end are within the same amplicon block
                continue
            else:
                # there is at least one gap
                blocks = flatref[startbis:endbis]
                newext = [(blocks[i], blocks[i+1]) for i in range(0, 2 * (len(blocks) // 2), 2)]
                if len(blocks) % 2:
                    newext.append((blocks[-1], item[1]))
                external.extend(newext)
        else:
            # start is outside a block
            blocks = flatref[startbis:endbis]
            if not blocks:
                external.append(tuple(item))
                continue
            newext = [(item[0], blocks.pop(0))]
            if blocks and len(blocks) > 1:  # Sometimes it's uncovered-covered-uncovered
                newext.extend([(blocks[i], blocks[i+1]) for i in range(2 * (len(blocks) // 2), 2)])
            if len(blocks) % 2:
                newext.append((blocks[-1], item[1]))
            external.extend(newext)
    return external


def main():

    opts = parse_cmd()
    # hotspot, indel, mnp, snp, custom = 0, 1, 2, 3, 4
    suffix = opts.sample_name
    outdir = opts.outdir
    platform = opts.platform
    filter_quality = opts.quality
    min_mapping_qv = opts.min_mapping_qv
    min_failed_size = opts.min_failed_size # ignore 1-bp dels
    over = opts.over # add this many nucleotides on each side of exons and CDS
    bam = opts.bamfile
    bed = opts.target
    wannovar = opts.wannovar
    annotation = opts.annotation
    genes = opts.genelist
    select = opts.select

    createdir(outdir)
    # Output file names
    reportfile = "{}/selector_report_{}.tsv".format(outdir, suffix)
    selectfile = "{}/selector_output_{}.csv".format(outdir, suffix)
    genesbed = "{}/selector_genes_{}.bed".format(outdir, suffix)
    failedbed = "{}/selector_failed_{}.bed".format(outdir, suffix)

    custom_min_coverage = 6
    custom_min_cov_each_strand = 3
    custom_strand_bias = 0.95

    if platform == "Proton":
        # Yeah against PEP8, shoot me
        min_coverage =        [6,   10,   5,    5,    custom_min_coverage][filter_quality]
        min_cov_each_strand = [3,    5,   0,    0,    custom_min_cov_each_strand][filter_quality]
        strand_bias =         [0.95, 0.9, 0.98, 0.98, custom_strand_bias][filter_quality]

    elif platform == "PGM":
        # Yeah against PEP8, shoot me
        min_coverage =        [6,   15,    6,    6,    custom_min_coverage][filter_quality]
        min_cov_each_strand = [3,    5,    0,    0,    custom_min_cov_each_strand][filter_quality]
        strand_bias =         [0.98, 0.85, 0.95, 0.95, custom_strand_bias][filter_quality]
    else:
        print("Unknown platform", platform)
        return 1

    run_params = [min_mapping_qv, min_coverage, min_cov_each_strand,
                  strand_bias, min_failed_size, over]

    if not 'oldgenes' in globals():
        genes = Genedata(genes, bed, annotation, over, genesbed)
    else:
        genes = oldgenes
        print("Using previously calculated `genes` data.")

    data = NGSData(wannovar, bam, genes)

    report(bam, run_params, genes, data, reportfile, failedbed) # TODO add "over"
    if select:
        data.select(output=selectfile)

    print("End.")
    return genes


if __name__ == '__main__':
    oldgenes = main()
