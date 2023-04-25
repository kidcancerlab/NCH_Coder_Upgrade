#!/gpfs0/export/opt/anaconda-2.3.0/bin/python
# 
# Author  : Ali Snedden
# Date    : 5/27/16
# License : MIT
# 
# PURPOSE :
#   Read in a GTF file and create a set of simulated fastq data that 
#   can be aligned by tophat to GRCH38. The purpose is to be able to 
#   control directly the relative ratio of different transcripts 
#   (ie isoforms) expressed. This will give us a way to test 
#   MISO, rMATS and other such programs.
#   
# VERSIONS:
#
#
#
#
# TESTING:
#   1. GTF is read in correctly and parsed by GTF_ENTRY.__init__
#      -> Tested by diffing output w/ original gtf file
#   2. Scales like crap ...O(N**2). This needs sped up b/c it 
#      cannot complete the whole genome and gtf in 12+ hours, which is 
#      unacceptable. Get's bogged down in finding exons for transcripts.
#      --> in link_exons_trans_and_genes() I fixed this issue.
#          After reading in exons, transcripts and genes, it does one final
#          pass through the list of GTF_ENTRYs and _assuming_ that 
#          all the transcripts are directly _after_ the gene but _before_
#          the next gene. Same goes for exons between transcripts.
#
# FUTURE: 
#   1. If I was a clever object oriented programmer, I'd appropriately use
#      inheritence in the parts of the classes, GENE, EXON and TRANSCRIPT
#
#   2. Improve algorithm so that I don't need to loop through Gtf 3 times
#      to read in all the exons, transcripts and genes.
#
#   3. Handle CDS, stop_codon, etc. and other features.  Currently I only
#      handle genes, transcripts and exons.
#
#   4. Add strandedness options. fr-firststrand, fr-secondstrand, fr-unstraned
#      See https://dbrg77.wordpress.com/2015/03/20/library-type-option-in-the-tuxedo-suite/
#      for discussion.
#      Use IGV for seeing if I got it correct.
#
#   5. Paired end reads.
#
#   6. Fix bug, for instance with 100 nt reads starting at 930330 for trans
#      ENST00000616125, it spans 3 exons, but I only report 2
#
import sys
import re
import time
import random
import copy
import operator
import numpy
from functions import read_gtf
from functions import get_exon_list
from functions import get_transcript_list
from functions import get_gene_list
from functions import read_genome
from functions import get_exon_seq
from functions import get_trans_seq
from functions import reverse_complement
from functions import get_list_of_unique_gtf_features
from functions import link_exons_trans_and_genes
from functions import create_gene_and_trans_lookup_dict
from functions import print_gtf_statistics
from functions import find_trans_that_differ_by_1_exon
from functions import read_config
from functions import print_transcripts_with_seqs
from functions import create_insert
from functions import create_fastq_file
from error import exit_with_error

def print_help(Exit=1):
    """This prints the Help"""
    print("\nUSAGE: ./simulate_fastq.py pathToGtf pathToSequenc pathToConfig numReads pathToFastq readType\n"
            "pathToGtf      = Path to the _ensembl_ GTF file\n"
            "pathToSequence = Path to the _ensembl_ whole genome file, fasta format\n"
            "pathToConfig   = Configuration file. Comment lines begin with '#'.\n"
            "                 contains list of transcripts and their relative abundance \n"
            "                 (an integer) and a line with ReadLength\n"
            "numReads       = the number of reads generated\n"
            "pathToFastq    = the output file name prefix\n"
            "readType       = either : single, paired-fr-first, paired-fr-second"
    )
    sys.exit(Exit)





#********************************************************************************
#********************************************************************************
#********************************  MAIN  ****************************************
#********************************************************************************
#********************************************************************************
def main():
    timeBegin = time.time()
    if(len(sys.argv) != 7):
        if(len(sys.argv) > 1 and (sys.argv[1] == "--help" or sys.argv[1] == "-h")):
            print_help(0)
        else:
            print_help(1)
    #pathToGtf = "/reference/homo_sapiens/GRCh38/ensembl/Annotation/Genes/gtf/Homo_sapiens.GRCh38.83.gtf"
    #pathToSeq = "/reference/homo_sapiens/GRCh38/ensembl/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    random.seed(42)
    pathToGtf = sys.argv[1]
    pathToSeq = sys.argv[2]
    pathToConfig = sys.argv[3]
    pathToFastq = sys.argv[5]
    readType  = sys.argv[6]

    gtfList   = read_gtf(pathToGtf)
    exonList  = get_exon_list(gtfList)
    transList = get_transcript_list(gtfList, exonList)
    geneList  = get_gene_list(gtfList, transList)
    chrmList  = read_genome(pathToSeq)
    uniqueFeatureList = get_list_of_unique_gtf_features(gtfList)
    get_exon_seq(exonList, chrmList)
    link_exons_trans_and_genes(gtfList, exonList, transList, geneList)
    # print_transcripts_with_seqs(transList)      # Debug link_exons_trans_and_genes()

    geneDict, transDict = create_gene_and_trans_lookup_dict(geneList, transList)
    print_gtf_statistics(exonList, transList, geneList)
    # find_trans_that_differ_by_1_exon(geneList, transList) # Uncomment for complete list
    readLength, desiredTransList, abundanceList, numOfReads = read_config(pathToConfig)

    numOfReads = int(sys.argv[4])

    if(readType != 'single' and readType != 'paired-fr-first' and
       readType != 'paired-fr-second'
    ):
        exit_with_error("ERROR!!! Incorrect value for {}".format(readType))
    else:
        ### Paired end reads are not working yet ###
        if(readType == 'paired-fr-first' or readType == 'paired-fr-second'):
            exit_with_error("ERROR!!! paired-fr-first and paired-fr-second \n"
                            "not yet implemented. \n\n"
                            "NOTE:: Both reads are tentatively found in the \n"
                            "       INSERT class. The second read is not used.\n"
                            "       The second read should definitely needs checked.\n")
        create_fastq_file(pathToFastq, desiredTransList, abundanceList, numOfReads,
                          readLength, transDict, transList, exonList, readType)
       

    print("Unique features in Gtf : ")
    for feature in uniqueFeatureList:
        print("\t%s"%(feature))
    timeEnd = time.time()
    print("Run time : %s"%(timeEnd - timeBegin))
    sys.exit(0)


if __name__ == "__main__":
    main()

