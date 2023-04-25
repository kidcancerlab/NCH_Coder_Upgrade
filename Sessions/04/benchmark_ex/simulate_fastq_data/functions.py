# Author : Ali Snedden
# Date   : 5/31/19
# License: MIT
# Purpose: 
#
# Notes :
#
# Questions:
#
# References :  
#
import sys
import re
import datetime
import random
import copy
import operator
import numpy
from classes import GTF_ENTRY
from classes import FASTQ_READ
from classes import CHROMOSOME
from classes import GENE
from classes import TRANSCRIPT
from classes import INSERT
from classes import EXON
from classes import reverse_complement
from error import exit_with_error


#********************************************************************************
#********************************************************************************
#******************************  FUNCTIONS **************************************
#********************************************************************************
#********************************************************************************
def read_gtf(PathToGtf = None):
    """
    ARGS:
        PathToGtf = path to gene transfer file

    RETURN:
        A list of GTF_ENTRYs

    DESCRIPTION:
        Reads in gtf file.

    DEBUG:
        Can reproduce input gtf file. Works as expected.

    FUTURE: 
    """
    if(PathToGtf[-3:] != "gtf"):
        exit_with_error("ERROR! You did not pass a file with the .gtf extention\n")

    gtfFile = open(PathToGtf, 'r')
    gtfList = []
    gtfEntry = None
    timeBegin = datetime.datetime.now()

    for line in gtfFile:
        if(line[0] == "#"):
            continue
        line = line.split("\t")
        # Format check
        if(len(line) != 9):
            exit_with_error("ERROR! There should be 9 tab separated columns in a"
                            " GTF file\nYou only have %i\n"%(len(line)))
        gtfEntry = GTF_ENTRY(Chromosome = line[0], Source = line[1], EntryType = line[2],
                             Start = line[3], Stop = line[4], Score = line[5], Strand = line[6],
                             Frame = line[7], Attribute = line[8])
        gtfList.append(gtfEntry)

    gtfFile.close()
    timeEnd = datetime.datetime.now()
    # Debug Code
    #for gtfEntry in gtfList:
    #    gtfEntry.print_entry()
    print("read_gtf()      run time = %s"%(timeEnd - timeBegin))
    return gtfList


def get_exon_list(gtfList):
    """
    ARGS:
        PathToGtf = path to gene transfer file

    RETURN:
        A list of GTF_ENTRYs

    DESCRIPTION:
        Note : this list is not unique. I.e. there are exons that will be duplicated.
        I make no effort to correct / eliminate duplicates b/c it would complicate
        mapping trans.exonIdxList for all the TRANSCRIPTs

    DEBUG:

    FUTURE: 
    """
    exonList = []
    timeBegin = datetime.datetime.now()
    for gtfEntry in gtfList:
        if(gtfEntry.etype == "exon"):
            exon = EXON(gtfEntry)
            exonList.append(exon)
    timeEnd = datetime.datetime.now()
    print("get_exon_list() run time = %s"%(timeEnd - timeBegin))
    return exonList


def get_transcript_list(gtfList, allExonList):
    """
    ARGS:
        GtfList = list of GTF_ENTRYs
        AllExonList = list of EXON entries

    RETURN:
        A list of TRANSCRIPTS

    DESCRIPTION:
        Scales like crap. O(N**2), which is too slow for the whole genome + gtf

    DEBUG: 
        Appears to correctly read in transcripts and associated exons.
        I checked by eye on first 87 lines of grch38.83.gtf

    FUTURE: 
    """
    transList = []
    timeBegin = datetime.datetime.now()
    for gtfEntry in gtfList:
        if(gtfEntry.etype == "transcript"):
            trans = TRANSCRIPT(gtfEntry)
            transList.append(trans)
    timeEnd = datetime.datetime.now()
    print("get_transcript_list() run time = %s"%(timeEnd - timeBegin))
    return transList



def get_gene_list(gtfList, allTransList):
    """
    ARGS:
        GtfList = list of GTF_ENTRYs
        AllTransList = list of TRANSCRIPT entries

    RETURN:
        A list of GENES

    DESCRIPTION:
        Scales like crap. O(N**2), which is too slow for the whole genome + gtf

    DEBUG: 
        Appears to correctly read in genes and associated transcripts.
        I checked by eye on first 87 lines of grch38.83.gtf

    FUTURE: 
    """
    geneList = []
    timeBegin = datetime.datetime.now()
    for gtfEntry in gtfList:
        if(gtfEntry.etype == "gene"):
            gene = GENE(gtfEntry)
            geneList.append(gene)
    timeEnd = datetime.datetime.now()
    print("get_gene_list() run time = %s"%(timeEnd - timeBegin))
    return geneList





def read_genome(seqPath):
    """
    ARGS:

    RETURN:
        A list of Chromosomes

    DESCRIPTION:
        Reads in a _whole_ genome fasta file and returns a list of Chromosomes. It
        follows the ensemble file format. Read in all sequences as upper case.

    DEBUG: 
        For a shortened whole genome fasta file, it correctly reads it in.

    FUTURE: 
    """
    seqFile = open(seqPath, "r")    # 1st fasta file 
    chrmList = []
    chrm = None
    seq = ""
    timeBegin = datetime.datetime.now()

    # Read file
    for line in seqFile:
        if(line[0] == "#"):                  # Skip ali's comments in shortened file
            continue
        if(line[0] == '>' or line[0] == '\n'):
            # Create chromosome, reset seq
            if(chrm is not None):
                #sys.stdout.write("chromosome : %s\n"%(chrm))
                chrmList.append(CHROMOSOME(chrm, seq))
                seq = ""
            #chrm = line[1:]            # Skip '>' for chrm name
            line = re.split(' ', line)      # only want chrm name
            chrm = (line[0])[1:]        # Skip '>' for chrm name
            continue
        else:
            line = re.split('\n', line)     # line is now a list of strings
            seq += line[0].upper()
    chrmList.append(CHROMOSOME(chrm, seq))  # pick up final chromosome.
    #sys.stdout.write("chromosome : %s\n"%(chrm))

    # Short debugging section, check that genome.fa duplicates original
    #for chrm in chrmList:
    #    sys.stdout.write("\n%s\n"%chrm.chrm)
    #    for i in range(0,len(chrm.seq)):
    #        sys.stdout.write("%c"%(chrm.seq[i]))
    #        if((i+1) % 60 == 0 and i != 0):
    #            sys.stdout.write("\n")

    seqFile.close()
    timeEnd = datetime.datetime.now()
    print("read_genome()   run time = %s"%(timeEnd - timeBegin))
    return chrmList




def get_exon_seq(exonList, chrmList):
    """
    ARGS:

    RETURN:
        None

    DESCRIPTION:
        Gets sequences for exons from chromosome list

    DEBUG: 
        Spot checked 3 exons, are all ok. More testing needed, however it is challenging
        to get a list of all the exons (incl. seqs) in a single file

    FUTURE: 
    """
    timeBegin = datetime.datetime.now()
    for exon in exonList:
        for chrm in chrmList:
            chrmLen = len(chrm.seq)

            if(chrm.chrm == exon.chrm):
                start = exon.start - 1        # -1 b/c python is 0 indexed, gtf file isnot
                end   = exon.stop
                if(start >= chrmLen or end >= chrmLen):
                    exit_with_error("ERROR!! start (%i) or stop (%i) Position > "
                                    "chromosome length (%i)\n"%(start, end, chrmLen))
                if(exon.strand == '+'):
                    exon.seq = chrm.seq[start:end]
                elif(exon.strand == '-'):
                    exon.seq = reverse_complement(chrm.seq[start:end])
                    tmp = exon.start
                    exon.start = exon.stop
                    exon.stop  = tmp
                else:
                    exit_with_error("ERROR! strand char = %s is invalid", exon.strand)
    timeEnd = datetime.datetime.now()
    print("get_exon_seq()  run time = %s"%(timeEnd - timeBegin))
    


def get_trans_seq(transList, exonList):
    """
    ARGS:

    RETURN:
        None

    DESCRIPTION:
        Gets sequences for transcripts from chromosome list

    DEBUG: 
        Tested on 2 transcripts, more testing required. Getting a transcript file
        with the transcripts and sequences is challenging though

    FUTURE: 
    """
    timeBegin = datetime.datetime.now()
    for trans in transList:
        exonNum = 0         # use to check that indexs are loaded in order
        prevExonNum = 0
        for exon in exonList:
            if(exon.transID == trans.transID):
                exonNum = int(exon.exonNum)
                if(exonNum - prevExonNum != 1):
                    exit_with_error("ERROR! exon numbers for %s are loaded out of "
                                    "order!\n"%(trans.transID))
                if(trans.seq is None):
                    trans.seq = exon.seq
                else:
                    trans.seq += exon.seq
                prevExonNum = exonNum
    timeEnd = datetime.datetime.now()
    print("get_trans_seq() run time = %s"%(timeEnd - timeBegin))




def get_list_of_unique_gtf_features(gtfList):
    """
    ARGS:
        gtfList : list of all GTF_ENTRYs

    RETURN:
        uniqueFeatureList : list of unique features

    DESCRIPTION:
        Finds all the unique features (ie column 3) of the GTF file

    DEBUG: 

    FUTURE: 
    """ 
    featureList = []
    uniqueFeatureList = []
    for gtfEntry in gtfList:
        featureList.append(gtfEntry.etype)
    for feature in featureList:
        if feature not in uniqueFeatureList:
            uniqueFeatureList.append(feature)
    return uniqueFeatureList



def link_exons_trans_and_genes(gtfList, exonList, transList, geneList):
    """
    ARGS:
        gtfList  : list of all GTF_ENTRYs
        exonList : list of all EXONS 
        transList: list of all TRANSCRIPTS
        geneList : list of all GENES

    RETURN:

    DESCRIPTION:
        Loops through gtfList and captures the indices of exons in exonList and passes
        it to the transcripts in transList.  Also captures indices of transcripts in
        transList and passes it to genes in geneList.
        
        Does this in one pass through gtfList and scales roughly O(N). Should be faster
        than previous versions.

    DEBUG: 
        1. I validated by using print_transcripts_with_seqs() and comparing against the
           biomart download for chromosome 1. My data file was _identical_ to biomart's. 

           For how this was done, see the debug comment in print_transcripts_with_seqs() 
        2. Checked Transcript.seq for reverse strand ('-') transcript. Used 
           ENST00000488147 it is correct.
        
    FUTURE: 
    """ 
    gIdx = 0            # Gene index, for geneList
    tIdx = 0            # Transcript index, for transList
    eIdx = 0            # Exon index, for exonList
    gtfListLen = len(gtfList)
    timeBegin = datetime.datetime.now()
    
    # Ugly / non-pythonic b/c cant find cleaner way of accessing the next gtfEntry in the list
    for i in range(len(gtfList)):
        if(gtfList[i].etype == "gene"):
            # Check that genes in geneList are same order as gtfList 
            if(gtfList[i].geneID != geneList[gIdx].geneID):
                exit_with_error(
                    "ERROR! gtfList[%i].geneID = %s and geneList[%i].geneID = %s"%(
                    i, gtfEntry.geneID, gIdx, geneList[gIdx].geneID)
                )
            j = i + 1

            # Get constituent transcripts between gene entries
            while(gtfList[j].etype != "gene"):
                if(gtfList[j].etype == "transcript"):

                    # Check that transcripts in transList are same order as gtfList 
                    # Checking transcripts after gene in gtf _actually_ are members of the gene
                    # Add trans info to appropriate geneList[]
                    if(gtfList[j].transID == transList[tIdx].transID and
                      gtfList[i].geneID == transList[tIdx].geneID    and
                      gtfList[i].geneID == geneList[gIdx].geneID):
                        geneList[gIdx].transList.append(transList[tIdx].transID)
                        geneList[gIdx].transIdxList.append(tIdx)
                        k = j + 1

                        # Get constituent exons between transcript entries
                        while(gtfList[k].etype != "transcript"):
                            if(gtfList[k].etype == "exon"):
                                # Check exons in exonList are same order as gtfList 
                                # Checking exons after trans in gtf are members trans
                                # Add exon info to appropriate transList[]
                                if(gtfList[k].transID == exonList[eIdx].transID and
                                  gtfList[i].geneID == exonList[eIdx].geneID   and
                                  gtfList[i].geneID == geneList[gIdx].geneID):
                                    transList[tIdx].exonList.append(exonList[eIdx].exonID)
                                    transList[tIdx].exonIdxList.append(eIdx)
                                    eIdx += 1
                                else:
                                    exit_with_error(
                                     "ERROR! gtfList[%i].transID = %s and exonList[%i]."
                                     "transID = %s\n\tgtfList[%i].geneID = %s and "
                                     "transList[%i].geneID = "
                                     "%s\n\tand geneList[%i].geneID = %s\n"%
                                     (k, gtfList[k].transID,eIdx, exonList[eIdx].transID,
                                     k, gtfList[k].geneID, tIdx, transList[tIdx].geneID,
                                     gIdx, geneList[gIdx].geneID)
                                    )
                            k += 1
                            if(k == gtfListLen):
                                break
                        tIdx += 1   
                    else:
                        exit_with_error(
                            "ERROR! gtfList[%i].transID= %s and transList[%i].transID = "
                            "%s\n\tgtfList[%i].geneID = %s and transList[%i].geneID = "
                            "%s\n\tand geneList[%i].geneID = %s\n"%
                            (j, gtfList[j].transID, tIdx, transList[tIdx].transID,
                            j, gtfList[j].geneID, tIdx, transList[tIdx].geneID, gIdx,
                            geneList[gIdx].geneID))
                j += 1
                if(j == gtfListLen):
                    break
            gIdx += 1
            
    # Now get transcript sequences.
    for trans in transList:
        trans.seq = ""
        for eIdx in trans.exonIdxList:
            trans.seq += exonList[eIdx].seq

    timeEnd = datetime.datetime.now()
    print("link_exons_trans_and_genes() run time = %s"%(timeEnd - timeBegin))

    

def create_gene_and_trans_lookup_dict(geneList, transList):
    """
    ARGS:
        geneList  : list of GENEs
        transList : list of TRANSCRIPTS

    RETURN:
        geneDict  : keys = geneID, value geneList indices
        transDict : keys = transID, value transList indices

    DESCRIPTION:
        Dictionaries which are associative arrays with the
        geneID and transID as the keys and the associated geneList and transList
        indices as the values.

    DEBUG: 
        Spot checked the resulting dictionaries, appears to be correct

    FUTURE: 
    """ 
    geneDict = {}
    transDict = {}
    timeBegin = datetime.datetime.now()
    
    for gIdx in range(len(geneList)):
        geneDict[(geneList[gIdx].geneID)[1:-1]] = gIdx      # [1:-1] to strip leading and trailing "
    for tIdx in range(len(transList)):
        transDict[(transList[tIdx].transID)[1:-1]] = tIdx

    timeEnd = datetime.datetime.now()
    print("create_gene_and_trans_lookup_dict() run time = %s"%(timeEnd - timeBegin))
    return geneDict, transDict



def print_gtf_statistics(exonList, transList, geneList):
    """
    ARGS:
        exonList : list of EXONs
        transList : list of TRANSCRIPTs
        geneList  : list of GENEs

    RETURN:

    DESCRIPTION:
        Prints interesting statistics
        

    DEBUG: 
        I spot checked the results of transPairsDiffFile.txt and it appears to 
        return sets of transcripts that actually differ by 1 exon.  I checked 
        by searching in ensembl.

        Using a short.gtf, I checked the using Excel (identical):
            meanExonLen
            sigmaExonLen
            meanTransLen
            sigmaTransLen
            maxTransLen
            minTransLen
            minExonLen
            maxExonLen
            exonWMaxLen
            exonWMinLen
            transWMaxLen
            transWMinLen
        
    FUTURE: 
    """ 
    # exons
    meanExonLen = 0
    sigmaExonLen= 0
    maxExonLen  = 0
    minExonLen  = 10000
    exonWMaxLen = 0
    exonWMinLen = 10000
    meanExonNum = 0
    sigmaExonLen = 0
    # transcripts
    meanTransLen= 0
    maxTransLen = 0
    minTransLen = 10000
    transWMaxLen = ""
    transWMinLen = ""
    sigmaTransLen = 0
    # gene
    meanTransNum = 0
    timeBegin = datetime.datetime.now()

    # genes 
    for gene in geneList:
        meanTransNum += len(gene.transIdxList)
    
    # transcripts
    for trans in transList:
        transLen = len(trans.seq)
        meanTransLen += transLen
        # Try max and min later...
        if(transLen < minTransLen):
            minTransLen = transLen
            transWMinLen = trans.transID
        if(transLen > maxTransLen):
            maxTransLen = transLen
            transWMaxLen = trans.transID
        meanExonNum += len(trans.exonIdxList)

    # Exons
    for exon in exonList:
        exonLen     = len(exon.seq)
        meanExonLen += exonLen
        #Find max and min
        if(exonLen < minExonLen):
            minExonLen = exonLen
            exonWMinLen = exon.exonID
        if(exonLen > maxExonLen):
            maxExonLen = exonLen
            exonWMaxLen = exon.exonID

    # Averages
    meanTransLen  = meanTransLen / float(len(transList))
    meanTransNum  = meanTransNum / float(len(geneList))
    meanExonNum   = meanExonNum  / float(len(transList))
    meanExonLen   = meanExonLen  / float(len(exonList))
    
    # Standard deviations 
    # Transcripts
    for trans in transList:
        transLen = len(trans.seq)
        sigmaTransLen += (transLen - meanTransLen)**2
    sigmaTransLen = (sigmaTransLen / float(len(transList) - 1))**0.5
    # Exons
    for exon in exonList:
        sigmaExonLen += (len(exon.seq) - meanExonLen)**2
    sigmaExonLen = (sigmaExonLen / float(len(exonList) - 1))**0.5
    

    print("\nMean exon length : %f"%(meanExonLen))
    print("Std. Dev. exon length : %f"%(sigmaExonLen))
    print("Exon w/ Min length : %s, len : %i"%(exonWMinLen, minExonLen))
    print("Exon w/ Max length : %s, len : %i\n"%(exonWMaxLen, maxExonLen))
    print("Mean trans length : %f"%(meanTransLen))
    print("Std. Dev. trans Length : %s"%(sigmaTransLen))
    print("Trans w/ Min length : %s, len : %i"%(transWMinLen, minTransLen))
    print("Trans w/ Max length : %s, len : %i"%(transWMaxLen, maxTransLen))
    print("Mean num of exons per trans : %f\n"%(meanExonNum))

    timeEnd = datetime.datetime.now()
    print("print_gtf_statistics() run time = %s"%(timeEnd - timeBegin))
                        
                    

def find_trans_that_differ_by_1_exon(geneList, transList):
    """
    ARGS:
        transList : list of TRANSCRIPTs
        geneList  : list of GENEs

    RETURN:

    DESCRIPTION:
        Prints interesting statistics
        

    DEBUG: 
        I spot checked the results of transPairsDiffFile.txt and it appears to 
        return sets of transcripts that actually differ by 1 exon.  I checked 
        by searching in ensembl.
    FUTURE: 
    """ 
    transPairsDiffFile = open("transPairsDiffFile.txt", "w+")
    transPairsDiffFile.write("# This file contains a list of transcripts that differ only by 1 exon\n")
    transPairsDiffFile.write("# exon transcript1 transcript2 exonDiff\n")

    # Find transcripts that vary by 1 exon. Write to transPairsDiffFile.txt
    # This for loop should probably be its own separate program, but since
    # I already have all the functions in this program
    # This loop could be commented out in the final version of the code.
    #
    # Find transcripts that differ by 1 exon
    for gene in geneList:
        for tidx1 in range(len(gene.transIdxList)):              # transcriptIdx1
            trans1 = transList[gene.transIdxList[tidx1]]
            tidx2 = tidx1 + 1
            while tidx2 < len(gene.transIdxList):
                trans2 = transList[gene.transIdxList[tidx2]]
                
                # Follow Differ by one exon only
                if((len(trans1.exonList) > len(trans2.exonList) + 1)   or 
                   (len(trans1.exonList) < len(trans2.exonList) - 1)):
                    tidx2 += 1              # without this, it will loop infinitely
                    continue    
                else:
                    # below two operations return the elements in argument 1 that are not in
                    # argument 2. Hence why we do it twice, exonDiff and exonDiff2 won't be the
                    # same when there is a difference of 1 exon.
                    exonDiff = list(set(trans1.exonList) - set(trans2.exonList))
                    exonDiff2 = list(set(trans2.exonList) - set(trans1.exonList))
                    #exonDiff = [exn for exn in trans1.exonList if exn not in trans2.exonList]
                    if(len(exonDiff) == 1 and len(exonDiff2) == 0):
                        transPairsDiffFile.write("%s %s %s %s\n"%(gene.geneID, trans1.transID,
                                                 trans2.transID, exonDiff))
                    if(len(exonDiff) == 0 and len(exonDiff2) == 1):
                        transPairsDiffFile.write("%s %s %s %s\n"%(gene.geneID, trans1.transID,
                                                 trans2.transID, exonDiff2))
                tidx2 += 1
    transPairsDiffFile.close()





def read_config(pathToConfig):
    """
    ARGS:
        pathToConfig : str, path to configuration file

    RETURN:
        readLength      = readlength desired
        desiredTransList= List of transcripts to use
        abundanceList   = list of relative abundances of transcripts

    DESCRIPTION:
        Config file format :
            1. Comment lines begin with '#'
            2. first non-header line begins with 'ReadLength'
            3. All subsequent lines must be transcripts with relative abundance
               The relative abundance can be any integer. Scaling is done 
               automatically. 
               E.g.
                    ENST00000488147 10
                    ENST00000473358 5
    DEBUG: 
        For small config files it reads in all the fields correctly.

    FUTURE: 
    """ 
    desiredTransList = []
    abundanceList = []              # integers used to get relative abundance of transcripts
    readLength = 0
    numOfReads = 0

    configFile = open(pathToConfig, 'r')
    for line in configFile:
        if(line[0] == "#"):
            continue
        line = (line.split("\n"))[0]        # remove trailing \n

        # Check for tabs, only spaces permitted
        if(re.search('\t',line)):
            exit_with_error("ERROR! Tabs not permitted in config file!\n")
        line = line.split(" ")

        # ReadLength
        if(line[0] == "ReadLength"):
            if(readLength == 0):
                readLength = int(line[1])
                continue
            else:
                exit_with_error("ERROR! multiple instances of ReadLength in config "
                                "file\n")

        # NumberOfReads
        if(line[0] == "NumberOfReads"):
            if(numOfReads == 0):
                numOfReads = int(line[1])
                continue
            else:
                exit_with_error("ERROR! multiple instances of ReadLength in config "
                                "file\n")
        
        # Transcripts
        if(re.search('ENST', line[0])):
            desiredTransList.append(line[0])
            abundanceList.append(int(line[1]))
        else:
            exit_with_error("ERROR! Incorrect transcript entry : %s\n"
                            " All entries should begin with 'ENST'\n"%(line))
            
    if(readLength == 0 or numOfReads == 0):
        exit_with_error("ERROR! ReadLength or NumberOfReads not specified in "
                        "config.txt\n")

    print("Config File Parameters : \nReadLength : %i\nNumberOfReads : %i"%(readLength,
          numOfReads))
    i = 0
    for trans in desiredTransList:
        print("%s %i"%(trans,abundanceList[i]))
        i += 1
    print("\n")

    return readLength, desiredTransList, abundanceList, numOfReads



def print_transcripts_with_seqs(transList):
    """
    ARGS:
        transList

    RETURN:
        None

    DESCRIPTION:
        This takes in a transList, makes a deep copy b/c I sort it.  I do not want 
        to mess up the actual order of transList, since the indices stored by 
        GENEs in transIdxList[] depends on the maintaining the order in transList.
    
        This function sorts all the transcripts numerically.
        This function is primarily used to debug link_exons_trans_and_genes()

    DEBUG: 
        1. Created chr1.fa and chr1.gtf, passed as args to this program
        2. Ran this function to print out _sorted_ transcript fasta file, called
           trans_list_and_seq.fa
        3. Went to http://useast.ensembl.org/biomart/martview, selected 
           chr1, sequences->cdna and downloaded the file. Downloaded as 
           biomart_export_cdna.fa 
        4. Created a program, validate/sort_biomart.py to read in biomart_export_cdna.fa 
           and generated sorted_biomart_export_cdna.fa
        5. ran 'diff trans_list_and_seq.fa validate/sorted_biomart_export_cdna.fa', it
           was identical.
           
           CONCLUSION:
               I correctly handle mapping the sequences to exons, and exons to 
               transcripts. If these weren't done correctly, (e.g. not handling 
               reverse complements), there is no way that this would work correctly.
        
               The following functions _MUST_ be correctly working (at least w/r/t
               transcripts):
                        GTF_ENTRY.__init__() 
                        CHROMOSOME.__init__()
                        TRANSCRIPT.__init__()
                        EXON.__init__()
                        read_gtf()  
                        get_exon_list()
                        get_transcript_list()
                        read_genome()
                        get_exon_seq()
                        get_trans_seq()
                        reverse_complement()
                        ---->  link_exons_trans_and_genes() <----
    FUTURE: 
    """ 
    transCopyList = copy.deepcopy(transList)
    outFile = open("trans_list_and_seq.fa", "w+")
    transCopyList.sort(key=operator.attrgetter('transNum'))
    
    for trans in transCopyList:
        outFile.write(">%s|%s\n"%(trans.geneID[1:-1], trans.transID[1:-1]))
        i = 0 

        for i in range(len(trans.seq)):
            outFile.write("%s"%(trans.seq[i]))
            if((i + 1) % 100 == 0):
                outFile.write("\n")

        outFile.write("\n")

    outFile.close()
    

def create_insert(Transcript = None, ReadLength = None, Mu = None, Sigma = None,
                  ExonList = None):
    """
    ARGS:
        Transcript : a TRANSCRIPT instance
        ReadLength : length of reads. Different from insert length
        Mu         : the mean of fragment length distribution
        Sigma      : the standard deviation of fragment length distribution
        ExonList   : 

    RETURN:
        AN INSERT of length n, where n fall in a distribution of rnorm(Mu,sigma)

    DESCRIPTION:

    DEBUG: 

    FUTURE: 
        1. Implement Proper solution where insert going into the Illumina adapter when
           stop - start < ReadLength
    """
    start = 0
    stop = 0
    timeBegin = datetime.datetime.now()
    transLength = len(Transcript.seq)

    # type check
    if(not isinstance(Transcript, TRANSCRIPT)):
        exit_with_error("ERROR! Transcript is not of class type TRANSCRIPT 1\n")

    insertLength = 0

    # Ensure inserts are at least as long as the readlength
    while(insertLength < ReadLength):
        start = random.randint (0, transLength -1 )    
        stop = start + int(numpy.random.normal(Mu, Sigma))
        # Avoid unrealistically short inserts
        if(stop - start < ReadLength):
            continue
        # Avoid inserts that are past end of transcripts.
        if(stop > transLength - 1):
            # Proper solution here would have insert going into the Illumina adapter
            stop = transLength - 1  
            if(stop - start < ReadLength): # Insert must be at least as large as a read
                continue
        insert = INSERT(Transcript = Transcript, StartWrtTrans = start,
                        StopWrtTrans = stop, ReadLength = ReadLength, 
                        ExonList = ExonList)

        insertLength = len(insert.seq)        

    timeEnd = datetime.datetime.now()
    # print("get_insert_list() run time = %s"%(timeEnd - timeBegin))
    return insert


def create_fastq_file(pathToFastq, desiredTransList, abundanceList, nReads,
                      readLength, transDict, transList, exonList, readType):
    """
    ARGS:
        pathToFastq      : Path to output fastq file
        desiredTransList : transcripts read from config file
        abundanceList    : list of integers that sum to a value used to normalize 
                           the number of reads. 
                                E.g. trans1 has 5 and trans2 has 10, 
                                     the ratio of trans1 : trans2 = 1 : 2
        nReads           : Rough number of desired reads, the ratios from abundanceList
                           is maintained at the expense of nReads. 
                                E.g from the above example if nReads = 10, 
                                    the actual number of reads would be 
                                    3 for trans1, 6 for trans2
        readLength       : length of reads
        transDict        : Dictionary used map transID quickly to the correct
                           transcript in transList
        transList        : List of TRANSCRIPTs. Contains sequence to pass to instance of 
                           FASTQ_READ()
        exonList         : List of EXONs. Passed to FASTQ_READ() to get metadata for each
                           fastq read. E.g. the start position and exons a read spans.
        readType         : either : single, paired-fr, paired-rf"
    RETURN:
        None. File written

    DESCRIPTION:
        Writes fastq file.

    DEBUG: 
        1. Blasted against ensembl database, spot checked a couple of transcripts.
           Need further debugging. 

           Took synthetic_sample.fastq and operated in vim on transcript ENST00000473358:
            Exons are ordered as : ENSE00001947070 ENSE00001922571 ENSE00001827679
            Copied synthetic_sample.fastq to poop.fq

            ****** in vim ******
            %s/^@Read.\{1,1000\}:start:\([0-9]\{1,100\}\):exons:\(.\{1,100\}\)\n\(^[A-Z]\{50\}\)\n^+\n\(^.\{50\}\)/\1\t\3\t\2/gc
            %s/:/\t/g   # remove colon sep exon names
            %s/"//g     # remove " around exon names

            ****** in shell, want exon reads start at (see order above) ******
            ****** Avoid grepping enties with start positions on the exon prior ******
            grep ENSE00001947070 poop.fq &> ENSE00001947070.txt
            grep ENSE00001922571 poop.fq | grep -v ENSE00001947070 &> ENSE00001922571.txt
            grep ENSE00001827679 poop.fq | grep -v ENSE00001922571 &> ENSE00001827679.txt
            awk '{print $1 "\t" $2}' ENSE00001947070.txt &> ENSE00001947070_1and2.txt
            awk '{print $1 "\t" $2}' ENSE00001922571.txt &> ENSE00001922571_1and2.txt
            awk '{print $1 "\t" $2}' ENSE00001827679.txt &> ENSE00001827679_1and2.txt
            awk '{print $2}' ENSE00001947070.txt | xargs -I{} grep -aob {} ENST00000473358.txt | 
                    awk 'BEGIN{FS=":"}{start = $1 + 29554; print start "\t" $2}' &> awk_out.txt
            awk '{print $2}' ENSE00001922571.txt | xargs -I{} grep -aob {} ENST00000473358.txt | 
                    awk 'BEGIN{FS=":"}{start = $1 + 30564 - 486; print start "\t" $2}' &> awk_out.txt
            awk '{print $2}' ENSE00001827679.txt | xargs -I{} grep -aob {} ENST00000473358.txt | 
                    awk 'BEGIN{FS=":"}{start = $1 + 30976 - 486 - 104; print start "\t" $2}' &> awk_out.txt
            Used diff to compare all the awk_out.txt to ENSE*_1and2.txt files.
                    CONCLUSION : they are identical. Therefor I get the correct start position from the
                                 correct sequences.
                    THEREFOR : I believe that create_fastq_file and FASTQ_READ() are working as expected.

        2. See debug comments of INSERT class.
           CONCLUSION : single end reads of transcripts/inserts on the '+' strand
                        in the sense direction work.

    FUTURE: 
        Include more error checking for goofy parameters, e.g. not enough reads for
        the ratios, etc.
    """
    abundanceSum = 0
    transIdx = 0
    readIdx = 0
    
    for abundance in abundanceList:
        abundanceSum += abundance
    #abundanceNormalization = abundanceNormalization / len(abundanceList)    # integer division
    if(abundanceSum < 1):
        exit_with_error("ERROR! abundanceSum = {}\nPlease enter abundance "
                        "values > 1\n".format(abundanceNormalization))


    if(readType == 'single'):
        pathToFastqR1 = pathToFastq + ".fq"
        fastqFileR1 = open(pathToFastqR1, "w+")
        fastqListR1 = []
        
    elif(readType == 'paired-fr-first' or readType == 'paired-fr-second'):
        pathToFastqR1 = pathToFastq + "-R1.fq"
        pathToFastqR2 = pathToFastq + "-R2.fq"
        fastqFileR1 = open(pathToFastqR1, "w+")
        fastqFileR2 = open(pathToFastqR2, "w+")
        fastqListR1 = []
        fastqListR2 = []
    else:
        exit_with_error("ERROR!!! Incorrect value for {}".format(readType))

    
    for transName in desiredTransList:
        try:
            trans = transList[transDict[transName]]
        except KeyError:
            exit_with_error("ERROR! {} is not a transcript annotated in your "
                            "gtf file\n".format(transName))

        for i in range(int(float(abundanceList[transIdx])/float(abundanceSum) * nReads)):
            insert = create_insert(trans, readLength, 150, 15, exonList)


            if(readType == 'single'):
                fastqEntry = FASTQ_READ(Insert = insert, ReadLength = readLength,
                                        MetaData = "@Read_num:%i"%(readIdx),
                                        ExonList=exonList, Direction = "forward")
                fastqListR1.append(fastqEntry)

            elif(readType == 'paired-fr-first'):
                fastqEntry = FASTQ_READ(Insert = insert, ReadLength = readLength,
                                        MetaData   = "@Read_num:%i"%(readIdx),
                                        ExonList   = exonList, Direction = "reverse")
                fastqListR1.append(fastqEntry)
                fastqEntry = FASTQ_READ(Insert = insert, ReadLength = readLength,
                                        MetaData   = "@Read_num:%i"%(readIdx),
                                        ExonList   = exonList, Direction = "forward")
                fastqListR2.append(fastqEntry)
            elif(readType == 'paired-fr-second'):
                fastqEntry = FASTQ_READ(Insert = insert, ReadLength = readLength,
                                        MetaData   = "@Read_num:%i"%(readIdx),
                                        ExonList   = exonList, Direction = "forward")
                fastqListR1.append(fastqEntry)
                fastqEntry = FASTQ_READ(Insert = insert, ReadLength = readLength,
                                        MetaData   = "@Read_num:%i"%(readIdx),
                                        ExonList   = exonList, Direction = "reverse")
                fastqListR2.append(fastqEntry)

            readIdx += 1
        transIdx += 1

    if(readType == 'single'):
        for fastqEntry in fastqListR1:
            fastqFileR1.write("%s\n%s\n+\n%s\n"%(fastqEntry.metadata,
                              fastqEntry.seq, fastqEntry.qual))
        fastqFileR1.close()
    else:
        for fastqEntry in fastqListR1:
            fastqFileR1.write("%s\n%s\n+\n%s\n"%(fastqEntry.metadata,
                              fastqEntry.seq, fastqEntry.qual))
        
        for fastqEntry in fastqListR2:
            fastqFileR2.write("%s\n%s\n+\n%s\n"%(fastqEntry.metadata,
                              fastqEntry.seq, fastqEntry.qual))

        fastqFileR1.close()
        fastqFileR2.close()



