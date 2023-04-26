# Author : Ali Snedden
# Date   : 5/15/19
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
import numpy as np
from error import exit_with_error
#********************************************************************************
#********************************************************************************
#*******************************  CLASSES  **************************************
#********************************************************************************
#********************************************************************************
class GTF_ENTRY:
    """
    Container to hold each GTF entry. 

    See Example : 
    /reference/homo_sapiens/GRCh38/ensembl/Annotation/Genes/gtf/Homo_sapiens.GRCh38.83.gtf

    See Specification:
    http://useast.ensembl.org/info/website/upload/gff.html
    http://www.gencodegenes.org/gencodeformat.html
    """
    def __init__(self, Chromosome = None, Source = None, EntryType= None, Start = None,
                 Stop = None, Score = None, Strand = None, Frame = None,
                 Attribute = None):
        """
        ARGS:
            Chromosome  : Chromosome, can only be 1,2,...,X,Y
            Source      : Database or Project name for entry
            EntryType   : exon, transcript, CDS, gene, etc
            Start       : start position on chromosome
            Stop        : stop position on chromosome
            Score       : UNKNOWN
            Strand      : '+' (forward) or '-' (reverse)
            Frame       : 0 indexed position of first base of codon
            Attribute   : semicolon sep list of tag-value pairs
            
        RETURN:
            NONE : Initializes GTF_ENTRY

        DESCRIPTION:
    
        DEBUG:
            Tested by reading in gtf file and printing out list of GTF_ENTRYs.
            compared using full Homo_sapiens.GRCh38.83.gtf and the output 
            was _identical_.

        FUTURE: 
        """
        self.chrm    = None        # str, Chromosome
        self.src     = None        # str, Source of data
        self.etype   = None        # str, exon, transcript, CDS, gene
        self.start   = None        # int, Start position on chrm
        self.stop    = None        # int, End position on chrm
        self.score   = None        # str, expects empty field, ie : '.'
        self.strand  = None        # str, '+' (forward) or '-' (reverse)
        self.frame   = None        # str, 0 indexed position of first base of codon
        self.attribute = None      # str, semicolon sep list of tag-value pairs
        # below vars are parsed from attribute list
        self.geneID  = None        # str, Gene ID ... starts with ENSG..
        self.geneName= None        # str, common gene name
        self.transID = None        # str, transcript ID ... starts with ENST...
        self.transName= None       # str, transcript Name.. = gene common name + num
        self.exonID  = None        # str, exon ID ... starts with ENSE
        self.exonNum = None        # int, exon number
        self.biotype = None        # str, gene_biotype

        if(Chromosome is not None):
            self.chrm = Chromosome
        else:
            exit_with_error("ERROR! Chromosome not specified!\n")
        if(Source is not None):
            self.src = Source
        else:
            exit_with_error("ERROR! Source not specified!\n")
        if(EntryType is not None):
            self.etype = EntryType
        else:
            exit_with_error("ERROR! EntryType not specified!\n")
        if(Start is not None):
            self.start = int(Start)
        else:
            exit_with_error("ERROR! Start not specified!\n")
        if(Stop is not None):
            self.stop = int(Stop)
        else:
            exit_with_error("ERROR! Stop not specified!\n")
        if(Score is not None):
            self.score = Score
            # Check for empty field
            if(self.score != '.'):
                exit_with_error("ERROR! Score = %s. Expected empty field\n"%(self.score))
        else:
            exit_with_error("ERROR! Score not specified!\n")
        if(Strand is not None):
            self.strand = Strand
            if(self.strand != '+' and self.strand != '-'):
                exit_with_error("ERROR! Strand = %s, wrong format\n"%(self.strand))
        else:
            exit_with_error("ERROR! Strand not specified!\n")
        if(Frame is not None):
            self.frame = Frame
        else:
            exit_with_error("ERROR! Frame not specified!\n")
        if(Attribute is not None):
            self.attribute = Attribute.split("\n")[0]

            # Now parse attribute list for relevant fields. Create dictionary for easy lookup
            attributeDict = {}
            attributeSplit = self.attribute.split(";")[:-1]     # lop off last due to ; at end
            for attrEntry in attributeSplit:
                attrEntry = attrEntry.strip()                   # Eliminate white space on ends
                key   = attrEntry.split(" ")[0]
                value = attrEntry.split(" ")[1]
                attributeDict[key] = value
            # Every key is not necessarily in the attributeDict, handle accordingly
            try:  self.geneID   = attributeDict["gene_id"]
            except KeyError:  pass
            try:  self.geneName = attributeDict["gene_name"]
            except KeyError:  pass
            try:  self.transID  = attributeDict["transcript_id"]
            except KeyError:  pass
            try:  self.transName= attributeDict["transcript_name"]
            except KeyError:  pass
            try:  self.exonID   = attributeDict["exon_id"]
            except KeyError:  pass
            try:  self.exonNum  = int((attributeDict["exon_number"])[1:-1])
            except KeyError:  pass
            try:  self.biotype  = attributeDict["gene_biotype"]
            except KeyError:  pass
        else:
            exit_with_error("ERROR! Attribute not specified!\n")


    def print_entry(self):
        print("%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s"%(self.chrm, self.src, self.etype,
              self.start, self.stop, self.score, self.strand, self.frame,
              self.attribute))


class FASTQ_READ:
    def __init__(self, Insert = None, ReadLength = None, MetaData = None,
                 ExonList = None, Direction = None):
        """
        ARGS:
            Insert     = an INSERT instance
            ReadLength = length of desired read.
            MetaData   = Read number
            ExonList   = 
            Direction  = either 'forward' or 'reverse'
                         'forward' : matches mRNA starting from 5' -> 3'
                                     (the way ribosome trascribes mRNA)
                         'reverse' : matches complement mRNA starting from 3' -> 5' 

                         E.g. 
                         5' =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=> 3'
                            -----> 'forward'
                                                     'reverse' <-------- 
        RETURN:
            NONE : Initializes GENE

        DESCRIPTION:
            Currently, synthetic reads are forced to be completely contained within
            the transcript. i.e, I am not permitting them to read into the
            adapters, we'll save that for future work

        NOTES :
            1. http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html

        DEBUG:
            1. See debugging comments from INSERT class

        FUTURE: 
            1. Test both types of paired end reads  
            2. Test reads/transcrips on the +/- strands
        """
        self.seq  = None
        self.qual = None
        self.readlen = None
        self.metadata = None
        self.readDirection  = Direction   # whether the fastq read is forward / reverse comp
        self.start  = 0      # coord wrt to chromosome start (gtf start def this way)
        self.stop   = 0      # coord wrt to chromosome stop (gtf stop def this way)

        # type check 
        if(not isinstance(Insert, INSERT)):
            exit_with_error("ERROR! Insert is not of class type INSERT\n")

        if(ReadLength is not None):
            # ReadLength = int(ReadLength)
            self.readlen = ReadLength
        else:
            exit_with_error("ERROR! ReadLength not specified!\n")

        self.get_qual()
        # ALSO : recall that Transcript.seq is always in the direction that
        # transcripts are transcribed.
        # See : http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html
        if(self.readDirection == "forward"):
            self.seq = Insert.seq[0:ReadLength]
            # start / stop are _inclusive_
            self.start = Insert.r1Start
            self.stop  = Insert.r1Stop  
            exonsSpannedL = Insert.r1ExonL
        elif(self.readDirection == "reverse"):
            seqLen = len (Insert.seq)
            self.seq = reverse_complement(Insert.seq[seqLen-ReadLength:seqLen])
            # start / stop are _inclusive_
            self.start = Insert.r2Start
            self.stop  = Insert.r2Stop
            exonsSpannedL = Insert.r2ExonL
        else:
            exit_with_error("ERROR!!! {} is invalid".format(self.readDirection))


        # Add useful information to reads
        # shouldn't be so complicated. will need to revise Transcript.ExonList to
        # simplify this part
        if(MetaData is None):
            exit_with_error("ERROR! MetaData not specified!\n")

        #exonsSpannedL = list(set(exonsSpannedL))
        self.metadata = "%s:trans:%s:start:%i:exons"%(MetaData,Insert.transcript.transID,
                        self.start)
        for exon in exonsSpannedL:
            self.metadata = "{}:{}:{}:{}".format(self.metadata, exon.exonID,
                            exon.start, exon.stop)




    def get_qual(self):
        """
        ARGS:
            NONE 
        RETURN:
            quality scores based on Illumina, Phread+33 between chars '!' (33)
            and 'J' (74), inclusive

        DESCRIPTION:
            Initializes quality information

        DEBUG:
            Should assign only G's. It does.

        FUTURE: 
            Perhaps add some randomness to this later or distribution to this later.
        """
        qual = ""
        for i in range(self.readlen):
            qual += 'G'
        self.qual = qual




class CHROMOSOME:
    def __init__(self, Chromosome = None, Sequence = None):
        """
        ARGS:
            Chromosome = chromosome string
            Sequence   = chromosome sequence

        RETURN:
            NONE : Initializes CHROMOSOME

        DESCRIPTION:
            Gets chromosome sequence so that we can use it to get the exon
            sequences.

        DEBUG:
            For a shortened whole genome fasta file, it correctly reads it in.
            See read_genome() for debugging code and futher details.

        FUTURE: 
        """
        self.chrm= None         # str, Chromosome
        self.seq     = []           # str, sequence

        if(Chromosome is not None):
            self.chrm = Chromosome
        else:
            exit_with_error("ERROR! Chromosome is not specified!\n")
        if(Sequence is not None):
            self.seq = Sequence
        else:
            exit_with_error("ERROR! Sequence is not specified!\n")




class GENE:
    def __init__(self, GtfEntry):
        """
        ARGS:
            GtfEntry = a single GTF_ENTRY class element

        RETURN:
            NONE : Initializes GENE

        DESCRIPTION:

        DEBUG:

        FUTURE: 
        """
        self.chrm   = None        # str, Chromosome
        self.start  = None        # int, Start position on chrm
        self.stop   = None        # int, End position on chrm
        self.strand = None        # str, '+' (forward) or '-' (reverse)
        self.geneID = None        # str, nominal geneID, but may belong to multiple genes
        self.geneName= None        # str, common gene name
        self.transList    = []    # transcript names that are part of this gene
        self.transIdxList = []    # list, use to quickly map to AllTransList
        transIdx    = 0
        
        # type check
        if(not isinstance(GtfEntry, GTF_ENTRY)):
            exit_with_error("ERROR! GtfEntry is not of class type GTF_ENTRY\n")

        if(GtfEntry.chrm is not None):
            self.chrm = GtfEntry.chrm
        else:
            exit_with_error("ERROR! GtfEntry.chrm is None\n")
        if(GtfEntry.start is not None):
            self.start = int(GtfEntry.start)
        else:
            exit_with_error("ERROR! GtfEntry.start is None\n")
        if(GtfEntry.stop is not None):
            self.stop = int(GtfEntry.stop)
        else:
            exit_with_error("ERROR! GtfEntry.stop is None\n")
        if(GtfEntry.strand is not None):
            self.strand = GtfEntry.strand
        else:
            exit_with_error("ERROR! GtfEntry.strand is None\n")
        if(GtfEntry.geneID is not None):
            self.geneID = GtfEntry.geneID
        else:
            exit_with_error("ERROR! GtfEntry.geneID is None\n")
        if(GtfEntry.geneName is not None):
            self.geneName = GtfEntry.geneName
        else:
            exit_with_error("ERROR! GtfEntry.geneName is None\n")

        # Get all the exons and their indices associated with this transcript
        # This is _HIGHLY_ inefficient. Scales as O(N**2)
        #for trans in AllTransList:
        #    if(trans.geneID == self.geneID):
        #        self.transList.append(trans.transID)
        #        self.transIdxList.append(transIdx)
        #    transIdx = transIdx + 1


class TRANSCRIPT:
    def __init__(self, GtfEntry):
        """
        ARGS:
            GtfEntry = a single GTF_ENTRY class element

        RETURN:
            NONE : Initializes TRANSCRIPT

        DESCRIPTION:

        FUTURE: 
        """
        self.seq    = None        # str, sequence
        self.chrm   = None        # str, Chromosome
        self.start  = None        # int, Start position on chrm
        self.stop   = None        # int, End position on chrm
        self.strand = None        # str, '+' (forward) or '-' (reverse)
        self.geneID = None        # str, nominal geneID, but may belong to multiple genes
        self.transID = None       # str, nominal transcript ID, may belong to mult. trans.
        self.transNum= None       # str, only number portion of transID for sorting by transcript num
        self.exonList= []         # exon names that are part of this transcript
        self.exonIdxList = []     # list, use to quickly map to AllExonList
        exonIdx     = 0

        self.copy = 1               # copy number of this particular transcript in the transcriptome
                                    # It is set to 1 at the moment.

        # type check
        if(not isinstance(GtfEntry, GTF_ENTRY)):
            exit_with_error("ERROR! GtfEntry is not of class type GTF_ENTRY\n")

        if(GtfEntry.chrm is not None):
            self.chrm = GtfEntry.chrm
        else:
            exit_with_error("ERROR! GtfEntry.chrm is None\n")
        if(GtfEntry.start is not None):
            self.start = int(GtfEntry.start)
        else:
            exit_with_error("ERROR! GtfEntry.start is None\n")
        if(GtfEntry.stop is not None):
            self.stop = int(GtfEntry.stop)
        else:
            exit_with_error("ERROR! GtfEntry.stop is None\n")
        if(GtfEntry.strand is not None):
            self.strand = GtfEntry.strand
        else:
            exit_with_error("ERROR! GtfEntry.strand is None\n")
        if(GtfEntry.geneID is not None):
            self.geneID = GtfEntry.geneID
        else:
            exit_with_error("ERROR! GtfEntry.geneID is None\n")
        if(GtfEntry.transID is not None):
            self.transID = GtfEntry.transID
            self.transNum = self.transID[4:]
        else:
            exit_with_error("ERROR! GtfEntry.transcriptID is None\n")


class INSERT:
    def __init__(self, Transcript = None, StartWrtTrans = 0, StopWrtTrans = 0,
                 ReadLength = 0, ExonList = None):
        """
        ARGS:
            Transcript = a TRANSCRIPT instance
            start = the start wrt to the transcript
            stop  = the stop wrt to the transcript
            sequence = the sequence of the insert

        RETURN:
            NONE : Initializes INSERT

        DESCRIPTION:
    
        DEBUG:
            1. Created a crappy shell script, extract_read.sh
               Looking at the reads in the file generated by 
               create_fastq_file() -> create_insert(), you'll see metadata look
               something like :

               @Read_num:0:trans:"ENST00000618181":start:935894:exons:"ENSE00002703998":935772:935896:"ENSE00002686739":939040:939129

               To check to see if the metadata correctly describes the read 
               (and by extension the read itself is likely correct), run :
    
               grep `bash src/simulate_fastq_data/extract_read.sh 935894 935772 \
                     935896 939040 939129` tmp.fq

               CONCLUSION : single end reads of transcripts/inserts on the '+'
                            strand in the sense direction work.

            2. Tested ENST00000488147 for the reverse strand of DNA
               --> Spot checked on ensembl reads with 1 or 2 exons. 
                   Appears to work, but there is a bug (see below)

        FUTURE: 
            1. Check that the Read 2 is correctly handled.
            2. BUG!!! The start and stop values put in the metadata when using
                      a transcript that is on the reverse strand is wrong. This 
                      needs fixed!
            3. Possible BUG : Added conditional to exclude duplicate exons 
               when exons are completely spanned. If I thought about the equalities
               harder, the additional conditional would likely be unnecessary
        """
        self.seq    = None        # str, sequence
        self.chrm   = None        # str, Chromosome
        self.start  = -1          # int, start position w/r/t the chromosome
        self.stop   = -1          # int, start position w/r/t the chromosome
        self.strand = None        # str, '+' (forward) or '-' (reverse)
        self.geneID = None        # str, nominal geneID, but may belong to multiple genes
        self.transID = None       # str, nominal transID, may belong to mult. trans.
        self.transNum= None       # str, only number portion of transID for sorting by transcript num
        if(StopWrtTrans - StartWrtTrans < ReadLength):
            exit_with_error("ERROR!! StopWrtTrans - StartWrtTrans ({}) < "
                            "ReadLength ({})\n".format(StopWrtTrans - StartWrtTrans,
                            ReadLength))
        # type check
        if(not isinstance(Transcript, TRANSCRIPT)):
            exit_with_error("ERROR! Transcript is not of class type TRANSCRIPT 2\n")
        else:
            self.transcript = Transcript

        if(Transcript.chrm is not None):
            self.chrm = Transcript.chrm
        else:
            exit_with_error("ERROR! Transcript.chrm is None\n")
        if(Transcript.strand is not None):
            self.strand = Transcript.strand
        else:
            exit_with_error("ERROR! Transcript.strand is None\n")
        if(Transcript.geneID is not None):
            self.geneID = Transcript.geneID
        else:
            exit_with_error("ERROR! Transcript.geneID is None\n")
        if(Transcript.transID is not None):
            self.transID = Transcript.transID
            self.transNum = self.transID[4:]
        else:
            exit_with_error("ERROR! Transcript.transcriptID is None\n")

        # get the sequence of the insert
        self.seq = Transcript.seq[StartWrtTrans:StopWrtTrans]

        # get the start and stop position of insert and associated reads (1 and 2)
        # relative to the chromosome. May discard read 2 if single end
        transPos= 0
        insertStart = 0
        insertStop  = 0
        exonL = [ExonList[idx] for idx in self.transcript.exonIdxList]
        if(Transcript.strand == '+'):     # Forward DNA strand
            exonL = sorted(exonL, key=operator.attrgetter('start'))
        elif(Transcript.strand == '-'):   # Reverse DNA strand
            exonL = sorted(exonL, key=operator.attrgetter('start'), reverse=True)
        else:
            exit_with_error("ERROR!! invalid value for Transcript.strand : {}".format(
                            Transcript.strand))
        exonSpanL= [] # List of exons spanned by insert
        # Read 1
        r1StartWrtTrans = StartWrtTrans 
        r1StopWrtTrans  = StartWrtTrans + ReadLength - 1 
        r1ExonSpanL = []
        r1Start = None
        r1Stop  = None
        # Read 2
        r2StartWrtTrans = StopWrtTrans
        r2StopWrtTrans  = StopWrtTrans - ReadLength + 1 
        r2ExonSpanL = []
        r2Start = None
        r2Stop  = None

        # Get all exons spanned by insert and reads 1 & 2.
        for exon in exonL:
            exonStart = transPos                    ## In transcript coords
            exonStop  = transPos + len(exon.seq)    ## In transcript coords
            ##### Insert #####
            ## Insert starts in exon
            if(StartWrtTrans >= exonStart and StartWrtTrans <= exonStop):
                exonSpanL.append(exon)
                insertStart = exon.start + (StartWrtTrans - exonStart)
            ## Insert spans exon
            if(StartWrtTrans <= exonStart and StopWrtTrans >= exonStop):
                exonSpanL.append(exon)
            ## Insert ends in exon
            if(StopWrtTrans >= exonStart and StopWrtTrans <= exonStop):
                exonSpanL.append(exon)
                insertStop = exon.start + (StopWrtTrans - exonStart)

            ##### Read 1 #####
            ## Insert starts in exon
            if(r1StartWrtTrans >= exonStart and r1StartWrtTrans <= exonStop):
                r1ExonSpanL.append(exon)
                r1Start = exon.start + (r1StartWrtTrans - exonStart)
            ## Insert spans exon
            if(r1StartWrtTrans <= exonStart and r1StopWrtTrans >= exonStop):
                if(exon not in r1ExonSpanL):
                    r1ExonSpanL.append(exon)
            ## Insert ends in exon
            if(r1StopWrtTrans >= exonStart and r1StopWrtTrans <= exonStop):
                # Prevent duplicates
                if(exon not in r1ExonSpanL):
                    r1ExonSpanL.append(exon)
                r1Stop = exon.start + (r1StopWrtTrans - exonStart)

            ##### Read 2 #####
            ## Insert starts in exon
            if(r2StopWrtTrans >= exonStart and r2StopWrtTrans <= exonStop):
                r2ExonSpanL.append(exon)
                r2Stop = exon.start + (r2StopWrtTrans - exonStart)
            ## Insert spans exon
            if(r2StopWrtTrans <= exonStart and r2StartWrtTrans >= exonStop):
                if(exon not in r2ExonSpanL):
                    r2ExonSpanL.append(exon)
            ## Insert ends in exon
            if(r2StartWrtTrans >= exonStart and r2StartWrtTrans <= exonStop):
                if(exon not in r2ExonSpanL):
                    r2ExonSpanL.append(exon)
                r2Start = exon.start + (r2StartWrtTrans - exonStart)

            transPos = transPos + len(exon.seq)

        ## Error Check ##
        if(len(r1ExonSpanL) == 0):
            exit_with_error("ERROR! Read 1 does _not_ span any exons!\n")
        if(len(r2ExonSpanL) == 0):
            exit_with_error("ERROR! Read 2 does _not_ span any exons!\n")
        if(r1Start is None or r1Stop is None or r2Start is None or r2Stop is None):
            exit_with_error("ERROR!! Invalid trans={}, r1Start,r1Stop,r2Start,r2Stop = "
                            "{},{},{},{}\n".format(Transcript.transID,r1Start,r1Stop,
                            r2Start,r2Stop))
        # Check for duplicate exons (bad!)
        if(len(r1ExonSpanL) != len(set(r1ExonSpanL))):
            exit_with_error("ERROR!! {} : len(r1ExonSpanL) {} != "
                            "len(set(r1ExonSpanL)) {}\n".format(Transcript.transID,
                            len(r1ExonSpanL),len(set(r1ExonSpanL))))
        if(len(r2ExonSpanL) != len(set(r2ExonSpanL))):
            exit_with_error("ERROR!! {} : len(r2ExonSpanL) {} != "
                            "len(set(r2ExonSpanL)) {}\n".format(Transcript.transID,
                            len(r2ExonSpanL),len(set(r2ExonSpanL))))
            
        self.start = insertStart    ## In chromosome coords
        self.stop  = insertStop     ## In chromosome coords
        self.exonL = exonSpanL      ## Exons spanned

        self.r1Start = r1Start      ## In chromosome coords
        self.r1Stop  = r1Stop       ## In chromosome coords
        self.r1ExonL = r1ExonSpanL  ## Exons spanned
        
        self.r2Start = r2Start      ## In chromosome coords
        self.r2Stop  = r2Stop       ## In chromosome coords
        self.r2ExonL = r2ExonSpanL  ## Exons spanned
        


class EXON:
    def __init__(self, GtfEntry):
        """
        ARGS:
            GtfEntry = a single GTF_ENTRY class element

        RETURN:
            NONE : Initializes EXON 

        DESCRIPTION:

        DEBUG:

        FUTURE: 
        """
        self.seq    = None        # str, sequence
        self.chrm   = None        # str, Chromosome
        self.start  = None        # int, Start position on chrm
        self.stop   = None        # int, End position on chrm
        self.strand = None        # str, '+' (forward) or '-' (reverse)
        self.exonNum= None        # int, '1' indexed exon position w/r/t other exons
        self.exonID = None        # str, exon ID ... starts with ENSE
        self.geneID = None        # str, nominal geneID, but may belong to multiple genes
        self.transID = None  # str, nominal transcript ID, may belong to mult. trans.

        # type check
        if(not isinstance(GtfEntry, GTF_ENTRY)):
            exit_with_error("ERROR! GtfEntry is not of class type GTF_ENTRY\n")
            
        if(GtfEntry.chrm is not None):
            self.chrm = GtfEntry.chrm
        else:
            exit_with_error("ERROR! GtfEntry.chrm is None\n")
        if(GtfEntry.start is not None):
            self.start = int(GtfEntry.start)
        else:
            exit_with_error("ERROR! GtfEntry.start is None\n")
        if(GtfEntry.stop is not None):
            self.stop = int(GtfEntry.stop)
        else:
            exit_with_error("ERROR! GtfEntry.stop is None\n")
        if(GtfEntry.strand is not None):
            self.strand = GtfEntry.strand
        else:
            exit_with_error("ERROR! GtfEntry.strand is None\n")
        if(GtfEntry.exonNum is not None):
            self.exonNum = GtfEntry.exonNum
        else:
            exit_with_error("ERROR! GtfEntry.exonNum is None\n")
        if(GtfEntry.exonID is not None):
            self.exonID = GtfEntry.exonID
        else:
            exit_with_error("ERROR! GtfEntry.exonID is None\n")
        if(GtfEntry.geneID is not None):
            self.geneID = GtfEntry.geneID
        else:
            exit_with_error("ERROR! GtfEntry.geneID is None\n")
        if(GtfEntry.transID is not None):
            self.transID = GtfEntry.transID
        else:
            exit_with_error("ERROR! GtfEntry.transcriptID is None\n")
            


def reverse_complement(seq):
    """
    ARGS:
        seq : sequence with _only_ A, T, C or G (case sensitive)

    RETURN:
        rcSeq : reverse complement of sequenced passed to it.

    DESCRIPTION:

    DEBUG: 
        Compared several sequences.  Is working.

    FUTURE: 
    """ 
    rcSeq = ""            # Reverse Complement sequence
    # Complement
    for char in seq:
        if(char == 'A' ):
            rcSeq += 'T' 
            continue
        if(char == 'T' ):
            rcSeq += 'A' 
            continue
        if(char == 'G' ):
            rcSeq += 'C' 
            continue
        if(char == 'C' ):
            rcSeq += 'G' 
            continue
        if(char == 'N' ):
            rcSeq += 'N' 
            continue

        if(char not in "ATCGN"):
            exit_with_error("ERROR! char %s is not a valid sequencing character!\n"%(char))
    # Revese
    rcSeq = rcSeq[::-1]
    return rcSeq

            



