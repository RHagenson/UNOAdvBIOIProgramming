#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Assignment #: Final Project
#
# Due date: Friday, May 8th
#
# Honor Pledge: On my honor as a student of the University of Nebraska at
#               Omaha, I have neither given nor received unauthorized help on
#               this programming assignment.
#
#    NAME: Ryan Hagenson
#    NUID: 972
#    EMAIL: rhagenso@nebrwesleyan.edu
#    Partners: None
#
# This program takes a human gene symbol as input and outputs a multisequence
#   FASTA file containing only orthologs of this gene.
# The input is used to find a matching FASTA record from a given multirecord
#   FASTA file (preferably a curated "NM_" record).
# This FASTA record is then BLASTed against a database (default: nt) outputing
#   matching records that are below a specified e-value to a temporary file.
# The final step is iterating through this temporary file and copying only the
#   orthologs to the originally specificed output name.

import os
import getopt
import sys
import re
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

# Function name: formatSeq
# Parameters: sequence/string, width/integer 
# Return value: string formated to the given width
# Description: formatSeq() takes two inputs, a sequence and a width, and 
# outputs the formatted sequence set to this line width.
def formatSeq(sequence, width):
    listSeq = []  # List for building sequence lines
    
    # Sequences under the length 'width' need no loop iteration
    if (len(sequence) < width):
        listSeq.append(sequence[0:])
    # Otherwise the sequence is chunked by width via a loop
    else:
        # Insert newline character after every character at width
        i = 0
        while ((i + width) < len(sequence)):
            listSeq.append(sequence[i:i + width]+'\n')
            i = i + width
        
        #Account for final line in sequence record
        if (i < len(sequence)):
            listSeq.append(sequence[i:])
        
    # Return the formatted sequence as a string
    return ''.join(listSeq)



def main():
##### USER-DEFINED FILES AND VALUES #####
    # The multirecord FASTA file that we extract our record from
    FASTA = open("/home/rhagenson/storage/human.rna.fna", 'r')

    messages = False                         # The verbosity boolean
    INPUT = ""                               # Input string 
    OUTPUT = sys.stdout                      # Output file handle; 
                                             # Default: stdout
    Seq = ""                                 # String to contain sequence
    FASTARecHead = ""                        # Allow printing of header from
                                             # the multirecord FASTA file
    FASTARecSeq = []                         # List to build the sequence from
                                             # the multirecord FASTA file
    CUTOFF = 70                              # The sequence line width, 
                                             # Default: 70 characters wide
    includeSeq = False                       # The boolean that controls the 
                                             # printing of sequence
    E_Value = 0.001                          # The e-value cutoff 
                                             # for the BLASTN query
    MaxSeqs = 100                            # The maximum number of records to
                                             # return from BLASTN 
    
    # Regex to match fastaHeader from BLASTN results
    fastaHeaderMatch = None
    fastaHeaderPattern = re.compile('gi\|[0-9]+\|ref|.+$')
    
    # The BLAST query temporary file
    FastaName = "/home/rhagenson/storage/FPFasta.%.5s.fna" % (os.getpid())
    FASTAFile = open(FastaName, 'w')
    
    # The temporary BLAST results file
    TempName = "/tmp/FPOrtholog.%.5s.xml" % (os.getpid())
    TEMP_OUT = open(TempName, 'w')
   

    # Directions for usage of this script as a series of print statements
    def usage():
        print('Usage: TabToFASTA.py [I STRING] [-O FILE] [-L width]\n')
        print('\t-V:\tverbose flag, prints messages about what is\n' , \
              '\t\tcurrently being done by the program')
        print('\t-I:\tinput file (tab-delimited); STDIN if not used.')
        print('\t-O:\toutput file (FASTA-formatted); STDOUT if not used.')
        print('\t-S:\tif used, the sequence is included in the output\n', \
              '\t\t(defaut is to not include the sequence)')
        print('\t-L:\tsequence width (default is 70 characters wide)')
        print('\t-E:\tsets the e-value cutoff for the BLASTN results')
        return
    
    # Enables command-line arguments via getopt and sys packages
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'VI:O:SL:')
    except getopt.GetoptError as err:
        # Redirect STDERR to STDOUT (ensures screen display)
        sys.stdout = sys.stderr
        # Print help information
        print(str(err))
        # Print usage information; usage() is the name of a
        # function (declared elsewhere) that displays usage
        # information (a series of print statements)
        usage()
        # Exit
        sys.exit(2)

    # The activity of each command line option
    for (opt, arg) in opts:
        if (opt == '-V'):
            messages = True
        if (opt == '-I'):
            INPUT = arg
        if (opt == '-O'):
            OUTPUT = open(arg, 'w')
        if (opt == '-L'):
            CUTOFF = int(arg)
        if (opt == '-S'):
            includeSeq = True
        if (opt == '-E'):
            E_Value = int(arg)


################ BELOW THIS POINT DO NOT CHANGE ###############
###############################################################
    
    FirstSeq = True  # After the first record is found switched off 
    FoundIt = False  # If a preferred curated (NM_) record is found stops
                     # looking through records and uses this record
    
    # Look for official symbol within multirecord FASTA file    
    # Check for NM_ record or take first record available
    line = FASTA.readline()
    while (line):
        if (FirstSeq) and (messages):
            print("Looking within the multirecord FASTA file", \
                  "for the provided gene symbol")
        
        # If a preferred curated (NM_) record use it
        if (INPUT in line) and ("NM_" in line):
            if (messages):
                print("A curated record has been found")
            
            # Rewrite the file at FastaName if a curated record is found 
            FASTAFile.close()
            BestFile = open(FastaName, 'w')
           
            # Empty FASTARecSeq and FASTARecHead 
            # prior to building the curated record
            FASTARecSeq = []
            FASTARecHead = ""

            # Write the header to file
            line = line.rstrip('\n')
            print(line, file=BestFile)
            
            # Capture the FASTA header for later printing
            FASTARecHead = line

            # Indicate a curated record has been found
            FoundIt = True

            # Read in and write the sequence lines to file
            line = FASTA.readline()
            while ('>gi' not in line) and (line != None):
                line = line.rstrip('\n')
                # Build the FASTA sequence for later printing
                FASTARecSeq.append(line)

                # Build the file for BLASTing
                print(line, file=BestFile)
                line = FASTA.readline()

        # End the search if a curated record is found
        if (FoundIt):
            BestFile.close()
            break 

        # If not a preferred record, build the record, but continue looking        
        if (INPUT in line) and (FirstSeq):
            if (messages):
                print("The first non-curated record has been found")

            # Indicate the first record has been found
            FirstSeq = False

            # Write the header to file
            line = line.rstrip('\n')
            print(line, file=FASTAFile)
            
            # Capture the FASTA header for later printing
            FASTARecHead = line

            # Read in and write the sequence lines to file
            line = FASTA.readline()
            while ('>gi' not in line) and (line != None):
                line = line.rstrip('\n')
                # Build the FASTA sequence for later printing
                FASTARecSeq.append(line)

                # Build the file for BLASTing
                print(line, file=FASTAFile)
                line = FASTA.readline()
        
        # Read in the next line if no conditions are matched
        line = FASTA.readline()

    if (messages):
        print("Continuing with the analysis by conducting a BLASTN", \
              "against the provided database with the choosen record")
        print("File location of choosen record:", FastaName)

    # Close the FASTA record prior to BLASTing
    FASTAFile.close()

    # Define the BLASTN
    blastnCmd = NcbiblastnCommandline( \
                query=str(FastaName), \
                db="/course-nfs/mpauley/shared/blastdb/nt", \
                max_target_seqs=MaxSeqs, \
                evalue=E_Value, \
                outfmt=5, \
                out=str(TempName))
    
    if (messages):
        print("The BLAST command that is running is:\n", blastnCmd)

    # Run the above command, BLASTing the record against the choosen database
    stdout, stderr = blastnCmd()

    # Close the BLAST results temp file prior to parsing the results
    TEMP_OUT.close()

    # Prior to parsing the BLAST results, print the query record to user OUTPUT
    print(FASTARecHead, file=OUTPUT)
    print(formatSeq(''.join(FASTARecSeq), CUTOFF), file=OUTPUT)
    
    # Parse results printing only those that are not human
    BLAST_RESULTS = open(TempName, 'r')
    blastRecords = NCBIXML.parse(BLAST_RESULTS)
    for blastRecord in blastRecords:
        for alignment in blastRecord.alignments:
            for hsp in alignment.hsps:
                if ("Homo sapiens" not in alignment.title):
                    # Match regex to extract FASTA header
                    fastaHeaderMatch = re.search(fastaHeaderPattern, \
                                                 alignment.title)
                    
                    # Always print the FASTA header to user OUTPUT
                    print(">", fastaHeaderMatch.group(), sep='', file=OUTPUT)  
                    
                    # If sequence flag is on, print the sequence for hsp,
                    # formatted according to the specified length, 
                    # to user OUTPUT
                    if (includeSeq):
                        print(formatSeq(hsp.sbjct, CUTOFF), file=OUTPUT)
    

    # Ensure all file handles are closed
    FASTAFile.close()
    BestFile.close()
    TEMP_OUT.close()
    OUTPUT.close()
    # Clean up BLAST results temp file
    os.remove(TempName)

if (__name__ == '__main__'):
    main()
