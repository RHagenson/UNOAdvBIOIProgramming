#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Assignment #: 4
#
# Due date: Thursday, February 19th
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
# This program takes a multi-record GenBank file and converts it into a
# multisequence FASTA-formatted file. Records in the output come right after
# one another, meaning there is no blank line separating them.

import sys
import getopt
import re


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
    INPUT = sys.stdin       # Input file handle; Default: stdin
    OUTPUT = sys.stdout     # Output file handle; Default: stdout
    CUTOFF = 70             # The sequence line width, 
                            # Default: 70 characters wide
    Seq = []                # List for containing sequence of record
    listRecord = [0, 1, 2]  # List for containing the desired parts of the
                            # record. 0 will be replaced by gi-number, 1 will
                            # be replaced by accession.version number, and 2
                            # will be replaced by the locus

    # Directions for usage of this script as a series of print statements
    def usage():
        print('Usage: GenBankToFASTA.py [-i FILE] [-o FILE] [-l width]')
        print('\t-i:\tinput file (tab-delimited); STDIN if not used.')
        print('\t-o:\toutput file (FASTA-formatted); STDOUT if not used.')
        print('\t-l:\tsequence width (defaults to 70 if not specified)')
        return
    
    # Enables command-line arguments via getopt and sys packages
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:o:l:')
    except getopt.GetoptError as err:
        # Redirect STDERR to STDOUT (insures screen display)
        sys.stdout = sys.stderr
        # Print help information
        print(str(err))
        # Print usage information; usage() is the name of a
        # function (declared elsewhere) that displays usage
        # information (a series of print statements)
        usage()
        # Exit
        sys.exit(2)

    # Define the actions of each command-line argument from above
    for (opt, arg) in opts:
        if (opt == '-i'):
            INPUT = open(arg, 'r')
        if (opt == '-o'):
            OUTPUT = open(arg, 'w')
        if (opt == '-l'):
            CUTOFF = int(arg)
    
    # Read in first line of INPUT
    line = INPUT.readline()
    # Continue reading in lines until EOF
    while (line):
        # If the line contains the locus, copy it to listRecord[2]
        if ("DEFINITION" in line):
            line = line.rstrip('\n')
            locusMatch = re.search(r'(?<=DEFINITION\s).*', line)
            listRecord[2] = locusMatch.group(0)

            # Append any extra lines to locus
            line = INPUT.readline()
            while ("ACCESSION" not in line):
                line = line.strip('\n')
                line = line.lstrip(' ')
                listRecord[2] = listRecord[2] + ' ' + line
                line = INPUT.readline()
               
        # If the line contains the GI and Accession numbers copy those to
        # listRecord[0] and listRecord[1], respectively
        if ("VERSION" in line):
            # Pattern match for GI number then replaces listRecord[0] with this
            # string
            GI_Match = re.search(r'(?<=GI:)\d+', line)
            if (GI_Match != None):
                listRecord[0] = GI_Match.group(0)
            
            # Pattern match for accession.version number then replaced
            # listRecord[1] with this string
            accessionMatch = re.search(r'[A-Z]{1,2}_?\d+\.\d+', line)
            if (accessionMatch != None):
                listRecord[1] = str(accessionMatch.group(0))
        
        # If the sequence follows begin copying it by reading in the next line
        if ("ORIGIN" in line):
            line = INPUT.readline()

            # Continue reading lines until end of record, denoted by '//'
            # As sequence lines are read they are added to list Seq
            while ("//" not in line):
                line = line.rstrip('\n')
                # Remove leading number
                line = line[10:]
                # Remove internal whitespace
                line = line.replace(' ', '')
                
                Seq.append(line)
                line = INPUT.readline()
        
        if ("//\n" in line):
            if ('NM_' in listRecord[1] or 
                'NR_' in listRecord[1] or
                'XM_' in listRecord[1] or
                'XR_' in listRecord[1]):
                # Print ref record header to file
                print(">gi", listRecord[0], "ref", listRecord[1], \
                      listRecord[2], sep='|', file=OUTPUT)
            else: 
                # Print gb record header to file
                print(">gi", listRecord[0], "gb", listRecord[1], \
                      listRecord[2], sep='|', file=OUTPUT)
            
            # If record has a sequence, print it to file
            if (Seq):
                print(formatSeq((''.join(Seq)).upper(), CUTOFF), file=OUTPUT)
                
                # Remove content of previous record
                Seq = []

        # Read in next line if none of the above conditions are met, most
        # likely due to another record being in the file
        line = INPUT.readline()

if __name__ == "__main__":
    main()
