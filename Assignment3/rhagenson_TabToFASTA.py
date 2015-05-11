#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Assignment #: 3
#
# Due date: Tuesday, February 10th 
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
# This program takes a tab-delimited file, assumming each line in the input
# file is a single record and contains the following four fields: GI number,
# accession.version number, gene name and symbol (the latter in parentheses),
# and sequence, all separated by a single tab character each. It converts this
# input into a multisequence GenBank FASTA-formatted file.

import getopt
import sys

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
    listRecord = []         # List for spliting each record (ie line)
    
    # Directions for usage of this script as a series of print statements
    def usage():
        print('Usage: TabToFASTA.py [-i FILE] [-o FILE] [-l width]')
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

    # Read in the first line from INPUT
    line = INPUT.readline() 
    while (line):
        # Remove newline character
        line = line.rstrip('\n')
        
        # Break each record into its parts
        listRecord = line.split('\t')
        
        # Print each record to OUTPUT, using formatSeq() to print the sequence
        # Note: There should be a single space before listRecord[2]
        print('>gi', listRecord[0], 'ref', listRecord[1], \
              ' '+listRecord[2], file=OUTPUT, sep='|')
        print(formatSeq(listRecord[3], CUTOFF), file=OUTPUT)
        
        # Read in next line from INPUT
        line = INPUT.readline()

if (__name__ == '__main__'):
    main()
