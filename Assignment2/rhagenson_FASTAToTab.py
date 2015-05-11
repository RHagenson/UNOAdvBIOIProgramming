#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Assignment #: 2
#
# Due date: Thursday, February 5th
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
# This program converts a multisequence human reference sequence GenBank FASTA
# file into a tab-delimited format. For a given FASTA record, output is a
# single tab-delimited line that includes GI number, the accession.version
# number, the gene name and symbol, and (optionally) the sequence itself.
# The organism name and preceeding text along with 
# the sequence type are all dropped.

import sys
import getopt

def main():
    # Declaration of variables
    INPUT = sys.stdin       # Input file handle; Default: stdin
    OUTPUT = sys.stdout     # Output file handle; Default: stdout
    includeSeq = False      # Boolean for including sequence; Default: False
    listHeader = []         # List for spliting header
    listSeq = []            # List for building sequence by line

    # Enables command-line options via getopt and sys packages
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:o:s')
    except getaopt.GetoptError as err:
        # Redirect STDERR to STDOUT (insures screen display)
        sys.stdout = sys.stderr
        # Print help information
        print(str(err))
        # Print usage information
        # usage()
        # Exit
        sys.exit(2)
    
    # Defines the action of each command-line option
    for (opt, arg) in opts:
        if (opt == '-i'):
            INPUT = open(arg, 'r')
        if (opt == '-o'):
            OUTPUT = open(arg, 'w')
        if (opt == '-s'):
            includeSeq = True

    # Reads INPUT line-by-line, writing to OUTPUT 
    line = INPUT.readline()
    while (line):
        if ('>' in line):
            # Prints previous record then clears variables for current record.
            # Placement here allows for the complete Seq to be built.
            if (listHeader):
                print(listHeader[1], listHeader[3], listHeader[4], \
                      ''.join(listSeq), sep='\t', file=OUTPUT)
                listHeader = []
                listSeq = []
            # Ensure one line output per record
            line = line.rstrip('\n')
            # Split by '|' for selective printing and description modification
            listHeader = line.split('|')
            # Remove organism name, trailing whitespace, and preceeding text
            listHeader[4] = \
                    listHeader[4][listHeader[4].find('Homo sapiens')+13:]
            # Remove sequence type by locating final comma and splicing
            # everything before it since sequence type 
            # always follows the final comma
            listHeader[4] = listHeader[4][:listHeader[4].rfind(',')]
        # For lines containing sequence, the newline is removed, 
        # and each line is added to listSeq
        if (includeSeq and '>' not in line):
            line = line.rstrip('\n')
            listSeq.append(line)
        # Read in next line from INPUT
        line = INPUT.readline()
   
    # Print final record
    print(listHeader[1], listHeader[3], listHeader[4], \
          ''.join(listSeq), sep='\t', file=OUTPUT)
    
    # Close files to release resources
    INPUT.close()
    OUTPUT.close()

if (__name__ == '__main__'):
    main()
