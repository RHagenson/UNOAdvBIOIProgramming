#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Assignment #: 8
#
# Due date: Tuesday, March 10th (Extended due to my trip)
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
# This program scans a multi-record FASTA file for records of the human fragile
# X mental retardation 1 gene (using the gene symbol and aliases as
# identifiers). For each record the program displays the accession number, GI,
# number of CGG/AGG repeats and  an appropriate message (“normal”, 
# “near-normal”, “permutation” or “full mutation”) depending on the repeat
# number.

import sys
import getopt
import re

def main(): 
    INPUT = sys.stdin               # Input file handle; Default: stdin
    OUTPUT = sys.stdout             # Output file handle; Default: stdout
    aliases = ['(FMR1)', '(FMRP)', 
               '(FRAXA)', '(POF)', 
               '(POF1)']            # Human fragile X mental retardation 1 gene
                                    # symbol and aliases
    sequence = []                   # Build sequence for whole analysis
    repeatNumber = 0                # Counts the number of CGG or AGG repeats
    message = ''                    # Is replaced with 'normal', 'near-normal',
                                    # 'permutation', or 'full mutation'
                                    # depending on repeat number
    GI_Match = None                 # Default, no match 
    accessionMatch = None           # Default, no match
    repeatMatch = None              # Default, no match


    # Enables command-line arguments via getopt and sys packages
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:o:')
    except getaopt.GetoptError as err:
        # Redirect STDERR to STDOUT (insures screen display)
        sys.stdout = sys.stderr
        # Print help information
        print(str(err))
        # Print usage information
        # usage()
        # Exit
        sys.exit(2)

    # Set behavior of command-line arguments
    for (opt, arg) in opts:
        if (opt == '-i'):
            INPUT = open(arg, 'r')
        if (opt == '-o'):
            OUTPUT = open(arg, 'w')
    
    # Compile the regex patterns to extract GI, accession number, and count
    # repeats
    GI_Pattern = re.compile('(?:>gi\|)(\d+)')
    accessionPattern = re.compile('[A-Z]{1,2}\_?\d+\.\d+')
    repeatPattern = re.compile('(?:GCG)((CGG|AGG)+)(?:CTG)')

    # Begin reading in data from INPUT
    line = INPUT.readline()
    while (line):
        if ('>' in line):
            # With a sequence find repeats and set message
            if (sequence):
                # Analyze sequence for human X fragile repeats
                repeatMatch = re.finditer(repeatPattern, \
                                          ''.join(sequence))
                if (repeatMatch):
                    for repeat in repeatMatch:
                        repeatNumber += int((len(repeat.group(1))/3))

                # Determine appropriate message
                if (200 < repeatNumber):
                    message = 'full mutation'
                if (55 <= repeatNumber <= 200):
                    message = 'permutation'
                if (45 <= repeatNumber <= 54):
                    message = 'near-normal'
                if (5 <= repeatNumber <= 44):
                    message = 'normal'

            # Print previous record to file
            if (message):
                print(GI_Match.group(1), accessionMatch.group(), \
                      repeatNumber, message, sep=', ', file=OUTPUT)

                # Reinitialize variables
                sequence = []
                repeatNumber = 0
                message = ''
                GI_Match = None
                accessionMatch = None
                repeatMatch = None

            for alias in aliases:
                if (alias in line):
                    # Extract both GI and accession number
                    GI_Match = re.search(GI_Pattern, line)
                    accessionMatch = re.search(accessionPattern, line)

        # Build sequence if valid record and in sequence lines
        elif (GI_Match != None) and ('>' not in line):
            line = line.rstrip('\n')
            # Build entire sequence into list
            sequence.append(line)
            
        # Read in next line
        line = INPUT.readline()

if (__name__ == '__main__'):
    main()
