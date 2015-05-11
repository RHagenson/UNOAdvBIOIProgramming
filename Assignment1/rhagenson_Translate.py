#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Assignment #: 1
#
# Due date: Thursday, January 29 
#
# Honor Pledge: On my honor as a student of the University of Nebraska at
#               Omaha, I have neither given nor received unauthorized help on
#               this programming assignment.
#
#    NAME: Ryan Hagenson
#    NUID: 972
#    EMAIL: rhagenso@nebrwesleyan.edu
#    Partners: Dr. William McClung (mcclung@nebrwesleyan.edu), caught an error 
#              in my thinking on line 82, I had forgotten the equals sign 
#              in the while statement so the last protein was sometimes
#              not translated
#              As well, he and I had a conversation
#              about what is happening "under the hood" of my code; ie: 
#              why computer scientist begin indexes at zero, 
#              how the getopt module works to pull in command line input,
#              what STOUT and STERR are and why line 46 is needed,
#              and, lastly, what the syntax ''.join(list) does to the list to
#              make it a string. 
#
# This program takes a user-input DNA sequence via the command-line argument -s
# and translates it into its protein equivalent in a given reading frame -f
# (default reading frame: 1). The program converts the user-input to uppercase
# prior to analysis.

import getopt
import sys

def main():
    sequence = ''     # User-input DNA sequence
    frame = 1         # The reading frame to translate 
                      # defaults to frame 1, if not specified
    translated = []   # List for translated protein codes, joined before print
    translation = {}  # Dictonary for the three-letter DNA to one-letter
                      # protein codes
    
    # Enables command-line arguments via getopt and sys packages
    try:
        opts, args = getopt.getopt(sys.argv[1:], 's:f:')
    except getaopt.GetoptError as err:
        # Redirect STDERR to STDOUT (insures screen display)
        sys.stdout = sys.stderr
        # Print help information
        print(str(err))
        # Print usage information
        # usage()
        # Exit
        sys.exit(2)

    for (opt, arg) in opts:
        if (opt == '-s'):
            sequence = str(arg).upper()
        if (opt == '-f'):
            if (1 <= int(arg) <= 3):
                frame = int(arg)

    translation = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
                   'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                   'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
                   'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                   'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                   'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                   'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                   'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                   'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
                   'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
                   'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
                   'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
                   'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                   'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                   'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                   'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

    # (frame - 1) corrects for index starting at zero
    i = (frame - 1)
    # Runs as long as there are at least three more position in the sequence
    # therefore does not attempt to translate incomplete codons
    while ((i + 3) <= len(sequence)):
        # current_codon becomes the next three letters in sequence
        current_codon = sequence[i:i+3]
        # If current_codon is a valid codon translate it to its protein code
        if (current_codon in translation):
            translated.append(translation[current_codon])
        # Otherwise indicate non-valid codon (ie non-base characters) with 'X'
        else:
            translated.append('X')
        # Shift i to the start of the next codon
        i += 3

    # Stringify translated list via join to empty list for printing
    print(''.join(translated))

if (__name__ == '__main__'):
    main()
