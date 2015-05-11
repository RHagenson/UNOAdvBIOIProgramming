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
#    Partners: None
#
# This program prompts the user for a DNA sequence (i.e. letters A, T, G, or C)
# then calculates and displays the GC-content of that sequence as a percentage.
# The input sequence is converted to uppercase prior to analysis

def main():
    sequence = ""   # User-input sequence
    GC_number = 0   # Integer for counting number of G/C bases
    baseChrs = 0    # Integer for counting total number of bases

    # Get the sequence from the user; case-insenstive
    sequence = str(input("Enter a sequence: ")).upper()

    # Determine the number of G's and C's, as well as total number of bases
    for i in sequence:
        if ((i == 'G') or (i == 'C')):
            GC_number += 1
            baseChrs += 1
        elif((i == 'A') or (i == 'T')):
            baseChrs += 1

    # Report the GC-percentage back to the user
    print("%%GC content: %.2f%%" % (GC_number / baseChrs * 100.00))

if (__name__ == '__main__'):
    main()
