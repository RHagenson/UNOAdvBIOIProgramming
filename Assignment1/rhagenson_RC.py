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
# and returns the Watson-Crick complement (i.e. A:T, T:A, G:C, C:G translation)
# as well as the reverse complement. The input sequence is converted to
# uppercase prior to analysis.

def main():
    sequence = ""       # User input sequence
    complement = []     # List for building the complement
    rComplement = []   # List for storing the reverse complement
    translation = {}    # Dictionary for Watson-Crick complements
 
    # Promt user for the sequence
    sequence = str(input("Enter a sequence: ")).upper()

    # The DNA base complement dictionary
    translation = {'A':'T', \
                   'T':'A', \
                   'G':'C', \
                   'C':'G'}
    
    # Interate through sequence to build complement
    for i in sequence:
        # For real base characters
        if (i == 'A' or i == 'T' or i == 'G' or i == 'C'):
            complement.append(translation[i])
        # For any other non-base character
        else:
            complement.append('X')
    
    # Build r_complement via negative step through complement
    rComplement = complement[::-1]

    # Join the lists into their string equivalents
    complement = ''.join(complement)
    rComplement = ''.join(rComplement)

    # Return the complement and reverse complement to the user
    print("Complement:", complement)
    print("Reverse complement:", rComplement)

if (__name__ == '__main__'):
    main()
