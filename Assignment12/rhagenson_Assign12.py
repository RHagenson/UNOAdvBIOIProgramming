#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Assignment #: 12
#
# Due date: April, 7nd
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
# This program computes and displays the frequency distribution of the amino
# acids from all the sequences in two given files. The first being of type
# fasta the second being of type swiss (uniprot). It computes the frequency for
# each file separately. It outputs the amino acids symbol, fasta file frequency
# of that amino acid, swiss file frequency of that amino acid, then the
# difference between fasta file and swiss file. The last two lines are the
# number of records in the fasta file followed by the number of records in the 
# swiss file.

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def main():
    ### USER-SPECIFIED VARIABLES  ###
    # FASTA file handle
    FASTA = open('/home/rhagenson/storage/human.protein.faa', 'r')
    # UniProt/Swiss file handle
    SWISS = open('/home/rhagenson/storage/uniprot_sprot_human.dat', 'r')

    ### DO NOT CHANGE BEYOND THIS POINT ###
    # FASTA variables
    fastaDict = {}          # Dictionary to hold amino acid:count for FASTA
    fastaTotalAA = 0        # Used to calculate percentage
    fastaNumRecords = 0     # Used to count the number records

    # SWISS variables
    swissDict = {}          # Dictionary to hold amino acid:count for SWISS
    swissTotalAA = 0        # Used to calculate percentage
    swissNumRecords = 0     # Used to count the number records
    
    # Compute frequencies for FASTA file
    for record in SeqIO.parse(FASTA, "fasta"):
        fastaTotalAA += len(record.seq)
        fastaNumRecords += 1

        for AA in record.seq:
            # For each new amino acid create a dictionary entry
            if (AA not in fastaDict):
                fastaDict[AA] = 1
            # For each existing amino acid interate count by 1
            else:
                fastaDict[AA] += 1

    # Compute frequencies for SWISS file
    for record in SeqIO.parse(SWISS, "swiss"):
        swissTotalAA += len(record.seq)
        swissNumRecords += 1
    
        for AA in record.seq:
            # For each new amino acid create a dictionary entry
            if (AA not in swissDict):
                swissDict[AA] = 1
            # For each existing amino acid iterate count by 1
            else:
                swissDict[AA] += 1

    # print header
    print('Results:')

    # Print processed results
    for AA in sorted(fastaDict):
        fastaPercent = round(((fastaDict[AA]/fastaTotalAA)*100), 2)
        swissPercent = round(((swissDict[AA]/swissTotalAA)*100), 2)
        difference = fastaPercent - swissPercent
        
        print(AA, \
              '{0:.2f}'.format(fastaPercent, 8) + '%', \
              '{0:.2f}'.format(swissPercent, 8) + '%', \
              '{0:.2f}'.format(difference, 8) + '%')

    # print number of records
    print('NCBI: ', fastaNumRecords, ' records')
    print('UniProt ', swissNumRecords, ' records')

if __name__ == '__main__':
    main()
