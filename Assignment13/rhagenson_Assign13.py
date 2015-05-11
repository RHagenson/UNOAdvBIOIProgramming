#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Assignment #: 13
#
# Due date: April 16th
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
# This program takes an XML blastn results file, for each record in the query
# database the FASTA sequence of the corresponding genomic sequence from the
# genomic hit is displayed. Information for the FASTA header is taken from a
# local database humanGenomic. An additional piece of information is included
# in the header: the smallest value of sbjct_start and largest value of 
# sbjct_end from the HSP class are appended to the header in brackets. These
# values are also used to trim the target sequence so it matches the query
# sequence.

import psycopg2
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast import Record
import re

# formatSeq() takes two inputs, a sequence and a width, and outputs the
# formatted sequence set to this line width.
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
    ### Steps
    # 1. parse BLAST_RESULTS for the accession number
    # 2. build FASTA header by SQL execute with accession number
    # 3. use sbjct_start and sbjct_end to trim target sequence to match query
    # 4. print header using SQL, then formatted sequence using formatSeq()
    
    ### User-defined variables, point this path to your XML input
    BLAST_RESULTS = \
    open('/home/mpauley/public/BIOI3500/geneListShort.vs.humanGenomic.blastXML', \
         'r')

    WIDTH = 70  # The desired width to format the sequence lines to

    ### Program varibales, DO NOT CHANGE BELOW THIS POINT
    accessionMatch = None                                   # Check for match
    accessionPattern = re.compile('[A-Z]{1,2}\_?\d+\.\d+')  # Regex to match
                                                            # accession number
    

    # Open connnection to SQL database
    # Details about the SQL database and login credentials have been removed #
    conn = psycopg2.connect("dbname="" user="" password="" host=""")
    
    # Create cursor for interaction with SQL database
    cursor = conn.cursor()

    # Parse the XML results file
    blastRecords = NCBIXML.parse(BLAST_RESULTS)
    for blastRecord in blastRecords:
        smallest_start = float('inf')  # Local variable for finding 
                                       # smallest sbjct_start
        largest_end = 0                # Local varibale for finding 
                                       # largest sbjct_end
        for description in blastRecord.descriptions:
            accessionMatch = re.search(accessionPattern, description.title)
            # Gather data from SQL database based on gene symbol
            cursor.execute("select gi, accno, description, seq \
                           from humanGenomic \
                           where accno='%.15s'" % (accessionMatch.group()))
            SQL_line = cursor.fetchone()

            # Iterate over alignments to find the record of interest
            for alignment in blastRecord.alignments:
                if (alignment.title == description.title):
                    for hsp in alignment.hsps:
                        if (hsp.sbjct_start < smallest_start):
                            smallest_start = hsp.sbjct_start
                        if (hsp.sbjct_end > largest_end):
                            largest_end = hsp.sbjct_end
        
            # Print the formatted output
                print('>gi', SQL_line[0], SQL_line[1],' ' + SQL_line[2] + ' ' + \
                     '[%.10s - %.10s]' % (smallest_start, largest_end), \
                      sep='|')
                if (SQL_line[3][smallest_start:largest_end+1]):
                    print(formatSeq(SQL_line[3][smallest_start:largest_end+1], \
                          WIDTH))

    # Store the first result of the above execute as a tuple
    symbolInfo = cursor.fetchone()
    
    # Close statments
    BLAST_RESULTS.close()
    cursor.close()
if (__name__ == '__main__'):
    main()
