#! /usr/bin/python3
# This program takes an XML blastn input and displays the title of the
# alignment, starting position, ending position, and the expect value for all
# HSPs that have an expect value <= 0.1

from Bio import SeqIO
from Bio.Blast import NCBIXML

def main():
    # Read in XML file of BLASTN results
    INPUT = open("../storage/blast/blastResults/twoRecords.v.humanGPlusT", 'r')
    
    # Parse the XML input file
    blastRecords = NCBIXML.parse(INPUT)
    
    # Iterate over each record for analysis
    for blastRecord in blastRecords:
        for alignment in blastRecord.alignments:
                firstRecord = True      # Boolean to print title only once
                
                # Iterate over high scoring pairs and print if expect <= 0.1
                hspNo = 0
                for hsp in alignment.hsps:
                        if (hsp.expect <= 0.1):
                            if (firstRecord):
                                print(alignment.title)
                                firstRecord = False
                            hspNo += 1
                            
                            # Final print statement for the tab-separated
                            # fields
                            print("HSP " + str(hspNo) + ":\t" + str(hsp.query_start) + \
                                  "\t" + str(hsp.query_end) + "\t" + str(hsp.expect))

    # Close file handle
    INPUT.close()

if (__name__ == '__main__'):
    main()
