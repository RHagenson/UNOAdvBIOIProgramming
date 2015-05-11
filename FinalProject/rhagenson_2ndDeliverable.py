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
# This program is a short testing script for working with the NCBI standalone BLAST utility 

import os
import getopt
import sys
from Bio.Blast.Applications import NcbiblastnCommandline

def main():
    ProcessId = os.getpid()
    print(ProcessId)
    Filename = "/tmp/FPOrtholog." + str(ProcessId)
    print(Filename)
    # The name of the temporary BLAST results file
    TEMP_OUT = open(Filename, 'w')  
    
    blastnCmd = NcbiblastnCommandline( \
                query="/home/rhagenson/storage/NM_000492.fna", \
                db="/course-nfs/mpauley/shared/blastdb/nt", \
                max_target_seqs=50, \
                evalue=0.01, \
                outfmt=5, \
                out="ter.t.it")

    stdout, stderr = blastnCmd()



if (__name__ == '__main__'):
    main()
