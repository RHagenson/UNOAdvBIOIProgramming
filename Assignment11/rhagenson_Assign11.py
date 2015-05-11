#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Assignment #: 11
#
# Due date: April, 2nd
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
# This program reads an input file of gene symbols, with one symbol per line,
# querying the database tables geneinfo, gene2pubmed, and gene2refseq writing 
# the gene ID, taxonomy ID, accession numbers (RNA, protein, and genomic) and 
# PubMed IDs for a gene with a given symbol (where symbol can be either the
# offical symbol or a synonym). Output is one line per gene ID with fields as
# follows separated by a tab (\t) character: gene symbol, gene ID, taxonomy ID,
# RNA accession numbers, protein accession numbers, genomic accession numbers,
# PubMed IDs. Separate multiple entries in a given field with pipes (|). 
# Indicate empty fields with a dash (-). For symbols not in the database, 
# all of the fields except the gene symbol should be dashes. All entries in a 
# field must be unique and in ascending order. Genes in multiple organisms 
# should be displayed in order of ascending taxonomy id.


import sys
import getopt
import psycopg2

def main(): 
    messages = True                 # Outputs messages to user on what is being
                                    # currently processed, set to False to
                                    # disable
    INPUT = sys.stdin               # Input file handle; Default: stdin
    OUTPUT = sys.stdout             # Output file handle; Default: stdout
    
    
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
   
    
    # Open connnection to SQL database
    # Details of the database location and login credentials have been removed #
    conn = psycopg2.connect("dbname="" user="" password="" host=""")
    
    # Create cursor for interaction with SQL database
    cursor = conn.cursor()

    for line in INPUT:
        if (messages):
            print('Now processing results for %s' % (line))
 
        line = line.rstrip('\n')
        
        # Gather data from SQL database based on gene symbol
        cursor.execute("select distinct geneinfo.symbol, geneinfo.gene_id, \
                        geneinfo.tax_id, rna_accession, pro_accession, \
                        gen_accession, pubmed_id from geneinfo, gene2pubmed, \
                        gene2refseq where geneinfo.gene_id = gene2pubmed.gene_id \
                        and geneinfo.gene_id = gene2refseq.gene_id \
                        and geneinfo.symbol = '%s' \
                        order by tax_id asc, geneinfo.gene_id asc;" % (line))
        

        sql_line = cursor.fetchone()  # Fetch the first result of the SQL query
       
        print(sql_line)
        # if symbol returns no results print it with NULL placeholder: '-'
        if (sql_line == None):
            print(line, '-', '-', '-', '-', '-', '-', sep='\t', file=OUTPUT)
            continue

        rnaList = []    # To hold multiple rna_accession prior to print
        proList = []    # To hold multiple pro_accession prior to print
        genList = []    # To hold multiple gen_accession prior to print
        pubList = []    # To hold multiple pubmed_id prior to print

        pastList = []   # To hold the symbol, gene_id, tax_id of past record
        # Store the first record's symbol, gene_id, and tax_id
        i = 0
        while (i < 3):
            pastList.append(sql_line[i])
            i += 1

        while(sql_line):
            if (messages):
                print('Now processing a new SQL record')
                print('Current line is:', sql_line)

            if (sql_line[1] != pastList[1]):
                # When and pastGene_ID do not match print the
                # previous record before moving on
                print(pastList[0], pastList[1], pastList[2], \
                      '|'.join(rnaList), '|'.join(proList), \
                      '|'.join(genList), '|'.join(pubList), \
                      file=OUTPUT, sep='\t')

                # Empty all lists for the next record
                pastList = []
                rnaList = []
                proList = []
                genList = []
                pubList = []

            # sql_line comes back as a tuple with
            # 0 as symbol, 1 as gene_id, 2 as tax_id
            # 3 as rna_accession, 4 as pro_accession
            # 5 as gen_accession, 6 as pubmed_id
                
            # Store the records symbol, gene_id, and tax_id
            i = 0
            while (i < 3):
                pastList.append(sql_line[i])
                i += 1


            # Build each of the lists for the current record
            if (sql_line[1] == pastList[1]):
                # Build rnaList
                if (sql_line[3] is not None):
                    if (str(sql_line[3]) not in rnaList):
                        rnaList.append(str(sql_line[3]))

                # Build proList
                if (sql_line[4] is not None):
                    if (str(sql_line[4]) not in proList):
                        proList.append(str(sql_line[4]))

                # Build genList
                if (sql_line[5] is not None):
                    if (str(sql_line[5]) not in genList):
                        genList.append(str(sql_line[5]))

                # Build pubList
                if (sql_line[6] is not None):
                    if (str(sql_line[6]) not in pubList):
                        pubList.append(str(sql_line[6]))

            # Get the next sql_line
            sql_line = cursor.fetchone()
       
        # Read in the next symbol from INPUT
        line = INPUT.readline()
    
        # Print final record of previous symbol 
        print(pastList[0], pastList[1], pastList[2], \
              '|'.join(rnaList), '|'.join(proList), \
              '|'.join(genList), '|'.join(pubList), \
              file=OUTPUT, sep='\t')
    
    # Free the world...or just the resources on the computer, whichever it is
    conn.close()
    cursor.close()
    INPUT.close()
    OUTPUT.close()

if (__name__ == '__main__'):
    main()
