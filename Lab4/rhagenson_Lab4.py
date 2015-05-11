#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Assignment #: 0
#
# Due date: Tuesday, March 17 (Extended since I lacked the lecture audio)
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
# This program takes a gene symbol, the displays the corresponding gene ID,
# taxonomy ID, locus tag, database references, chromosome, map location,
# decription, type of gene, nomenclature symbol, full name, nomenclature
# status, other designation, and modified date. All on a new line preceded by
# its description (e.g. Gene ID: 1080).

import psycopg2

def main():
    labels = ('Gene ID: ', \
              'Taxonomy ID: ', \
              'Locus Tab: ', \
              'Database references: ', \
              'Chromosome: ', \
              'Map location: ', \
              'Description: ', \
              'Type of gene: ', \
              'Nomenclature symbol: ', \
              'Full name: ', \
              'Nomenclature status: ', \
              'Other designation: ', \
              'Modified date: ')
    
    # Open connnection to SQL database
    # Database details have been removed #
    conn = psycopg2.connect("dbname="" user="" password="" host=""")
    
    # Create cursor for interaction with SQL database
    cursor = conn.cursor()

    # Gather data from SQL database based on gene symbol
    cursor.execute("select gene_id, tax_id, locustag, dbxrefs, \
                   chromosome, map_location, description, \
                   type_of_gene, nom_symbol, fullname, \
                   nom_status, other_desig, mod_date \
                   from geneinfo where symbol='CFTR';")
    
    # Store the first result of the above execute as a tuple
    symbolInfo = cursor.fetchone()
    
    # Print each label and its information
    i = 0
    while (i < len(labels)):
        print(labels[i], symbolInfo[i])
        i += 1

if (__name__ == '__main__'):
    main()
