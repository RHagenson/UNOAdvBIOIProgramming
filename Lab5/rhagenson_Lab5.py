#! /usr/bin/python3

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def main():
    INPUT = open('../storage/human.rna.fna', 'r')   # The given file to be 
                                                    # read for 
                                                    # this assignment

    number_of_sequences = 0
    number_of_nucleotides = 0
    shortest_sequence = float('inf')
    longest_sequence = 0
    average_sequence = 0

    for record in SeqIO.parse(INPUT, "fasta"):
        if (record.seq):
            number_of_sequences += 1
        
        number_of_nucleotides += len(record)
        
        if (len(record) < shortest_sequence):
            shortest_sequence = len(record)

        if (longest_sequence < len(record)):
            longest_sequence = len(record)
    
    average_sequence = number_of_nucleotides / number_of_sequences
    
    # Print the 
    print('Number of sequences:', number_of_sequences)
    print('Number of nucleotides:', number_of_nucleotides)
    print('Shortest sequence:', shortest_sequence)
    print('Longest sequence:', longest_sequence)
    print('Average sequence: %4.2f' % average_sequence)


if __name__ == '__main__':
    main()
