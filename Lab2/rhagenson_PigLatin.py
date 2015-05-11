#! /usr/bin/python3
# Name: Ryan Hagenson
# Class: BIOI 3500
# Lab #: 2
#
# Due date: Tuesday, February 24th
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
# This program reads in a single word from the user and concerts the word using
# regular expression to Pig Latin. The input is case insensitive, however the
# output will be displayed in all lowercase.

import re

def main():
    userInput = ''                                        # The string
                                                          # containing the
                                                          # user-input word
    pattern = re.compile('([^aeiouAEIOU]+)([A-Za-z]+)$')  # Matches if the 
                                                          # first letter is 
                                                          # a consonant, 
                                                          # followed by zero
    # prompts the user for the word to be converted
    userInput = input('Enter an English word:\t')

    # if the pattern is matched, print the consonant-based Pig Latin using
    # re.sub()
    if (pattern.match(userInput)):
        print(re.sub(pattern, '\\2\\1'+'ay', userInput).lower())
    # else print the vowel-based Pig Latin
    else:
        print(userInput.lower() + 'way')

if (__name__ == '__main__'):
    main()
