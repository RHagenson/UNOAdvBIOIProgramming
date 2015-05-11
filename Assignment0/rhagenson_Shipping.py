#! /usr/bin/python3
# Name: Ryan
# Class: BIOI 3500
# Assignment #: 0
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
# This program prompts the user to enter the distance a package is to be 
# shipped and the weight of the package.  It then calculates and displays
# the shipping charges for that package.

import math

def main():
    weight = 0           # Integer storing weight of package, 
                         # value limit 1 - 60
    miles = 0            # Integer storing distance of package shipment,
                         # must be positive value
    distanceUnits = 0.0  # Float storing distance / 100 for
                         # multiplication with rate
    rate = 0.0           # Float which holds rate based on package
                         # weight

    # Prints welcome message, and asks user for the weight of package and
    # distance of trip
    print("Welcome to the You Package It We Savage It " + \
          "Shipping Company. \n")
    
    # Continue to ask for the weight until it is an acceptable value
    while not (1 <= weight <= 60):
        weight = int(input("How heavy is your package in pounds (1-60)? "))
    
    # Continue to ask for the miles until it is an acceptable value
    while not (1 <= miles):
        miles = int(input("How far will you be " + \
                          "shipping the package in miles? "))

    # Determines rate based on weight
    if (weight <= 10): 
        rate = 5.01
    elif (weight <= 20):
        rate = 7.02
    elif (weight <= 30):
        rate = 9.03
    elif (weight <= 40):
        rate = 11.04
    else:  # weight <= 60
        rate = 13.05

    # Calculates number to multiply by rate for total cost
    distanceUnits = math.ceil(miles / 100)

    # Prints total shipping cost
    if (miles == 1):
        print("\nYour total shipping cost for 1 mile is", 
              rate * distanceUnits, ".\n")
    else:
        print("\nYour total shipping cost for", miles, "miles is", 
              "$%.2f.\n" % (rate * distanceUnits))

if (__name__ == '__main__'):
    main()
