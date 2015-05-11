// Name        : Mark Pauley
// Class       : CIST 1400
// Program #   : 9 
// Due Date    : 8 October 2013

// This program prompts the user to enter the distance a package is to be 
// shipped and the weight of the package.  It then colculates and displays
// the shipping charges for that package.

import java.util.Scanner;

public class Shipping
{
    public static void main(String []args)
    {
        Scanner input =             // Scanner object for user input
                        new Scanner(System.in);

        int miles = 0,              // Variable for user input (miles)  
            weight = 0;             // Variable for user input (weight)
        double rate = 0,            // Variable for rate per 100 miles
               distanceUnits = 0;   // Variable for calculating shipping cost
       
        // Prints welcome message
        System.out.printf("Welcome to the You Package It We Savage It " +
                          "Shipping Company.\n\n");

        // Prompts user for an integer between 0 and 60 for weight
        do 
        {
            System.out.print("How heavy is your package in pounds (1 - 60)?  ");
            weight = input.nextInt();
        } while ((weight <= 0) || (weight >= 60));

        // Prompts user for an integer greater than zero for miles
        do 
        {
            System.out.print("How far will you be shipping the package in " +
                             "miles?  ");
            miles = input.nextInt();
        } while (miles <= 0);
        

        // Finds shipping rate per 100 miles based on weight
        if (weight <= 10)
        {    
            rate = 5.01;
        }
        else if (weight <= 20)
        {
            rate = 7.02;
        }
        else if (weight <= 30)
        {
            rate = 9.03;
        }
        else if (weight <= 40)
        {
            rate = 11.04;
        }
        else                        // weight <= 60
        {
            rate = 13.05;
        }

        // Calculates number to multiply by rate for total cost
        distanceUnits = miles / 100;
        if ((miles % 100) != 0)
        {
            distanceUnits++;
        }

        // Prints total shipping cost
        System.out.printf("\nYour total shipping cost for %d miles is $%.2f.\n",
                          miles, rate * distanceUnits);
    }
}
