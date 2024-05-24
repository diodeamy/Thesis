#include <iostream>
#include <string>


#include "miscellaneous.h"
using namespace std;




/* This function reads user input and checks for an integer input. Otherwise gives error promting the user to insert a new number.
*/
int getIntegerFromUser()
{
    bool test = false;
    int choice;
    
    while (not test)
    {
        string input;
        cin >> input;
    
        if ( isStringInteger( input, &choice) )
            test = true;
        else cout << "Please insert your choice as a integer number:   ";
    }
    
    return choice;
}




/* This function reads user input and checks for an integer input. Otherwise gives error promting the user to insert a new number. The integer has to be larger or equal than 'minValue'.
*/
int getIntegerFromUser(int minValue)
{
    bool test = false;
    int choice;
    
    while (not test)
    {
        string input;
        cin >> input;
    
        if ( isStringInteger( input, &choice) and choice>=minValue )
            test = true;
        else cout << "Please insert your choice as a number larger or equal than '" << minValue << "' :   ";
    }
    
    return choice;
}



/* This function reads user input and checks for an integer input. Otherwise gives error promting the user to insert a new number. The integer has to be in the range 'minValue' to 'maxValue'.
*/
int getIntegerFromUser(int minValue,
		      int maxValue)
{
	bool test = false;
	int choice;
	
	while (not test)
	{
		string input;
		cin >> input;
		
		if ( isStringInteger( input, &choice) and choice>=minValue and choice<=maxValue )
			test = true;
		else cout << "Please insert your choice as a number in the range '" << minValue << " to " << maxValue << "' :   ";
	}
	
	return choice;
}




/* This function reads user input and checks for a double. Otherwise gives error promting the user to insert a new number. The real number can take any value.
*/
double getDoubleFromUser()
{
	double choice;
	bool test = false;
	
	while (not test)
	{
		string input;
		cin >> input;
		
		if ( isStringDouble( input, &choice) )
			test = true;
		else cout << "Please insert a valid real number:   ";
	}
	
	return choice;
}



/* This function reads user input and checks for a double. Otherwise gives error promting the user to insert a new number. The real number must be larger than minValue.
*/
double getDoubleFromUser(double minValue)
{
	double choice;
	bool test = false;
	
	while (not test)
	{
		string input;
		cin >> input;
		
		if ( isStringDouble( input, &choice) and choice>minValue )
			test = true;
		else cout << "Please insert a valid real number which is larger than " << minValue << ":   ";
	}
	
	return choice;
}



/* This function reads user input and checks for a double. Otherwise gives error promting the user to insert a new number. The real number must be larger than minValue.
*/
double getDoubleFromUser(double minValue, double maxValue)
{
    double choice;
    bool test = false;
	
    while (not test)
    {
        string input;
        cin >> input;
		
        if ( isStringDouble( input, &choice) and choice>minValue and choice<maxValue)
            test = true;
        else cout << "Please insert a valid real number which is between " << minValue << " and " << maxValue <<":   ";
    }
	
    return choice;
}




/* This function takes as input an array of character and checks (starting with the first element) how many of its elements in a row are digits. If at least the first element is a digit, than it returns true and also it returns the continuous sequence of digits that it found.
*/
bool isStringInteger(string const str,
		  int *number)
{
	bool test = false;
	
	if ( (int)str[0]>=(int)'0' and (int)str[0]<=(int)'9' )
	{	// than the first character in the char array is a digit, so read the sequence of digits
		*number = atoi( str.c_str() );	// atoi() reads the number until if find other chracters, than it stops
		test = true;
	}
	
	return test;
}



/* This function takes as input an array of character and checks if the input is a double. If at least the first element is a digit (or '.' or '-' or '-.' followed by a digit), it returns true and also it returns the double precission number that it found.
*/
bool isStringDouble(string const str,
		   double *result)
{
	bool test = false;
	
		// if the first character in the char array is a digit, read the double
	if ( (int)str[0]>=(int)'0' and (int)str[0]<=(int)'9' )
	{
		*result = atof( str.c_str() );	// atof() reads the number until if find other chracters, than it stops
		test = true;
	}
		// if first character is a '.' but the 2nd is a number, so the input is still a double
	else if ( (int)str[0]>=(int)'.' and (int)str[1]>=(int)'0' and (int)str[1]<=(int)'9' )
	{
		*result = atof( str.c_str() );
		test = true;
	}
		// if first character is a '-' but the 2nd is a number, so the input is still a double
	else if ( (int)str[0]>=(int)'-' and (int)str[1]>=(int)'0' and (int)str[1]<=(int)'9' )
	{
		*result = atof( str.c_str() );
		test = true;
	}
		// if first 2 characters are a '-.' but the 3rd is a number, so the input is still a double
	else if ( (int)str[0]>=(int)'-' and (int)str[1]>=(int)'.' and (int)str[2]>=(int)'0' and (int)str[2]<=(int)'9' )
	{
		*result = atof( str.c_str() );
		test = true;
	}
	
	return test;
}



