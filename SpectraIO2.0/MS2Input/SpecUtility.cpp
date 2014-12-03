/*
 * SpecUtility.cpp
 *
 *  Created on: 2009-7-13
 *      Author: fan
 */

#include <string>
#include "../../include/predefine.h"
#include "SpecUtility.h"

using namespace std;

bool InvalidChar ( char c )
{
	return c > 57 || c < 43;
}

bool Mass_Less(const SPEC_SORTER_INDEX_INFO & a, const SPEC_SORTER_INDEX_INFO & b)
{
	return a.lfMass < b.lfMass;
}

bool Mass_Greater(const SPEC_SORTER_INDEX_INFO & a, const SPEC_SORTER_INDEX_INFO & b)
{
	return a.lfMass > b.lfMass;
}

string string_toupper(string  str)   
{   
	int i, j;   
            
	j = str.length();   
                    
	for(i = 0; i < j;i++)   
		str[i]   =   toupper((int)str[i]);   
          
	return str;
} 
