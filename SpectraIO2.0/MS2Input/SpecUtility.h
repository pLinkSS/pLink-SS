/*
 * SpecUtility.h
 *
 *  Created on: 2009-7-13
 *      Author: fan
 */
#include <string>

#ifndef SPECUTILITY_H_
#define SPECUTILITY_H_

using namespace std;

bool InvalidChar ( char c );

bool Mass_Less(const SPEC_SORTER_INDEX_INFO & a, const SPEC_SORTER_INDEX_INFO & b);

bool Mass_Greater(const SPEC_SORTER_INDEX_INFO & a, const SPEC_SORTER_INDEX_INFO & b);

string string_toupper(string  str);

#endif /* SPECUTILITY_H_ */
