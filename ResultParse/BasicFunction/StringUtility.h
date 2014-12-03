#ifndef STRINGUTILITY_H_
#define STRINGUTILITY_H_

#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "../../include/sdk.h"
#include "../bio_analysis_sdk.h"

using namespace std;

void FetchLetter(string & strTemp);

void ReturnAfternEqual(string & strTemp);

void DeleteFrontBlank(string & strTemp);
void GetBeforeDot(string & strTemp);


bool GetLine(ifstream & fin, string & strRet);
bool GetLine(ifstream & fin, string & strBuf, string & strRet);

void ReadLineC(FILE * fin, string & strRet);

bool CheckFileValid(const string & strFile);

void StringToHex(string & strValue);

void GetBetweenColon(string & strTemp);

void ReadUntilString(ifstream & fin, string & strBuf, const string & strVal);
void ReadUntilString(ifstream & fin, const string & strVal);

void ThrowFirstOfRemain(ifstream & fin, string & strBuf, const string & strItem, string & strVal);
void
ThrowFirstOfRemain(ifstream & fin, const string & strItem, string & strVal);

void RemainOfString(ifstream & fin, string & strBuf, const string & strItem, string & strVal);

//string GetScanNum(const string & strTemp);
//string SplitTitle(const string & title);

void RemainOfString(ifstream & fin, const string & strItem, string & strVal);
void DeleteRedundancy(vector<string> & vStrTemp);

void GetNSplitStringByComma(string & strTemp, vector<string> &vstrTemp);
void SplitStringByComma(const string & strTemp, vector<string> & vstrTemp);

string GetSpectraFileName(const string & strFilePath, const int & FirstScan, const int & LastScan,
		const int & ChargeState);

bool GetBeforeEqual(string & strTemp, string & strResult);

string stringToUpper(const string & strVal);

int setFindInt(const set<string> & setVal, const string & strVal);

#endif /* FUNCTIONS_H_ */
