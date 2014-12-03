#include "StringUtility.h"

//using namespace proteomics_sdk;

//const size_t MAX_CHAR_IN = 200 * 1024 * 1024;

const int MAX_BUFFER_SIZE = 819200;

char ss[MAX_BUFFER_SIZE] =
{ 0 };

bool CheckFileValid(const string & strFile)
{
	return 0 == access(strFile.c_str(), 0x04);
}

void ReturnAfternEqual(string & strTemp)
{
	size_t pos = strTemp.find_first_of('=');
	strTemp = strTemp.substr(pos + 1);
}

void ReadUntilString(ifstream & fin, string & strBuf, const string & strVal)
{
	string strTemp;
	while (1)
	{
		GetLine(fin, strBuf, strTemp);
		if (strTemp == strVal || strTemp[0] == 0) //the last line of the file
			break;
	}
}

void ReadUntilString(ifstream & fin, const string & strVal)
{
	string strTemp;
	while (1)
	{
		GetLine(fin, strTemp);
		if (strTemp == strVal || strTemp[0] == 0) //the last line of the file
			break;
	}
}

void GetNSplitStringByComma(string & strTemp, vector<string> & vstrTemp)
{
	//todo 存在错误当strTemp.size() == 0 时
	size_t nProts = 0;
	sscanf(strTemp.c_str(), PRI_SIZE_T, &nProts);
	strTemp += ',';
	string strV;
	size_t npos = strTemp.find_first_of(',');

	vstrTemp.clear();
	for (size_t j = 1; j <= nProts; ++j)
	{
		size_t temp = npos + 1;
		npos = strTemp.find_first_of(',', temp);
		if (string::npos == npos)
		{
			CErrInfo info("ResultParser", "SplitStringByComma",
					"Cannot parse the proteins list: " + strTemp);
			throw runtime_error(info.Get());
		}

		strV = strTemp.substr(temp, npos - temp);
		vstrTemp.push_back(strV);
	}
}

void SplitStringByComma(const string & strTemp, vector<string> & vstrTemp)
{
	//todo 存在错误当strTemp.size() == 0 时
	string strVal;
	vstrTemp.clear();

	//	if (strTemp.size() == 0)
	//		return;

	for (size_t i = 0; i < strTemp.size(); i++)
	{
		if (strTemp[i] == ',')
		{
			vstrTemp.push_back(strVal);
			strVal.clear();
			continue;
		}

		strVal += strTemp[i];
	}
	vstrTemp.push_back(strVal);
}

void FetchLetter(string & strTemp)
{
	string strVal;

	for (size_t i = 0; i < strTemp.size(); i++)
		if (isalpha(strTemp[i]))
			strVal += strTemp[i];
	strTemp = strVal;
}

void ThrowFirstOfRemain(ifstream & fin, string & strBuf, const string & strItem, string & strVal)
{
	RemainOfString(fin, strBuf, strItem, strVal);
	strVal.erase(strVal.begin());
}

void ThrowFirstOfRemain(ifstream & fin, const string & strItem, string & strVal)
{
	RemainOfString(fin, strItem, strVal);
	strVal.erase(strVal.begin());
}

void RemainOfString(ifstream & fin, string & strBuf, const string & strItem, string & strVal)
{
	string strLine;
	while (1)
	{
		GetLine(fin, strBuf, strLine);

		if (strLine.substr(0, strItem.length()) == strItem)
		{
			strVal = strLine.substr(strItem.length());
			break;
		}
	}
	//cout << strLine << endl;
}

void RemainOfString(ifstream & fin, const string & strItem, string & strVal)
{
	string strLine;
//	cout << "Item=" << strItem << endl;
	while (1)
	{
		GetLine(fin, strLine);

		//cout << "Line=" << strLine << endl;
		if (strLine.substr(0, strItem.length()) == strItem)
		{
			strVal = strLine.substr(strItem.length());
			break;
		}
	}
}

const size_t MAX_CHAR_IN = 50000000;
char buf[MAX_CHAR_IN + 1];
size_t strlcnt = 0;

void _GetStream(ifstream & fin, string & strBuf)
{
	fin.read(buf, MAX_CHAR_IN);

	int size = fin.gcount();
	buf[size] = 0;

	strBuf.clear();
	strBuf = buf;
	strlcnt = 0;
}

bool GetLine(ifstream & fin, string & strBuf, string & strRet)
{
	char chr;
	strRet.clear();

	while (1)
	{
		if (strlcnt >= strBuf.length() && fin.eof())
			return false;

		if (strlcnt >= strBuf.length())
			_GetStream(fin, strBuf);

		while (1)
		{
			chr = strBuf[strlcnt++];

			if (strlcnt >= strBuf.length())
				_GetStream(fin, strBuf);

			if (0 == chr || 0xa == chr || 0xd == chr)
				break;

			strRet += chr;
		}
		size_t len = strRet.length();
		while (len > 0 && (0xa == strRet[len - 1] || 0xd == strRet[len - 1] || 32
				== strRet[len - 1] || 0 == strRet[len - 1] || 9 == strRet[len - 1]))
		{
			strRet[len - 1] = 0;
			--len;
		}
		if (len == 0)
			continue;//read black line continue
		else
			break;
	}
	return true;
}

bool GetLine(ifstream & fin, string & strRet)
{
	strRet.clear();
	if (fin.eof())
		return false;
	try
	{
		while (!fin.eof())
		{
			string is;
			getline(fin, is);
			strRet = is;
			if (strRet.size() != 0)
				break;
		}
	} catch (ifstream::failure e)
	{
		cout << "Exception opening/reading file";
		return false;
	}
	if (strRet.size() == 0)
		return false;
	return true;

	/*// 2013.10.28 old backup
	memset(ss, '\0', sizeof(ss));

	strRet = "";

	while (1)
	{
		if (fin.eof())
		{
			if (0 == strRet.length())
				return false;
			else
				return true;
		}

		fin.getline(ss, MAX_BUFFER_SIZE);
		size_t len = strlen(ss);
		//cout << ss << endl;

		//	printf("ss = %s\n", ss);
		while (len > 0 && (0xa == ss[len - 1] || 0xd == ss[len - 1]))
		{
			ss[len - 1] = 0;
			--len;
		}

		if (len == 0)
			continue;//read black line continue

		strRet += ss;
		if (len == MAX_BUFFER_SIZE - 1)
		{
			strRet = strRet.substr(0, strRet.length() - 1);
		}
		if (len != MAX_BUFFER_SIZE - 1)
			break;
	}
	return true;*/
}

void ReadLineC(FILE * fin, string & strRet)
{
	const int MAX_BUFFER_SIZE = 81920;

	char buffer[MAX_BUFFER_SIZE] =
	{ 0 };
	char chr;
	int i = 0;
	do
	{
		chr = getchar();
		buffer[i++] = chr;
	} while (chr != '\n');
	buffer[i - 1] = '\0';
	strRet = buffer;
	//cout << strRet;
}

void DeleteFrontBlank(string & strTemp)
{
	while (1)
	{
		if (32 != strTemp[0])
			break;

		strTemp.erase(strTemp.begin());
	}
}

void StringToHex(string & strValue)
{
	size_t nPos = 0;

	while ((nPos = strValue.find("%")) != string::npos)
	{
		string strHexValue = strValue.substr(nPos + 1, 2);

		int nVal1 = 0, nVal2 = 0;

		if ('9' >= strHexValue[1])
			nVal1 = strHexValue[1] - '0';
		else
		{
			if (isupper(strHexValue[1]))
				nVal1 = strHexValue[1] - 'A' + 10;
			else
				nVal1 = strHexValue[1] - 'a' + 10;
		}
		if ('9' >= strHexValue[0])
			nVal2 = strHexValue[0] - '0';
		else
		{
			if (isupper(strHexValue[0]))
				nVal2 = strHexValue[0] - 'A' + 10;
			else
				nVal2 = strHexValue[0] - 'a' + 10;
		}

		nVal1 += nVal2 * 16;
		string strHexString;
		strHexString += (char) (nVal1);

		strValue.replace(nPos, 3, strHexString);
	}
}

void GetBetweenColon(string & strTemp)
{
	size_t nPos = strTemp.find_first_of('"');
	strTemp = strTemp.substr(nPos + 1);
	nPos = strTemp.find_first_of('"');
	strTemp = strTemp.substr(0, nPos);
}

void GetBeforeDot(string & strTemp)
{
	return;
	size_t pos = strTemp.find_first_of('.');
	if (pos != string::npos)
		strTemp = strTemp.substr(0, pos);
}

void DeleteRedundancy(vector<string> & vStrTemp)
{
	set<string> setStrTemp;

	for (size_t t = 0; t < vStrTemp.size(); t++)
	{
		GetBeforeDot(vStrTemp[t]);
		setStrTemp.insert(vStrTemp[t]);
	}

	vStrTemp.clear();
	for (set<string>::iterator it = setStrTemp.begin(); it != setStrTemp.end(); it++)
	{
		vStrTemp.push_back(*it);
	}
}

string SplitTitle(const string & title)
{

	string strTemp, str1, str2, ans;

	size_t pos = title.find_last_of('.');
	strTemp = title.substr(0, pos);

	pos = strTemp.find_last_of('.');
	strTemp = strTemp.substr(0, pos);

	pos = strTemp.find_last_of('.');
	str1 = strTemp.substr(pos + 1);

	strTemp = strTemp.substr(0, pos);
	pos = strTemp.find_last_of('.');
	str2 = strTemp.substr(pos + 1);

	ans += strTemp.substr(0, pos);
	ans += ", ";
	ans += str2;

	if (str1 != str2)
	{
		ans += " - ";
		ans += str1;
	}

	return ans;
}

string GetScanNum(const string & strTemp)
{
	string strVal = strTemp;

	size_t pos = strVal.find_last_of('.');

	strVal = strVal.substr(0, pos);
	pos = strVal.find_last_of('.');

	strVal = strVal.substr(0, pos);

	return strVal;
}

string GetSpectraFileName(const string & strFilePath, const int & FirstScan, const int & LastScan,
		const int & ChargeState)
{
	char strTmp[PATH_MAX] =
	{ '\0' };

	string strVal;
	strVal.clear();

	strVal = strFilePath;
	sprintf(strTmp, ".%d.%d.%d.dta", FirstScan, LastScan, ChargeState);

	strVal += strTmp;

	return strVal;
}

bool GetBeforeEqual(string & strTemp, string & strResult)
{
	strResult.clear();

	size_t pos = strTemp.find_first_of("=");

	if (pos != string::npos)
	{
		strResult = strTemp.substr(0, pos);
		strTemp = strTemp.substr(pos + 1);

		return true;
	}
	return false;
}

string stringToUpper(const string & strVal)
{
	string strRet = "";

	for (size_t t = 0; t < strVal.size(); t++)
	{
		strRet += toupper(strVal[t]);
	}

	return strRet;
}

int setFindInt(const set<string> & setVal, const string & strVal)
{
	int cnt = 0;

	for (set<string>::const_iterator it = setVal.begin(); it != setVal.end(); it++, cnt++)
	{
		if (*it == strVal)
			return cnt;
	}

	return -1;
}
