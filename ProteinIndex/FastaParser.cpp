#include <string>
#include <fstream>
#include <iostream>
#include "FastaParser.h"

using namespace std;
using namespace ProteinIndex;

CFastaParser::CFastaParser(void)
{
}

CFastaParser::~CFastaParser(void)
{
}

// Get the number of all protein entries in the sequence database.
size_t CFastaParser::GetEntryNum(ifstream& ifDB) const
{
	if (!ifDB.is_open())
	{		
		proteomics_sdk::CErrInfo err("CFastaParser", "GetEntryNum","DB can't open!");
		throw runtime_error(err.Get().c_str());
		return 0;
	}
	ifDB.clear();
	ifDB.seekg(0,std::ios::beg);

	size_t tNum = 0;
//	std::streamsize tLength = DATABASE_FILE_BUF_SIZE;
	//char* pszBuf = new char[DATABASE_FILE_BUF_SIZE];

	while (!ifDB.eof()) 
	{
		string str;
		getline(ifDB, str);
		//ifDB.getline(pszBuf, tLength - 1);
		if ( str[0] == '>') 
		{	
			tNum ++ ;
			//cout <<tNum<<endl;
//			if(tNum==24472)
//			{
//				cout << pszBuf;
//			}
		}
	}

	// Important! clear and pointer to the starting position.
	ifDB.clear();
	ifDB.seekg(0,std::ios::beg);

	//delete pszBuf;
	cout <<"protein number is\t"<<tNum<<endl;
	return tNum;
}

// Read one Protein Entry at current postion  in file and return the CProtein object.//llq��ȡһ��ĵ��׵�AC��DE��SQ
bool  CFastaParser::ReadOnePrteinEntry(proteomics_sdk::CProtein& proEntry, ifstream& ifDB)
{
	proEntry.m_strAC.clear();
	proEntry.m_strDE.clear();
	proEntry.m_strSQ.clear();

	if(ifDB==NULL || !ifDB.is_open())
	{
		proteomics_sdk::CErrInfo err_info("CFastaParser", "ReadOnePrteinEntry","Open fasta file error!");
		throw runtime_error(err_info.Get().c_str());
	}
	if (ifDB.eof()) 
		return false;

//	std::streamsize tLength = DATABASE_FILE_BUF_SIZE;
//	char* pszBuf = new char[DATABASE_FILE_BUF_SIZE];
	//string  strSQ, strAC, strDE;
	string str;
	getline(ifDB,str);
	str = str.substr(0, str.find_last_not_of('\r')+1);
	if ( str.empty() || str[0] != '>' ) 
	{
		proteomics_sdk::CErrInfo err_info("CFastaParser", "ReadOnePrteinEntry","FASTA Format error!One Protein Entry must start with char '>'!");
		throw runtime_error(err_info.Get().c_str());
		return false;
	}

	
	size_t	tPos = str.find(' ', 0);

	proEntry.m_strAC = str.substr(1, tPos - 1);
	proEntry.m_strDE = str.substr(tPos + 1);

	ChangeInvalidChar(proEntry.m_strDE);
	_ReadEntrySQ(proEntry.m_strSQ,ifDB);

//	delete pszBuf;
//	return true;
	return true;
}

// Read protein entry's ACCESSION ID ( AC )by specified file position
bool CFastaParser::ReadEntryAC(string& strAC, long   lPosEntryStart, ifstream& ifDB)
{
	std::streamsize tLength = MAX_AC_LENGTH_PRO;
	char* pszBuf = new char[MAX_AC_LENGTH_PRO];


	ifDB.seekg(lPosEntryStart, ios::beg);

	if (ifDB.peek() == '>') 					//skip the char '>'
		ifDB.get();

	ifDB.getline(pszBuf, tLength - 1, ' ');

	strAC = pszBuf ;
	strAC = strAC.substr(0, strAC.find_last_not_of('\r')+1);
	delete pszBuf;
	return true;
}

//Read protein entry's Description ( DE ) by specified file position
bool CFastaParser::ReadEntryDE(string& strDE, long   lPosEntryStart, ifstream& ifDB)
{
	std::streamsize tLength = MAX_DE_LENGTH_PRO;
	char* pszBuf = new char[MAX_DE_LENGTH_PRO];

	ifDB.seekg(lPosEntryStart, ios::beg);

	ifDB.getline(pszBuf, tLength - 1, ' ');	// skip the AC name

//	ifDB.getline(pszBuf, tLength - 1);			// start read DE
	getline(ifDB, strDE);
//	strDE = pszBuf ;
	strDE = strDE.substr(0, strDE.find_last_not_of('\r')+1);
	ChangeInvalidChar(strDE);

	delete pszBuf;
	return true;
}

//Read protein entry's Sequence ( SQ ) by specified file position
bool CFastaParser::ReadEntrySQ(string& strSQ, long   lPosEntryStart, ifstream& ifDB)
{
//	std::streamsize tLength = DATABASE_FILE_BUF_SIZE;
//	char* pszBuf = new char[DATABASE_FILE_BUF_SIZE];

	ifDB.seekg(lPosEntryStart, ios::beg);

	//skip the AC and DE
	while(!ifDB.eof() 
		&& ifDB.get() != '\n' ) ;

	// start read SQ
	strSQ = "";
	string str;
	while(!ifDB.eof() 
		&& ifDB.peek() != '>' ) 
	{
//		ifDB.getline(pszBuf, tLength - 1);
		getline(ifDB, str);
		strSQ += str;
		strSQ = strSQ.substr(0, strSQ.find_last_not_of('\r')+1);
	}
//	delete pszBuf;

	return true;
}

void CFastaParser::ChangeInvalidChar(string& strContent)
{
	char chInvalid = 0x01;
	_ChangeChar(strContent, chInvalid);
	chInvalid = '<';
	_ChangeChar(strContent, chInvalid);
	chInvalid = '&';
	_ChangeChar(strContent, chInvalid);

}

void CFastaParser::_ChangeChar(string& strContent, char chInvalidToChange)
{
	string::size_type idx = 0;
	idx = strContent.find(chInvalidToChange, idx);
	while (idx != string::npos ) {
		strContent[idx] = ' ';
		idx = strContent.find(chInvalidToChange, idx);
	}
}
// //read entry for SQ
bool CFastaParser::_ReadEntrySQ(string& strSQ, ifstream& ifDB)
{
//	std::streamsize tLength = DATABASE_FILE_BUF_SIZE;
	string  strLine;

	
	
	while(!ifDB.eof() 
		&& ifDB.peek() != '>' ) 
	{
		string str;

		getline(ifDB, str);
		
		size_t tlen = str.size();
		for(size_t t=0; t < tlen; ++t)
		{
			if(str[t] >='A' && str[t] <='Z') strSQ += str[t]; 
		}
//		strSQ += str;
//		strSQ = strSQ.substr(0, strSQ.find_last_not_of('\r')+1);
	}

	return true;
}

