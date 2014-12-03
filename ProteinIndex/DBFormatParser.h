#pragma once

#include <fstream>
#include "../include/sdk.h"
#include "ProteinIndex.h"
using namespace proteomics_sdk;
namespace ProteinIndex{
class CDBFormatParser
{
public:

	CDBFormatParser(void)
	{
	}

	virtual ~CDBFormatParser(void)
	{
	}

	// Get the number of all protein entries in the sequence database.
	virtual size_t GetEntryNum(ifstream& ifDB) const  = 0;	

	// Read one Protein Entry at current postion  in file and return the CProtein object.
	virtual	bool ReadOnePrteinEntry(CProtein& proEntry, ifstream& ifDB) = 0;

	// Read protein entry's ACCESSION ID ( AC )by specified file position
	virtual bool ReadEntryAC(string& strAC, long   lPosEntryStart, ifstream& ifDB) = 0;

	//Read protein entry's Description ( DE ) by specified file position
	virtual bool ReadEntryDE(string& strDE, long   lPosEntryStart, ifstream& ifDB) = 0;

	//Read protein entry's Sequence ( SQ ) by specified file position
	virtual bool ReadEntrySQ(string& strSQ, long   lPosEntryStart, ifstream& ifDB) = 0;
	
};
}
