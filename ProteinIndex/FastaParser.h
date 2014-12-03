#pragma once
#include <string>
#include <fstream>
#include "../include/sdk.h"
#include "ProteinIndex.h"
#include "DBFormatParser.h"

using namespace std;
namespace ProteinIndex{
class CFastaParser :
	public CDBFormatParser
{

public:
	CFastaParser(void);
	~CFastaParser(void);

	/// Get the number of all the protein entries in the sequence database. Note: Before using this function, must call OpenFile() earlier, and after this function must call CloseFile().
	size_t GetEntryNum(ifstream& ifDB) const;	

	// Read one Protein Entry at current postion  in file and return the CProtein object.
	bool ReadOnePrteinEntry(CProtein& proEntry, ifstream& ifDB);

	// Read protein entry's ACCESSION ID ( AC )by specified file position
	bool ReadEntryAC(string& strAC, long   lPosEntryStart, ifstream& ifDB);

	//Read protein entry's Description ( DE ) by specified file position
	bool ReadEntryDE(string& strDE, long   lPosEntryStart, ifstream& ifDB);

	//Read protein entry's Sequence ( SQ ) by specified file position
	bool ReadEntrySQ(string& strSQ, long   lPosEntryStart, ifstream& ifDB);

	// Change some Invalid chars such as '&','>'in DE to space char ' '. These invalid chars can interrupt the display of the protein information. 
	void ChangeInvalidChar(string& strContent);

protected:

	/// Transform the char var 'chInvalidToChange' into space char ' '.
	void _ChangeChar(string& strContent, char chInvalidToChange);

	/// read entry for sequence
	bool _ReadEntrySQ(string& strSQ, ifstream& ifDB);
};
}
