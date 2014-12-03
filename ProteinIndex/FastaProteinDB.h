#pragma once

#include <string>
#include "../include/sdk.h"
#include "ProteinIndex.h"
#include "FastaParser.h"

using namespace std;
using namespace proteomics_sdk;
namespace ProteinIndex{

class CFastaParser;
	
class CFastaProteinDB
{
public:
	CFastaProteinDB(void);
	~CFastaProteinDB(void);

	string&	GetDBPath(void);
	string& GetDBName(void);

	// Use the No.of protein in the database---proteinID to get the information of one protein.
	CProtein& GetProteinByID(size_t proteinID);


	CProtein& ReadOnePrteinEntry();

	// Read one Protein Entry at current postion  in FASTA file and return the CProtein object
	bool ReadOnePrteinEntry(proteomics_sdk::CProtein& rPro);

	inline CFastaParser& GetParser(void){
		return m_Parser;
	};

	void SetDBPath(const string& filepath);
	void SetDBName(const string& dbname);

	inline bool	Empty(){ return GetDBName().empty() || GetDBPath().empty(); }

	/// Open file stream
	void OpenFile(const char* szFilePath);
	void OpenFile();
	/// Close file stream
	void CloseFile(void);

	inline bool IsEof(){
		return m_ifIn.eof();
	};

	inline bool IsOpen()
	{
		return m_ifIn.is_open();
	}

	// To the start position of the file.
	void SeekBegin(void);

	// Get the number of all the protein entries in the sequence database
	size_t GetEntryNum(void);

	// Get current protein entry that read from sequence database.
	CProtein& GetProEntry(void);

protected:

	CFastaParser m_Parser ;

	CProtein m_proEntry ;

	// The name of protein database. such as "SWISSPROT","nrdb".
	string m_strDBName ;

	// The path of the DB. such as "c:\\database\\nr.fasta"
	string m_strPath;

private:
	///	file stream
	ifstream m_ifIn ;

};
}
