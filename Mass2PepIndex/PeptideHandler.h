#ifndef PEPTIDEHANDLER_H_
#define PEPTIDEHANDLER_H_

#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "MetaReader.h"
#include "Mass2PepIndex.h"
//#include "ace/OS.h"
//#include "ace/Mem_Map.h"

#include "../include/sdk.h"
#include "../include/predefine.h"

using namespace std;
typedef unsigned char uchar;

namespace Mass2PepIndex
{

class CPeptideHandler
{
public:
	CPeptideHandler();
	virtual ~CPeptideHandler();
	
	void Init(const string & strWorkDir, const string &strMetaName);
	void Close();

	size_t GetMaxPepNum();
	string GetDBName();	
	size_t GetMultiplier();
	
	string GenerFileName(uchar uFileIdx,string strApp);
//	size_t GetIDByMass(size_t);
	void   GetCurFileID(size_t tPepID, uchar &cDatNum, size_t &tInvertedFilesPos);
	uchar  GetFileIDByMass(size_t tMass, long &ts,long &te);
	size_t GetCurrentPepFileEnd();
	size_t GetCurrentPepFileHead(); // added at 2014.3.18
	
	size_t 	CountIdxStructSize();
	void 	ReadIdxFromFile(size_t &tMass, long &lPos);
	
	size_t 	CountPepStructSize();

	void OpenInvertedFiles(FILE ** aDat, string strDat, size_t m_tDatNum);
	void CloseInvertedFiles(FILE ** aDat, string strDat, size_t m_tDatNum);
	
	bool InitDiskDatFileByPepID(size_t tPepID, size_t tStartPos, string &strIDX, string &strDat);
	void InitDiskDatFileByFileID(uchar uFileID, size_t tStartPos, string &strIDX, string &strDat);

	void SeekToPlaceByPepID(size_t tPepID); //added at 2014.3.18, change file pointer to the right place in current file.
//	void CloseAllFILE(ACE_Mem_Map * aDat, string strDat, size_t m_tDatNum);

//	bool InitMMapDatFileByPepID(size_t tProID, size_t tStartPos,string &strIDX, string &strDat);
//	void InitMMapDatFileByFileID(uchar uFileID, size_t tStartPos,string &strIDX, string &strDat);

	FILE * m_fpIDX;
	FILE * m_fpDAT;
	
//	ACE_Mem_Map m_MapIDX;
//	ACE_Mem_Map m_MapDAT;
	
	CMetaReader m_metaReader;	
protected:
	size_t 	m_tPriorPepNum;
	
	short	m_uCurFileID;
	
	string	_GenerFileName();
	size_t m_tHeadPepNum; // added at 2014.3.18
};

}

#endif /*PEPTIDEHANDLER_H_*/
