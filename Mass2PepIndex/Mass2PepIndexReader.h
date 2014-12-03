#ifndef MASS2PEPINDEXREADER_H_
#define MASS2PEPINDEXREADER_H_

#include "../ProteomicsSDK/ProteomicsSDK.h"
#include "../ProteinIndex/ProteinReader.h"
#include "../ProteinIndex/ProteinReaderFactory.h"
#include "Mass2PepIndex.h"

using namespace std;
using namespace ProteinIndex;
using namespace proteomics_sdk;

//class CProteinRandomReaderFactory;

namespace Mass2PepIndex
{

class CMetaReader;

class CMass2PepIndexReader
{
public:
	virtual ~CMass2PepIndexReader()
	{
		
	}
	
	virtual void Open(const std::string strWorkPath, const std::string strMetaName, ProteinRandomReaderType eType=ProteinIndex::Protein_Random_Index_Disk, bool bSQEX = false) = 0;
	virtual void Close() = 0;
	virtual bool GetNext(PEP_SQ & stPepSQ) = 0;
	virtual bool GetByID(PEP_SQ & stPepSQ, size_t tSQID) = 0;	
	virtual bool GetProIDByPos(unsigned char cDatNum, size_t tPos, vector<size_t> & vtProID) = 0;
	virtual size_t GetIDByMass(double tMass) = 0;//czhou
	virtual void MoveTo(double dMass) = 0;
	virtual size_t GetSize() = 0; // add by yjwu
};

}

#endif /*MASS2PEPINDEXREADER_H_*/

