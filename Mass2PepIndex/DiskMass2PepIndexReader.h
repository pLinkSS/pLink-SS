#ifndef DISKMASS2PEPINDEXREADER_H_
#define DISKMASS2PEPINDEXREADER_H_

#include "Mass2PepIndexReader.h"
#include "PeptideHandler.h"

using namespace std;
namespace Mass2PepIndex{
class CDiskMass2PepIndexReader : public CMass2PepIndexReader
{
public:
	CDiskMass2PepIndexReader();
	virtual ~CDiskMass2PepIndexReader();
	
	virtual void Open(const std::string strWorkPath, const std::string strMetaName, ProteinRandomReaderType eType=ProteinIndex::Protein_Random_Index_Disk, bool bSQEX = false);
	virtual void Close();
	virtual bool GetNext(PEP_SQ & stPepSQ);
	virtual bool GetByID(PEP_SQ & stPepSQ, size_t tSQID);
	virtual bool GetProIDByPos(unsigned char cDatNum, size_t tPos, vector<size_t> & vtProID);//{return true;};
	virtual size_t GetIDByMass(double tMass);
	virtual void MoveTo(double dMass);
	bool	GetPepInfoByID(PEP_INFO_EX_TAG & stPep, size_t tSQID);
	size_t GetSize() 
	{
		return m_pepHandler.GetMaxPepNum();
	}
	
	CPeptideHandler m_pepHandler;
	
	size_t m_tDatNum;
	FILE ** aPos; 
	FILE ** aPid;
	
protected:
	
	
	ProteinIndex::CProteinRandomReaderFactory m_ProRandomFactory;	
	ProteinIndex::CRandomReader * m_proLoader;
	
	string strIDX;
	string strDAT;

	size_t m_tPepID;
	bool   m_bSQEX;	
	
	void 	_ReadPepFromFile(PEP_INFO_EX_TAG & sDatExTagBlock);

};

}

#endif /*DISKMASS2PEPINDEXREADER_H_*/
