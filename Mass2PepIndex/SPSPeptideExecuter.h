#ifndef SPSPEPTIDEEXCUTER_H_
#define SPSPEPTIDEEXCUTER_H_

#include <string>
#include <stdexcept>
#include <iostream>
#include <algorithm>

#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/option.h"

#include "../ProteinIndex/ProteinWriterFactory.h"
#include "../ProteinIndex/ProteinReader.h"
#include "../ProteinIndex/DiskRandomSQReader.h"
#include "../ProteinIndex/ProteinReaderFactory.h"
#include "../ProteinIndex/PepCalcFunc.h"

#include "PepFilter.h"
#include "MetaCreater.h"
#include "DigestSimulator.h"
#include "PeptideHandler.h"


using namespace std;
using namespace ProteinIndex;

namespace Mass2PepIndex{

typedef  PEP_INFO_EX_TAG SPS_PEP_STRUCT;


bool SPS_PEP_STRUCT_CMP(const SPS_PEP_STRUCT & PepFir, const SPS_PEP_STRUCT &PepSec);//used for sort
bool HEAP_STRUCT_CMP(const _HEAP_ITEM & PepFir, const _HEAP_ITEM &PepSec);


class CSPSPeptideExecuter
{
public:
	CSPSPeptideExecuter();
	virtual ~CSPSPeptideExecuter();
	
	void 	Close();	
	
	void	Init(const PINDEX_HEAD& pIndexHead, const PINDEX_ITEM & pIndexItem) ;
	void 	EmptyPepRecord(size_t tProID);
	
	void	GetPepInfor(size_t&, SPS_PEP_STRUCT&);
	size_t	GetUsedPepNum();
	size_t 	GetMaxPepNum();
	void 	AppendPepIdxBySQ(string &strProSQ, size_t tProID, size_t tPos);
	void	GetUsefulPep();
	
	void 	WritePepIndex(size_t tProID);
	void	WriteMeta();
	void 	MergePepIndex();
	void	MergePepIndexWithPro2Pep();
	vector<SPS_PEP_STRUCT>	m_vsDatExTagBlock;
protected:
	
	long		m_sPos,m_ePos;

	PINDEX_HEAD m_PIndexHead;
	PINDEX_ITEM m_PIndexItem;
	
	CDigestSimulator m_Digester;
	CPepCalcFunc	 m_PepFunc;
	CMetaCreater	 m_CMetaCreater;

	proteomics_sdk::CEnzyme m_Enzyme;
	CPepFilter* m_pFilter;
		
	
	
	size_t m_tBegProID;
	
	vector <string> m_vsTmpFileName;
	
	long m_tUniquePepSQNum;
	size_t m_tMultplier;
	size_t m_tCleaveWay;

	std::string m_strAAListPath;
	std::string m_strWorkDir;
	
	PEP_INFO_EX			m_stDatBlock;
	vector<PEP_INFO_EX>	m_vsDatBlock;
	
	vector<size_t> 		m_vpPro;
//	vector<size_t>		m_tDatID;
	PEP2PRO				m_sPep2Pro;
	vector<PEP2PRO>		m_vsPep2Pro;
	size_t				m_tPep2ProSize;
	
	size_t				m_tBitNum;
	vector<size_t>	    m_vtPepNum;
	vector<double>		m_vdfMassRange;

//	bool				m_bESort;
	size_t				m_tMaxMemorySize;
	double				m_dfMassRange;
//	
	size_t				m_tOutTxt;
	META_ITEM 			m_stMetaItem;
	
	FILE *				m_fPep2Pro;
	FILE *				m_fDebug;
	string				m_strDebugFileName;
	
	size_t				m_tDepth;
	long long int		m_lltUniquePepSQNum;
	long long int		m_lltPepSQNum;
	
//	size_t 	_CountPepStructSize();
//	size_t 	_CountIdxStructSize();
//	
	void 	_WriteIdx(FILE *fIdx, size_t tMass, long llPos);
	void	_ReadIdx(FILE *fIdx, size_t &tMass, long &llPos);

	void	_WriteTmpPep(FILE *fIdx, SPS_PEP_STRUCT & sDatExTagBlock);
	void 	_ReadTmpPep(FILE *fIdx, SPS_PEP_STRUCT & sDatExTagBlock);

	void	_WriteTargetPep(FILE *fIdx, SPS_PEP_STRUCT & sDatExTagBlock);
	void	_WriteTargetPro(FILE *fDAT, FILE *fPRO);
//	void 	_ReadTargetPep(FILE *fIdx, SPS_PEP_STRUCT & sDatExTagBlock,int nIsTargeFile );

//	string 	_GenerTargeFileName(string strType = "");
//	void 	_GenerTmpFileName(string &strPepDatName);

	string	_GenerFileName(size_t tFileID, string strType);

	bool 	_SPS_PEP_STRUCT_EQU(const SPS_PEP_STRUCT & PepFir, const SPS_PEP_STRUCT &PepSec);

	void    _InitMeta();
	void	_InitOtherPara(CPepFilter &filter);

	void 	_SetWorkDir(const char * szDestDir);
	
};
}
#endif /*SPSPEPTIDEEXECUTER_H_*/
