#include "../include/sdk.h"
#include "SPSPeptideExecuter.h"

namespace Mass2PepIndex{

struct SortPep 
{
	bool operator()(const PEP_INFO_EX & PepFir, const PEP_INFO_EX & PepSec)const
	{
		if(PepFir.tMass != PepSec.tMass) return PepFir.tMass < PepSec.tMass;
//		if(PepFir.lfTag != PepSec.lfTag) return PepFir.lfTag < PepSec.lfTag;


		std::string strSQ1, strSQ2;
		
//		proLoaderInPep2Pro->ReadPepSQ(strSQ1,PepFir.tPos,PepFir.cLen);
//		proLoaderInPep2Pro->ReadPepSQ(strSQ2,PepSec.tPos,PepSec.cLen);

		if(strSQ1!=strSQ2) return strSQ1<strSQ2;

		return PepFir.tPos < PepSec.tPos;
	}
};

struct SortPepTag 
{
	bool operator()(const PEP_INFO_EX_TAG & PepFir, const PEP_INFO_EX_TAG & PepSec)const
	{
		if(PepFir.tMass != PepSec.tMass) return PepFir.tMass < PepSec.tMass;
		if(PepFir.lfTag != PepSec.lfTag) return PepFir.lfTag < PepSec.lfTag;
		return PepFir.tPos<PepSec.tPos;
	}
};

bool SPS_PEP_STRUCT_CMP(const SPS_PEP_STRUCT & PepFir, const SPS_PEP_STRUCT &PepSec)
{
	if(PepFir.tMass != PepSec.tMass) return PepFir.tMass < PepSec.tMass;
	if(PepFir.lfTag != PepSec.lfTag) return PepFir.lfTag < PepSec.lfTag;
	if(PepFir.tProID != PepSec.tProID) return PepFir.tProID < PepSec.tProID;
	return PepFir.tPos<PepSec.tPos;
}

bool HEAP_STRUCT_CMP(const HEAP_ITEM & PepFir, const HEAP_ITEM &PepSec)
{
	return SPS_PEP_STRUCT_CMP(PepSec.sDatExTagBlock,PepFir.sDatExTagBlock);
}


CSPSPeptideExecuter::CSPSPeptideExecuter()
:m_tMultplier(DEFAULT_MULTPLIER)
,m_tPep2ProSize(1097152)
,m_tBitNum(DEFAULT_BIT_NUM)
,m_tMaxMemorySize(DEFAULT_MEMORY_SIZE)
,m_dfMassRange(DEFAULT_PEP_RANGE)
,m_tOutTxt(0)
,m_fPep2Pro(0)
,m_fDebug(0)
,m_strDebugFileName("CSPSingleScanIndexCreator.txt")
,m_tDepth(1)
{
	m_pFilter = NULL;
	m_tCleaveWay = 0;
}

CSPSPeptideExecuter::~CSPSPeptideExecuter()
{
}

void CSPSPeptideExecuter::Close()
{
	
}

void CSPSPeptideExecuter::Init(const PINDEX_HEAD& pIndexHead, const PINDEX_ITEM & pIndexItem) 
{
	try
	{
		m_PIndexHead = pIndexHead;
		m_PIndexItem = pIndexItem;
				
		m_vsTmpFileName.clear();
		
		CEnzyme	enzyme;
		
		CPepFilter filter;	
		filter.SetMaxMissedClvg(m_PIndexItem.nMaxMissSite);
		filter.SetMaxPepLength(m_PIndexItem.tMaxPepLength);
		filter.SetMinPepLength(m_PIndexItem.tMinPepLength);
		filter.SetMaxPepMass(m_PIndexItem.tMaxPepMass);
		filter.SetMinPepMass(m_PIndexItem.tMinPepMass);
		m_pFilter = new CPepFilter();
		*m_pFilter = filter;
		
		m_tMaxMemorySize = m_PIndexItem.tMaxMemSize;
		m_dfMassRange = (double)(m_PIndexItem.tMassRange);
		_SetWorkDir(m_PIndexHead.strOutpath.c_str());
		m_tOutTxt = m_PIndexItem.tOutTxt;
		
		string strEnzymeName = m_PIndexItem.vtEnzymeNames[0];
		CEnzymeConf enzyme_conf(m_PIndexHead.stEnzymeList);
		if(!enzyme_conf.GetEnzyme(strEnzymeName, enzyme))
		{
			proteomics_sdk::CErrInfo err_info("CIndexer", "CreatePepMassIndexes","Can't get enzyme!");
			throw runtime_error(err_info.Get().c_str());
		}
		
		m_tCleaveWay = m_PIndexItem.nCleaveWay;
		m_Enzyme = enzyme;
		m_strAAListPath = m_PIndexHead.strAAList;	
		m_PepFunc.Init(m_strAAListPath);
		
		_InitMeta();
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "Init");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "Init");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CSPSPeptideExecuter::EmptyPepRecord(size_t tProID)
{
	m_vsDatExTagBlock.clear();
	m_tBegProID = tProID;
}

size_t CSPSPeptideExecuter::GetMaxPepNum()
{
	return m_tMaxMemorySize / sizeof(SPS_PEP_STRUCT);
}

void CSPSPeptideExecuter::GetPepInfor(size_t& t,SPS_PEP_STRUCT& pepInfor)
{
	pepInfor = m_vsDatExTagBlock[t];
}

size_t CSPSPeptideExecuter::GetUsedPepNum()
{
	return m_vsDatExTagBlock.size();
}

void CSPSPeptideExecuter::AppendPepIdxBySQ(string &strProSQ, size_t tProID, size_t tPos)
{
	if(0 == m_tCleaveWay)
	{
		m_Digester.DigestPros2(strProSQ, m_Enzyme,m_PepFunc, m_CMetaCreater.GetIsMono(), m_pFilter->GetMaxMissedClvg(), m_pFilter, m_vsDatExTagBlock, tPos, tProID, m_CMetaCreater.m_vstMetaItems.size());	
	}
	else if(1 == m_tCleaveWay)
	{
		m_Digester.SemiDigestPros1(strProSQ, m_Enzyme,m_PepFunc, m_CMetaCreater.GetIsMono(), m_pFilter->GetMaxMissedClvg(), m_pFilter, m_vsDatExTagBlock, tPos, tProID, m_CMetaCreater.m_vstMetaItems.size());
	}
	else if(2 == m_tCleaveWay)
	{
		m_Digester.NoneEnzymeDigestPros1(strProSQ, m_PepFunc, m_CMetaCreater.GetIsMono(), m_pFilter->GetMaxMissedClvg(), m_pFilter, m_vsDatExTagBlock, tPos, tProID, m_CMetaCreater.m_vstMetaItems.size());
	}		
}

void CSPSPeptideExecuter::GetUsefulPep()
{
	try
	{
		if(0 == m_vsDatExTagBlock.size()) return;
		sort(m_vsDatExTagBlock.begin(), m_vsDatExTagBlock.end(), SPS_PEP_STRUCT_CMP);
		
		size_t tUsefulNum = 1;
		for(size_t t = 1; t < m_vsDatExTagBlock.size(); ++t)
		{
			if(false == _SPS_PEP_STRUCT_EQU(m_vsDatExTagBlock[t], m_vsDatExTagBlock[t -1]))
			{
				m_vsDatExTagBlock[tUsefulNum++] = m_vsDatExTagBlock[t];
			}
		}
		
		if(tUsefulNum != m_vsDatExTagBlock.size())
		{
			m_vsDatExTagBlock.erase(m_vsDatExTagBlock.begin()+tUsefulNum, m_vsDatExTagBlock.end());
		}
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "GetUsefulPep");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "GetUsefulPep");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

}

void CSPSPeptideExecuter::WritePepIndex(size_t tProID)
{
	try
	{
		sort(m_vsDatExTagBlock.begin(), m_vsDatExTagBlock.end(), SPS_PEP_STRUCT_CMP);

		string strPepDatName = _GenerFileName(m_vsTmpFileName.size(), PEP_INFO_TMP);
		m_vsTmpFileName.push_back(strPepDatName);
		
		FILE  * fDat = fopen(strPepDatName.c_str(), "wb");	
		if(!fDat)
		{
			CErrInfo info("CSPSPeptideExecuter", "WritePepIndex", "Can't Create Tmp File!");
			cerr << info.Get() << endl;
			throw runtime_error(info.Get().c_str());
		}
		fseek(fDat, 0, SEEK_SET);
		size_t tUniquePepSQNum = 0;
		if(0 != m_vsDatExTagBlock.size())
		{
			_WriteTmpPep(fDat,m_vsDatExTagBlock[0]);	
			tUniquePepSQNum = 1;
		}
		
		for(size_t t = 1; t < m_vsDatExTagBlock.size(); ++t)
		{
//			if(false == _SPS_PEP_STRUCT_EQU(m_vsDatExTagBlock[t], m_vsDatExTagBlock[t -1]))
//			{
				_WriteTmpPep(fDat,m_vsDatExTagBlock[t]);
				++tUniquePepSQNum;
//			}
		}

		fflush(fDat);
		fclose(fDat);
		
		cout << "------------------------------------" << endl;
		cout << "Generated peptide number:" << m_vsDatExTagBlock.size() << endl;
//		cout << "Unique Peptide Num:" << tUniquePepSQNum << endl;
		cout << "------------------------------------" << endl;
		
		m_CMetaCreater.m_stMetaHead.tPepSQNum += m_vsDatExTagBlock.size();
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "WritePepIndex");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "WritePepIndex");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CSPSPeptideExecuter::MergePepIndex()
{
	if(0 == m_vsTmpFileName.size()) return;
	 
	try
	{
		FILE **m_fTmpFile = new FILE*[m_vsTmpFileName.size()];
		
		HEAP_ITEM  *heap = new HEAP_ITEM[m_vsTmpFileName.size()];
		size_t tCnt = 0;
		
		for(size_t t = 0; t < m_vsTmpFileName.size(); ++t)
		{
			m_fTmpFile[t] = fopen(m_vsTmpFileName[t].c_str(), "rb");
			fseek(m_fTmpFile[t],0,SEEK_SET);
			
			heap[tCnt].tSer = t;
			_ReadTmpPep(m_fTmpFile[t], heap[tCnt].sDatExTagBlock);
			push_heap(heap,heap + (++tCnt), HEAP_STRUCT_CMP);
		}		
		
		size_t tFileID = 0;
		FILE *fpIDX = fopen(_GenerFileName(tFileID,PEP_INFO_IDX).c_str(), "wb");
		FILE *fpDAT = fopen(_GenerFileName(tFileID,PEP_INFO_DAT).c_str(), "wb");

		m_tUniquePepSQNum = 0;
		size_t tLastAllMassPepNum = 0;
		size_t tLastUniMassPepNum = 0;
		size_t tMaxPepNum = GetMaxPepNum();
		size_t tLastMass = heap[0].sDatExTagBlock.tMass;
		HEAP_ITEM lHeap;
		lHeap.sDatExTagBlock.tMass = 0;//record of last sDatExTagBlock, and there are no value of 0;**********
		lHeap.sDatExTagBlock.lfTag = -1;
		do
		{
			if(false == _SPS_PEP_STRUCT_EQU(lHeap.sDatExTagBlock, heap[0].sDatExTagBlock))
			{
				if(tLastUniMassPepNum + DEFAULT_MAX_SQ2PEP_TIME >= tMaxPepNum)//**
				{						
					fclose(fpIDX);
					fclose(fpDAT);
					++tFileID;
					fpIDX = fopen(_GenerFileName(tFileID,PEP_INFO_IDX).c_str(), "wb");
					fpDAT = fopen(_GenerFileName(tFileID,PEP_INFO_DAT).c_str(), "wb");

					m_stMetaItem.tMinMass = tLastMass;
					m_stMetaItem.tMaxMass = heap[0].sDatExTagBlock.tMass;
					m_stMetaItem.tPepSQNum = tLastAllMassPepNum;
					m_stMetaItem.tUniquePepSQNum = tLastUniMassPepNum;
					m_stMetaItem.tUniqueMassNum = 0;
					
					m_CMetaCreater.m_vstMetaItems.push_back(m_stMetaItem);			
					
					tLastMass = heap[0].sDatExTagBlock.tMass;
					
					tLastAllMassPepNum = 0;
					tLastUniMassPepNum = 0;
				}
				
				if(heap[0].sDatExTagBlock.tMass != lHeap.sDatExTagBlock.tMass)
				{
					_WriteIdx(fpIDX, heap[0].sDatExTagBlock.tMass, ftell(fpDAT));
				}
				
				_WriteTargetPep(fpDAT, heap[0].sDatExTagBlock);	
				lHeap = heap[0];
				
				++tLastUniMassPepNum;
				++m_tUniquePepSQNum;
			}
			++tLastAllMassPepNum;
			
			HEAP_ITEM tHeap = heap[0];
			pop_heap(heap, heap + (tCnt--), HEAP_STRUCT_CMP);
			
			if(0 == feof(m_fTmpFile[tHeap.tSer]))
			{
				heap[tCnt].tSer = tHeap.tSer;
				_ReadTmpPep(m_fTmpFile[tHeap.tSer], heap[tCnt].sDatExTagBlock);
				push_heap(heap,heap + (++tCnt), HEAP_STRUCT_CMP);
			}
			else
			{
//				fclose(m_fTmpFile[tCurPep]);
//				remove(m_vsTmpFileName[tCurPep].c_str());
			}
			
			if(!tCnt)
			{
				m_stMetaItem.tMinMass = tLastMass;
				m_stMetaItem.tMaxMass = tHeap.sDatExTagBlock.tMass;
				m_stMetaItem.tPepSQNum = tLastAllMassPepNum;
				m_stMetaItem.tUniquePepSQNum = tLastUniMassPepNum;
				m_stMetaItem.tUniqueMassNum = 0;
				
				m_CMetaCreater.m_vstMetaItems.push_back(m_stMetaItem);
			}
		}while(tCnt);	
		
		m_CMetaCreater.m_stMetaHead.tUniquePepSQNum = m_tUniquePepSQNum;
		m_CMetaCreater.m_stMetaHead.tIdxNum = m_CMetaCreater.m_vstMetaItems.size();
		
		for(size_t t = 0; t < m_vsTmpFileName.size(); ++t)
		{
			fclose(m_fTmpFile[t]);
			remove(m_vsTmpFileName[t].c_str());
		}
	
		fclose(fpIDX);
		fclose(fpDAT);

		delete[] heap;
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "MergePepIndex");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "MergePepIndex");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CSPSPeptideExecuter::MergePepIndexWithPro2Pep()
{
	if(0 == m_vsTmpFileName.size()) return;
	 
	try
	{
		FILE **m_fTmpFile = new FILE*[m_vsTmpFileName.size()];
		
		HEAP_ITEM  *heap = new HEAP_ITEM[m_vsTmpFileName.size()];
		size_t tCnt = 0;
		
		for(size_t t = 0; t < m_vsTmpFileName.size(); ++t)
		{
			m_fTmpFile[t] = fopen(m_vsTmpFileName[t].c_str(), "rb");
			fseek(m_fTmpFile[t],0,SEEK_SET);
			
			heap[tCnt].tSer = t;
			_ReadTmpPep(m_fTmpFile[t], heap[tCnt].sDatExTagBlock);
			push_heap(heap,heap + (++tCnt), HEAP_STRUCT_CMP);
		}		
		
		size_t tFileID = 0;
		FILE *fpIDX = fopen(_GenerFileName(tFileID,PEP_INFO_IDX).c_str(), "wb");
		FILE *fpDAT = fopen(_GenerFileName(tFileID,PEP_INFO_DAT).c_str(), "wb");
		FILE *fpPID = fopen(_GenerFileName(tFileID,PEP_INFO_PID).c_str(), "wb");
		FILE *fpPOS = fopen(_GenerFileName(tFileID,PEP_INFO_POS).c_str(), "wb");
		m_tUniquePepSQNum = 0;
		size_t tLastAllMassPepNum = 0;
		size_t tLastUniMassPepNum = 0;
		size_t tMaxPepNum = GetMaxPepNum();
		size_t tLastMass = heap[0].sDatExTagBlock.tMass;
		HEAP_ITEM lHeap;
		lHeap.sDatExTagBlock.tMass = 0;//record of last sDatExTagBlock, and there are no value of 0;**********
		lHeap.sDatExTagBlock.lfTag = -1;
		m_vpPro.clear();
		
		do
		{
			if(false == _SPS_PEP_STRUCT_EQU(lHeap.sDatExTagBlock, heap[0].sDatExTagBlock))
			{
				if(tLastUniMassPepNum + DEFAULT_MAX_SQ2PEP_TIME >= tMaxPepNum)//**
				{						
					fclose(fpIDX);
					fclose(fpDAT);
					_WriteTargetPro(fpPOS, fpPID);
					fclose(fpPID);
					fclose(fpPOS);
					++tFileID;
					fpIDX = fopen(_GenerFileName(tFileID,PEP_INFO_IDX).c_str(), "wb");
					fpDAT = fopen(_GenerFileName(tFileID,PEP_INFO_DAT).c_str(), "wb");
					fpPID = fopen(_GenerFileName(tFileID,PEP_INFO_PID).c_str(), "wb");
					fpPOS = fopen(_GenerFileName(tFileID,PEP_INFO_POS).c_str(), "wb");					
					m_stMetaItem.tMinMass = tLastMass;
					m_stMetaItem.tMaxMass = heap[0].sDatExTagBlock.tMass;
					m_stMetaItem.tPepSQNum = tLastAllMassPepNum;
					m_stMetaItem.tUniquePepSQNum = tLastUniMassPepNum;
					m_stMetaItem.tUniqueMassNum = 0;
					
					m_CMetaCreater.m_vstMetaItems.push_back(m_stMetaItem);			
					
					tLastMass = heap[0].sDatExTagBlock.tMass;
					
					tLastAllMassPepNum = 0;
					tLastUniMassPepNum = 0;
				}
				
				if(heap[0].sDatExTagBlock.tMass != lHeap.sDatExTagBlock.tMass)
				{
					_WriteIdx(fpIDX, heap[0].sDatExTagBlock.tMass, ftell(fpDAT));
				}
				
				if(tLastUniMassPepNum )
					_WriteTargetPro(fpPOS, fpPID);
				
				m_vpPro.clear();
				m_vpPro.push_back(heap[0].sDatExTagBlock.tProID);

				_WriteTargetPep(fpDAT, heap[0].sDatExTagBlock);
								
				lHeap = heap[0];
				
				++tLastUniMassPepNum;
				++m_tUniquePepSQNum;
			}
			else 
			{
				if(heap[0].sDatExTagBlock.tProID != m_vpPro[m_vpPro.size()-1]) m_vpPro.push_back(heap[0].sDatExTagBlock.tProID);				
			}
			
			++tLastAllMassPepNum;
			
			HEAP_ITEM tHeap = heap[0];
			pop_heap(heap, heap + (tCnt--), HEAP_STRUCT_CMP);
			
			if(0 == feof(m_fTmpFile[tHeap.tSer]))
			{
				heap[tCnt].tSer = tHeap.tSer;
				_ReadTmpPep(m_fTmpFile[tHeap.tSer], heap[tCnt].sDatExTagBlock);
				push_heap(heap,heap + (++tCnt), HEAP_STRUCT_CMP);
			}
			else
			{
//				fclose(m_fTmpFile[tCurPep]);
//				remove(m_vsTmpFileName[tCurPep].c_str());
			}
			
			if(!tCnt)
			{
				m_stMetaItem.tMinMass = tLastMass;
				m_stMetaItem.tMaxMass = tHeap.sDatExTagBlock.tMass;
				m_stMetaItem.tPepSQNum = tLastAllMassPepNum;
				m_stMetaItem.tUniquePepSQNum = tLastUniMassPepNum;
				m_stMetaItem.tUniqueMassNum = 0;
				
				m_CMetaCreater.m_vstMetaItems.push_back(m_stMetaItem);
			}
		}while(tCnt);	
		
		_WriteTargetPro(fpPOS, fpPID);
		m_CMetaCreater.m_stMetaHead.tUniquePepSQNum = m_tUniquePepSQNum;
		m_CMetaCreater.m_stMetaHead.tIdxNum = m_CMetaCreater.m_vstMetaItems.size();
		
		for(size_t t = 0; t < m_vsTmpFileName.size(); ++t)
		{
			fclose(m_fTmpFile[t]);
			remove(m_vsTmpFileName[t].c_str());
		}
	
		fclose(fpIDX);
		fclose(fpDAT);
		fclose(fpPID);
		fclose(fpPOS);		
		
		delete[] heap;
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "MergePepIndex");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "MergePepIndex");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CSPSPeptideExecuter::WriteMeta()
{
	m_CMetaCreater.WriteMeta();
}

void CSPSPeptideExecuter::_WriteIdx(FILE *fIdx, size_t tMass, long llPos)
{	
	if(!fIdx)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteIdx", "File pointer should not be none!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}	
	
	try
	{
		fwrite((char *)&tMass, sizeof(size_t), 1, fIdx);
		fwrite((char *)&llPos, sizeof(long), 1, fIdx);	
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteIdx");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteIdx");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CSPSPeptideExecuter::_ReadIdx(FILE *fIdx, size_t &tMass, long &llPos)
{	
	if(!fIdx)
	{
		CErrInfo info("CSPSPeptideExecuter", "_ReadIdx", "File pointer should not be none!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}	
	
	try{
		fread((char *)&tMass, sizeof(size_t), 1, fIdx);
		fread((char *)&llPos, sizeof(long), 1, fIdx);	
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "_ReadIdx");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "_ReadIdx");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CSPSPeptideExecuter::_WriteTmpPep(FILE *fIdx, SPS_PEP_STRUCT &  sDatExTagBlock)
{
	if(!fIdx)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteTmpPep", "File pointer should not be none!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}	
	
	try
	{
		fwrite((char *)&sDatExTagBlock.tMass, sizeof(size_t), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.lfTag, sizeof(double), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.tProID, sizeof(size_t), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.tPos, sizeof(size_t), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.cLen, sizeof(unsigned char), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.cMiss, sizeof(unsigned char), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.cEnd, sizeof(unsigned char), 1, fIdx);
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteTmpPep");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteTmpPep");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CSPSPeptideExecuter::_ReadTmpPep(FILE *fIdx, SPS_PEP_STRUCT & sDatExTagBlock )
{
	if(!fIdx)
	{
		CErrInfo info("CSPSPeptideExecuter", "_ReadTmpPep", "File pointer should not be none!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

	try
	{
		fread((char *)&sDatExTagBlock.tMass, sizeof(size_t), 1, fIdx);
		fread((char *)&sDatExTagBlock.lfTag, sizeof(double), 1, fIdx);
		fread((char *)&sDatExTagBlock.tProID, sizeof(size_t), 1, fIdx);
		fread((char *)&sDatExTagBlock.tPos, sizeof(size_t), 1, fIdx);
		fread((char *)&sDatExTagBlock.cLen, sizeof(unsigned char), 1, fIdx);
		fread((char *)&sDatExTagBlock.cMiss, sizeof(unsigned char), 1, fIdx);
		fread((char *)&sDatExTagBlock.cEnd, sizeof(unsigned char), 1, fIdx);
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "_ReadTmpPep");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "_ReadTmpPep");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CSPSPeptideExecuter::_WriteTargetPep(FILE *fIdx, SPS_PEP_STRUCT &  sDatExTagBlock)
{
	if(!fIdx)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteTargetPep", "File pointer should not be none!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}	
	
	try
	{
		fwrite((char *)&sDatExTagBlock.tMass, sizeof(size_t), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.tProID, sizeof(size_t), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.tPos, sizeof(size_t), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.cLen, sizeof(unsigned char), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.cMiss, sizeof(unsigned char), 1, fIdx);
		fwrite((char *)&sDatExTagBlock.cEnd, sizeof(unsigned char), 1, fIdx);
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteTargetPep");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteTargetPep");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CSPSPeptideExecuter::_WriteTargetPro(FILE *fPOS, FILE *fPRO)
{
	if(!fPOS || !fPRO)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteTargetPro", "File pointer should not be none!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}	
	
	try
	{
		size_t tmp = ftell(fPRO);
		fwrite((char *)&tmp, sizeof(size_t), 1, fPOS);
		tmp = m_vpPro.size();
		fwrite((char *)&tmp, sizeof(size_t), 1, fPRO);
		for(size_t t = 0; t < m_vpPro.size(); ++t)
		{
			fwrite((char *)&m_vpPro[t], sizeof(size_t), 1, fPRO);
		}
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteTargetPro");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "_WriteTargetPro");		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

string CSPSPeptideExecuter::_GenerFileName(size_t tFileID, string strType)
{
	string strMOrA = ".mono";
	if(!m_PIndexItem.bMono) strMOrA = ".avrg";
	string strFileName = m_strWorkDir + m_CMetaCreater.GetDBName() + "." + m_CMetaCreater.GetEnzymeName() + strMOrA;
	char szProID[10];
	sprintf(szProID,".%d",tFileID);	
	return strFileName + string(szProID) + strType;
}

bool CSPSPeptideExecuter::_SPS_PEP_STRUCT_EQU(const SPS_PEP_STRUCT & PepFir, const SPS_PEP_STRUCT &PepSec)
{
	return (PepFir.tMass == PepSec.tMass &&	PepFir.lfTag == PepSec.lfTag);
}


void CSPSPeptideExecuter::_InitMeta()
{			
	try
	{
		m_CMetaCreater.Load(m_PIndexHead.strOutpath.c_str(), m_PIndexItem.strDBName.c_str(), m_PIndexItem.vtEnzymeNames[0].c_str());
		m_CMetaCreater.m_stMetaHead.nCleaveWay = (long)m_PIndexItem.nCleaveWay;
		m_CMetaCreater.m_stMetaHead.tMultiplier = DEFAULT_MULTPLIER;
		m_CMetaCreater.m_stMetaHead.tPepSQNum = 0;
		m_CMetaCreater.m_stMetaHead.tUniquePepSQNum = 0;
		m_CMetaCreater.m_stMetaHead.tUniqueMassNum = 0;
		m_CMetaCreater.m_stMetaHead.tMinMass = m_PIndexItem.tMinPepMass * m_CMetaCreater.m_stMetaHead.tMultiplier;
		m_CMetaCreater.m_stMetaHead.tMaxMass = m_PIndexItem.tMaxPepMass * m_CMetaCreater.m_stMetaHead.tMultiplier;
		m_CMetaCreater.m_stMetaHead.tMinLength = m_PIndexItem.tMinPepLength;
		m_CMetaCreater.m_stMetaHead.tMaxLength = m_PIndexItem.tMaxPepLength;
		m_CMetaCreater.m_stMetaHead.tIdxNum = 0;//*************
		m_CMetaCreater.m_stMetaHead.tMissCleaveSites = m_PIndexItem.nMaxMissSite;
		m_CMetaCreater.m_stMetaHead.bMono = m_PIndexItem.bMono;
		m_CMetaCreater.m_stMetaHead.bPep2Pro = m_PIndexItem.bPep2Pro;
		
		m_CMetaCreater.m_vstMetaItems.clear();
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "_InitMeta");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "_InitMeta");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void 	CSPSPeptideExecuter::_SetWorkDir(const char * szDestDir)
{
	try
	{
		string strDestDir( szDestDir );
		string::size_type idx = 0;
		idx = strDestDir.find('\\', idx);
		while (idx != string::npos ) {
			strDestDir[idx] = '/';
			idx = strDestDir.find('\\', idx);
		}
	
		idx = strDestDir.find_last_not_of(" ");
		if (strDestDir[idx] != '/' ) 
		{
			m_strWorkDir = strDestDir.substr(0,idx + 1) + "/";
		}
		else
			m_strWorkDir = strDestDir;
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSPeptideExecuter", "SetWorkDir");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSPeptideExecuter", "SetWorkDir");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}
}
