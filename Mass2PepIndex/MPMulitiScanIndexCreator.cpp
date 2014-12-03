#include "Mass2PepIndex.h"
#include "MassRangeCounter.h"
#include "MassRangeWriter.h"


#include "MPMulitiScanIndexCreator.h"

namespace Mass2PepIndex{
CMPMulitiScanIndexCreator::CMPMulitiScanIndexCreator()
{
}

CMPMulitiScanIndexCreator::~CMPMulitiScanIndexCreator()
{
}

void CMPMulitiScanIndexCreator::Close()
{
	
}

void CMPMulitiScanIndexCreator::Init(const PINDEX_HEAD &pIndexHead, const PINDEX_ITEM & pIndexItem) 
{
	m_PIndexHead = pIndexHead;
	m_PIndexItem = pIndexItem;
}

void CMPMulitiScanIndexCreator::WriteIndex()
{
	_ScanMassRange();
	
	_Dofork();
	
	_WriteMeta();
}

void CMPMulitiScanIndexCreator::_ScanMassRange()
{	
	CMassRangeCounter massRangeCounter;
	m_vdfMassRange = massRangeCounter.GetMassRange();
}

void CMPMulitiScanIndexCreator::_Excute(double lfLowMass, double lfHighMass)
{
	CMassRangeWriter massRangeWriter;
	massRangeWriter.WriteIndex( lfLowMass,  lfHighMass);
//	META_ITEM metaItem;
	massRangeWriter.GetMetaITEM();
}

void CMPMulitiScanIndexCreator::_Dofork()
{
//	int nProNum = m_vdfMassRange.size() -1;
//	int pid[nProNum];
//	
//	string strTmpFileName = "TmpFile";
//	char strNum[10];//restore the kth of the tmp meta file	
	
//	for(size_t tp=0; tp<nProNum; ++tp)
//	{		
//		sprintf(strNum,"%d", tp);
//		string strFileName = strTmpFileName + * new string(strNum);		//临时文件名，文件路径是否要指定？
//		
////		pid[tp]	= fork();	//*********************
//		
//		if(0 == pid[tp])
//		{
//			META_ITEM metaItem; 
//			_Excute(m_vdfMassRange[tp], m_vdfMassRange[tp + 1]);
////			 peptide2ProDatCreater.GetMetaITEM();
//		
//			FILE *finf;
//			finf = fopen( strFileName.c_str(), "w" );			
//			fprintf(finf,"%d %d %d %d %d",metaItem.tPepSQNum, metaItem.tUniquePepSQNum, metaItem.tMinMass,  metaItem.tMaxMass, metaItem.tUniqueMassNum);
//			
//			exit(0);//**************************************	
//		}

//		MassRangeWriter::WriteIndex(m_vdfMassRange, );
//		PrintToFile();
//	}
}

void CMPMulitiScanIndexCreator::_WriteMeta()
{
	
}

}

