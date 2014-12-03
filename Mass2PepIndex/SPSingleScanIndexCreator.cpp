
#include "SPSingleScanIndexCreator.h"

namespace Mass2PepIndex{
CSPSingleScanIndexCreator::CSPSingleScanIndexCreator()
{
}

CSPSingleScanIndexCreator::~CSPSingleScanIndexCreator()
{
}

void CSPSingleScanIndexCreator::Close()
{
	
}

void CSPSingleScanIndexCreator::Init(const PINDEX_HEAD &pIndexHead, const PINDEX_ITEM & pIndexItem, bool bMono) 
{
	try
	{
		m_PIndexHead = pIndexHead;
		m_PIndexItem = pIndexItem;

		if(bMono)
		{
			m_PIndexItem.bMono = true;
			m_PIndexItem.bAvrg = false;
		}
		else
		{
			m_PIndexItem.bMono = false;
			m_PIndexItem.bAvrg = true;		
		}
		
		m_peptideExecuter.Init(m_PIndexHead, m_PIndexItem);

		ProteinIndex::ProteinRandomReaderType eType = (ProteinRandomReaderType)(m_PIndexItem.bProIndexMap);
		
		ProteinIndex::CProteinRandomReaderFactory ProRandomFactory;		
		m_proteinExecuter = ProRandomFactory.GetRandomSQReader(eType);
		m_proteinExecuter->Open(m_PIndexHead.strOutpath,m_PIndexItem.strDBName);
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSingleScanIndexCreator", "Init");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSingleScanIndexCreator", "Init");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CSPSingleScanIndexCreator::WriteIndex()
{
	try
	{
		long llMaxPepNum = m_peptideExecuter.GetMaxPepNum();
		size_t tProNum = m_proteinExecuter->GetProNum();
		string strProSQ = "";
				
		m_peptideExecuter.EmptyPepRecord(0);

		for(size_t tProID = 0; tProID < tProNum; ++tProID)
		{
			size_t tPos = m_proteinExecuter->GetByID( strProSQ, tProID );

			long llTmpNum = DEFAULT_MAX_SQ2PEP_TIME * strProSQ.size();
			long llUsedNum = m_peptideExecuter.GetUsedPepNum();

			if(llUsedNum && llUsedNum + llTmpNum > llMaxPepNum)	
			{
				m_peptideExecuter.WritePepIndex(tProID);
				m_peptideExecuter.EmptyPepRecord(tProID);//the range is [...)
				
				cout << "------------------------------------" << endl;
				cout << "Finished:" << tProID << " at all of " << tProNum << endl;
				cout << "------------------------------------" << endl;
			}
			long tLast = m_peptideExecuter.GetUsedPepNum();
			m_peptideExecuter.AppendPepIdxBySQ(strProSQ, tProID, tPos);
	
			for(long tp = m_peptideExecuter.m_vsDatExTagBlock.size() -1; tp >= tLast; --tp)
			{
				m_peptideExecuter.m_vsDatExTagBlock[tp].cEnd = 0;
				size_t pL = m_peptideExecuter.m_vsDatExTagBlock[tp].tPos - tPos;
				if(0 == pL || (1 == pL && 'M' == strProSQ[0])) 
				{
					m_peptideExecuter.m_vsDatExTagBlock[tp].cEnd |=1;
				}
				
				if(m_peptideExecuter.m_vsDatExTagBlock[tp].tPos + m_peptideExecuter.m_vsDatExTagBlock[tp].cLen == tPos + strProSQ.size())
				{
					m_peptideExecuter.m_vsDatExTagBlock[tp].cEnd |=2;
				}
			}
			
			llUsedNum += llTmpNum;
			
		}
	
		if(0 != m_peptideExecuter.GetUsedPepNum())//any Peptide left
		{
			m_peptideExecuter.WritePepIndex(tProNum);
			m_peptideExecuter.EmptyPepRecord(tProNum);//the range is [...)		

			cout << "------------------------------------" << endl;
			cout << "Finished:" << tProNum << " at all of " << tProNum << endl;
			cout << "------------------------------------" << endl;
		}
	
		cout << "************************************************************************" << endl;
		cout << "It is merging and may take long time for large database, please wait.. "  << endl;
		cout << "************************************************************************" << endl;
		
		if(m_PIndexItem.bPep2Pro)
			m_peptideExecuter.MergePepIndexWithPro2Pep();
		else
			m_peptideExecuter.MergePepIndex();
		m_peptideExecuter.WriteMeta();
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CSPSingleScanIndexCreator", "WriteIndex");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CSPSingleScanIndexCreator", "WriteIndex");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

}

}
