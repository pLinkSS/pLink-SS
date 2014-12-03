#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
//#include "ACAlgorithmClass.h"
#include "test.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "../ProteinIndex/ProteinIndex.h"
#include "../ProteinIndex/ProteinReaderFactory.h"
#include "../ProteinIndex/ProteinWriterFactory.h"
#include "../ProteinIndex/ProteinReader.h"
#include "../ProteinIndex/ProteinWriter.h"
#include "XLinkACInfer.h"
using namespace ProteinIndex;
using namespace proteomics_sdk;
using namespace std;

static bool EVLesser( const pair<double,bool>& elem1, const pair<double,bool>& elem2 )
{
	return elem1.first < elem2.first;
}

bool CXLinkACInfer::Less(const CAssignedProtein & a, const CAssignedProtein & b)
{
	
	if(a.m_nPeptide == b.m_nPeptide)
	{
		//TODO: add the coverage compare function here;


		if(a.m_vpContainPep.size() == b.m_vpContainPep.size())
		{

			if(fabs(a.GetPepAverRank() - b.GetPepAverRank()) <= EPS)//float equal
			{
				return true;//todo return the shorter one
			}
			return a.GetPepAverRank() < b.GetPepAverRank();//a.GetTotalEValue() < b.GetTotalEValue();
		}
		return a.m_vpContainPep.size() > b.m_vpContainPep.size();
	}
	return a.m_nPeptide > b.m_nPeptide;
}

bool CXLinkACInfer::PeptideLess(const pair<const CXLinkPepResult *, size_t> & a, const pair<const CXLinkPepResult *, size_t> & b)
{
	return a.first->m_lfEvalue < b.first->m_lfEvalue;
}


vector<CAssignedProtein> & CXLinkACInfer::Infer(vector<CXLinkMatchResult> & vResults, const vector<CSpectrum> & vSpectra,
		FILTER_CRITERIA_INFO & cs)
{
	m_pTrace->Debug("in CXLinkACInfer::Infer");

	size_t tTotalPep = 0;
	for(size_t i = 0;i < vResults.size();++i)
	{
		tTotalPep += 2*vResults[i].m_vPeptideResults.size();
	}
	
	vector<string> vstrSrc;
	vstrSrc.reserve(tTotalPep);
	vector<size_t> vId1;
	vId1.reserve(tTotalPep);
	vector<size_t> vId2;
	vId2.reserve(tTotalPep);
	vector<size_t> vId3;
	vId3.reserve(tTotalPep);
	
	//fill in the vector
	for(size_t i = 0;i < vResults.size();++i)
	{
		for(size_t j = 0;j < vResults[i].m_vPeptideResults.size();++j)
		{
			vstrSrc.push_back(vResults[i].m_vPeptideResults[j].m_AlphaPeptide.m_szSequence);
			vId1.push_back(i);
			vId2.push_back(j);
			vId3.push_back(0);
			
			if(vResults[i].m_vPeptideResults[j].m_bPair)
			{
				vstrSrc.push_back(vResults[i].m_vPeptideResults[j].m_BetaPeptide.m_szSequence);	
				vId1.push_back(i);
				vId2.push_back(j);
				vId3.push_back(1);
			}
			
			vResults[i].m_vPeptideResults[j].m_vAlphaProteinID.clear();
			vResults[i].m_vPeptideResults[j].m_vAlphaProteinAC.clear();
			vResults[i].m_vPeptideResults[j].m_vBetaProteinID.clear();
			vResults[i].m_vPeptideResults[j].m_vBetaProteinAC.clear();
			vResults[i].m_vPeptideResults[j].m_vAlphaProteinSite.clear();
			vResults[i].m_vPeptideResults[j].m_vBetaProteinSite.clear();
		}
	}
	
	AC ac;
	ac.Init_(m_pCond);
	
//	for (int i = 0; i < vstrSrc.size(); i++)
//		cout<<vstrSrc[i]<<endl;

	ac.Create_(vstrSrc);
	
	//READE index.AC FILE
	CDBConf db(m_pCond->m_strDBConfPath.c_str());
	size_t tProNum = 0;

	for(int nDBid = 0 ;nDBid < m_pCond->m_vSelectedDBName.size();++ nDBid)
	{
		string strDBPath = db.GetPath(m_pCond->m_vSelectedDBName[nDBid].c_str(), m_pCond->m_SelectedEnzyme.m_strName.c_str());
		CProteinRandomReaderFactory factory;
		// add by emily : protein index
		CRandomReader *pIndexACLoader = factory.GetRandomACReader(Protein_Random_Index_Map);
		CRandomReader *pIndexSQLoader = factory.GetRandomSQReader(Protein_Random_Index_Map);
	
	//	CRandomReader *pIndexACLoader = factory.GetRandomACReader(Protein_Random_Index_Disk);//czhou
	//	CRandomReader *pIndexSQLoader = factory.GetRandomSQReader(Protein_Random_Index_Disk);
		try
		{
			// add by emily : protein index
			
			// todo : modify for PEPTIDE_PAIR_DOUBLE_ONLINE
			pIndexACLoader->Open(strDBPath, m_pCond->m_vSelectedDBName[nDBid]);
			pIndexSQLoader->Open(strDBPath, m_pCond->m_vSelectedDBName[nDBid]);
		}
		catch(runtime_error & e)
		{
			CErrInfo info("CXLinkACferor", "Run", "in the function ProteinIndex::Open.");
			info.Append("strDBPath=" + strDBPath);
			info.Append("DBName=" + m_pCond->m_vSelectedDBName[nDBid]);
			throw runtime_error(info.Get(e).c_str());
		}
		catch(...)
		{
			CErrInfo info("CXLinkACferor", "Run", "caught an unknown exception in the function ProteinIndex::Open.");
			info.Append("strDBPath=" + strDBPath);
			info.Append("DBName=" + m_pCond->m_vSelectedDBName[nDBid]);
			throw runtime_error(info.Get().c_str());		
		}
		tProNum = pIndexACLoader->GetProNum();
		
	char szDebug[50];
	sprintf(szDebug, "pro number: %d", tProNum);
	m_pTrace->Debug(szDebug);

	for(size_t i = 0;i < tProNum;++i)
	{
	
//		string bs = "begin";
//		DDebug(bs);
		
		string strAC;
		// add by emily : protein index
		pIndexACLoader->GetByID(strAC, i);
//		DDebug(strAC);
		
		string strPro;
		// add by emily : protein index
		pIndexSQLoader->GetByID(strPro, i);
//		DDebug(strPro);
		
		
		set<pair<size_t,size_t> > setIDs;
		setIDs.clear();
		//AC::Query(strPro, setIDs);
		
		ac.Query_(strPro, setIDs);
		
		for(set<pair<size_t,size_t> >::iterator it = setIDs.begin();
			it != setIDs.end();
			++it)
		{
			//vResults[vId1[*it]].m_vPeptideResults[vId2[*it]].m_vProteinID.push_back(i);
			size_t tSeqId;
			size_t tProSite;
			tSeqId = (*it).first;
			tProSite = (*it).second;
			tProSite = tProSite + 1 - vstrSrc[tSeqId].length();
			
			if(vId3[tSeqId] == 0)
			{
				vResults[vId1[tSeqId]].m_vPeptideResults[vId2[tSeqId]].m_vAlphaProteinID.push_back(i);
				vResults[vId1[tSeqId]].m_vPeptideResults[vId2[tSeqId]].m_vAlphaProteinAC.push_back(strAC);
				vResults[vId1[tSeqId]].m_vPeptideResults[vId2[tSeqId]].m_vAlphaProteinSite.push_back(tProSite);
			}
			else
			{
				vResults[vId1[tSeqId]].m_vPeptideResults[vId2[tSeqId]].m_vBetaProteinID.push_back(i);
				vResults[vId1[tSeqId]].m_vPeptideResults[vId2[tSeqId]].m_vBetaProteinAC.push_back(strAC);
				vResults[vId1[tSeqId]].m_vPeptideResults[vId2[tSeqId]].m_vBetaProteinSite.push_back(tProSite);
			}
		}
	}
	
//	string bs = "end";
//	DDebug(bs);
	
	// add by emily : protein index
	pIndexACLoader->Close();
	pIndexSQLoader->Close();
	delete pIndexACLoader;
	delete pIndexSQLoader;
	
	}
	m_pTrace->Debug("no infer...");
	
	m_vFoundProteins.clear();
	return m_vFoundProteins;
	
}

double CXLinkACInfer::Filter(const vector<CXLinkMatchResult> & vResults, FILTER_CRITERIA_INFO & cs)
{

	return 0.0;
}
