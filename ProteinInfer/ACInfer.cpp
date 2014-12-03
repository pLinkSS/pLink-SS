#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "ACAlgorithm.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "../ProteinIndex/ProteinIndex.h"
#include "../ProteinIndex/ProteinReaderFactory.h"
#include "../ProteinIndex/ProteinWriterFactory.h"
#include "../ProteinIndex/ProteinReader.h"
#include "../ProteinIndex/ProteinWriter.h"
#include "ACInfer.h"
using namespace ProteinIndex;
using namespace proteomics_sdk;
using namespace std;

static bool EVLesser( const pair<double,bool>& elem1, const pair<double,bool>& elem2 )
{
	return elem1.first < elem2.first;
}

//sameset and subset
int CACInfer::_GetRoot(int  root[], int nOrd)
{
	int r = root[nOrd];
	if(r == nOrd || r == -1)
	{
		root[nOrd] = -1;
		return nOrd;
	}
	else
	{
		return root[nOrd] = _GetRoot(root, r);
	}
}

bool CACInfer::_IsSame(int x, int y, vector<string> * vAssignPep)
{
	return vAssignPep[x].size() == vAssignPep[y].size();
}

bool CACInfer::Less(const CAssignedProtein & a, const CAssignedProtein & b)
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

bool CACInfer::PeptideLess(const pair<const CPeptideResult *, size_t> & a, const pair<const CPeptideResult *, size_t> & b)
{
	return a.first->m_lfEvalue < b.first->m_lfEvalue;
}

//void DDebug(string &ds)
//{
//	return;
//	ofstream os(".\\inferdebug.txt",ios::app);
//	os << ds << endl;
//	os.close();
//}
vector<CAssignedProtein> & CACInfer::Infer(vector<CSimpleMatchResult> & vResults, const vector<CSpectrum> & vSpectra,
		FILTER_CRITERIA_INFO & cs)
{
	size_t tTotalPep = 0;
	for(size_t i = 0;i < vResults.size();++i)
	{
		tTotalPep += vResults[i].m_vPeptideResults.size();
	}
	vector<string> vstrSrc;
	vstrSrc.reserve(tTotalPep);
	vector<size_t> vId1;
	vId1.reserve(tTotalPep);
	vector<size_t> vId2;
	vId2.reserve(tTotalPep);
	
	//fill in the vector
	for(size_t i = 0;i < vResults.size();++i)
	{
		for(size_t j = 0;j < vResults[i].m_vPeptideResults.size();++j)
		{
			vstrSrc.push_back(vResults[i].m_vPeptideResults[j].m_peptide.m_szSequence);
			vResults[i].m_vPeptideResults[j].m_vProteinID.clear();
			vResults[i].m_vPeptideResults[j].m_vProteinAC.clear();
			vId1.push_back(i);
			vId2.push_back(j);
		}
	}
	CTrace * pTrace = CTrace::GetInstance(*m_pCond);
	pTrace->Debug("inferring algorithm initializing..");
	AC::Init(m_pCond);
	pTrace->Debug("creating the tree...");
	AC::Create(vstrSrc);
	pTrace->Debug("creating the tree completed.");
	
	//READE index.AC FILE
	CDBConf db(m_pCond->m_strDBConfPath.c_str());
	string strDBPath = db.GetPath(m_pCond->m_vSelectedDBName[0].c_str(), m_pCond->m_SelectedEnzyme.m_strName.c_str());
	CProteinRandomReaderFactory factory;
	// add by emily : protein index
	CRandomReader *pIndexACLoader = factory.GetRandomACReader(Protein_Random_Index_Map);
	CRandomReader *pIndexSQLoader = factory.GetRandomSQReader(Protein_Random_Index_Map);

//	CRandomReader *pIndexACLoader = factory.GetRandomACReader(Protein_Random_Index_Disk);//czhou
//	CRandomReader *pIndexSQLoader = factory.GetRandomSQReader(Protein_Random_Index_Disk);
	try
	{
		// add by emily : protein index
		pIndexACLoader->Open(strDBPath, m_pCond->m_vSelectedDBName[0]);
		pIndexSQLoader->Open(strDBPath, m_pCond->m_vSelectedDBName[0]);
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CACferor", "Run", "in the function ProteinIndex::Open.");
		info.Append("strDBPath=" + strDBPath);
		info.Append("DBName=" + m_pCond->m_vSelectedDBName[0]);
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CACferor", "Run", "caught an unknown exception in the function ProteinIndex::Open.");
		info.Append("strDBPath=" + strDBPath);
		info.Append("DBName=" + m_pCond->m_vSelectedDBName[0]);
		throw runtime_error(info.Get().c_str());		
	}
	size_t tProNum = pIndexACLoader->GetProNum();
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
		
		
		set<size_t> setIDs;
		AC::Query(strPro, setIDs);

		for(set<size_t>::iterator it = setIDs.begin();
			it != setIDs.end();
			++it)
		{
			//vResults[vId1[*it]].m_vPeptideResults[vId2[*it]].m_vProteinID.push_back(i);
			vResults[vId1[*it]].m_vPeptideResults[vId2[*it]].m_vProteinID.push_back(i);
			vResults[vId1[*it]].m_vPeptideResults[vId2[*it]].m_vProteinAC.push_back(strAC);
			
		}
	}
	
//	string bs = "end";
//	DDebug(bs);
	
	// add by emily : protein index
	pIndexACLoader->Close();
	pIndexSQLoader->Close();
	delete pIndexACLoader;
	delete pIndexSQLoader;
	//the same as simpleinfer
	m_vFoundProteins.clear();
	
	double lfThreshold = Filter(vResults, cs);
	size_t tSpecOrder=0;

	map<size_t, int> mapProteinRecord;


	string * strProteinName = new string[tProNum];
	size_t * tProteinID = new size_t[tProNum];
	vector< pair<const CPeptideResult * ,size_t> > * vpContainPep = 
		new vector< pair<const CPeptideResult * ,size_t> >[tProNum];
	double * lfEvaluate = new double[tProNum];
	double * lfFPR = new double[tProNum];
	int * root = new int[tProNum];
	
	vector<string> * vAssignPep = new vector<string>[tProNum];
	fill(lfEvaluate, lfEvaluate + tProNum, 1);
	fill(lfFPR, lfFPR + tProNum, 1);
	fill(root, root + tProNum, -1);
	int nTotalPtn = 0;

	for(vector<CSimpleMatchResult>::const_iterator it_result = vResults.begin(); it_result != vResults.end(); ++it_result,++tSpecOrder)
	{
		if(it_result->m_vPeptideResults.empty())
			continue;

		for(vector<CPeptideResult>::const_iterator it_list =  it_result->m_vPeptideResults.begin();
			it_list != it_result->m_vPeptideResults.end()/* && !bOnce*/;
			++ it_list)
		{
			if(!it_list->m_bEV || (it_list->m_lfEvalue > lfThreshold))
				continue;

			if((int)it_list->m_peptide.m_tLength < cs.nPepLength || it_list->m_peptide.m_lfMass < cs.lfPepMass)
				continue;

			size_t sizePro = it_list->m_vProteinID.size(), i, j, k;
			vector<size_t> vFirstAppearProtein;

			vector<int> vnDifferOrd;
			int nR = -1;
			for(i = 0;i < sizePro;++i)
			{
				if(it_list->m_vProteinAC[i].substr(0, m_strFPSign.size()) == m_strFPSign)
					continue;
				int nOrd = mapProteinRecord[it_list->m_vProteinID[i]];
				if(0 == nOrd)
				{
					nOrd = mapProteinRecord[it_list->m_vProteinID[i]] = (++nTotalPtn);
					vFirstAppearProtein.push_back(it_list->m_vProteinID[i]);
					strProteinName[nTotalPtn] = it_list->m_vProteinAC[i];
					tProteinID[nTotalPtn] = it_list->m_vProteinID[i];
				}
				else
				{
					int nRoot = _GetRoot(root, nOrd);
					nR = nRoot;
					if(nRoot != nOrd)
					{
						if(find(it_list->m_vProteinID.begin(), 
							it_list->m_vProteinID.end(), 
							tProteinID[nRoot]) != it_list->m_vProteinID.end())
						{
							root[nOrd] = nRoot;
							nR = nRoot;
						}
						else
						{
							if(_IsSame(nOrd, nRoot, vAssignPep))
							{
								root[nOrd] = -1;
								root[nRoot] = nOrd;
								nR = nOrd;
							}
							else
							{
								//create a new root
								
								//root[nOrd] = -1;
								vnDifferOrd.push_back(nOrd);
								nR = nOrd;
							}
						}
					}
				}
				vAssignPep[nOrd].push_back(it_list->m_peptide.m_szSequence);
				vpContainPep[nOrd].push_back( make_pair(&(*it_list),tSpecOrder));
				lfEvaluate[nOrd] *= it_list->m_lfEvalue;
				lfFPR[nOrd] *= it_list->m_lfFPR;


			}//end for(i = 0;i < sizePro;++i)


			if(!vnDifferOrd.empty())
			{
				bool *visit = new bool[tProNum];
				memset(visit, 0, sizeof(bool) * tProNum);
				size_t sizeOrd = vnDifferOrd.size();
				for(i = 0;i < sizeOrd;++i)
				{
					if(visit[i])
						continue;
					visit[i] = 1;
					vector<int> vTempOrd;
					vTempOrd.push_back(vnDifferOrd[i]);
					size_t tSize = vAssignPep[vnDifferOrd[i]].size();
					int temproot = _GetRoot(root, vnDifferOrd[i]);
					int newroot = vnDifferOrd[i];
					for(j = i + 1;j < sizeOrd;++j)
					{
						if(!visit[j] && _GetRoot(root, vnDifferOrd[j]) == temproot)
						{
							visit[j] = 1;
							vTempOrd.push_back(vnDifferOrd[j]);
							if(vAssignPep[vnDifferOrd[j]].size() > tSize)
							{
								newroot = vnDifferOrd[j];
								tSize = vAssignPep[vnDifferOrd[j]].size();
							}
						}
					}
					for(j = 0;j < vTempOrd.size();++j)
						root[vTempOrd[j]] = newroot;
					for(j = 0;j < sizePro;++j)
					{
						int testOrd = mapProteinRecord[it_list->m_vProteinID[i]];
						if(find(vTempOrd.begin(), vTempOrd.end(), testOrd) == vTempOrd.end())
						{
							bool f = true;
							for(k = 0;k < vAssignPep[vnDifferOrd[i]].size();++k)
							{
								if(find(vAssignPep[testOrd].begin(), vAssignPep[testOrd].end(), vAssignPep[vnDifferOrd[i]][k]) == vAssignPep[testOrd].end())
								{
									f = false;
									break;
								}
							}
							if(f)
							{
								root[newroot] = testOrd;
								break;
							}
						}
					}
					root[newroot] = -1;
				}
				delete [] visit;
			}

			if(-1 == nR)
			{
				vector<size_t>::iterator it = vFirstAppearProtein.begin();
				if(it != vFirstAppearProtein.end())
				{
					nR = mapProteinRecord[*it];
					++it;
					while(it != vFirstAppearProtein.end())
					{
						root[mapProteinRecord[*it]] = nR;
						++it;
					}
				}
			}
			else
			{
				vector<size_t>::iterator it = vFirstAppearProtein.begin();
				while(it != vFirstAppearProtein.end())
				{
					root[mapProteinRecord[*it]] = nR;
					++it;
				}
			}
		}

	}//end for(vector<CMatchResult>::iterator it_result = vResults.begin(); it_result != vResults.end(); ++it_result,++tSpecOrder)

	//form the return value
	int i, j = 0;
	int * pNewOrd = new int [nTotalPtn + 1];
	memset(pNewOrd, 0, sizeof(int) * (1 + nTotalPtn));
	for(i = 1;i <= nTotalPtn;++i)
	{
		if(root[i] != -1)
			continue;
		CAssignedProtein APro;
		APro.EmptyContents();
		APro.m_vpContainPep = vpContainPep[i];
		APro.m_strAC = strProteinName[i];
		APro.m_nProID = tProteinID[i];
		APro.m_nPeptide = 0;
		vector<string>::iterator it_pep = vAssignPep[i].begin();
		for(;it_pep != vAssignPep[i].end(); ++it_pep)
		{
			//if cannot find it in previous seq, ++
			if(it_pep == find(vAssignPep[i].begin(), it_pep, *it_pep))
				++APro.m_nPeptide;
		}

		sort(APro.m_vpContainPep.begin(), APro.m_vpContainPep.end(), PeptideLess);
		m_vFoundProteins.push_back(APro);
		pNewOrd[i] = j++;
	}
	for(i = 1;i <= nTotalPtn;++i)
	{
		if(i == _GetRoot(root,i))
			continue;
		if(_IsSame(i, root[i], vAssignPep))
		{
			m_vFoundProteins[pNewOrd[root[i]]].m_vSameSet.push_back(strProteinName[i]);
			m_vFoundProteins[pNewOrd[root[i]]].m_vSameSetID.push_back(tProteinID[i]);
		}
		else
		{
			m_vFoundProteins[pNewOrd[root[i]]].m_vSubSet.push_back(strProteinName[i]);
			m_vFoundProteins[pNewOrd[root[i]]].m_vSubSetID.push_back(tProteinID[i]);
		}
	}
	delete[] pNewOrd;

	sort(m_vFoundProteins.begin(), m_vFoundProteins.end(), Less);
	// end form the return value

	delete [] strProteinName;
	delete [] tProteinID;
	delete [] vpContainPep;
	delete [] lfEvaluate;
	delete [] lfFPR;
	delete [] root;
	
	delete [] vAssignPep;
	return m_vFoundProteins;
	
}

double CACInfer::Filter(const vector<CSimpleMatchResult> & vResults, FILTER_CRITERIA_INFO & cs)
{
	vector<pair<double,bool> > EVHistogram;

	double lfFPRThreshold = cs.lfFPR;
	for(vector<CSimpleMatchResult>::const_iterator it_result = vResults.begin(); it_result != vResults.end(); ++it_result)
	{
		if(!it_result->m_vPeptideResults.empty())
		{
			vector<CPeptideResult>::const_iterator itCandidate = it_result->m_vPeptideResults.begin();
			if( (itCandidate->m_peptide.m_lfMass < cs.lfPepMass) || (itCandidate->m_peptide.m_tLength < cs.nPepLength))
				continue;
			double lfEV = itCandidate->m_lfEvalue;
			bool bReverse = false;
			for(size_t i=0 ; i < itCandidate->m_vProteinAC.size(); i++)
			{
				if(m_strFPSign == itCandidate->m_vProteinAC[i].substr(0, m_strFPSign.size()))
				{
					bReverse = true;
					break;
				}
			}
			EVHistogram.push_back(make_pair(lfEV, bReverse));
		}
	}

	double lfThreshold = 1;

	if(!EVHistogram.empty())
	{
		sort(EVHistogram.begin(), EVHistogram.end(), EVLesser);
	
		lfThreshold = EVHistogram[0].first;
		
		
		for(size_t i=0, nReverse= 0 ; i<EVHistogram.size(); i++)
		{
			if(EVHistogram[i].second)
				nReverse ++;

			double lfFPR =  2 * (double) nReverse / (double)i;
	
			if(lfFPRThreshold > lfFPR && !EVHistogram[i].second && lfThreshold < EVHistogram[i].first)
				lfThreshold = EVHistogram[i].first;
			
			//cout << i << ' ' << lfThreshold << ' ' << lfFPR << endl;
		}

		if(0 >= lfThreshold)
			lfThreshold = 1;
	}
	return lfThreshold;
}
