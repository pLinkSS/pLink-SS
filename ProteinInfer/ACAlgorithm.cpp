#include<iostream>
#include<cstdio>
#include<vector>
#include<string>
#include<map>
#include<fstream>
#include<queue>
#include<set>
#include "../include/sdk.h"
#include "ACAlgorithm.h"
using namespace std;
using namespace AC;
using namespace proteomics_sdk;
#define MAX_STATE 5000000

Node g_states[MAX_STATE];
map<size_t, vector<size_t> > mapOutput;
vector<string> vTotalPep;
size_t tCurrentState;
size_t tTotQuery;
CCondition * pCond = NULL;
bool AC::Empty()
{
	return 1 == tCurrentState;
}
void AC::Init(CCondition * pC)
{
	memset(g_states, 0, sizeof(g_states));
	tCurrentState = 1;
	vTotalPep.clear();
	mapOutput.clear();
	pCond = pC;
}
size_t AC::Find(size_t tCurr, char cVal)
{
	if(0 == g_states[tCurr].m_tLeft)
	{
		return 0;
	}
	if(g_states[tCurr].m_tLeftEdge == cVal)
	{
		return g_states[tCurr].GetChild();
	}
	size_t tProc = g_states[tCurr].GetChild();
	while(1)
	{
		if(0 == g_states[tProc].GetBrother())
			break;
		if(g_states[tProc].m_tRightEdge == cVal)
			return g_states[tProc].GetBrother();
		tProc = g_states[tProc].GetBrother();
	}
	if(1 == tCurr)
		return 1;
	return 0;
}
void AC::MergeOutput(size_t tDstState, size_t tSrcState)
{
	if(mapOutput.find(tSrcState) == mapOutput.end())
		return;
	mapOutput[tDstState].insert(mapOutput[tDstState].begin(), mapOutput[tSrcState].begin(), mapOutput[tSrcState].end());
}
void AC::CreateFailFunction()
{
	g_states[1].SetF(1);
	queue<size_t> q;
	if(0 != g_states[1].GetChild())
	{
		q.push(g_states[1].GetChild());
		g_states[g_states[1].GetChild()].SetF(1);
		size_t tProc = g_states[1].GetChild();
		while(0 != g_states[tProc].GetBrother())
		{
			q.push(g_states[tProc].GetBrother());
			tProc = g_states[tProc].GetBrother();
			g_states[tProc].SetF(1);
		}
	}
	while(!q.empty())
	{
		size_t tCurr = q.front();
		q.pop();
		size_t tProc = g_states[tCurr].GetChild();
		if(tProc == 0)
		{
			continue;
		}
		q.push(tProc);

		size_t nF = g_states[tCurr].GetF();
		size_t tRet = 0;
		while(!(tRet = Find(nF, g_states[tCurr].m_tLeftEdge)))
		{
			nF = g_states[nF].GetF();
		}
		g_states[tProc].SetF(tRet);
		MergeOutput(tProc, tRet);
		while(1)
		{
			size_t tTemp = tProc;
			tProc = g_states[tProc].GetBrother();
			if(0 == tProc)
				break;
			q.push(tProc);
			size_t nF = g_states[tCurr].GetF();
			while(!(tRet = Find(nF, g_states[tTemp].m_tRightEdge)))
			{
				nF = g_states[nF].GetF();
			}
			g_states[tProc].SetF(tRet);
			MergeOutput(tProc, tRet);
		}
	}
}
void AC::Create(const vector<string> & vSrc)
{
	vTotalPep = vSrc;
	tCurrentState = 1;
	//for(size_t i = 0;i < vSrc[0].length();++i)
	//{
	//	g_states[tCurrentState].m_tLeft = tCurrentState + 1;
	//	g_states[tCurrentState].m_tLeftEdge = vSrc[0][i];
	//	++tCurrentState;
	//}
	for(size_t i = 0;i < vSrc.size();++i)
	{
		size_t tCursor = 1;
		size_t tLen = 0;
		for(tLen = 0;tLen < vSrc[i].length();++tLen)
		{
			//find whether vSrc[i][j] existed
			if(0 == g_states[tCursor].GetChild())//no child node
			{
				break;
			}
			else
			{
				if(g_states[tCursor].m_tLeftEdge == vSrc[i][tLen])
				{
					tCursor = g_states[tCursor].GetChild();
					continue;
				}
				else
				{
					size_t tBrother = g_states[tCursor].GetChild();
					while(1)
					{
						if(tBrother == 0)
							break;
						if(g_states[tBrother].m_tRightEdge == vSrc[i][tLen])
						{
							break;
						}
						tBrother = g_states[tBrother].GetBrother();
					}
					if(tBrother == 0)
					{
						break;
					}
					else
					{
						tCursor = g_states[tBrother].GetBrother();
						continue;
					}
				}
			}

		}
		if(tLen < vSrc[i].length())
		{
			//create new states from the tCursor node
			if(0 == g_states[tCursor].GetChild())
			{
				g_states[tCursor].SetChild(tCurrentState + 1, vSrc[i][tLen]);
				++tCurrentState;
				for(size_t j = tLen + 1;j < vSrc[i].length();++j)
				{
					g_states[tCurrentState].SetChild(tCurrentState + 1, vSrc[i][j]);
					++tCurrentState;
					if(tCurrentState == MAX_STATE)
					{
						CErrInfo info("SearchEngine", "ACAlgorithm", "in the function Create");
						info.Append("too many states");
						throw invalid_argument(info.Get().c_str());
					}
				}
			}
			else
			{
				size_t tBrother = g_states[tCursor].GetChild();
				while(g_states[tBrother].GetBrother() != 0)
				{
					tBrother = g_states[tBrother].GetBrother();
				}
				g_states[tBrother].SetBrother(tCurrentState + 1, vSrc[i][tLen]);
				++tCurrentState;
				for(size_t j = tLen + 1;j < vSrc[i].length();++j)
				{
					g_states[tCurrentState].SetChild(tCurrentState + 1, vSrc[i][j]);
					++tCurrentState;
					if(tCurrentState == MAX_STATE)
					{
						CErrInfo info("SearchEngine", "ACAlgorithm", "in the function Create");
						info.Append("too many states");
						throw invalid_argument(info.Get().c_str());
					}
				}


			}
			mapOutput[tCurrentState].push_back(i);
		}
		else
		{
			mapOutput[tCursor].push_back(i);
		}
	}
	CTrace * pTrace = CTrace::GetInstance();
	pTrace->Debug("creating fail function..");
	CreateFailFunction();
	pTrace->Debug("creating fail function completed.");

}

bool AC::Transfer(size_t & tCurr, char cVal)
{
	if(0 == g_states[tCurr].m_tLeft)
	{
		tCurr = g_states[tCurr].GetF();
		return false;
	}
	if(g_states[tCurr].m_tLeftEdge == cVal)
	{
		tCurr = g_states[tCurr].GetChild();
		return true;
	}
	size_t tProc = g_states[tCurr].GetChild();
	while(1)
	{
		if(0 == g_states[tProc].GetBrother())
		{
			size_t tTemp = tCurr;
			tCurr = g_states[tCurr].GetF();
			return 1 == tTemp ? true : false;
		}
		if(g_states[tProc].m_tRightEdge == cVal)
		{
			tCurr = g_states[tProc].m_tRight;
			return true;
		}
		tProc = g_states[tProc].GetBrother();
	}
	return false;
}
/*
bool AC::Transfer(size_t & tCurr, char cVal)
{
	 if(0 == g_states[tCurr].m_tLeft)
	 {
		  tCurr = g_states[tCurr].GetF();
		  return (tCurr == 0 ? true : false);
	 }
	 //if(EqualChar(g_states[tCurr].m_tLeftEdge, cVal))
	 if(g_states[tCurr].m_tLeftEdge == cVal)
	 {
		  tCurr = g_states[tCurr].GetChild();
		  return true;
	 }
	 size_t tProc = g_states[tCurr].GetChild();
	 while(1)
	 {
		  if(0 == g_states[tProc].GetBrother())
		  {
			   size_t tTemp = tCurr;
			   tCurr = g_states[tCurr].GetF();
			   return 0 == tCurr ? true : false;
		  }
		  //if(EqualChar(g_states[tCurr].m_tRightEdge, cVal))
		  if(g_states[tProc].m_tRightEdge == cVal)
		  {
			   tCurr = g_states[tProc].m_tRight;
			   return true;
		  }
		  tProc = g_states[tProc].GetBrother();
	 }
	 return false;
}
*/
/*

bool AC::Transfer(size_t & tCurr, char cVal)
{
	if(0 == g_states[tCurr].m_tLeft)
	{
		tCurr = g_states[tCurr].GetF();
		if(tCurr == 0)
		{
			tCurr = 1;
			return true;
		}
		else
		{
			return false;
		}
	}
//	if((g_states[tCurr].m_tLeftEdge == cVal) != EqualChar(g_states[tCurr].m_tLeftEdge, cVal))
//	{
//		cout << "diff: " << cVal << ' ' << g_states[tCurr].m_tLeftEdge << endl;
//	}
	//if(EqualChar(g_states[tCurr].m_tLeftEdge, cVal))
	if(g_states[tCurr].m_tLeftEdge == cVal)
	{
		tCurr = g_states[tCurr].GetChild();
		return true;
	}
	size_t tProc = g_states[tCurr].GetChild();
	while(1)
	{
		if(0 == g_states[tProc].GetBrother())
		{
			size_t tTemp = tCurr;
			tCurr = g_states[tCurr].GetF();
			if(tCurr == 0)
			{
				tCurr = 1;
				return true;
			}
			else
			{
				return false;
			}
		}
		//if(EqualChar(g_states[tCurr].m_tRightEdge, cVal))
		if(g_states[tProc].m_tRightEdge == cVal)
		{
			tCurr = g_states[tProc].m_tRight;
			return true;
		}
		tProc = g_states[tProc].GetBrother();
	}
	return false;
}
*/
int AC::Output(const string &strSrc, size_t tPepEnd, size_t tCurr, set<size_t> & setIDs)
{
	if(mapOutput.find(tCurr) == mapOutput.end())
		return 0;
	for(size_t i = 0;i < mapOutput[tCurr].size();++i)
	{
		//是否是特异性酶切
		string & str = vTotalPep[mapOutput[tCurr][i]];
//		if(str == "GFHIHEFGDATNGCVSAGPHFNPFK")
//		{
//			cout << tCurr<< endl << strSrc << endl << tPepEnd << endl << strSrc[tPepEnd - str.length()] << endl;
//			getchar();
//		}
		const CEnzyme & enzyme = pCond->m_SelectedEnzyme;
		const string & strSites = enzyme.m_strCleave;
		const string & strNotCleave = enzyme.m_strNotCleave;
		switch(pCond->m_nCleaveWay)
		{
			case 0:
				if(enzyme.GetIsNTerm())
				{
					if((tPepEnd == str.length() - 1
							|| (str.length() == tPepEnd && strSrc[0] == 'M') 
							|| (strSites.find(str[0]) != string::npos && strNotCleave.find(strSrc[tPepEnd - str.length()]) == string::npos))
						&&(tPepEnd == strSrc.length() - 1 
							|| (strSites.find(strSrc[tPepEnd + 1]) != string::npos && strNotCleave.find(str[str.length() - 1]) == string::npos)))
					{
						setIDs.insert(mapOutput[tCurr][i]);
					}
				}
				else
				{
					if((tPepEnd == str.length() - 1 
							|| (str.length() == tPepEnd && strSrc[0] == 'M') 
							|| (strSites.find(strSrc[tPepEnd - str.length()]) != string::npos && strNotCleave.find(str[0]) == string::npos))
						&&(tPepEnd == strSrc.length() - 1 
							|| (strSites.find(str[str.length() - 1]) != string::npos) && strNotCleave.find(strSrc[tPepEnd + 1]) == string::npos))
					{
						setIDs.insert(mapOutput[tCurr][i]);
					}
				}
				break;
			case 1:
				if(enzyme.GetIsNTerm())
				{
					if((tPepEnd == str.length() - 1
							|| (str.length() == tPepEnd && strSrc[0] == 'M') 
							|| (strSites.find(str[0]) != string::npos && strNotCleave.find(strSrc[tPepEnd - str.length()]) == string::npos))
						||(tPepEnd == strSrc.length() - 1 
							|| (strSites.find(strSrc[tPepEnd + 1]) != string::npos && strNotCleave.find(str[str.length() - 1]) == string::npos)))
					{
						setIDs.insert(mapOutput[tCurr][i]);
					}
				}
				else
				{
					if((tPepEnd == str.length() - 1 
							|| (str.length() == tPepEnd && strSrc[0] == 'M') 
							|| (strSites.find(strSrc[tPepEnd - str.length()]) != string::npos && strNotCleave.find(str[0]) == string::npos))
						||(tPepEnd == strSrc.length() - 1 
							|| (strSites.find(str[str.length() - 1]) != string::npos) && strNotCleave.find(strSrc[tPepEnd + 1]) == string::npos))
					{
						setIDs.insert(mapOutput[tCurr][i]);
					}
				}
				break;
			case 2:
				setIDs.insert(mapOutput[tCurr][i]);
				break;
			default:
				setIDs.insert(mapOutput[tCurr][i]);
		}
//		switch(tCleaveWay)
//		{
//			case 0:		
//				if((tPepEnd == str.length() - 1 || (str.length() == tPepEnd && strSrc[0] == 'M') || strSrc[tPepEnd - str.length()] == 'K' || strSrc[tPepEnd - str.length()] == 'R')
//					&&(tPepEnd == strSrc.length() - 1 || str[str.length() - 1] == 'K' || str[str.length() - 1] == 'R'))
//				{
//					setIDs.insert(mapOutput[tCurr][i]);
//				}
//				break;
//			case 1:
//				if((tPepEnd == str.length() - 1 || (str.length() == tPepEnd && strSrc[0] == 'M') || strSrc[tPepEnd - str.length()] == 'K' || strSrc[tPepEnd - str.length()] == 'R')
//					||(tPepEnd == strSrc.length() - 1 || str[str.length() - 1] == 'K' || str[str.length() - 1] == 'R'))
//				{
//					setIDs.insert(mapOutput[tCurr][i]);
//				}
//				break;
//			case 2:
//				setIDs.insert(mapOutput[tCurr][i]);
//				break;
//			default:
//				setIDs.insert(mapOutput[tCurr][i]);
//		}

		
	}
	return mapOutput[tCurr].size();
}

void AC::Query(const string & strSrc, set<size_t> & setIDs)
{
	if(!Empty())
	{
		size_t tCurr = 1;
		for(size_t i = 0;i < strSrc.size();++i)
		{
			if(Transfer(tCurr, strSrc[i]))
			{
				Output(strSrc, i, tCurr, setIDs);
			}
			else
				--i;
		}
	}
}

void AC::Query(const string & strSrc, set< pair<size_t, size_t> > & setIDs)
{
	if(!Empty())
	{
		size_t tCurr = 1;
		for(size_t i = 0;i < strSrc.size();++i)
		{
			if(Transfer(tCurr, strSrc[i]))
			{
				Output(strSrc, i, tCurr, setIDs);
			}
			else
				--i;
		}
	}
}

int AC::Output(const string &strSrc, size_t tPepEnd, size_t tCurr, set< pair<size_t, size_t> > & setIDs)
{
	if(mapOutput.find(tCurr) == mapOutput.end())
		return 0;
//	cout<<"-------------Pro: "<<strSrc<<endl;
	for(size_t i = 0;i < mapOutput[tCurr].size();++i)
	{
		//是否是特异性酶切
		string & str = vTotalPep[mapOutput[tCurr][i]];
//		cout<<str<<" ";
		const CEnzyme & enzyme = pCond->m_SelectedEnzyme;
		const string & strSites = enzyme.m_strCleave;
		const string & strNotCleave = enzyme.m_strNotCleave;
		switch(pCond->m_nCleaveWay)
		{
			case 0:
				if(enzyme.GetIsNTerm())
				{
					if((tPepEnd == str.length() - 1
							|| (str.length() == tPepEnd && strSrc[0] == 'M')
							|| (strSites.find(str[0]) != string::npos && strNotCleave.find(strSrc[tPepEnd - str.length()]) == string::npos))
						&&(tPepEnd == strSrc.length() - 1
							|| (strSites.find(strSrc[tPepEnd + 1]) != string::npos && strNotCleave.find(str[str.length() - 1]) == string::npos)))
					{
						setIDs.insert( make_pair(mapOutput[tCurr][i], tPepEnd) );
					}
				}
				else
				{
					if((tPepEnd == str.length() - 1
							|| (str.length() == tPepEnd && strSrc[0] == 'M')
							|| (strSites.find(strSrc[tPepEnd - str.length()]) != string::npos && strNotCleave.find(str[0]) == string::npos))
						&&(tPepEnd == strSrc.length() - 1
							|| (strSites.find(str[str.length() - 1]) != string::npos) && strNotCleave.find(strSrc[tPepEnd + 1]) == string::npos))
					{
						setIDs.insert( make_pair(mapOutput[tCurr][i], tPepEnd) );
					}
				}
				break;
			case 1:
				if(enzyme.GetIsNTerm())
				{
					if((tPepEnd == str.length() - 1
							|| (str.length() == tPepEnd && strSrc[0] == 'M')
							|| (strSites.find(str[0]) != string::npos && strNotCleave.find(strSrc[tPepEnd - str.length()]) == string::npos))
						||(tPepEnd == strSrc.length() - 1
							|| (strSites.find(strSrc[tPepEnd + 1]) != string::npos && strNotCleave.find(str[str.length() - 1]) == string::npos)))
					{
						setIDs.insert( make_pair(mapOutput[tCurr][i], tPepEnd) );
					}
				}
				else
				{
					if((tPepEnd == str.length() - 1
							|| (str.length() == tPepEnd && strSrc[0] == 'M')
							|| (strSites.find(strSrc[tPepEnd - str.length()]) != string::npos && strNotCleave.find(str[0]) == string::npos))
						||(tPepEnd == strSrc.length() - 1
							|| ((strSites.find(str[str.length() - 1]) != string::npos) && strNotCleave.find(strSrc[tPepEnd + 1]) == string::npos)))
					{
						setIDs.insert( make_pair(mapOutput[tCurr][i], tPepEnd) );
					}
				}
				break;
			case 2:
				setIDs.insert( make_pair(mapOutput[tCurr][i], tPepEnd) );
				break;
			default:
				setIDs.insert( make_pair(mapOutput[tCurr][i], tPepEnd) );
		}
	}
//	cout<<endl;
	return mapOutput[tCurr].size();
}


