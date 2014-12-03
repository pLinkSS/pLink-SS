#include <string>
#include <vector>
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "../include/predefine.h"
#include "common.h"
#include "SimpleGenerator.h"

using namespace std;
using namespace proteomics_sdk;

void CSimpleGenerator::Init(const CAAConf & aaconf, bool bMono)
{
	int nDeny = 0x2904202;
	//srand(time(0));
	srand(AAC_TIME_SEED);
	CMapAAMass mapAAMass = aaconf.GetMapAAMass();
	for(size_t i = 0;i < 26;++i)
	{
		if(nDeny & (1 << i))
			m_lfMass[i] = 0.0;
		else
		if(bMono)
			m_lfMass[i] = mapAAMass.m_mapAAMass[i + 'A'].m_lfMonoMass;
		else
			m_lfMass[i] = mapAAMass.m_mapAAMass[i + 'A'].m_lfAvrgMass;
	}	
	/*
	//int nDeny = 0x2904202;
	int nDeny = 0x2904A02;
	//srand(time(0));
	srand(AAC_TIME_SEED);
	
	CMapAAMass mapAAMass = aaconf.GetMapAAMass();
	size_t nId = 0;
	for(size_t i = 0;i < 26;++i)
	{
		if(nDeny & (1 << i))
		{
			// do nothing
			
			// exclude B J L O U X Z
			// therefore AA_NUM = 19
		}
		else
		{
			if(bMono)
			{
				m_lfMass[nId] = mapAAMass.m_mapAAMass[i + 'A'].m_lfMonoMass;
				m_cAA[nId] = char(i + 'A');
				++nId;
			}
			else
			{
				m_lfMass[nId] = mapAAMass.m_mapAAMass[i + 'A'].m_lfAvrgMass;
				m_cAA[nId] = char(i + 'A');
				++nId;
			}
		}
	}	
	
	*/
};
//void GetPeptide(double lfMass, const string & strIn, string & strOut)
//{
//	strOut = "";
//	double mass = 0.0;
//	size_t tPos = 0;
//	while(mass < lfMass && tPos < strIn.length())
//	{
//		strOut += strIn[tPos];
//		mass += m_lfMass[strIn[tPos] - 'A'];
//		++tPos;
//	}
//}
void CSimpleGenerator::GetPeptide(double lfMass, string & strOut)
{
	strOut = "";
	double mass = 0.0;

	
	/*
	while(mass < lfMass)
	{
		int x = abs(rand()) % AA_NUM;
		strOut += m_cAA[x];
		mass += m_lfMass[x];
	}
	*/
	
	while(mass < lfMass)
	{
		int x = abs(rand()) % 26;
		while(m_lfMass[x] < 0.01)
		{
			x = abs(rand()) % 26;
		}
		strOut += (char)(x + 'A');
		mass += m_lfMass[x];
	}
	
	
};

void CSimpleGenerator::Batch(double lfMass, size_t tNum, vector<string> & vstrOut)
{
	srand(AAC_TIME_SEED);
	vstrOut.clear();
	vstrOut.reserve(tNum);
	string strTemp;
	for(size_t i = 0;i < tNum;++i)
	{
		GetPeptide(lfMass, strTemp);
		vstrOut.push_back(strTemp);
	}
};

void CSimpleGenerator::Close(void)
{
};
