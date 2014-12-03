#include<vector>
#include "../include/sdk.h"
#include "SpectraSearcher.h"
using namespace std;
namespace proteomics_search
{
void CSpectraSearcher::Construct(const vector<SPEC_SORTER_INDEX_INFO> & vInfo)
{
	if(vInfo.empty())
		return;
	
	Close();
	m_nStart = (int)vInfo[0].lfMax;
	m_pData = new int [(int)vInfo[vInfo.size() - 1].lfMax - m_nStart + 1];
	memset(m_pData, 0, sizeof(int) * ((int)vInfo[vInfo.size() - 1].lfMax - m_nStart + 1));
	int nTemp = m_nStart - 1;
	for(int i = 0;i < (int)vInfo.size();++i)
	{
		if((int)vInfo[i].lfMax > nTemp)
		{
			for(int j = (int)vInfo[i].lfMax;j > nTemp;--j)
			{
				m_pData[j - m_nStart] = i;
			}
			nTemp = (int)vInfo[i].lfMax;
		}
	}
}
void CSpectraSearcher::Close()
{
	if(m_pData != NULL)
	{
		delete [] m_pData;
		m_pData = NULL;
	}
	m_nStart = 0;
}
int CSpectraSearcher::FindLowerBound(double & lfMass)
{
	return m_pData[(int)lfMass];
}
CSpectraSearcher::CSpectraSearcher(): m_pData(NULL), m_nStart(0){}
CSpectraSearcher::~CSpectraSearcher()
{
	if(m_pData != NULL)
	{
		delete [] m_pData;
		m_pData = NULL;
	}
}
}
