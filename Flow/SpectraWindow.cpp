#include<vector>
#include "../include/sdk.h"
#include "SpectraWindow.h"
using namespace std;
namespace proteomics_search
{

CSpectraWindow::CSpectraWindow()
{
	vSpectraWindows.clear();
	m_tCurrentWindow = 0;
}

CSpectraWindow::~CSpectraWindow()
{
	vSpectraWindows.clear();
	m_tCurrentWindow = 0;
}

bool CSpectraWindow::Construct(const vector<SPEC_SORTER_INDEX_INFO> & vInfo)
{
	vSpectraWindows.clear();
	// suppose vInfo has been sorted 
	size_t tWindowStart = 0,tWindowEnd = 0;
	struct SPECTRA_WINDOW_INFO stWindowInfo;
	
	while(tWindowStart < vInfo.size() || tWindowEnd < vInfo.size())
	{
		if(tWindowEnd >= vInfo.size())
		{
			return false;
		}

		if(tWindowStart < vInfo.size() && vInfo[tWindowStart].lfMin < vInfo[tWindowEnd].lfMax)
		{
			if(vSpectraWindows.size())
			{
				vSpectraWindows[vSpectraWindows.size()-1].lfMax = vInfo[tWindowStart].lfMin; 
			}
			stWindowInfo.lfMin = vInfo[tWindowStart].lfMin;
			stWindowInfo.lfMax = vInfo[tWindowStart].lfMin;
			stWindowInfo.nStartId = -1;
			stWindowInfo.tSpectraNum = 0;
			vSpectraWindows.push_back(stWindowInfo);
			tWindowStart++;
		}
		else
		{
			if(vSpectraWindows.size())
			{
				vSpectraWindows[vSpectraWindows.size()-1].lfMax = vInfo[tWindowEnd].lfMax; 
			}
			
			double lfLB = vInfo[tWindowEnd].lfMin;
			for(size_t i=vSpectraWindows.size()-1;i>=0;--i)
			{
				if(vSpectraWindows[i].lfMin >= lfLB)
				{
					vSpectraWindows[i].tSpectraNum ++;
					if(vSpectraWindows[i].nStartId == -1)
					{
						vSpectraWindows[i].nStartId = tWindowEnd;
					}
				}
				else
					break;
			}
			
			stWindowInfo.lfMin = vInfo[tWindowEnd].lfMax;
			stWindowInfo.lfMax = vInfo[tWindowEnd].lfMax;
			stWindowInfo.nStartId = -1;
			stWindowInfo.tSpectraNum = 0;
			vSpectraWindows.push_back(stWindowInfo);
			
			tWindowEnd++;
			if(vInfo[tWindowStart].lfMin == vInfo[tWindowEnd].lfMax)
			{
				tWindowStart++;
			}
		}
	}
	
	/*
	for(int i=0;i<vSpectraWindows.size();++i)
	{
		cout << i << "	:" << vSpectraWindows[i].tSpectraNum << "	" 
			<< vSpectraWindows[i].nStartId << "	"
			<< vSpectraWindows[i].lfMin << "	"
			<< vSpectraWindows[i].lfMax << endl << endl;
	}
	
	system("pause");
	*/
	
	return true;
}

void CSpectraWindow::Begin()
{
	m_tCurrentWindow = 0;
}

bool CSpectraWindow::GetCurrentWindowInfo(struct SPECTRA_WINDOW_INFO & stInfo)
{
	if(m_tCurrentWindow == vSpectraWindows.size())
		return false;
	stInfo = vSpectraWindows[m_tCurrentWindow];
	m_tCurrentWindow++;
	return true;
}

}


