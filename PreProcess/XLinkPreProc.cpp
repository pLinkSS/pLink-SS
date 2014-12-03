#include <algorithm>
#include <iostream>
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "XLinkPreProc.h"

using namespace std;
using namespace proteomics_sdk;



CXLinkPreProc::CXLinkPreProc()
{
}

CXLinkPreProc::~CXLinkPreProc()
{
}

void CXLinkPreProc::Init(CCondition & condition)
{
	CTrace * pTrace = CTrace::GetInstance();
	pTrace->Debug("In CXLinkPreProc::Init(CCondition & condition)");
	
	if(!condition.ValidateAll())
	{
		CErrInfo info("CXLinkPreProc", "Init", "The validation of inputing condition is failed.");
		throw invalid_argument(info.Get().c_str());
	}

	m_Condition = condition;
	
	pTrace->Debug("Out CXLinkPreProc::Init(CCondition & condition)");
}

bool CXLinkPreProc::Instrument(InstrumentType eType) const
{
	return true;
}

bool CXLinkPreProc::_IsMatched(double lfCalMz,double lfExpMz)
{
	double lfTol = m_Condition.m_lfFragmentTol;
	if(m_Condition.m_strFragmentTolType == "ppm")
	{
		lfTol *= 0.000001*lfExpMz;
	}


	if(lfCalMz <= lfExpMz + lfTol && lfCalMz >= lfExpMz - lfTol)
		return true;
	else
		return false;
}

void CXLinkPreProc::Run(CSpectrum Srp, CSpectrum &Scp) 
{	
	try {
		if(!m_Condition.ValidateAll())
		{
			CErrInfo info("CXLinkPreProc", "Run", "The validation of inputing condition is failed.");
			throw invalid_argument(info.Get().c_str());
		}
		if(0 == Srp.m_tPeaksNum)
		{
			CErrInfo info("CXLinkPreProc", "Run", "The spectrum is empty.");
			throw invalid_argument(info.Get().c_str());
		}

		vector<bool > Yxf;
		vector<struct PeakInfo> Fx;
		for(size_t i = 0 ;i < Srp.m_tPeaksNum ; ++i)
		{
			Yxf.push_back(true);
			Fx.push_back(PeakInfo());
		}

		sort(Srp.m_pPeaks, Srp.m_pPeaks + Srp.m_tPeaksNum, proteomics_sdk::CSpectrum::mz_lesser);

		vector<double > Zzh;
		Zzh.push_back(108);
		Zzh.push_back(153);
		Zzh.push_back(200);
		int Zh = 0;
		for(size_t i = 0 ;i < Srp.m_tPeaksNum && Zh < Zzh.size(); ++ i)
		{
			if(Srp.m_pPeaks[i].lfMz < Zzh[Zh])
				continue;
			else if(Srp.m_pPeaks[i].lfMz > Zzh[Zh]+1)
			{
				while(Srp.m_pPeaks[i].lfMz > Zzh[Zh]+1)
				{
					Zh++;
					if(Zh >= Zzh.size())
						break;
				}
			}
			else
			{
				if(int(Srp.m_pPeaks[i].lfMz) == Zzh[Zh])
				{
					Fx[i].nPeakType = 5;
					Fx[i].nCharge = 0;
				}
			}
		}
	
		vector<double > Zl;
		for(int i = Srp.m_nCharge ;i >= 1 ; -- i)
		{
			Zl.push_back((Srp.m_lfMH + (i-1)*IonMass_Proton)/i);
		}

		int Dh = 0;
		for(size_t i = 0 ;i < Srp.m_tPeaksNum ; ++ i)
		{
			if(Srp.m_pPeaks[i].lfMz < Zl[Dh]-1)
				continue;
			else if(Srp.m_pPeaks[i].lfMz > Zl[Dh]+1)
			{
				while(Srp.m_pPeaks[i].lfMz > Zl[Dh]+1 && Dh < Srp.m_nCharge-1)
				{
					Dh++;
				}
			}
			else
			{
				if(_IsMatched(Zl[Dh],Srp.m_pPeaks[i].lfMz))
				{
					Fx[i].nPeakType = 4;
					Fx[i].nCharge = Srp.m_nCharge - Dh;
				}
			}
		}

		for(size_t i = 0 ;i < Srp.m_tPeaksNum ; ++ i)
		{
			for(size_t j = i + 1 ; j < Srp.m_tPeaksNum ; ++ j)
			{
				double Pc = Srp.m_pPeaks[j].lfMz - Srp.m_pPeaks[i].lfMz;
				if(Pc > 1.5)
					break;
				for(int k = Srp.m_nCharge ; k >= 1 ; --k)
				{
					if(_IsMatched(Srp.m_pPeaks[i].lfMz+IonMass_Proton/k,Srp.m_pPeaks[j].lfMz))
					{
						double Bl = Srp.m_pPeaks[j].lfIntensity / Srp.m_pPeaks[i].lfIntensity;
						if(Bl <= 2)
						{
							if(Fx[j].nPeakType != 0)
								continue;

							if(Fx[i].nPeakType == 2 || Fx[i].nPeakType == 3 || Fx[i].nPeakType == 5)
								continue;

							if(Fx[j].nCharge != 0 && Fx[j].nCharge != k)
								continue;

							if(Fx[i].nCharge != 0 && Fx[i].nCharge != k)
								continue;

							Fx[j].nPeakType = 1;
							Fx[j].nMainPeakId = i;
							Fx[j].nCharge = k;
							Fx[i].nCharge = k;
							break;
						}
					}
				}
			}
		}
		for(int i = Srp.m_tPeaksNum - 1 ;i >= 0 ; --i)
		{
			for(int j = i - 1 ; j >= 0 ; -- j)
			{
				double Pc = Srp.m_pPeaks[i].lfMz - Srp.m_pPeaks[j].lfMz;

				if(Pc > H2OMass_Aver + 1.5)
					break;

				for(int k = Srp.m_nCharge ; k >= 1 ; --k)
				{
					if(_IsMatched(Srp.m_pPeaks[j].lfMz+NH3Mass_Mono/k,Srp.m_pPeaks[i].lfMz))
					{
						double Bl = Srp.m_pPeaks[j].lfIntensity / Srp.m_pPeaks[i].lfIntensity;
						if(Bl < 1.5)
						{
							if(Fx[j].nPeakType != 0)
								continue;

							if(Fx[i].nPeakType != 0 && Fx[i].nPeakType != 4)
								continue;

							if(Fx[j].nCharge != 0 && Fx[j].nCharge != k)
								continue;

							if(Fx[i].nCharge != 0 && Fx[i].nCharge != k)
								continue;

							Fx[j].nPeakType = 2;
							Fx[j].nMainPeakId = i;
							Fx[j].nCharge = k;
							Fx[i].nCharge = k;
							break;
						}
					}
					else if(_IsMatched(Srp.m_pPeaks[j].lfMz+H2OMass_Mono/k,Srp.m_pPeaks[i].lfMz))
					{
						double Bl = Srp.m_pPeaks[j].lfIntensity / Srp.m_pPeaks[i].lfIntensity;
						if(Bl < 1.5)
						{
							if(Fx[j].nPeakType != 0)
								continue;

							if(Fx[i].nPeakType != 0 && Fx[i].nPeakType != 4)
								continue;

							if(Fx[j].nCharge != 0 && Fx[j].nCharge != k)
								continue;

							if(Fx[i].nCharge != 0 && Fx[i].nCharge != k)
								continue;

							Fx[j].nPeakType = 3;
							Fx[j].nMainPeakId = i;
							Fx[j].nCharge = k;
							Fx[i].nCharge = k;
							break;
						}
					}
				}
			}
		}

		int Fs = Srp.m_tPeaksNum;

		double Zq = 0.0;
		for(size_t i = 0 ;i < Srp.m_tPeaksNum ; ++ i)
		{
			if(Fx[i].nPeakType == 0 )
			{
				if(Srp.m_pPeaks[i].lfIntensity > Zq)
					Zq = Srp.m_pPeaks[i].lfIntensity;
				continue;
			}

			Yxf[i] = false;
			Fs -- ;
		}

	#ifdef FOR_PSM
		for(size_t i = 0;i < Srp.m_tPeaksNum;++i)
		{
			if(Yxf[i])
			{
				if(Srp.m_pPeaks[i].lfIntensity/Zq < 0.01)
				{
					Yxf[i] = false;
					Fs --;
				}
			}
		}
	#endif

		Scp.clear();
		Scp.CopyBasicItems(Srp);
		Scp.m_tPeaksNum = Fs;
		Scp.m_pPeaks = new CPeak[Scp.m_tPeaksNum];
		Scp.m_lfSqtMaxInten = 0.0;

		size_t i = 0 ,  j = 0;
		for(i = 0 , j = 0 ;i < Srp.m_tPeaksNum;++i)
		{
			if(Yxf[i])
			{
				Scp.m_pPeaks[j].lfMz = Srp.m_pPeaks[i].lfMz;
				Scp.m_pPeaks[j].nMz = (int)(0.5 + Srp.m_pPeaks[i].lfMz * MZMULTIPLIER);
				Scp.m_pPeaks[j].nCharge = Fx[i].nCharge;

				if(Srp.m_pPeaks[i].lfIntensity < 0.0)
				{
					CErrInfo info("CXLinkPreProc", "Run", "intensity<0");
					info.Append("file_path=" + Scp.m_strFilePath);
					throw runtime_error(info.Get().c_str());
				}
	#ifdef FOR_PSM
				Scp.m_pPeaks[j].lfIntensity = Srp.m_pPeaks[i].lfIntensity/Zq;
	#else
				Scp.m_pPeaks[j].lfIntensity = sqrt(Srp.m_pPeaks[i].lfIntensity);

				if(Scp.m_pPeaks[j].lfIntensity > Scp.m_lfSqtMaxInten)
					Scp.m_lfSqtMaxInten = Scp.m_pPeaks[j].lfIntensity;
	#endif
				j++;
			}
		}

		Scp.Create_Hash_Table();
	} catch(exception & e) {
		CErrInfo info("XLinkPreProc", "Run",
							"in the function CPreProcess::Run.", e);
		throw runtime_error(info.Get(e).c_str());
	} catch(...) {
		CErrInfo info("XLinkPreProc", "Run",
							"caught an unknown exception in the function CPreProcess::Run.");
		throw runtime_error(info.Get().c_str());
	}
}

void CXLinkPreProc::Close(void)
{
	m_Condition.clear();
}
