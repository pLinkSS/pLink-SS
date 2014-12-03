#include "../include/sdk.h"
#include "../include/predefine.h"
#include <fstream>
#include <iostream>
using namespace std;
using namespace proteomics_sdk;


#include "XLinkPFDFileIO.h"
#ifdef WIN32
const char cSlash = '\\';
#else
const char cSlash = '/';
#endif
CXLinkPFDFileIO::CXLinkPFDFileIO()
{
	m_pTrace = CTrace::GetInstance();
}

CXLinkPFDFileIO::~CXLinkPFDFileIO()
{
}

bool CXLinkPFDFileIO::WriteFile(string strOutputFile, const vector<CXLinkMatchResult> & vResults,
			const vector<CSpectrum> & vSpectra)
{
	const int MAX_BUFFER_SIZE = 81920;
	char szBuf[MAX_BUFFER_SIZE] = {0};
	string strTXT;
	for(size_t i=0; i<vResults.size(); ++i)
	{
		sprintf(szBuf, "%s\n", vSpectra[i].m_strFilePath.c_str());
		strTXT += szBuf;
		
		// modify by emily
		// save ev-coef too for the convenience of refined - search
		sprintf(szBuf, "%d %.5lf %.5lf %.5lf %.5lf %.5lf %d %u %d %u\n", vSpectra[i].m_nCharge,
				vSpectra[i].m_lfIntensity, vSpectra[i].m_lfMH, vSpectra[i].m_lfMZ,
				vSpectra[i].m_stEVCoef.lfCoef0,vSpectra[i].m_stEVCoef.lfCoef1,vSpectra[i].m_nEvalueNo,
				vResults[i].m_tScore, vResults[i].m_vPeptideResults.size(),vResults[i].m_tRealCandidate);
		strTXT += szBuf;
		
		if(vResults[i].m_vPeptideResults.empty())
			continue;

		for(size_t j = 0;j < vResults[i].m_vPeptideResults.size();++j)
		{
			const CXLinkPepResult & pep_res = vResults[i].m_vPeptideResults[j];

			// score and e-value
			sprintf(szBuf, "%.6lf %E ", pep_res.m_lfScore, pep_res.m_lfEvalue);
			strTXT += szBuf;
			// alpha E-value and beta E-value added at 2014.9.10
			sprintf(szBuf, "%E %E ", pep_res.m_lfAlphaEvalue, pep_res.m_lfBetaEvalue);
			strTXT += szBuf;
			// linker // add link site 1. 2014.7.16.
			sprintf(szBuf,"%d %d %d %d ",pep_res.m_XLink.m_eXLinkType,pep_res.m_XLink.m_tAlphaSite + 1,pep_res.m_XLink.m_tBetaSite + 1,pep_res.m_XLink.m_nLinkerId);
			strTXT += szBuf;
			// pep1
			sprintf(szBuf,"%c %s %c %.6lf %d %d ",
					pep_res.m_AlphaPeptide.m_cPrev, pep_res.m_AlphaPeptide.m_szSequence, pep_res.m_AlphaPeptide.m_cNext, pep_res.m_AlphaPeptide.m_lfMass, 
					pep_res.m_AlphaPeptide.m_tLength, pep_res.m_vAlphaProteinID.size());
			strTXT += szBuf;
			
			// proteinId for pep1
			for(size_t k = 0; k < pep_res.m_vAlphaProteinID.size(); k++)
			{
				sprintf(szBuf, "%d ",pep_res.m_vAlphaProteinID[k]);
				strTXT += szBuf;		
			}
			// mod for pep1
			sprintf(szBuf, "%d ", pep_res.m_AlphaPeptide.m_tModCnt);
			strTXT += szBuf;
			for(size_t k = 0; k < pep_res.m_AlphaPeptide.m_tModCnt; k++)
			{
				sprintf(szBuf, "%d %d ",pep_res.m_AlphaPeptide.m_tModSites[k][0] + 1,  // add mod site 1. 2014.1.14.
						pep_res.m_AlphaPeptide.m_tModSites[k][1]);
				strTXT += szBuf;			
			}
			
			if(pep_res.m_bPair)
			{
				// pep2
				sprintf(szBuf,"%c %s %c %.6lf %d %d ",
						pep_res.m_BetaPeptide.m_cPrev, pep_res.m_BetaPeptide.m_szSequence, pep_res.m_BetaPeptide.m_cNext, pep_res.m_BetaPeptide.m_lfMass, 
						pep_res.m_BetaPeptide.m_tLength, pep_res.m_vBetaProteinID.size());
				strTXT += szBuf;
				
				// proteinId for pep2
				for(size_t k = 0; k < pep_res.m_vBetaProteinID.size(); k++)
				{
					sprintf(szBuf, "%d ",pep_res.m_vBetaProteinID[k]);
					strTXT += szBuf;		
				}
				// mod for pep2
				sprintf(szBuf, "%d ", pep_res.m_BetaPeptide.m_tModCnt);
				strTXT += szBuf;
				for(size_t k = 0; k < pep_res.m_BetaPeptide.m_tModCnt; k++)
				{
					sprintf(szBuf, "%d %d ",pep_res.m_BetaPeptide.m_tModSites[k][0] + 1, // add mod site 1. 2014.1.14.
							pep_res.m_BetaPeptide.m_tModSites[k][1]);
					strTXT += szBuf;			
				}
			}
			
			// add match info
			for(size_t k=0;k<pep_res.m_AlphaPeptide.m_tLength;++k)
			{
				sprintf(szBuf,"%d ",pep_res.m_stMatchInfo.aPepConf[k]);
				strTXT += szBuf;
			}

			if(pep_res.m_bPair)
			{
				for(size_t k=0;k<pep_res.m_BetaPeptide.m_tLength;++k)
				{
					sprintf(szBuf,"%d ",pep_res.m_stMatchInfo.aPepConf[pep_res.m_AlphaPeptide.m_tLength + k]);
					strTXT += szBuf;
				}
			}
			
			sprintf(szBuf,"%f %f",pep_res.m_stMatchInfo.lfMatchedSpecInt,pep_res.m_stMatchInfo.lfUnMatchedSpecInt);
			strTXT += szBuf;

			strTXT += "\n";
		}
	}
	
	ofstream of(strOutputFile.c_str());
	of<<strTXT;
	of.close();

	return true;

}
bool CXLinkPFDFileIO::LoadFile(string strOutputPath, vector<CXLinkMatchResult> & vResults,
		vector<CSpectrum> & vSpectra, size_t & tLoad)
{	
	ifstream ifIn(strOutputPath.c_str());
	string strTemp;
	size_t tTemp;
	CXLinkMatchResult mr;
	CXLinkPepResult pep_res;
	pep_res.m_bEV = true;
	CSpectrum spec;
	while(getline(ifIn, strTemp))
	{
		if(!ifIn)
			break;
		vSpectra.push_back(spec);
		vResults.push_back(mr);
		vSpectra[tLoad].m_strFilePath = strTemp;
		ifIn >> vSpectra[tLoad].m_nCharge >> vSpectra[tLoad].m_lfIntensity >> vSpectra[tLoad].m_lfMH
		>> vSpectra[tLoad].m_lfMZ >> vSpectra[tLoad].m_stEVCoef.lfCoef0 >> vSpectra[tLoad].m_stEVCoef.lfCoef1
		>> vSpectra[tLoad].m_nEvalueNo >> vResults[tLoad].m_tScore >> tTemp >> vResults[tLoad].m_tRealCandidate;
		
		if(0 == tTemp)
		{
			++tLoad;
			getline(ifIn, strTemp);
			continue;
		}
		for(size_t j = 0;j < tTemp;++j)
		{
			size_t tTempSize, tTempID;
			ifIn >> pep_res.m_lfScore >> pep_res.m_lfEvalue ;
			ifIn >> pep_res.m_lfAlphaEvalue >> pep_res.m_lfBetaEvalue; //added at 2014.9.10
			int nTmpType;
			ifIn >> nTmpType >> pep_res.m_XLink.m_tAlphaSite >> pep_res.m_XLink.m_tBetaSite >> pep_res.m_XLink.m_nLinkerId;
			pep_res.m_XLink.m_eXLinkType = (XLinkType)nTmpType;
			
			if(pep_res.m_XLink.m_eXLinkType <= 2)
			{
				pep_res.m_bPair = false;
			}
			else
			{
				pep_res.m_bPair = true;
			}
			
			ifIn >> pep_res.m_AlphaPeptide.m_cPrev >> pep_res.m_AlphaPeptide.m_szSequence
						 >> pep_res.m_AlphaPeptide.m_cNext >> pep_res.m_AlphaPeptide.m_lfMass >> pep_res.m_AlphaPeptide.m_tLength >> tTempSize;
			
			pep_res.m_vAlphaProteinID.clear();
			
			for(size_t k = 0; k < tTempSize; k++)
			{
				ifIn >> tTempID;
				pep_res.m_vAlphaProteinID.push_back(tTempID);	
			}
			
			ifIn >> pep_res.m_AlphaPeptide.m_tModCnt;
			
			for(size_t k = 0; k < pep_res.m_AlphaPeptide.m_tModCnt; k++)
			{
				ifIn >> pep_res.m_AlphaPeptide.m_tModSites[k][0] >> pep_res.m_AlphaPeptide.m_tModSites[k][1];		
			}
			
			if(pep_res.m_bPair)
			{
				ifIn >> pep_res.m_BetaPeptide.m_cPrev >> pep_res.m_BetaPeptide.m_szSequence
									 >> pep_res.m_BetaPeptide.m_cNext >> pep_res.m_BetaPeptide.m_lfMass >> pep_res.m_BetaPeptide.m_tLength >> tTempSize;
						
				pep_res.m_vBetaProteinID.clear();
				
				for(size_t k = 0; k < tTempSize; k++)
				{
					ifIn >> tTempID;
					pep_res.m_vBetaProteinID.push_back(tTempID);	
				}
				
				ifIn >> pep_res.m_BetaPeptide.m_tModCnt;
				
				for(size_t k = 0; k < pep_res.m_BetaPeptide.m_tModCnt; k++)
				{
					ifIn >> pep_res.m_BetaPeptide.m_tModSites[k][0] >> pep_res.m_BetaPeptide.m_tModSites[k][1];		
				}
			}
			
			// add match info

			for(size_t k=0;k<pep_res.m_AlphaPeptide.m_tLength;++k)
			{
				ifIn >> pep_res.m_stMatchInfo.aPepConf[k];
			}

			if(pep_res.m_bPair)
			{
				for(size_t k=0;k<pep_res.m_BetaPeptide.m_tLength;++k)
				{
					ifIn >> pep_res.m_stMatchInfo.aPepConf[pep_res.m_AlphaPeptide.m_tLength + k];
				}
			}
			
			ifIn >> pep_res.m_stMatchInfo.lfMatchedSpecInt >> pep_res.m_stMatchInfo.lfUnMatchedSpecInt;
			
			
			vResults[tLoad].m_vPeptideResults.push_back(pep_res);
		}
		getline(ifIn, strTemp);
		++tLoad;
		
	}
	ifIn.close();
	return true;
}

bool CXLinkPFDFileIO::LoadAll(string strOutput, int nSpectraTotal, int nTotalFile, vector<CXLinkMatchResult> & vResults,
		vector<CSpectrum> & vSpectra, size_t & tTotalSize, const string & strIdentifier)
{
	size_t tLoad = 0;
	vSpectra.clear();
	vResults.clear();
	vSpectra.reserve(nSpectraTotal);
	vResults.reserve(nSpectraTotal);
	if(strOutput[strOutput.length() - 1] != cSlash)
	{
		strOutput += cSlash;
	}
	char szBuf[20] = {0};
	for(int i = 0;i < nTotalFile;++i)
	{
		sprintf(szBuf, "load file: %d", i+1);
		m_pTrace->Debug(szBuf);
		string strFile = strOutput;
		sprintf(szBuf, "%d.%s.pfd", i, strIdentifier.c_str());
		strFile += szBuf;
		LoadFile(strFile, vResults, vSpectra, tLoad);
	}
	
	tTotalSize = tLoad;
	return true;
}
