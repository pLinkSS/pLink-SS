/*
 * DTAOutput.cpp
 *
 *  Created on: 2010-3-22
 *      Author: fan
 */
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include "DTAOutput.h"

CDTAOutput::CDTAOutput() : m_bMono(true), m_strOutPath(""), m_nOutCount(0), m_tMaxPeaksNum(MAX_PEAKS_NUM)
{
}

CDTAOutput::~CDTAOutput() {
	CloseOutFile();
}

MS2FormatType CDTAOutput::GetType(void)
{
	return PFF_DTA;
}

void CDTAOutput::WriteAll(string strPath, vector<CSpectrum> & S)
{
	int nTot = 0;
	try
	{
		StartWrite(strPath.c_str(),nTot);
	}
	catch(exception & e)
	{
		CErrInfo info("CDTAOutput", "WriteAll", "StartWrite failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDTAOutput", "WriteAll", "caught an unknown exception from StartWrite().");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}

	CSpectrum spec;
	int idx = 0;
	for(size_t i = 0;i < S.size();++i)
	{
		spec = S[i];
		WriteNext(spec, idx);
	}
	EndWrite();
}

void CDTAOutput::StartWrite(string strPath, int & nTotal)
{
	string strDir = strPath;
	if(strDir[strPath.length() - 1] != SLASH)
	{
		strDir += SLASH;
	}
	
	m_strOutPath = strDir;		
	m_nOutCount = 0;
#ifdef WIN32
	mkdir(m_strOutPath.c_str());
#else
	mkdir(m_strOutPath.c_str(), 777);
#endif
}

void CDTAOutput::WriteNext(CSpectrum & spec, int &idx)
{
	WriteDTA(m_strOutPath, spec);
	++m_nOutCount;
	idx = m_nOutCount;
}


void CDTAOutput::WriteDTA(string strPath, CSpectrum & spec)
{
	string strEntireName = strPath + spec.m_strFilePath + ".dta";
	try
	{
		OpenOutFile(strEntireName);
	}
	catch(exception & e)
	{
		CErrInfo info("CDTAOutput", "WriteDTA", "OpenOutFile failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDTAOutput", "WriteDTA", "caught an unknown exception from OpenOutFile().");
		info.Append("strEntireName=" + strEntireName);
		throw runtime_error(info.Get().c_str());
		
	}
	
	WriteDTAHead(spec);
	WriteMZAndItensity(spec);
	CloseOutFile();
}

void CDTAOutput::WriteDTAHead(CSpectrum& spec)
{
	fprintf(m_ofOut, "%lf %d\n", spec.m_lfMH, spec.m_nCharge);
}

void CDTAOutput::WriteMZAndItensity(CSpectrum& spec, bool bAdmitBlank)
{
	for (size_t i = 0; i < spec.m_tPeaksNum; i++)
	{
		fprintf(m_ofOut, "%lf %lf\n", spec.m_pPeaks[i].lfMz, spec.m_pPeaks[i].lfIntensity);
	}
}
void CDTAOutput::WriteNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx)
{
	CSpectrum spec;
	
	for(int i = 0; i < nNumber; ++i)
	{
		spec = S[i];
		WriteNext(spec, idx);	
	}
}
void CDTAOutput::EndWrite(void)
{
	++m_nOutCount;
}

void CDTAOutput::OpenOutFile(string strFilePath)
{
	if (( m_ofOut = fopen(strFilePath.c_str(), "wt+") )== NULL)
	{
		CErrInfo info("CDTAOutput", "OpenOutFile", "m_ofOut.fail() when open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());		
	}
}

void CDTAOutput::CloseOutFile(void)
{
	fclose(m_ofOut);
};

