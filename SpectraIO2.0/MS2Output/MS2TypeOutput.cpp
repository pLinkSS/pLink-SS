/*
 * MS2TypeOutput.cpp
 *
 *  Created on: 2010-3-11
 *      Author: fan
 */

#include "MS2TypeOutput.h"

//using namespace proteomics_sdk;
using namespace std;

namespace proteomics_sdk {

CMS2TypeOutput::CMS2TypeOutput():m_bMono(true),m_strOutPath(""), m_nOutCount(0),  m_tMaxPeaksNum(MAX_PEAKS_NUM)
{

}

CMS2TypeOutput::~CMS2TypeOutput() {
	CloseOutFile();
}

void CMS2TypeOutput::WriteAll(string strPath, vector<CSpectrum> & S)
{
	int nTot = 0;
	try
	{
		StartWrite(strPath.c_str(),nTot);
	}
	catch(exception & e)
	{
		CErrInfo info("CMGFOutput", "WriteAll", "StartLoad failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMGFOutput", "WriteAll", "caught an unknown exception from StartLoad().");
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

void CMS2TypeOutput::StartWrite(string strPath,int & nTotal)
{
	try
	{
		OpenOutFile(strPath);
	}
	catch(exception & e)
	{
		CErrInfo info("CMS2TypeOutput", "StartWrite", "OpenOutFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMS2TypeOutput", "StartWrite", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}
	
	WriteFileHead();
	nTotal = 0;	//Currently set it to 0
	m_nOutCount = 0;
}

void CMS2TypeOutput::OpenOutFile(string strFilePath)
{
	if (( m_ofOut = fopen(strFilePath.c_str(), "wt+") )== NULL)
	{
		CErrInfo info("CMS2TypeOutput", "OpenOutFile", "m_ofOut.fail() when open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());		
	}

	m_strOutPath = strFilePath;
}

void CMS2TypeOutput::WriteFileHead()
{
	fprintf(m_ofOut, "H	Extractor	pXtract\n");
	fprintf(m_ofOut, "H	ExtractorVersion	2.0\n");
	fprintf(m_ofOut, "H	Comments	pXtract written by pFind Studio, 2010\n");
}

void CMS2TypeOutput::WriteNext(CSpectrum & spec, int & idx)
{
	WriteMS2Head(spec);
	WriteMZAndItensity(spec, true);	
}

void CMS2TypeOutput::WriteMS2Head(CSpectrum& spec)
{
	fprintf(m_ofOut, "S\t%06d\t%06d\t%lf\n", spec.m_tScanNo, spec.m_tScanNo, spec.m_lfMZ);
	fprintf(m_ofOut, "I	RetTime	%lf\n", spec.m_lfRetentionTime);
	fprintf(m_ofOut, "I	ActivationType	%s\n", spec.m_strActivationType.c_str());
	fprintf(m_ofOut, "I	InstrumentType	%s\n", spec.m_strInstrumentType.c_str());
	fprintf(m_ofOut, "Z	%d	%lf\n", spec.m_nCharge, spec.m_lfMH);
}


void CMS2TypeOutput::WriteMZAndItensity(CSpectrum& spec, bool bAdmitBlank)
{
	for (size_t i = 0; i < spec.m_tPeaksNum; i++)
	{
		fprintf(m_ofOut, "%lf %lf\n", spec.m_pPeaks[i].lfMz, spec.m_pPeaks[i].lfIntensity);
	}
}

void CMS2TypeOutput::EndWrite(void)
{
	CloseOutFile();
}

void CMS2TypeOutput::CloseOutFile(void)
{

	fclose(m_ofOut);
};


void CMS2TypeOutput::WriteNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx)
{
	CSpectrum spec;
	
	for(int i = 0; i < nNumber; ++i)
	{
		spec = S[i];
		WriteNext(spec, idx);	
	}
}


}
