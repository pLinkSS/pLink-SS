/*
 * DTASInput.cpp
 *
 *  Created on: 2009-7-13
 *      Author: fan
 */

#include "DTASInput.h"

using namespace std;

namespace proteomics_sdk {

CDTASInput::CDTASInput() : m_bMono(true), m_strInPath(""), m_nInCount(0),  m_tMaxPeaksNum(MAX_PEAKS_NUM)
{

}

CDTASInput::~CDTASInput()
{
	CloseInFile();
}

void CDTASInput::LoadAll(string strPath, vector<CSpectrum> & S)
{
	int nTot = 0;
	try
	{
		StartLoad(strPath.c_str(),nTot);
	}
	catch(exception & e)
	{
		CErrInfo info("CDtasIO", "LoadAll", "StartLoad failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDtasIO", "LoadAll", "caught an unknown exception from StartLoad().");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}
	S.clear();
	CSpectrum spec;
	int idx = 0;
	for(int i = 0;i < nTot;++i)
	{
		LoadNext(spec, idx);
		S.push_back(spec);
	}
	EndLoad();
}


void CDTASInput::OpenInFile(string strFilePath)
{
	m_ifIn.open(strFilePath.c_str());
	if(m_ifIn.fail())
	{
		CErrInfo info("CDtasIO", "OpenInFile", "open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());
	}
	m_strInPath = strFilePath;
}

void CDTASInput::CloseInFile(void)
{
	if(m_ifIn.is_open())
		m_ifIn.clear();
	m_ifIn.close();
};

void CDTASInput::ReadDTAHead(CSpectrum& spec)
{
	string strLine;
	if(!GetLine(strLine,true))
	{
		CErrInfo info("CDtasIO", "ReadDTAHead", "GetLine failed.");
		throw runtime_error(info.Get().c_str());
	}

	sscanf(strLine.c_str(),"%lf %d", &spec.m_lfMH, &spec.m_nCharge);

	//calculate MZ
	if(m_bMono)
		spec.m_lfMZ = (spec.m_lfMH + (spec.m_nCharge - 1) * IonMass_Proton) / spec.m_nCharge;
	else
		spec.m_lfMZ = (spec.m_lfMH + (spec.m_nCharge - 1) * IonMass_Aver_H) / spec.m_nCharge;
	spec.m_lfIntensity = 0.0;
}


void CDTASInput::ReadMZAndItensity(CSpectrum& spec, bool bAdmitBlank)
{
	string str("init for not null");

	spec.clear();
	spec.ReAlloc(m_tMaxPeaksNum);


	CPeak peak;

	while(m_ifIn.good()
		&& !m_ifIn.eof()
		&& str != "")
	{
		if(!GetLine(str, bAdmitBlank))
		{
			return;
		}
		if("" == str)
			break;
		sscanf(str.c_str(),"%lf %lf", &peak.lfMz, &peak.lfIntensity);
		spec.m_pPeaks[spec.m_tPeaksNum++] = peak;
		if(spec.m_tPeaksNum == m_tMaxPeaksNum)
		{
			spec.ReAlloc(m_tMaxPeaksNum <<= 1);
		}
	}
	spec.ReAlloc(spec.m_tPeaksNum);
}

bool CDTASInput::GetLine(string& str, bool bAdmitBlank)
{
	while(!m_ifIn.eof())
	{
		if(!m_ifIn.good())
		{
			CErrInfo info("CDtasIO", "GetLine", "test m_ifIn.good() failed.");
			throw runtime_error(info.Get().c_str());
		}
		getline(m_ifIn, str);
		//remove 0a, 0d, ' '
		int tPos = (int)str.length() - 1;
		while(tPos >= 0 && (0xa == str[tPos] || 0xd == str[tPos] || '\n' == str[tPos]))
		{
			str.erase(str.begin() + tPos);
			--tPos;
		}

		if(bAdmitBlank && "" == str)
		{
			continue;
		}
		return true;
	}
	str.clear();
	return false;
}


bool CDTASInput::IsDTAFileName(string strFileName)
{
	while(strFileName.size() > 0 && strFileName[strFileName.size() - 1] == ' ')
	{
		strFileName.erase(strFileName.size() - 1, 1);
	}
	if(strFileName.size() <= 4)
		return false;

	string strExt = strFileName.substr(strFileName.size() - 4, 4);
	for(int i = 1;i < 4;++i)
	{
		tolower(strExt[i]);
	}
	return ".dta" == strExt;
}

bool CDTASInput::GetDtaFileName(string & strFileName)
{
	while(m_ifIn.good() && !m_ifIn.eof())
	{
		getline(m_ifIn, strFileName);
        if(IsDTAFileName(strFileName))
		{
			return true;
		}
	}
	strFileName.clear();
	return false;
}


void CDTASInput::StartLoad(string strPath,int & nTotal)
{
	try
	{
		OpenInFile(strPath);
	}
	catch(exception & e)
	{
		CErrInfo info("CDtasIO", "StartWrite", "OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDtasIO", "StartWrite", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}
	m_strInPath = strPath;
	nTotal = 0;
	string strFileName;
	while(GetDtaFileName(strFileName))
		nTotal++;
	
/*
	char szBuf[256] = {0};
	sprintf(szBuf, "Total Spectra: %d", nTotal);
	CTrace::GetInstance()->Info(szBuf, MODULE_SPECTRAIO);
	*/
	CloseInFile();
	OpenInFile(strPath);
	m_nInCount = 0;
}

void CDTASInput::LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int & idx)
{
	S.clear();
	CSpectrum spec;

	for(int i = 0; i < nNumber; ++i)
	{
		LoadNext(spec, idx);
		S.push_back(spec);
	}

}

void CDTASInput::LoadNext(CSpectrum & spec, int & idx)
{
	string strDtaFileName;
	if(!GetDtaFileName(strDtaFileName))
	{
		CErrInfo info("CDtasIO", "LoadNext", "GetDtaFileName failed.");
		throw runtime_error(info.Get().c_str());
	}

	ReadDTAHead(spec);
	ReadMZAndItensity(spec, false);
	spec.m_strFilePath = strDtaFileName.substr(0, strDtaFileName.length()-4);//暂时不添加.dta后缀
	idx = m_nInCount;
	++m_nInCount;
}

void CDTASInput::LoadPrev(CSpectrum & spec, int &idx)
{
	//TODO:Haven't done
}
void CDTASInput::LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo)
{
	//TODO:wait to be done.
}

void CDTASInput::EndLoad(void)
{
	CloseInFile();
}

vector<string> CDTASInput::GetAllSpecName(int begin, int step, int num)
{
	vector<string>vstrRet;
	return vstrRet;
}

MS2FormatType CDTASInput::GetType(void)
{
	return PFF_DTAS;
}
}
