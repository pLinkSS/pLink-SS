/*
 * RAWInput.cpp
 *
 *  Created on: 2009-7-21
 *      Author: fan
 */

#ifdef WIN32
	const char cSlash = '\\';
#else
	const char cSlash = '/';
#endif
	

#include "RAWInput.h"
#include <stdio.h>
#include <stdlib.h>

namespace proteomics_sdk {

CRAWInput::CRAWInput():m_bMono(true), m_strInPath(""), m_strFileName(""), m_nInCount(0),m_tMaxPeaksNum(MAX_PEAKS_NUM) {

}

CRAWInput::~CRAWInput() {
	CloseInFile();
}

void CRAWInput::OpenInFile(string strFilePath)
{
//	m_ifIn.open(strFilePath.c_str(), ios::binary);
	m_ifIn.open(strFilePath.c_str());
	if(m_ifIn.fail())
	{
		CErrInfo info("CMS2TypeInput", "OpenInFile", "m_ifIn.fail() when open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());
	}

//	OSVERSIONINFO OsVersionInfo;
	m_strInPath = strFilePath;
}

void CRAWInput::StartLoad(string strPath, int & nTotal)
{
	//TODO: 加入判断操作系统代码

	string command = "RawInput.exe "+ strPath;

	ReplaceSlash(command);

	cout<<"Reading raw, please wait...."<<endl;
	system(command.c_str());
	cout<<"Reading finish."<<endl;

	m_strFileName = strPath;
	size_t iPlace = m_strFileName.rfind(cSlash);
	m_strFileName.erase(0, iPlace+1);
	m_strFileName.replace(m_strFileName.length()-4, 4, "");
	
	strPath.replace(strPath.length()-3, 3, "pfd");
	


//	cout<<"Before the openinfile"<<endl;
	try
	{
		OpenInFile(strPath);
	}
	catch(exception & e)
	{
		CErrInfo info("CRAWTypeInput", "StartLoad", "OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CRAWTypeInput", "StartLoad", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}

//	cout<<"After the openinfile"<<endl;
	//计算共有谱图数
	nTotal = 0;

	bool b = true;

	string strLine = "";

//	cout<<"Before the caculate of spectrum number"<<endl;
	while(b)
	{
		b = GetLine(strLine, true);
		if(b && strLine.substr(0, 1) == "H")
		{
			++nTotal;
		}
	}

	if(m_ifIn.is_open())
		m_ifIn.clear();
	m_ifIn.close();

/*
	char szBuf[256] = {0};
	sprintf(szBuf, "Total Spectra: %d", nTotal);
	CTrace::GetInstance()->Info(szBuf, MODULE_SPECTRAIO);
	*/
//	关了重开
	try
	{
		OpenInFile(strPath);
	}
	catch(exception & e)
	{
		CErrInfo info("CRAWTypeInput", "StartLoad", "OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CRAWTypeInput", "StartLoad", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}

	//ReadMS2Head();	//暂时不用这里跳，在LoadNext里跳

	m_nInCount = 0;
	strBFmark = "";
}

bool CRAWInput::GetLine(string & str, bool bAdmitBlank)
{
	while(!m_ifIn.eof())
	{
		if(!m_ifIn.good())
		{
			CErrInfo info("CRAWInput", "GetLine", "test m_ifIn.good() failed.");
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


void CRAWInput::LoadNext(CSpectrum & spec, int &idx)
{
	ReadPFDHead(spec);
	ReadMZAndItensity(spec, false);


	idx = m_nInCount;
	++m_nInCount;
}

void CRAWInput::ReadPFDHead(CSpectrum & spec)
{
	bool b = true;	//用于记录获取下一行是否成功

	while (b && (strBFmark.substr(0, 1) != "H"))
	{
		b = GetLine(strBFmark, true);
	}

	GetLine(strBFmark, true);

	sscanf(strBFmark.c_str(), "%d\t%lf", &spec.m_tScanNo, &spec.m_lfRetentionTime);

	GetLine(strBFmark, true);

	char cActivationType[10000] = {0};
	char cInstrumentType[10000] = {0};
	sscanf(strBFmark.c_str(), "%s\t%s", cActivationType, cInstrumentType);
	spec.m_strActivationType = cActivationType;
	spec.m_strInstrumentType = cInstrumentType;

	GetLine(strBFmark, true);

	if (strBFmark[0] == '\0')
		strBFmark.erase(strBFmark.begin());
	sscanf(strBFmark.c_str(), "%d%d%lf", &spec.m_tPrecursorScanNo, &spec.m_nCharge, &spec.m_lfMZ);
	
	char c_strFilePath[100] = {0};
//	sprintf(c_strFilePath, ".%d.%d.%d.dta", spec.m_tScanNo, spec.m_tScanNo, spec.m_nCharge);//暂时不添加.dta后缀
	sprintf(c_strFilePath, ".%d.%d.%d", spec.m_tScanNo, spec.m_tScanNo, spec.m_nCharge);
	
	spec.m_strFilePath.assign(c_strFilePath);
	
	spec.m_strFilePath = m_strFileName + spec.m_strFilePath;

	if(m_bMono)
		spec.m_lfMH = spec.m_lfMZ * spec.m_nCharge - (spec.m_nCharge - 1) * IonMass_Proton;
	else
		spec.m_lfMH = spec.m_lfMZ * spec.m_nCharge - (spec.m_nCharge - 1) * IonMass_Aver_H;

	return;
}

void CRAWInput::ReadMZAndItensity(CSpectrum& spec, bool bAdmitBlank)
{
	bool b = true;	//用于记录获取下一行是否成功

	while (b && (strBFmark.substr(0, 1) != "M"))
	{
		b = GetLine(strBFmark, true);
	}

	while (!((b && isdigit(*strBFmark.substr(0, 1).c_str()) ) || strBFmark.substr(0, 1) == "H"))
		b = GetLine(strBFmark, true);

	CPeak peak;
	spec.clear();
	spec.ReAlloc(m_tMaxPeaksNum);

	while(m_ifIn.good()
		&& !m_ifIn.eof()
		&& strBFmark.substr(0, 1) != "H")
	{
		sscanf(strBFmark.c_str(),"%lf\t%lf", &peak.lfMz, &peak.lfIntensity);
		spec.m_pPeaks[spec.m_tPeaksNum++] = peak;
		if(spec.m_tPeaksNum == m_tMaxPeaksNum)
		{
			spec.ReAlloc(m_tMaxPeaksNum <<= 1);
		}
		GetLine(strBFmark, bAdmitBlank);
	}

	spec.ReAlloc(spec.m_tPeaksNum);
}

void CRAWInput::LoadAll(string strPath, vector<CSpectrum> & S)
{
	int nTot = 0;
	try
	{
		StartLoad(strPath.c_str(),nTot);
	}
	catch(exception & e)
	{
		CErrInfo info("CRAWInput", "LoadAll", "StartLoad failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CRAWInput", "LoadAll", "caught an unknown exception from StartLoad().");
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

void CRAWInput::LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx)
{
	S.clear();
	CSpectrum spec;

	for(int i = 0; i < nNumber; ++i)
	{
		LoadNext(spec, idx);
		S.push_back(spec);
	}
}

void CRAWInput::EndLoad(void)
{
	CloseInFile();
}

void CRAWInput::CloseInFile()
{
	if(m_ifIn.is_open())
		m_ifIn.clear();
	m_ifIn.close();
}

void CRAWInput::ReplaceSlash(string & command)
{
	int i = 0;

	while (1)
	{
		i = command.find("\\", i);
		if (i == -1)
			break;
		command.replace(i, 1, "\\\\");
		i+=2;
	}

	return;
}

void CRAWInput::LoadPrev(CSpectrum & spec, int &idx)
{

}

void CRAWInput::LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo)
{

}

vector<string> CRAWInput::GetAllSpecName(int begin, int step, int num)
{
	vector<string>vstrRet;
	return vstrRet;
}

}
