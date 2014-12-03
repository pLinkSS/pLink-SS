/*
 * DTAInput.cpp
 *
 *  Created on: 2009-7-7
 *      Author: fan
 */

#include "DTAInput.h"
#include "SpecUtility.h"
 
CDTAInput::CDTAInput() : m_bMono(true), m_strInPath(""), m_nInCount(0), m_nTotal(0), m_tMaxPeaksNum(MAX_PEAKS_NUM)
{
	m_vstrPartFilePath.empty();
	m_vstrTotalFilePath.empty();
	m_mapTitle2Loc.clear();
}

CDTAInput::~CDTAInput() {
	CloseInFile();
}

MS2FormatType CDTAInput::GetType(void)
{
	return PFF_DTA;
}

void CDTAInput::LoadAll(string strPath, vector<CSpectrum> & S)
{
	int nTot = 0;
	try
	{
		StartLoad(strPath.c_str(),nTot);
	}
	catch(exception & e)
	{
		CErrInfo info("CDtaIO", "LoadAll", "StartLoad failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDtaIO", "LoadAll", "caught an unknown exception from StartLoad().");
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

void CDTAInput::StartLoad(string strPath, int & nTotal)
{
	proteomics_sdk::CFileFind ff;
	if(!ff.Open(strPath))
	{
		CErrInfo info("CDtaIO", "StartLoad", "CFileFind::open failed when open " + strPath + ".");
		throw runtime_error(info.Get().c_str());
	}
	if(strPath[strPath.length() - 1] != SLASH)
	{
		strPath += SLASH;
	}
	
	nTotal = 0;
	while(ff.GetNextFile())
	{
		size_t len = strlen(ff.m_pent->d_name);

		if(len < 4 || strcmp(string_toupper(ff.m_pent->d_name + len - 4).c_str(), ".DTA") != 0)
			continue;
		
		m_vstrTotalFilePath.push_back(strPath + ff.m_pent->d_name);
		m_vstrPartFilePath.push_back(ff.m_pent->d_name);
		m_mapTitle2Loc[ff.m_pent->d_name] = nTotal;
		++nTotal;
	}
	m_strInPath = strPath;
	m_nInCount = -1;
	m_nTotal = nTotal;
/*
	char szBuf[256] = {0};
	sprintf(szBuf, "Total Spectra: %d", nTotal);
	CTrace::GetInstance()->Info(szBuf, MODULE_SPECTRAIO);*/
}
void CDTAInput::LoadNext(CSpectrum & spec, int &idx)
{
	string strNewFile, strFilePath;
	if(m_nInCount+1 != m_nTotal)
		strFilePath = m_vstrTotalFilePath[m_nInCount + 1];
	else
		return;

	try
	{
		OpenInFile(string(strFilePath));
	}
	catch(exception & e)
	{
		CErrInfo info("CDtaIO", "LoadNext", "OpenInFile: Cannot open " + strFilePath + ".");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDtaIO", "LoadNext", "caught an unknown exception from OpenInFile: Cannot open " + strFilePath + ".");
		throw runtime_error(info.Get().c_str());
	}

	string strTemp = strFilePath;
	size_t pos = strTemp.find_last_of(SLASH);
//	strNewFile = strTemp.substr(pos + 1, strTemp.size() - 1 - pos).c_str();//暂时不添加.dta后缀
	strNewFile = strTemp.substr(pos + 1, strTemp.size() - 5 - pos).c_str();
	spec.m_strFilePath = strNewFile;
	try
	{
		ReadDTAHead(spec);
		ReadMZAndItensity(spec, true);
	}
	catch(exception & e)
	{
		CErrInfo info("CDtaIO", "LoadNext", "ReadDTAHead & ReadMZandIntensity: Cannot read " + spec.m_strFilePath + ".");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDtaIO", "LoadNext", "caught an unknown exception from ReadDTAHead & ReadMZandIntensity: Cannot read " + spec.m_strFilePath + ".");
		throw runtime_error(info.Get().c_str());
	}
	CloseInFile();
	++m_nInCount;
	idx = m_nInCount;
}

void CDTAInput::LoadPrev(CSpectrum & spec, int &idx)
{
	string strNewFile, strFilePath;
	if(m_nInCount >= 1)
		strFilePath = m_vstrTotalFilePath[m_nInCount - 1];
	else return;
	try
	{
		OpenInFile(strFilePath);
	}
	catch(exception & e)
	{
		CErrInfo info("CDtaIO", "LoadNext", "OpenInFile: Cannot open " + strFilePath + ".");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDtaIO", "LoadNext", "caught an unknown exception from OpenInFile: Cannot open " + strFilePath + ".");
		throw runtime_error(info.Get().c_str());
	}

	string strTemp = strFilePath;
	size_t pos = strTemp.find_last_of(SLASH);
	//strNewFile = strTemp.substr(pos + 1, strTemp.size() - 1 - pos).c_str();//暂时不添加.dta后缀
	strNewFile = strTemp.substr(pos + 1, strTemp.size() - 5 - pos).c_str();	
	spec.m_strFilePath = strNewFile;
	try
	{
		ReadDTAHead(spec);
		ReadMZAndItensity(spec, true);
	}
	catch(exception & e)
	{
		CErrInfo info("CDtaIO", "LoadNext", "ReadDTAHead & ReadMZandIntensity: Cannot read " + spec.m_strFilePath + ".");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDtaIO", "LoadNext", "caught an unknown exception from ReadDTAHead & ReadMZandIntensity: Cannot read " + spec.m_strFilePath + ".");
		throw runtime_error(info.Get().c_str());
	}
	CloseInFile();
	m_nInCount--;
	idx = m_nInCount;
}

void CDTAInput::LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx)
{
	S.clear();
	S.reserve(nNumber);

	CSpectrum spec;
	for(int i = 0; i < nNumber; ++i)
	{
		LoadNext(spec, idx);
		S.push_back(spec);
	}
}
void CDTAInput::LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo)
{
	//wait to be done.
	string strNewFile;
	
	try
	{
		OpenInFile(m_strInPath+strTitle);
	}
	catch(exception & e)
	{
		CErrInfo info("CDtaIO", "LoadNext", "OpenInFile: Cannot open " + m_strInPath + strTitle + ".");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDtaIO", "LoadNext", "caught an unknown exception from OpenInFile: Cannot open " + m_strInPath + strTitle + ".");
		throw runtime_error(info.Get().c_str());
	}

	spec.m_strFilePath = strTitle;
	try
	{
		ReadDTAHead(spec);
		ReadMZAndItensity(spec, true);
	}
	catch(exception & e)
	{
		CErrInfo info("CDtaIO", "LoadNext", "ReadDTAHead & ReadMZandIntensity: Cannot read " + spec.m_strFilePath + ".");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDtaIO", "LoadNext", "caught an unknown exception from ReadDTAHead & ReadMZandIntensity: Cannot read " + spec.m_strFilePath + ".");
		throw runtime_error(info.Get().c_str());
	}
	CloseInFile();
	m_nInCount = m_mapTitle2Loc[strTitle];
	nSpecNo = m_nInCount;
}
void CDTAInput::EndLoad(void)
{
	++m_nInCount;
}

void CDTAInput::LoadDTA(const char* pszDTAFilePath, const char* pszFileOnly, vector<CSpectrum> & S)
{
	string strDTAFile(pszDTAFilePath);
	OpenInFile(strDTAFile);
	CSpectrum spec;
	spec.m_strFilePath = pszFileOnly;

	ReadDTAHead(spec);
	ReadMZAndItensity(spec, true);
	if(spec.m_tPeaksNum > 0)
	{
		S.push_back(spec);
	}

	CloseInFile();
}

void CDTAInput::OpenInFile(string strFilePath)
{
	m_ifIn.open(strFilePath.c_str());
	if(m_ifIn.fail())
	{
		CErrInfo info("CDtaIO", "OpenInFile", "open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());
	}
	//m_strInPath = strFilePath;
}

void CDTAInput::CloseInFile(void)
{
	if(m_ifIn.is_open())
		m_ifIn.clear();
	m_ifIn.close();
};

void CDTAInput::ReadDTAHead(CSpectrum& spec)
{
	string strLine;
	if(!GetLine(strLine,true))
	{
		CErrInfo info("CDtaIO", "ReadDTAHead", "GetLine failed.");
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

void CDTAInput::ReadMZAndItensity(CSpectrum& spec, bool bAdmitBlank)
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
		sscanf(str.c_str(),"%lf %lf", &peak.lfMz, &peak.lfIntensity);
		spec.m_pPeaks[spec.m_tPeaksNum++] = peak;
		
		if(spec.m_tPeaksNum == m_tMaxPeaksNum)
		{
			spec.ReAlloc(m_tMaxPeaksNum <<= 1);
		}
	}
	spec.ReAlloc(spec.m_tPeaksNum);
}

bool CDTAInput::GetLine(string& str, bool bAdmitBlank)
{
	while(!m_ifIn.eof())
	{
		if(!m_ifIn.good())
		{
			CErrInfo info("CDtaIO", "GetLine", "test m_ifIn.good() failed.");
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


vector<string> CDTAInput::GetAllSpecName(int begin, int step, int num)
{
	vector<string>vstrRet;
	for(size_t i = 0,  j = (size_t)begin ; j<m_vstrPartFilePath.size() && i<(size_t)num ; ++i, j+=step)
	{
		vstrRet.push_back(m_vstrPartFilePath[j]);
	}

	return vstrRet;
}
