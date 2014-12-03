/*
 * PKLInput.cpp
 *
 *  Created on: 2009-7-13
 *      Author: fan
 */


#include "PKLInput.h"

using namespace std;

namespace proteomics_sdk {

CPKLInput::CPKLInput():m_bMono(true), m_strInPath(""),  m_nInCount(0), m_tMaxPeaksNum(MAX_PEAKS_NUM)
{

}

CPKLInput::~CPKLInput()
{
	CloseInFile();
}

void CPKLInput::LoadAll(string strPath, vector<CSpectrum> & S)
{
	int nTot = 0;
	try
	{
		StartLoad(strPath.c_str(),nTot);
	}
	catch(exception & e)
	{
		CErrInfo info("CPKLInput", "LoadAll", "StartLoad failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPKLInput", "LoadAll", "caught an unknown exception from StartLoad().");
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

void CPKLInput::OpenInFile(string strFilePath)
{
	m_ifIn.open(strFilePath.c_str());
	if(m_ifIn.fail())
	{
		CErrInfo info("CPKLInput", "OpenInFile", "open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());
	}
	m_strInPath = strFilePath;
}

void CPKLInput::CloseInFile(void)
{
	if(m_ifIn.is_open())
		m_ifIn.clear();
	m_ifIn.close();
};


void CPKLInput::ReadPKLHead(CSpectrum& spec)
{

	string strLine;

	do
	{
		if(!GetLine(strLine,false))
		{
			CErrInfo info("CPKLInput", "ReadPKLHead", "GetLine failed.");
			throw runtime_error(info.Get().c_str());
		}
	}
	while(strLine=="");

	sscanf(strLine.c_str(),"%lf%lf%d", &spec.m_lfMZ, &spec.m_lfIntensity, &spec.m_nCharge);

	if(m_bMono)
		spec.m_lfMH = spec.m_lfMZ * spec.m_nCharge - (spec.m_nCharge - 1) * IonMass_Proton;
	else
		spec.m_lfMH = spec.m_lfMZ * spec.m_nCharge - (spec.m_nCharge - 1) * IonMass_Aver_H;


}


void CPKLInput::ReadMZAndItensity(CSpectrum& spec)
{
	string str = "init for not null";

	spec.clear();
	spec.ReAlloc(m_tMaxPeaksNum);
	CPeak peak;
	spec.m_tPeaksNum = 0;
	while(m_ifIn.good()
		&& !m_ifIn.eof()
		&& str != "")
	{
		if(!GetLine(str, false) || str=="")
		{
			break;
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

bool CPKLInput::GetLine(string& str, bool bAdmitBlank)
{

	while(!m_ifIn.eof())
	{
		if(!m_ifIn.good())
		{
			CErrInfo info("CPklIO", "GetLine", "test m_ifIn.good() failed.");
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
			/*
			cout << "continue : str = " << str << endl;
			char tmpch;
			cin >> tmpch;
			*/

			continue;
		}


		return true;
	}

	str.clear();

	return false;
}


void CPKLInput::StartLoad(string strPath,int & nTotal)
{
	try
	{
		OpenInFile(strPath);
	}
	catch(exception & e)
	{
		CErrInfo info("CPKLInput", "StartWrite", "OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPKLInput", "StartWrite", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}
	m_strInPath = strPath;
	nTotal = 0;
	string strBFmark;


	while(GetLine(strBFmark,true))
	{
		if(IsPKLHead(strBFmark))
			++nTotal;
	}
	
/*
	char szBuf[256] = {0};
	sprintf(szBuf, "Total Spectra: %d", nTotal);
	CTrace::GetInstance()->Info(szBuf, MODULE_SPECTRAIO);
*/
	CloseInFile();
	OpenInFile(strPath);
	m_nInCount = 0;
}
void CPKLInput::LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int & idx)
{
	S.clear();
	CSpectrum spec;

	for(int i = 0; i < nNumber; ++i)
	{
		LoadNext(spec, idx);
		S.push_back(spec);
	}
}

void CPKLInput::LoadNext(CSpectrum & spec, int & idx)
{
	ReadPKLHead(spec);

	ReadMZAndItensity(spec);



	//--------rename the single PKL file---------------//
	int pos = (int)m_strInPath.size();
	while(pos >= 0 && ('\\' != m_strInPath[pos] && '/' != m_strInPath[pos])) --pos;
	string strFileName = m_strInPath.substr(pos + 1,(int)m_strInPath.size() - pos - 1);
	char ptemp[PATH_MAX];
	sprintf(ptemp, "%d", m_nInCount);
	string nOrder(ptemp);
	//spec.m_strFilePath = strFileName + "." + nOrder + ".dta";//暂时不添加.dta后缀
	spec.m_strFilePath = strFileName + "." + nOrder ;
	//-------------end rename-------------------------//
	idx = m_nInCount;
	++m_nInCount;
}
void CPKLInput::EndLoad(void)
{
	CloseInFile();
}

bool CPKLInput::IsPKLHead(string strLine)
{
	istringstream isIn(strLine);
	string strSub1, strSub2, strSub3;
	isIn >> strSub1 >> strSub2 >> strSub3;
	//replace_if(strSub3.begin(), strSub3.end(), InvalidChar, "");
	return strSub3 != "";
}

void CPKLInput::LoadPrev(CSpectrum & spec, int &idx)
{

}

void CPKLInput::LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo)
{

}

vector<string> CPKLInput::GetAllSpecName(int begin, int step, int num)
{
	vector<string>vstrRet;
	return vstrRet;
}

}
