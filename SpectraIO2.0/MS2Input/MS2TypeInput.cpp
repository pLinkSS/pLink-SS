/*
 * MS2TypeInput.cpp
 *
 *  Created on: 2009-7-14
 *      Author: fan
 */


#define USEMEMORY 1000000

#ifdef WIN32
	const char cSlash = '\\';
#else
	const char cSlash = '/';
#endif

#include "MS2TypeInput.h"
#include "SpecUtility.h"

namespace proteomics_sdk {

CMS2TypeInput::CMS2TypeInput():m_bMono(true), m_strInPath(""), m_strFileName(""), m_nInCount(0), m_nTotal(0), m_tMaxPeaksNum(MAX_PEAKS_NUM), strBFmark("")
{
	m_mapTitle2Loc.clear();
	m_vSpecLoc.clear();
}

CMS2TypeInput::~CMS2TypeInput() {
	CloseInFile();
}

void CMS2TypeInput::OpenInFile(string strFilePath)
{
	m_ifIn.open(strFilePath.c_str());
	if(m_ifIn.fail())
	{
		CErrInfo info("CMS2TypeInput", "OpenInFile", "m_ifIn.fail() when open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());
	}
	m_strInPath = strFilePath;
}

void CMS2TypeInput::GetScanNo(char *p, int &iScanNoFrom, int &iScanNoTo)
{
	int i = 0;
	iScanNoFrom = 0, iScanNoTo = 0;
	while(p[i])
	{
		if(p[i]>='0'&&p[i]<='9')
			break;
		++i;
	}
	while(p[i]>='0'&&p[i]<='9')
	{
		iScanNoFrom *=10;
		iScanNoFrom += p[i]-'0';
		++i;
	}
	++i;
	while(p[i]>='0'&&p[i]<='9')
	{
		iScanNoTo *=10;
		iScanNoTo += p[i]-'0';
		++i;
	}
}

int CMS2TypeInput::GetCharge(char *p)
{
	int i = 2, ret = 0;
	while(p[i]>='0'&&p[i]<='9')
	{
		ret *= 10;
		ret += p[i]-'0';
		++i;
	}
	return ret;
}

void CMS2TypeInput::StartLoad(string strPath, int & nTotal)
{
	/*try
	{
		OpenInFile(strPath);
	}
	catch(exception & e)
	{
		CErrInfo info("CMS2TypeInput", "StartLoad", "OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMS2TypeInput", "StartLoad", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}
	
	m_strFileName = strPath;
	size_t iPlace = m_strFileName.rfind("\\");
	m_strFileName.erase(0, iPlace+1);
	m_strFileName.replace(m_strFileName.length()-4, 4, "");
	
//计算共有谱图数
	nTotal = 0;

	bool b = true;
	string strBFmark = "this is a string";
	string strLine = "";

	while(b)
	{
		b = GetLine(strLine, true);
		if(b && strLine.substr(0, 1) == "S")
		{
			++nTotal;
		}
	}

	if(m_ifIn.is_open())
		m_ifIn.clear();
	m_ifIn.close();*/
	
	nTotal = 0;
	
	char Buff[USEMEMORY];
	
	FILE * fp;
	
	try
	{
		fp = fopen(strPath.c_str(), "rb");
	}
	catch(exception & e)
	{
		CErrInfo info("CMS2TypeInput", "StartLoad", "OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMS2TypeInput", "StartLoad", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}
	
	m_strFileName = strPath;
	size_t iPlace = m_strFileName.rfind(cSlash);
	m_strFileName.erase(0, iPlace+1);
	m_strFileName.replace(m_strFileName.length()-4, 4, "");
	
	int ByteCnt = 0, iTmpByteCnt = 0, iTmpScanNoFrom = 0, iTmpScanNoTo = 0;
	m_vSpecLoc.clear();
	m_mapTitle2Loc.clear();
	while(fgets(Buff,USEMEMORY,fp))
	{
		if('S'==Buff[0] && '\t' == Buff[1])
		{
			iTmpByteCnt = ByteCnt;
			m_vSpecLoc.push_back(ByteCnt);
			GetScanNo(Buff, iTmpScanNoFrom, iTmpScanNoTo);
		}
		else if('Z'==Buff[0] && '\t' == Buff[1])
		{
			string strFilePath;
			char c_strFilePath[100] = {0};
			sprintf(c_strFilePath, ".%d.%d.%d", iTmpScanNoFrom, iTmpScanNoTo, GetCharge(Buff));

			strFilePath.assign(c_strFilePath);

			strFilePath = m_strFileName + strFilePath;
			strFilePath = string_toupper(strFilePath);
			m_mapTitle2Loc[strFilePath]=nTotal++;
		}
		ByteCnt += (int)(strlen(Buff));
	}

	fclose(fp);
	
/*
	char szBuf[256] = {0};
	sprintf(szBuf, "Total Spectra: %d", nTotal);
	CTrace::GetInstance()->Info(szBuf, MODULE_SPECTRAIO);
	*/
//关了重开
	
	try
	{
		OpenInFile(strPath);
	}
	catch(exception & e)
	{
		CErrInfo info("CMS2TypeInput", "StartLoad", "OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMS2TypeInput", "StartLoad", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}

	//ReadMS2Head();	//暂时不用这里跳，在LoadNext里跳

	m_nTotal = nTotal;
	m_nInCount = 0;
	
	/*
	try
	{
		m_ifFile = fopen(strPath.c_str(), "r");
	}
	catch(exception & e)
	{
		CErrInfo info("CMS2TypeInput", "StartLoad", "OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMS2TypeInput", "StartLoad", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}
	m_nInCount = 0;
	*/
}


//跳过所有H部分
void CMS2TypeInput::ReadFileHead()
{
	string strLine;
	bool b = GetLine(strLine, true);
	while(b && strLine.substr(0, 1) != "S")
	{
		b = GetLine(strLine, true);
	}
}

bool CMS2TypeInput::GetLine(string & str, bool bAdmitBlank)
{
	
	while(!m_ifIn.eof())
	{
		if(!m_ifIn.good())
		{
			CErrInfo info("CMGFInput", "GetLine", "test m_ifIn.good() failed.");
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
	
	/*
	while(!feof(m_ifFile))
	{
		char buff[USEMEMORY] = {0};
		fgets( buff, USEMEMORY, m_ifFile );
		str = buff;
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
	*/
}

void CMS2TypeInput::LoadNext(CSpectrum & spec, int &idx)
{
	if(m_nInCount >= m_nTotal)
		return ;

	ReadMS2Head(spec);

	ReadMZAndItensity(spec, false);


	idx = m_nInCount;
	++m_nInCount;

}

void CMS2TypeInput::LoadPrev(CSpectrum & spec, int &idx)
{
	if(m_nInCount-2>=0)
	{
		if(m_ifIn.eof())
		{
			m_ifIn.clear();
		}
		m_ifIn.seekg(m_vSpecLoc[m_nInCount-2]);//TODO:seekg
		GetLine(strBFmark, true);
		LoadNext(spec, idx);
		m_nInCount-=2;
	}
}

void CMS2TypeInput::ReadMS2Head(CSpectrum & spec)
{

//	cout<<"Enter ReadMS2Head"<<endl;

	bool b = true;	//用于记录获取下一行是否成功

	while (b && (strBFmark.substr(0, 1) != "S"))
	{
		b = GetLine(strBFmark, true);
	}

	spec.m_tScanNo =atoi(strBFmark.substr(2, 6).c_str());
	spec.m_lfMZ = atof(strBFmark.substr(16, strBFmark.length()-16).c_str());

	b = GetLine(strBFmark, true);

	while (b && "I" == strBFmark.substr(0, 1))
	{
		if ("RetTime" == strBFmark.substr(2, 7))
		{
			spec.m_lfRetentionTime =atof(strBFmark.substr(10, 8).c_str());
		}

		if ("ActivationType" == strBFmark.substr(2, 14))
		{
			spec.m_strActivationType = strBFmark.substr(17, 3);
		}

		if ("InstrumentType" == strBFmark.substr(2, 14))
		{
			spec.m_strInstrumentType = strBFmark.substr(17, 4);
		}

		if ("PrecursorScan" == strBFmark.substr(2, 12))
		{
			spec.m_tPrecursorScanNo = atoi(strBFmark.substr(15, strBFmark.length()-15).c_str());
		}

		b = GetLine(strBFmark, true);
	}


	if (b && "Z" == strBFmark.substr(0, 1))
	{
		spec.m_nCharge = atoi(strBFmark.substr(2, 1).c_str());
		if (strBFmark.substr(3, 1) != "\t")
		{
			//cout<<strBFmark.substr(2, 2)<<endl;
			spec.m_nCharge = atoi(strBFmark.substr(2, 2).c_str());
		}
		b =  GetLine(strBFmark, true);
	}

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

//	cout<<"Exit ReadMS2Head"<<endl;

}

void CMS2TypeInput::ReadMZAndItensity(CSpectrum& spec, bool bAdmitBlank)
{

//	cout<<"Enter ReadMZAndItensity"<<endl;
	bool b = true;

	while (!(b && isdigit(*strBFmark.substr(0, 1).c_str())))
		b = GetLine(strBFmark, true);

	CPeak peak;
	spec.clear();
	spec.ReAlloc(m_tMaxPeaksNum);

	while(m_ifIn.good()
		&& !m_ifIn.eof()
		&& strBFmark.substr(0, 1) != "S")
	{
		sscanf(strBFmark.c_str(),"%lf %lf", &peak.lfMz, &peak.lfIntensity);
		spec.m_pPeaks[spec.m_tPeaksNum++] = peak;
		if(spec.m_tPeaksNum == m_tMaxPeaksNum)
		{
			spec.ReAlloc(m_tMaxPeaksNum <<= 1);
		}
		GetLine(strBFmark, bAdmitBlank);
	}

	spec.ReAlloc(spec.m_tPeaksNum);
//	cout<<"Exit ReadMZAndItensity"<<endl;
}

void CMS2TypeInput::LoadAll(string strPath, vector<CSpectrum> & S)
{
	int nTot = 0;
	try
	{
		StartLoad(strPath.c_str(),nTot);
	}
	catch(exception & e)
	{
		CErrInfo info("CMGFInput", "LoadAll", "StartLoad failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMGFInput", "LoadAll", "caught an unknown exception from StartLoad().");
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

void CMS2TypeInput::LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx)
{
	S.clear();
	CSpectrum spec;

	for(int i = 0; i < nNumber; ++i)
	{
		LoadNext(spec, idx);
		S.push_back(spec);
	}
}

void CMS2TypeInput::LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo)
{
	int idx = 0;
	m_nInCount = m_mapTitle2Loc[strTitle];
	nSpecNo = m_nInCount;
	if(m_ifIn.eof())
		m_ifIn.clear();
	m_ifIn.seekg(m_vSpecLoc[m_nInCount]);//TODO:seekg
	GetLine(strBFmark, true);
	LoadNext(spec,idx);
}
vector<string> CMS2TypeInput::GetAllSpecName(int begin, int step, int num)
{
	vector<string>vstrRet;

	map<string, int>::iterator it;
	int i = 0;
	for(it = m_mapTitle2Loc.begin() ; i < begin && it!=m_mapTitle2Loc.end() ; ++it, ++i)
	{

	}
	for(int j = 0; j<num && it!=m_mapTitle2Loc.end() ; ++j)
	{
		vstrRet.push_back(it->first);
		int i = 0;
		while(i<step)
			++it, ++i;
	}
	return vstrRet;
}
void CMS2TypeInput::EndLoad(void)
{
	
	CloseInFile();
	
	/*
	fclose(m_ifFile);*/
}

void CMS2TypeInput::CloseInFile()
{
	if(m_ifIn.is_open())
		m_ifIn.clear();
	m_ifIn.close();
}

}
