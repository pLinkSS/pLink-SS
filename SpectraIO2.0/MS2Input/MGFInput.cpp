/*
 * MGFInput.cpp
 *
 *  Created on: 2009-7-13
 *      Author: fan
 */
#include "MGFInput.h"
#include "SpecUtility.h"

//using namespace proteomics_sdk;
using namespace std;

 
namespace proteomics_sdk {

CMGFInput::CMGFInput():m_bMono(true),m_strInPath(""), m_nInCount(0),  m_tMaxPeaksNum(MAX_PEAKS_NUM)
{
	m_bMissCharge = false;
	m_buffString = "";
	m_lfLineCount = 0;
}

CMGFInput::~CMGFInput() {
	CloseInFile();
}

void CMGFInput::LoadAll(string strPath, vector<CSpectrum> & S)
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

void CMGFInput::OpenInFile(string strFilePath)
{
	m_ifIn.open(strFilePath.c_str());
	if(m_ifIn.fail())
	{
		CErrInfo info("CMGFInput", "OpenInFile", "m_ifIn.fail() when open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());
	}
	m_strInPath = strFilePath;
}

void CMGFInput::CloseInFile(void)
{
	if(m_ifIn.is_open())
		m_ifIn.clear();
	m_ifIn.close();
};

void CMGFInput::ReadFileHead()
{
	string strLine;
	bool b = GetLine(strLine, true);
	while(b && strLine.substr(0, 5) != "BEGIN")
	{
		if(strLine.substr(0, 7) == "CHARGE=")
		{
			strLine = strLine.substr(7, strLine.size() - 7);
			strLine += ' ';
			string temp;
			for(size_t i = 0;i < strLine.size();++i)
			{
				if(!isalpha(strLine[i]) && !isdigit(strLine[i]) && strLine[i] != '+'
					&& strLine[i] != '-')
				{
					if(temp != "")
					{
						if(temp != "and")
						{
							char c = temp[temp.size() - 1];
							temp.substr(0, temp.size() - 1);
							m_vCharge.push_back(atoi(temp.c_str()) * (c == '+'?1:-1));
						}
						temp = "";
					}
				}
				else
				{
					temp += strLine[i];
				}
			}
		}
		b = GetLine(strLine, true);
	}
}

void CMGFInput::ReadMGFHead(CSpectrum& spec)
{
	string strLine;

	m_bMissCharge = true;
	
	m_buffString = "";
	
	GetLine(strLine, true);
	
	while((*(strLine.substr(0, 1).c_str())>'9' || *(strLine.substr(0,1).c_str()) < '0') && "END" != strLine.substr(0, 3))
	{
//		cout<<*(strLine.substr(0, 1).c_str())<<endl;
//		cout<<strLine.substr(0, 3)<<endl;
		if("TITLE=" == strLine.substr(0,6))
		{
			spec.m_strFilePath = strLine.substr(6,(int)strLine.size() - 6);
		}
		else if("CHARGE=" == strLine.substr(0,7))
		{
			spec.m_nCharge = atoi(strLine.substr(7,(int)strLine.size() - 7).c_str());
			m_bMissCharge = false;
		}
		else if("PEPMASS=" == strLine.substr(0,8))
		{
			spec.m_lfMZ = atof(strLine.substr(8,(int)strLine.size() - 8).c_str());
		}
		else if ("RETTIME=" == strLine.substr(0, 8))
		{
			spec.m_lfRetentionTime = atof(strLine.substr(8, (int)strLine.size() -8).c_str());
		}
		
		GetLine(strLine,true);
		
	}

	if ((*(strLine.substr(0, 1).c_str())<='9' && *(strLine.substr(0,1).c_str()) >= '0'))
		m_buffString = strLine;
	else
		m_buffString = "END";
	
	if(m_bMissCharge)
	{
		if(m_vCharge.size() == 0)
		{
			spec.m_nCharge = 2;	//todo :暂时将遗漏电荷的谱图电荷记为2电荷
			
			char szBuf[255]={0};
			sprintf(szBuf, "Missing charge value on line %.0lf", m_lfLineCount-1);
			
			CTrace::GetInstance()->Alert(szBuf, MODULE_SPECTRAIO);

			m_bMissCharge = false;
		}
		else
		{
			m_strMissCharge = strLine;
			m_vCurrCharge = m_vCharge;
		}
	}


	spec.m_lfIntensity = 0.0;
	if(m_bMono)
		spec.m_lfMH = spec.m_lfMZ * spec.m_nCharge - (spec.m_nCharge - 1) * IonMass_Proton;
	else
		spec.m_lfMH = spec.m_lfMZ * spec.m_nCharge - (spec.m_nCharge - 1) * IonMass_Aver_H;
	
//	cout<<"Now out of the MGFHead, the strLine now is "<<strLine<<endl;
}


void CMGFInput::ReadMZAndItensity(CSpectrum& spec, bool bAdmitBlank)
{
	string str = m_buffString;
//	cout<<"The str now is "<<str<<endl;

	if(m_bMissCharge)
		str = m_strMissCharge;
	else if (str == "END")
		return;


//	cout<<"Now out of the choice, the strLine now is "<<str<<endl;

	CPeak peak;
	spec.clear();
	spec.ReAlloc(m_tMaxPeaksNum);
	while(m_ifIn.good()
		&& !m_ifIn.eof()
		&& str.substr(0, 3) != "END")
	{
		sscanf(str.c_str(),"%lf %lf", &peak.lfMz, &peak.lfIntensity);
		
		spec.m_pPeaks[spec.m_tPeaksNum++] = peak;
		if(spec.m_tPeaksNum == m_tMaxPeaksNum)
		{
			spec.ReAlloc(m_tMaxPeaksNum <<= 1);
		}
		GetLine(str, bAdmitBlank);
	}
	spec.ReAlloc(spec.m_tPeaksNum);
}

bool CMGFInput::GetLine(string& str, bool bAdmitBlank)
{
	while(!m_ifIn.eof())
	{
		if(!m_ifIn.good())
		{
			CErrInfo info("CMGFInput", "GetLine", "test m_ifIn.good() failed.");
			throw runtime_error(info.Get().c_str());
		}
		getline(m_ifIn, str);
		m_lfLineCount = 1 + m_lfLineCount;
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

void CMGFInput::StartLoad(string strPath,int & nTotal)
{
	try
	{
		OpenInFile(strPath);
	}
	catch(exception & e)
	{
		CErrInfo info("CMGFInput", "StartLoad", "OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMGFInput", "StartWrite", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}
	
	m_vSpecLoc.clear();
	m_mapTitle2Loc.clear();
	m_vCharge.clear();
	ReadFileHead();
	CloseInFile();

	m_strInPath = strPath;
	nTotal = 0;
	string strBFmark = "this is a string";
	const int MAX_BUFFER_SIZE = 81920;
	char szTemp[MAX_BUFFER_SIZE] = {0};

	FILE * fp = fopen(strPath.c_str(), "rb");
	bool b = true;
	int ByteCnt = 0, istrLen, iTmpByteCnt = 0;
	
	//判断是否为空文件
	struct stat statbuf;
  fstat(fileno(fp), &statbuf);
	if (0 == statbuf.st_size)
	{
		fclose(fp);
		nTotal = 0;
		return;
	}
	
//	int iLineBegin = 0;
	
	bool bGetPrecursor = 0;//是否获得了电荷信息
	while(b)
	{
		
		if("BEGIN" == strBFmark.substr(0, 5))
		{
			
			while(fgets(szTemp, MAX_BUFFER_SIZE, fp))
			{
				istrLen = strlen(szTemp);
				strBFmark = szTemp;
				ByteCnt += istrLen ;
				if(strBFmark.substr(0, 3) == "END")
				{
					if (!bGetPrecursor)
					{
						nTotal += m_vCharge.size();
					}
					bGetPrecursor = 0;
					break;
				}
				else if(strBFmark.substr(0, 7) == "CHARGE=")
				{
					bGetPrecursor = 1;
					++nTotal;
				}
				else if("TITLE=" == strBFmark.substr(0, 6))
				{
					strBFmark = string_toupper(strBFmark);
					if(strBFmark[strBFmark.size()-1]=='\n')
						strBFmark.erase(strBFmark.size()-1, 1);
					if(strBFmark[strBFmark.size()-1]=='\r')
						strBFmark.erase(strBFmark.size()-1, 1);

					while(m_vSpecLoc.size()>(unsigned)nTotal)//no charge
						++nTotal;
					m_mapTitle2Loc[strBFmark.substr(6, strBFmark.size()-6)] = nTotal;
					while(m_vSpecLoc.size()<(unsigned)nTotal)//no charge
						m_vSpecLoc.push_back(0);

					m_vSpecLoc.push_back(iTmpByteCnt);
				}
			}
		}
		iTmpByteCnt = ByteCnt;
		b = (fgets(szTemp, MAX_BUFFER_SIZE, fp) != NULL);
		istrLen = strlen(szTemp);
		ByteCnt += istrLen ;
		strBFmark = szTemp;
	}

//	//remove the head line CHARGE=
//	--nTotal;
	fclose(fp);

	/*
	char szBuf[256] = {0};
	sprintf(szBuf, "Total Spectra: %d", nTotal);
	CTrace::GetInstance()->Info(szBuf, MODULE_SPECTRAIO);
	*/
	OpenInFile(strPath);
	m_nInCount = 0;
	m_lfLineCount = 0;
	m_nTotal = nTotal;
}

void CMGFInput::LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx)
{
	S.clear();
	CSpectrum spec;
	

	for(int i = 0; i < nNumber; ++i)
	{
		LoadNext(spec, idx);
		S.push_back(spec);
	}
}

void CMGFInput::LoadNext(CSpectrum & spec, int & idx)
{
	if(m_nInCount >= m_nTotal)
		return ;

	if(m_bMissCharge && m_vCurrCharge.size() != m_vCharge.size())
	{
		if(!m_vCurrCharge.empty())
		{
			spec = m_specTemp;
			spec.m_nCharge = m_vCurrCharge.back();
			m_vCurrCharge.pop_back();
			spec.m_lfIntensity = 0.0;
		if(m_bMono)
			spec.m_lfMH = spec.m_lfMZ * spec.m_nCharge - (spec.m_nCharge - 1) * IonMass_Proton;
		else
			spec.m_lfMH = spec.m_lfMZ * spec.m_nCharge - (spec.m_nCharge - 1) * IonMass_Aver_H;
			idx = m_nInCount;
			++m_nInCount;
			if(spec.m_strFilePath == "")
			{
				//--------rename the single file---------------//
				int pos = (int)m_strInPath.size();
				while(pos >= 0 && ('\\' != m_strInPath[pos] && '/' != m_strInPath[pos])) --pos;
				string strFileName = m_strInPath.substr(pos + 1,(int)m_strInPath.size() - pos - 1);
				char * ptemp = new char [PATH_MAX];
				sprintf(ptemp, "%d", m_nInCount);
				string nOrder(ptemp);
	//		spec.m_strFilePath = strFileName + "." + nOrder + ".dta";//暂时不添加.dta后缀
				spec.m_strFilePath = strFileName + "." + nOrder;
				delete[] ptemp;
			}
			else
			{
				if(".0.dta" == spec.m_strFilePath.substr(spec.m_strFilePath.size() - 6, 6))
					spec.m_strFilePath = spec.m_strFilePath.substr(0, spec.m_strFilePath.size() - 6);
				else if(".dta" == spec.m_strFilePath.substr(spec.m_strFilePath.size() - 4, 4))
					spec.m_strFilePath = spec.m_strFilePath.substr(0, spec.m_strFilePath.size() - 4);
				char p[205];
				sprintf(p, "%d", spec.m_nCharge);
				string temp(p);
				spec.m_strFilePath += ".";
				spec.m_strFilePath += p;
	//		spec.m_strFilePath += ".dta";//暂时不添加.dta后缀
			}

			//-------------end rename-------------------------//
			return;
		}
		else
			m_bMissCharge = false;
	}

	string strBFmark;
	while(GetLine(strBFmark, true))
		if("BEGIN" == strBFmark.substr(0, 5))
			break;
			
	ReadMGFHead(spec);
	ReadMZAndItensity(spec, true);

	idx = m_nInCount;
	++m_nInCount;

	if(m_bMissCharge)
	{
		m_vCurrCharge = m_vCharge;
		spec.m_nCharge = m_vCurrCharge.back();
		m_vCurrCharge.pop_back();
		spec.m_lfIntensity = 0.0;
		if(m_bMono)
			spec.m_lfMH = spec.m_lfMZ * spec.m_nCharge - (spec.m_nCharge - 1) * IonMass_Proton;
		else
			spec.m_lfMH = spec.m_lfMZ * spec.m_nCharge - (spec.m_nCharge - 1) * IonMass_Aver_H;
		m_specTemp = spec;
	}

	if(spec.m_strFilePath == "")
	{
		//--------rename the single file---------------//
		int pos = (int)m_strInPath.size();
		while(pos >= 0 && ('\\' != m_strInPath[pos] && '/' != m_strInPath[pos])) --pos;
		string strFileName = m_strInPath.substr(pos + 1,(int)m_strInPath.size() - pos - 1);
		char * ptemp = new char [PATH_MAX];
		sprintf(ptemp, "%d", m_nInCount);
		string nOrder(ptemp);
//	spec.m_strFilePath = strFileName + "." + nOrder + ".dta";//暂时不添加.dta后缀
		spec.m_strFilePath = strFileName + "." + nOrder;
		delete[] ptemp;
	}
	else
	{
		////暂时不添加.dta后缀
		
		if(m_bMissCharge)
		{
			if(".0.dta" == spec.m_strFilePath.substr(spec.m_strFilePath.size() - 6, 6))
				spec.m_strFilePath = spec.m_strFilePath.substr(0, spec.m_strFilePath.size() - 6);
			else if(".dta" == spec.m_strFilePath.substr(spec.m_strFilePath.size() - 4, 4))
				spec.m_strFilePath = spec.m_strFilePath.substr(0, spec.m_strFilePath.size() - 4);
			char p[205];
			sprintf(p, "%d", spec.m_nCharge);
			string temp(p);
			spec.m_strFilePath += ".";
			spec.m_strFilePath += p;
//		spec.m_strFilePath += ".dta";
		}
/*	else
			if(".dta" != spec.m_strFilePath.substr(spec.m_strFilePath.size() - 4, 4))
				spec.m_strFilePath += ".dta";*/
		
	}
	//-------------end rename-------------------------//

}
void CMGFInput::LoadPrev(CSpectrum & spec, int &idx)
{
	if(m_nInCount-2>=0)
	{
		if(m_ifIn.eof())
			m_ifIn.clear();
		m_ifIn.seekg(m_vSpecLoc[m_nInCount-2]);//TODO:seekg may have the size limit
		LoadNext(spec, idx);
		m_nInCount-=2;
	}
	else
		return;
}
void CMGFInput::LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo)
{
	//wait to be done.
	int idx = 0;
	m_nInCount = m_mapTitle2Loc[strTitle];
	nSpecNo = m_nInCount;
	if(m_ifIn.eof())
		m_ifIn.clear();
	m_ifIn.seekg(m_vSpecLoc[m_nInCount]);
	LoadNext(spec,idx);
}
vector<string> CMGFInput::GetAllSpecName(int begin, int step, int num)
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
void CMGFInput::EndLoad(void)
{
	CloseInFile();
}

}
