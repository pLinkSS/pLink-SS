/*
 * MGFOutput.cpp
 *
 *  Created on: 2010-2-25
 *      Author: fan
 */

#include "MGFOutput.h"

//using namespace proteomics_sdk;
using namespace std;

namespace proteomics_sdk {

CMGFOutput::CMGFOutput():m_bMono(true),m_strOutPath(""), m_nOutCount(0),  m_tMaxPeaksNum(MAX_PEAKS_NUM)
{
	m_bMissCharge = false;
	m_buffString = "";
}

CMGFOutput::~CMGFOutput() {
	CloseOutFile();
}

void CMGFOutput::WriteAll(string strPath, vector<CSpectrum> & S)
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
//	S.clear();
	CSpectrum spec;
	int idx = 0;
	for(size_t i = 0;i < S.size();++i)
	{
//		S.pop_back(spec);
		spec = S[i];
		WriteNext(spec, idx);
	}
	EndWrite();
}

void CMGFOutput::StartWrite(string strPath,int & nTotal)
{
	try
	{
		OpenOutFile(strPath);
	}
	catch(exception & e)
	{
		CErrInfo info("CMGFOutput", "StartWrite", "OpenOutFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMGFOutput", "StartWrite", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}
	
	m_vCharge.clear();
	WriteFileHead();/*
	CloseInFile();

	m_strInPath = strPath;
	nTotal = 0;
	string strBFmark = "this is a string";
	const int MAX_BUFFER_SIZE = 81920;
	char szTemp[MAX_BUFFER_SIZE] = {0};

	FILE * fp = fopen(strPath.c_str(), "r");
	bool b = true;
	while(b)
	{
		if("BEGIN" == strBFmark.substr(0, 5))
		{
			while(fgets(szTemp, MAX_BUFFER_SIZE, fp) != NULL)
			{
				strBFmark = szTemp;
				if(strBFmark.substr(0, 3) == "END")
				{
					nTotal += m_vCharge.size();
					break;
				}
				else
				{
					if(strBFmark.substr(0, 7) == "CHARGE=")
					{
						++nTotal;
						break;
					}
				}
			}
		}
		b = (fgets(szTemp, MAX_BUFFER_SIZE, fp) != NULL);
		strBFmark = szTemp;
	}
//	//remove the head line CHARGE=
//	--nTotal;
	fclose(fp);
	OpenInFile(strPath);*/
	nTotal = 0;	//Currently set it to 0
	m_nOutCount = 0;
}

void CMGFOutput::OpenOutFile(string strFilePath)
{
	if (( m_ofOut = fopen(strFilePath.c_str(), "wt+") )== NULL)
	{
		CErrInfo info("CMGFOutput", "OpenOutFile", "m_ofOut.fail() when open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());		
	}
	/*m_ifIn.open(strFilePath.c_str());
	if(m_ifIn.fail())
	{
		CErrInfo info("CMGFInput", "OpenInFile", "m_ifIn.fail() when open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());
	}*/
	m_strOutPath = strFilePath;
}

void CMGFOutput::WriteFileHead()
{
	/*
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
	*/
	//Currently write nothin.
}

void CMGFOutput::WriteNext(CSpectrum & spec, int & idx)
{
	WriteMGFHead(spec);
	WriteMZAndItensity(spec, true);
	fprintf(m_ofOut, "%s\n\n", "END IONS");
	
}

void CMGFOutput::WriteMGFHead(CSpectrum& spec)
{
	fprintf(m_ofOut, "BEGIN IONS\n");
	fprintf(m_ofOut, "TITLE=%s\n", GetTitle(spec).c_str());
	fprintf(m_ofOut, "CHARGE=%d+\n", spec.m_nCharge);
	fprintf(m_ofOut, "PEPMASS=%lf\n", spec.m_lfMZ);
	

}

string CMGFOutput::GetTitle(CSpectrum& spec)
{
	return spec.m_strFilePath;
}

void CMGFOutput::WriteMZAndItensity(CSpectrum& spec, bool bAdmitBlank)
{
	for (size_t i = 0; i < spec.m_tPeaksNum; i++)
	{
		fprintf(m_ofOut, "%lf %lf\n", spec.m_pPeaks[i].lfMz, spec.m_pPeaks[i].lfIntensity);
	}
	/*
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
	spec.ReAlloc(spec.m_tPeaksNum);*/
}

void CMGFOutput::EndWrite(void)
{
	CloseOutFile();
}

void CMGFOutput::CloseOutFile(void)
{
	/*
	if(m_ifIn.is_open())
		m_ifIn.clear();
	m_ifIn.close();*/
	fclose(m_ofOut);
};
/*
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

*/

void CMGFOutput::WriteNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx)
{
	CSpectrum spec;
	

	for(int i = 0; i < nNumber; ++i)
	{
//		S.pop_back(spec);
		spec = S[i];
		WriteNext(spec, idx);	
	}
}


}
