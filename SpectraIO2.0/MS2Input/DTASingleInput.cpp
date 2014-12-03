
#include "DTASingleInput.h"

CDTASingleInput::CDTASingleInput(void)
{
}

CDTASingleInput::~CDTASingleInput(void)
{
}
void CDTASingleInput::StartLoad(string strPath, int & nTotal)
{
	try
	{
		OpenInFile(strPath);
	}
	catch(exception & e)
	{
		CErrInfo info("CDTASingleInput", "StartWrite", "OpenInFile: " + strPath + " failed.");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDTASingleInput", "StartWrite", "caught an unknown exception from OpenInFile: " + strPath + " failed.");
		info.Append("strPath=" + strPath);
		throw runtime_error(info.Get().c_str());
	}
	m_strInPath = strPath;
	nTotal = 1;

	CloseInFile();
//	OpenInFile(strPath);
}

void CDTASingleInput::LoadAll(string strPath, vector<CSpectrum> & S)
{

	try
	{
		OpenInFile(string(strPath));
	}
	catch(exception & e)
	{
		CErrInfo info("CDTASingleInput", "LoadAll", "OpenInFile: Cannot open " + strPath + ".");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDTASingleInput", "LoadAll", "caught an unknown exception from OpenInFile: Cannot open " + strPath + ".");
		throw runtime_error(info.Get().c_str());
	}

	CSpectrum spec;
	try
	{
		ReadMHAndCharge(spec);
		ReadMZAndItensity(spec);
	}
	
	catch(exception & e)
	{
		CErrInfo info("CDTASingleInput", "LoadAll", "ReadMHAndCharge & ReadMZandIntensity: Cannot read " + spec.m_strFilePath + ".");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDTASingleInput", "LoadAll", "caught an unknown exception from ReadMHAndCharge & ReadMZandIntensity: Cannot read " + spec.m_strFilePath + ".");
		throw runtime_error(info.Get().c_str());
	}

	S.clear();
	S.push_back(spec);

	CloseInFile();
}
void CDTASingleInput::LoadNext(CSpectrum & spec, int & idx)
{

	try
	{
		OpenInFile(string(m_strInPath));
	}
	catch(exception & e)
	{
		CErrInfo info("CDTASingleInput", "LoadNext", "OpenInFile: Cannot open " + m_strInPath + ".");
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDTASingleInput", "LoadNext", "caught an unknown exception from OpenInFile: Cannot open " + m_strInPath + ".");
		throw runtime_error(info.Get().c_str());
	}

	ReadMHAndCharge(spec);
	ReadMZAndItensity(spec);

	CloseInFile();
}
void CDTASingleInput::LoadNextBatch(int nNumber, vector<CSpectrum> & S, int &idx)
{
	idx = 0;
	LoadAll(m_strInPath, S);
}
void CDTASingleInput::LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo)
{
	//wait to be done.
}
void CDTASingleInput::LoadPrev(CSpectrum & spec, int &idx)
{
	//wait to be done.
}
void CDTASingleInput::EndLoad(void)
{
	CloseInFile();
}
void CDTASingleInput::CloseInFile(void)
{
	if(m_ifIn.is_open())
		m_ifIn.clear();
	m_ifIn.close();
}
void CDTASingleInput::OpenInFile(string strFilePath)
{
	m_ifIn.open(strFilePath.c_str());
	if(m_ifIn.fail())
	{
		CErrInfo info("CDTASingleInput", "OpenInFile", "open file: " + strFilePath + " failed.");
		throw runtime_error(info.Get().c_str());
	}
	m_strInPath = strFilePath;
}
void CDTASingleInput::ReadMHAndCharge(CSpectrum& spec)
{
	string strLine;
	if(!GetLine(strLine,true))
	{
		CErrInfo info("CDTASingleInput", "ReadDTAHead", "GetLine failed.");
		throw runtime_error(info.Get().c_str());
	}
	sscanf(strLine.c_str(),"%lf %d", &spec.m_lfMH, &spec.m_nCharge);
}
void CDTASingleInput::ReadMZAndItensity(CSpectrum& spec)
{
	string str = "init for not null";

	spec.clear();
	spec.ReAlloc(MAX_PEAKS_NUM);
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

		size_t PeakNumMax = MAX_PEAKS_NUM;
		if((unsigned)MAX_PEAKS_NUM == spec.m_tPeaksNum)
		{
			spec.ReAlloc(PeakNumMax <<= 1);
		}

	}
	spec.ReAlloc(spec.m_tPeaksNum);
}
bool CDTASingleInput::GetLine(string& str, bool bAdmitBlank)
{
	while(!m_ifIn.eof())
	{
		if(!m_ifIn.good())
		{
			CErrInfo info("CDTASingleInput", "GetLine", "test m_ifIn.good() failed.");
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

vector<string> CDTASingleInput::GetAllSpecName(int begin, int step, int num)
{
	vector<string>vstrRet;
	vstrRet.push_back(m_strInPath);
	return vstrRet;
}