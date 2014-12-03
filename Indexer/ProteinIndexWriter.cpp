#include <string>
#include "../include/sdk.h"

#include "../Mass2PepIndex/Mass2PepIndex.h"
#include "ProteinIndexWriter.h"

#define DATABASE_FILE_BUF_SIZE 40960
#define SQ_WRITE_LINE_SIZE 60

#define AC_LENGTH_PRO 100
#define DE_LENGTH_PRO 80

using namespace std;
using namespace ProteinIndex;
using namespace Mass2PepIndex;
using namespace proteomics_sdk;

CProteinIndexWriter::CProteinIndexWriter()
:m_tDepth(0)
{
}

CProteinIndexWriter::~CProteinIndexWriter()
{
}

void CProteinIndexWriter::Init(const PINDEX_HEAD &pIndexHead, const PINDEX_ITEM& pIndexItem) 
{
	m_PIndexHead = pIndexHead;
	m_PIndexItem = pIndexItem;
}

void CProteinIndexWriter::Write(CProteinWriter* pproIdxWriter) 
{
	_ShowFBeginInfo("CProteinIndexWriter::Write");
	double lf = clock();
	

	if(1 == m_PIndexItem.nAutoReverse)
	{
		_WriteTargetDecoyData(m_PIndexItem.strFastaFile, "");
		m_PIndexItem.strFastaFile = m_PIndexItem.strFastaFile.substr(0, m_PIndexItem.strFastaFile.size()-6) + "_comb_rever.fasta";
	}
#ifdef CONSOLE
	cout << "\tOut path: " << m_PIndexHead.strOutpath << endl;
	cout << "\tDB name: " << m_PIndexItem.strDBName << endl;
	cout << "\tstrFastaFile:"<< (m_PIndexItem.strFastaFile )<< endl;
	cout << endl;
#endif
	
	pproIdxWriter->Init(m_PIndexHead.strOutpath.c_str(), m_PIndexItem.strDBName.c_str(), m_PIndexItem.strFastaFile.c_str(), m_PIndexItem.tMaxFileSize,m_PIndexHead,m_PIndexItem);
	try
	{
		if(m_PIndexItem.bPro) pproIdxWriter->WriteProteinIndex();
		if(m_PIndexItem.tDatFlow) pproIdxWriter->WriteProteinSuffix(m_PIndexItem.tDatFlow);
		_ShowFEndInfo("CProteinIndexWriter::Write", lf);
		
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinIndexWriter", "Write");
		err_info.Append("FASAT Name = " + (string)(m_PIndexItem.strFastaFile));
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinIndexWriter", "Write","Unknown Error!");
		err_info.Append("FASAT Name = " + (string)(m_PIndexItem.strFastaFile));
		throw runtime_error(err_info.Get().c_str());
	}
}

CProtein CProteinIndexWriter::_ReadNext(ifstream &m_ifIn)
{
	string str, strSQ, strAC, strDE;

	char* pszBuf = new char[DATABASE_FILE_BUF_SIZE];
	m_ifIn.getline(pszBuf, DATABASE_FILE_BUF_SIZE - 1);
	str = pszBuf;
	
	_ChangeInvalid(str);

	size_t sizeLen = str.length();

	size_t tPos = str.find(' ', 0);
	if ( 0 == sizeLen || pszBuf[0] != '>' )
		cout << "Read One Protein Entry procedure is error! strlen == 0 or first char != >!"<<endl;

	if(tPos > AC_LENGTH_PRO)
		strAC = str.substr(1, AC_LENGTH_PRO);

	else
		strAC = str.substr(1, tPos - 1);


	if(tPos + DE_LENGTH_PRO < sizeLen)
		strDE = str.substr(tPos + 1, DE_LENGTH_PRO);
	else
		strDE = str.substr(tPos + 1);
	_ChangeInvalid(strDE);
	if(0 == strDE.size()) strDE = strAC;
	if('>' == strDE[0]) strDE = strDE.substr(1, strDE.size() - 1);

	strSQ = _ReadEntrySQ(m_ifIn,pszBuf);
	_ChangeInvalid(strSQ);

	delete[] pszBuf;

	CProtein pro(strAC, strDE, strSQ);

	return pro;
}

//read entry for SQ
string CProteinIndexWriter::_ReadEntrySQ(ifstream &m_ifIn, char*  pszBuf)
{
	string strSQ, strLine;

	if (m_ifIn.fail()) {
		do {
			m_ifIn.clear();
			m_ifIn.getline(pszBuf, DATABASE_FILE_BUF_SIZE - 1);
		} while(!m_ifIn.eof() && m_ifIn.fail());
	}

	while(!m_ifIn.eof() 
		&& m_ifIn.peek() != '>' ) 
	{
		m_ifIn.getline(pszBuf, DATABASE_FILE_BUF_SIZE - 1);
		
		strLine = pszBuf;
		strLine.erase(0,strLine.find_first_not_of("\r\t\n "));// added at 2013.12.13
		strLine.erase(strLine.find_last_not_of("\r\t\n ")+1);
		strSQ += strLine;

		if (m_ifIn.eof())		//avoid end line is blank!
			continue;

		if (m_ifIn.fail()) {		//avoid large line!
			m_ifIn.clear();
		}


	}
	return strSQ;
}

void CProteinIndexWriter::_WriteProtein(ofstream &ofFile, const CProtein &protein, const size_t nSQSize, const size_t line)
{
	ofFile << '>' << protein.m_strAC << ' ' << protein.m_strDE << endl;

	for(size_t k = 0; k < line ; k++)
		ofFile << protein.m_strSQ.substr(k*SQ_WRITE_LINE_SIZE, SQ_WRITE_LINE_SIZE) << endl;
	ofFile << protein.m_strSQ.substr(line*SQ_WRITE_LINE_SIZE) << endl;
}

void CProteinIndexWriter::_WriteTargetDecoyData(string strSourcePath, string strNewPath)
{
	_ShowFBeginInfo("CProteinIndexWriter::WriteTargetDecoyData");
	double lf = clock();
	
	ifstream m_ifIn;
	m_ifIn.open(strSourcePath.c_str());
	if(m_ifIn.fail())
	{
		proteomics_sdk::CErrInfo err_info("CProteinIndexWriter","WriteTargetDecoyData");
		throw runtime_error(err_info.Get().c_str());
	}

	if(strNewPath.empty())
		strNewPath =  strSourcePath.substr(0, strSourcePath.size()-6) +"_comb_rever.fasta";
	
	ofstream ofFile(strNewPath.c_str());

	m_ifIn.clear(); 
	m_ifIn.seekg(0,std::ios::beg);
	
	for(size_t i=0; !m_ifIn.eof(); ++i)
	{

		CProtein protein = _ReadNext(m_ifIn);
		size_t nSQSize = protein.m_strSQ.size();
		size_t j = nSQSize/ SQ_WRITE_LINE_SIZE;
		_WriteProtein(ofFile, protein, nSQSize, j);

		CProtein reverse = protein.Reverse();
		_WriteProtein(ofFile, reverse, nSQSize, j);
		/////////////////////////////////////////////////////////////////////////////////////
	}

	ofFile.close();

	if(m_ifIn.is_open())
	{
		m_ifIn.clear();
		m_ifIn.close();
	}
	_ShowFEndInfo("CProteinIndexWriter::WriteTargetDecoyData", lf);
	return;
}

void CProteinIndexWriter::_ChangeInvalid(string& strContent)
{
	string chFace = "";
//	_ChangeChar(strContent, chFace);
	chFace = "<";
	_ChangeChar(strContent, chFace);
	chFace = "&";
	_ChangeChar(strContent, chFace);
	chFace = "  ";
	_ChangeChar(strContent, chFace);
	chFace = "\r";
	_ChangeChar(strContent, chFace);
}

void CProteinIndexWriter::_ChangeChar(string& strContent, string chChanged)
{
	string::size_type idx = 0;
	idx = strContent.find(chChanged, idx);
	while (idx != string::npos ) {
//		strContent[idx] = ' ';
		strContent.erase(idx,1);
		idx = strContent.find(chChanged, idx);
	}
}

void	CProteinIndexWriter::_ShowFBeginInfo(string strFName)
{
	try
	{
#ifdef CONSOLE
		for(size_t t=0; t<=m_tDepth; ++t)
		{
			cout << "--" ;
		}
		
		cout << "Begin " << strFName << " " ;
		
		for(size_t t=0; t<=m_tDepth; ++t)
		{
			cout << "--" ;
		}
		
		cout << endl << endl;
		++m_tDepth;
#endif
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CProteinIndexWriter", "_ShowFBeginInfo");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CProteinIndexWriter", "_ShowFBeginInfo");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}	
}

void	CProteinIndexWriter::_ShowFEndInfo(string strFName, double dfTime)
{
	try
	{
#ifdef CONSOLE
		cout << "Time spent on " << strFName << " is:" << (clock() - dfTime)/1000 << " seconds." << endl << endl;
			
		for(size_t t=0; t<m_tDepth; ++t)
		{
			cout << "--" ;
		}
		
		cout << "End of " << strFName << " " ;
		
		for(size_t t=0; t<m_tDepth; ++t)
		{
			cout << "--" ;
		}
		
		cout << endl << endl;
		--m_tDepth;
#endif
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CProteinIndexWriter", "_ShowFEndInfo");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CProteinIndexWriter", "_ShowFEndInfo");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}


