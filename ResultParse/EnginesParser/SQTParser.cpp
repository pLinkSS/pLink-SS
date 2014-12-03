#include "SQTParser.h"

//using namespace proteomics_sdk;
using namespace bio_analysis;

extern ostringstream osspBuildLog;

CSQTParser::CSQTParser()
{
}

CSQTParser::~CSQTParser()
{
}

void CSQTParser::_ParserFront(ifstream & fin, CParseredDataSetInfo & pbres)
{
	CConditionInfo cond;
	//ToDo:Add code to make up of CConditionInfo
	pbres.m_vConditions.push_back(cond);
	//pbres.m_vCondition = cond;
}

void CSQTParser::_GetSQTMod(string & strVal, vector<CModificationSiteInfo> & vModTemp)
{
	string strTemp;
	size_t k = 0;
	CModificationSiteInfo ModTemp;
	for (size_t t = 0; t < strVal.size(); t++)
	{
		if (strVal[t] >= 'A' && strVal[t] <= 'Z')
		{
			strTemp += strVal[t];
			k++;
		}
		else if (strVal[t] == '(')
		{
			ModTemp.m_strModName.clear();
			//ModTemp.m_tPos = k - 1;todo
			ModTemp.m_tPos = k;
			for (; t < strVal.size(); t++)
			{
				ModTemp.m_strModName += strVal[t];
				if (strVal[t] == ')')
					break;
			}
			vModTemp.push_back(ModTemp);
		}
	}
	strVal.clear();
	strVal = strTemp;
}

void CSQTParser::_ParseSpectra(ifstream & fin, CParseredDataSetInfo & pbres)
{
	CMatchSpectraInfo mr;
	CMatchPeptideInfo mp;
	int FisrtScan, LastScan;
	int FirstM = 1, FirstS = 1;
	char strTmp[100];
	double lfTmp, score, E_Value;
	int nTmp;
	mp.clear();
	mr.clear();
	int SpectraNum = 0;
	while (1)
	{
		string strTemp, strVal;
		bool eof = GetLine(fin, strTemp);
		if (eof == false)
		{
			sort(mp.m_vMod.begin(), mp.m_vMod.end(), Mod_Sort);
			DeleteRedundancy(mp.m_vProteinAC);
			mr.m_vPeptides.push_back(mp);
			mp.clear();
			pbres.m_mapResultData[mr.m_nCharge].push_back(mr);
			mr.clear();
			//cout << SpectraNum << endl; system
			osspBuildLog << SpectraNum << endl;
#ifdef DEBUG
			cout << SpectraNum << endl;
#endif
			break;
		}
		//cout << strTemp << endl;
		DeleteFrontBlank(strTemp);
		if ('S' == strTemp[0])
		{
			//cout << strTemp << endl;
			if (FirstS == 1)
				FirstS = 0;
			else
				pbres.m_mapResultData[mr.m_nCharge].push_back(mr);
			mr.clear();
			FirstM = 1;
			sscanf(strTemp.c_str(), "S %d %d %d %lf %s %lf", &FisrtScan, &LastScan, &mr.m_nCharge,
					&lfTmp, strTmp, &mr.m_lfMH);
			mr.m_nSpectraNum = ++SpectraNum;
			mr.m_strFileName = GetSpectraFileName(m_strFilePath, FisrtScan, LastScan, mr.m_nCharge);
		}
		else if ('M' == strTemp[0])
		{
			if (FirstM == 1)
				FirstM = 0;
			else
			{
				sort(mp.m_vMod.begin(), mp.m_vMod.end(), Mod_Sort);
				DeleteRedundancy(mp.m_vProteinAC);
				mr.m_vPeptides.push_back(mp);
			}
			mp.clear();
			char chr, SQtmp[1000];
			sscanf(strTemp.c_str(), "M %d %d %lf %lf %lf %lf %d %d %s %c", &mp.m_tRank, &nTmp,
					&mp.m_lfCalc_MH, &E_Value, &score, &lfTmp, &nTmp, &nTmp, SQtmp, &chr);
			mp.m_vlfScores.push_back(score);
			mp.m_vlfScores.push_back(E_Value);
			if (mr.m_vPeptides.size() == 1)
				mr.m_vPeptides[0].m_vlfScores[1] = E_Value;
			mp.m_lfDelta = mr.m_lfMH - mp.m_lfCalc_MH;
			mp.m_lfPPM = CBioMethods::TransDeltaMH2PPM(mp.m_lfDelta, mr.m_lfMH);
			//mp.m_lfCalc_M = mp.m_lfCalc_MH - protonH;
			mp.m_cNext = mp.m_cPrev = '-';
			strVal = SQtmp;
			//cout << strVal << endl;
			size_t posNP = strVal.find_first_of('.');
			if (posNP != string::npos)
				mp.m_cPrev = strVal[posNP - 1];
			strVal = strVal.substr(posNP + 1);
			posNP = strVal.find_last_of('.');
			if (posNP != strVal.length() - 1 && posNP != string::npos)
				mp.m_cNext = strVal[posNP + 1];
			strVal = strVal.substr(0, posNP);
			//cout << strVal << endl;
			vector<CModificationSiteInfo> vModTemp;
			_GetSQTMod(strVal, vModTemp);
			mp.m_vMod = vModTemp;
			mp.m_strSQ = strVal;
		}
		else if ('L' == strTemp[0])
		{
			char protmp[1000];
			string strProTemp;
			sscanf(strTemp.c_str(), "L %s", protmp);
			strProTemp = protmp;
			mp.m_vProteinAC.push_back(strProTemp);
		}
		//		else
		//		{
		//			CErrInfo info("SQTParser", "_ParseSpectra",
		//					"The file has some error!");
		//			throw runtime_error(info.Get());
		//		}
	}
}

void CSQTParser::Parse(const string & strFilePath, CParseredDataSetInfo & pbres)
{
	ifstream fin(strFilePath.c_str());
	if (!fin.good())
	{
		CErrInfo info("SQTParser", "Parse", "Cannot access the file: " + strFilePath);
		throw runtime_error(info.Get());
	}
	pbres.m_eEngineType = ST_SQT;
	size_t pos = strFilePath.find_last_of('/');
	m_strFilePath = strFilePath.substr(pos + 1);
	m_strFilePath = m_strFilePath.substr(0, m_strFilePath.length() - 4);
	//cout << m_strFilePath << endl;
	_ParserFront(fin, pbres);
	_ParseSpectra(fin, pbres);
	//cout << "SQT Run Over!" << endl;
	fin.close();
}
