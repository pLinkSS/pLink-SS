#include "MascotParser.h"

//using namespace proteomics_sdk;
using namespace bio_analysis;

extern ostringstream osspBuildLog;

CMascotParser::CMascotParser()
{
}

CMascotParser::~CMascotParser()
{
}

void CMascotParser::_ParseParameters(ifstream & fin, string & strBuf, CConditionInfo & cond)
{
	//cout << "Start to Parse Parameters..." << endl;
	string strTemp;
	strTemp = "Content-Type: application/x-Mascot; name="; //notice //x
	strTemp += 34;
	strTemp += "parameters";
	strTemp += 34; //here must  separate to four step,or you get a wrong answer
	ReadUntilString(fin, strBuf, strTemp);
	cond.m_bPepMono = true;
	cond.m_lfPepBase = 0.0;
	cond.m_lfFragBase = 0.0;
	ThrowFirstOfRemain(fin, strBuf, "TOL", strTemp);
	sscanf(strTemp.c_str(), "%lf", &cond.m_lfPepTol);
	ThrowFirstOfRemain(fin, strBuf, "TOLU", cond.m_strPepTolType);
	ThrowFirstOfRemain(fin, strBuf, "ITOL", strTemp);
	sscanf(strTemp.c_str(), "%lf", &cond.m_lfFragTol);
	ThrowFirstOfRemain(fin, strBuf, "ITOLU", cond.m_strFragTolType);
	ThrowFirstOfRemain(fin, strBuf, "PFA", strTemp);
	sscanf(strTemp.c_str(), "%d", &cond.m_nMaxMiss);
	ThrowFirstOfRemain(fin, strBuf, "DB", strTemp);
	cond.m_vDatabase.push_back(strTemp);

	ThrowFirstOfRemain(fin, strBuf, "MASS", strTemp);
	if (strTemp == "Monoisotopic")
		cond.m_bFragMono = true;
	else
		cond.m_bFragMono = false;

	ThrowFirstOfRemain(fin, strBuf, "CLE", cond.m_Enzyme);
	ThrowFirstOfRemain(fin, strBuf, "FILE", cond.m_strInputPath);
	ThrowFirstOfRemain(fin, strBuf, "INSTRUMENT", cond.m_strInstrument);
}

void CMascotParser::_ParseMasses(ifstream & fin, string & strBuf, CConditionInfo & cond) //todo
{
	//cout << "Start to Parse Masses..." << endl;
	string strTemp, strVal;
	strTemp = "Content-Type: application/x-Mascot; name="; //notice //x
	strTemp += 34;
	strTemp += "masses";
	strTemp += 34;
	ReadUntilString(fin, strBuf, strTemp);
	while (1)
	{
		GetLine(fin, strBuf, strTemp);
		if (strTemp.substr(2) == m_strBoundary)
			break; //end of the part
		if (strTemp.substr(0, 5) == "delta")
		{
			size_t pos = strTemp.find_first_of('=');
			size_t nNum = atoi(strTemp.substr(5, pos - 5).c_str());
			pos = strTemp.find_first_of(',');
			strVal = strTemp.substr(pos + 1);
			FetchLetter(strVal);
			if (nNum > cond.m_vVarMod.size())
				cond.m_vVarMod.resize(nNum);
			cond.m_vVarMod[nNum - 1] = strVal;
		}
		if (strTemp.substr(0, 16) == "FixedModResidues")
			continue;
		if (strTemp.substr(0, 8) == "FixedMod")
		{
			size_t pos = strTemp.find_first_of('=');
			size_t nNum = atoi(strTemp.substr(8, pos - 8).c_str());
			pos = strTemp.find_first_of(',');
			strVal = strTemp.substr(pos + 1);
			FetchLetter(strVal);
			if (nNum > cond.m_vFixMod.size())
				cond.m_vFixMod.resize(nNum);
			cond.m_vFixMod[nNum - 1] = strVal;
		}
	}
	//cout << "Masses Over..." << endl;
}

void CMascotParser::_ParseHeader(ifstream & fin, string & strBuf, CConditionInfo & cond,
		int & nTotal)
{
	//cout << "Start to Parse Header..." << endl;
	string strTemp;
	strTemp = "Content-Type: application/x-Mascot; name="; //notice //x
	strTemp += 34;
	strTemp += "header";
	strTemp += 34;
	ReadUntilString(fin, strBuf, strTemp);
	ThrowFirstOfRemain(fin, strBuf, "queries", strTemp);
	sscanf(strTemp.c_str(), "%d", &cond.m_nTotalSpec);
	nTotal = cond.m_nTotalSpec;
	//cout << "Header Over..." << endl;
}

void CMascotParser::_ParseSummary(ifstream & fin, string & strBuf, int *naCharge)
{
	//cout << "Start to Parse Summary..." << endl;
	string strTemp;
	strTemp = "Content-Type: application/x-Mascot; name="; //notice //x
	strTemp += 34;
	strTemp += "summary";
	strTemp += 34;
	ReadUntilString(fin, strBuf, strTemp);
	while (1)
	{
		GetLine(fin, strBuf, strTemp);
		if (strTemp.substr(2) == m_strBoundary)
			break;
		if (strTemp.substr(0, 4) == "qexp")
		{
			size_t pos = strTemp.find_first_of('=');
			int nNum = atoi(strTemp.substr(4, pos - 4).c_str());
			pos = strTemp.find_first_of(',');
			strTemp = strTemp.substr(pos + 1);
			istringstream isCharge(strTemp);
			char cTemp;
			isCharge >> naCharge[nNum - 1] >> cTemp;
			if (cTemp != '+')
				naCharge[nNum - 1] = -naCharge[nNum - 1];
		}
	}
	//cout << "Summary Over..." << endl;
}

void CMascotParser::_ParsePeptides(ifstream & fin, string & strBuf,
		vector<CMatchSpectraInfo> & vSpectraInfo, CConditionInfo & cond, const int * anCharge)
{
	//cout << "Start to Parse Peptides..." << endl;
	string strTemp, strVal;
	strTemp = "Content-Type: application/x-Mascot; name="; //notice //x
	strTemp += 34;
	strTemp += "peptides";
	strTemp += 34;
	ReadUntilString(fin, strBuf, strTemp);
	CMatchSpectraInfo mr;
	CMatchPeptideInfo mp;
	int indexTemp = -1, indexQ = -1, indexP = -1;
	size_t pos1, pos2;
	while (1)
	{
		GetLine(fin, strBuf, strTemp);
		if (strTemp.substr(2) == m_strBoundary)
		{
			mr.m_nSpectraNum = indexQ;
			vSpectraInfo.push_back(mr); //the last one
			mr.m_vPeptides.clear();
			mr.clear();
			break;
		}
		pos2 = strTemp.find_first_of('=');
		char szchr = strTemp.substr(pos2 - 1, 1).c_str()[0];
		if (szchr == 't' || szchr == 'l')
			continue;
		pos1 = strTemp.find_first_of('_');
		indexTemp = indexQ;
		indexQ = atoi(strTemp.substr(1, pos1 - 1).c_str());
		indexP = atoi(strTemp.substr(pos1 + 2, pos2 - pos1 - 2).c_str());
		//here you must push_back the -1 to the vSpectraInfo,or you can't valid add title
		if (indexQ != indexTemp && indexTemp != -1)
		{
			mr.m_nSpectraNum = indexTemp;
			vSpectraInfo.push_back(mr);
			mr.m_vPeptides.clear();
			mr.clear();
		}
		strVal = strTemp.substr(pos2 + 1);
		if (strVal == "-1")
			continue;

		mp.m_tRank = indexP;
		vector<string> vstrLater;
		SplitStringByComma(strVal, vstrLater);
		sscanf(vstrLater[1].c_str(), "%lf", &mp.m_lfCalc_MH);
		sscanf(vstrLater[2].c_str(), "%lf", &mp.m_lfDelta);
		if (mr.m_vPeptides.size() == 0)
			mr.m_lfMH = mp.m_lfCalc_MH + mp.m_lfDelta;
		mp.m_lfPPM = CBioMethods::TransDeltaMH2PPM(mp.m_lfDelta, mr.m_lfMH);
		mp.m_strSQ = vstrLater[4];
		CModificationSiteInfo ModTemp;
		strTemp = vstrLater[6];
		size_t tV = vstrLater[6].size();
		for (size_t i = 0; i < tV; i++)
		{
			if (strTemp[i] == '0')
				continue;
			//todo
//			if (i == 0)
//				ModTemp.m_tPos = 0;
//			else if (i == tV - 1)
//				ModTemp.m_tPos = i - 2;
//			else
//				ModTemp.m_tPos = i - 1;
			//todo
			ModTemp.m_tPos = i;
			ModTemp.m_strModName = cond.m_vVarMod[strTemp[i] - '0' - 1];
			mp.m_vMod.push_back(ModTemp);
		}
		double lfTemp;
		sscanf(vstrLater[7].c_str(), "%lf", &lfTemp);
		mp.m_vlfScores.push_back(lfTemp);
		tV = vstrLater.size();
		if (tV >= 11)
		{
			vstrLater[10].erase(vstrLater[10].begin());
			vstrLater[10].erase(vstrLater[10].begin());
			GetBetweenColon(vstrLater[10]);
			mp.m_vProteinAC.push_back(vstrLater[10]);
			for (size_t i = 11; i < tV; i++)
			{
				GetBetweenColon(vstrLater[i]);
				mp.m_vProteinAC.push_back(vstrLater[i]);
			}
		}
		DeleteRedundancy(mp.m_vProteinAC);
		//		for (size_t t = 0; t < mp.m_vProteinAC.size(); t++)
		//		{
		//			GetBeforeDot(mp.m_vProteinAC[t]);
		//		}
		////////////////
		mp.m_cNext = mp.m_cPrev = '-';
		GetLine(fin, strBuf, strTemp);
		pos1 = strTemp.find_first_of('=');
		mp.m_cPrev = strTemp[pos1 + 1];
		mp.m_cNext = strTemp[pos1 + 3];
		sort(mp.m_vMod.begin(), mp.m_vMod.end(), Mod_Sort);
		mr.m_vPeptides.push_back(mp);
		mp.m_vMod.clear();
		mp.m_vProteinAC.clear();
		mp.m_vlfScores.clear();
		mp.clear();
	}
	//cout << "Peptides Over..." << endl;
}

void CMascotParser::_ParseQuery(ifstream & fin, string & strBuf,
		vector<CMatchSpectraInfo> & vSpectraInfo, const int & nTotal)
{
	//cout << "Start to Parse Query..." << endl;
	string strTemp, strVal, strTitle = "title";
	char szchr[PATH_MAX];
	strVal = "Content-Type: application/x-Mascot; name="; //notice //x
	for (int i = 1; i <= nTotal; i++)
	{
		sprintf(szchr, "%d", i);
		strTemp = strVal;
		strTemp += 34;
		strTemp += "query";
		strTemp += szchr;
		strTemp += 34;
		ReadUntilString(fin, strBuf, strTemp);
		ThrowFirstOfRemain(fin, strBuf, strTitle, strTemp);
		StringToHex(strTemp);
		vSpectraInfo[i - 1].m_strFileName = strTemp;
	}
	//cout << "Query Over..." << endl;
}

void CMascotParser::Parse(const string & strFilePath, CParseredDataSetInfo & pbres)
{
	pbres.m_eEngineType = ST_MASCOT;
	//cout << "Start to Parse Mascot..." << endl;
	if (!CheckFileValid(strFilePath))
	{
		CErrInfo info("MacotParser", "Parse", "Cannot access the file: " + strFilePath);
		throw runtime_error(info.Get());
	}
	//FILE * fp = fopen(strFilePath.c_str(), "r"); //open the file
	ifstream fin(strFilePath.c_str());
	if (!fin.good())
	{
		CErrInfo info("MascotParser.cpp", "Parse", "Cannot open the file: " + strFilePath);
		throw runtime_error(info.Get());
	}
	//istringstream iss;
	//iss.clear();
	string strBuf;
	ThrowFirstOfRemain(fin, strBuf, "Content-Type: multipart/mixed; boundary", m_strBoundary);
	CConditionInfo cond;
	int nTotal;
	vector<CMatchSpectraInfo> vSpectraInfo;
	_ParseParameters(fin, strBuf, cond);
	_ParseMasses(fin, strBuf, cond);
	_ParseHeader(fin, strBuf, cond, nTotal);
	//cout << "nTotal = " << nTotal << endl;
	//cout << nTotal << endl;
	//flog << nTotal << endl;
	osspBuildLog << nTotal << endl;
#ifdef DEBUG
	cout << nTotal << endl;
#endif
	//pBuildLog(oss.str());
	int *anCharge = new int[nTotal];
	_ParseSummary(fin, strBuf, anCharge);
	_ParsePeptides(fin, strBuf, vSpectraInfo, cond, anCharge);
	_ParseQuery(fin, strBuf, vSpectraInfo, nTotal);
	pbres.m_vConditions.push_back(cond);
	//pbres.m_vCondition = cond;
	cond.m_vDatabase.clear();
	cond.m_vFixMod.clear();
	cond.m_vVarMod.clear();
	for (int i = 0; i < nTotal; i++)
	{
		vSpectraInfo[i].m_nCharge = anCharge[i];
		//if (vSpectraInfo[i].m_vPeptides.size() > 0)
		pbres.m_mapResultData[anCharge[i]].push_back(vSpectraInfo[i]);
		vSpectraInfo[i].m_vPeptides.clear();
	}
	vSpectraInfo.clear();
	delete[] anCharge;
	fin.close();
}
