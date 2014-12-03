#include "pFindParser.h"

//using namespace proteomics_sdk;
using namespace bio_analysis;

//extern ofstream flog;
extern ostringstream osspBuildLog;

CpFindParser::CpFindParser()
{
}

CpFindParser::~CpFindParser()
{
}

void CpFindParser::_ParseConditonInfo(ifstream & fin, CParseredDataSetInfo & pbres, size_t & nTotal)
{
	//cout << "hfchen1"<< endl;
	string strTemp, szTemp;
	CConditionInfo cond;
	ThrowFirstOfRemain(fin, "InputPath", cond.m_strInputPath);
	//m_vDatabase
	ThrowFirstOfRemain(fin, "Database", strTemp);
	//cout << "over1" << endl;
	GetNSplitStringByComma(strTemp, cond.m_vDatabase);
	//cout << "over2" << endl;
	//m_Enzyme
	ThrowFirstOfRemain(fin, "Enzyme", cond.m_Enzyme);
	//m_vFixMod
	//cout << "ove3r" << endl;
	ThrowFirstOfRemain(fin, "Fixed_modifications", strTemp);
	GetNSplitStringByComma(strTemp, cond.m_vFixMod);
	//cond.m_vVarMod
	//cout << "over4" << endl;
	ThrowFirstOfRemain(fin, "Variable_modifications", strTemp);
	GetNSplitStringByComma(strTemp, cond.m_vVarMod);
	//cout << "over" << endl;
	//bool m_bPepMono;
	ThrowFirstOfRemain(fin, "Peptide_Mass", strTemp);
	if (strTemp == "Monoisotopic")
		cond.m_bPepMono = true;
	else
		cond.m_bPepMono = false;
	//cout << "fasdfad" << endl;
	//m_m_lfPepBase
	//ThrowFirstOfRemain(fin, strBuf, "Peptide_Mass_Tolerance_Base", strTemp);
	//sscanf(strTemp.c_str(), "%lf", &cond.m_lfPepBase);
	cond.m_lfPepBase = 0.0;
	ThrowFirstOfRemain(fin, "Peptide_Mass_Tolerance", strTemp);
	sscanf(strTemp.c_str(), "%lf", &cond.m_lfPepTol);
	ThrowFirstOfRemain(fin, "Peptide_Mass_Tolerance_Type", cond.m_strPepTolType);

	//////////////////////
	//bool m_bPepMono;
	ThrowFirstOfRemain(fin, "Fragment_Mass", strTemp);
	if (strTemp == "Monoisotopic")
		cond.m_bFragMono = true;
	else
		cond.m_bFragMono = false;
	//m_m_lfPepBase
	//ThrowFirstOfRemain(fin, strBuf, "Fragment_Mass_Tolerance_Base", strTemp);
	//sscanf(strTemp.c_str(), "%lf", &cond.m_lfFragBase);
	cond.m_lfFragBase = 0.0;
	ThrowFirstOfRemain(fin, "Fragment_Mass_Tolerance", strTemp);
	sscanf(strTemp.c_str(), "%lf", &cond.m_lfFragTol);
	ThrowFirstOfRemain(fin, "Fragment_Mass_Tolerance_Type", cond.m_strFragTolType);
	////////////////////////
	ThrowFirstOfRemain(fin, "Max_Missed_Cleavages", strTemp);
	sscanf(strTemp.c_str(), "%d", &cond.m_nMaxMiss);
	ThrowFirstOfRemain(fin, "Instrument_type", cond.m_strInstrument);
	ThrowFirstOfRemain(fin, "Spectra", strTemp);
	nTotal = atoi(strTemp.c_str());
	cond.m_nTotalSpec = nTotal;
	pbres.m_vConditions.push_back(cond);
	//pbres.m_vCondition = cond;
	//cout << "over" << endl;
}

void CpFindParser::_ParsePeptideInfor(ifstream & fin, CMatchSpectraInfo & mr,
		CMatchPeptideInfo & mp, size_t tRank)
{
	string strTemp;
	double lfTemp;
	const size_t MAX_LEN = 100;
	char szTemp[MAX_LEN] =
	{ '\0' };

	sprintf(szTemp, "NO"PRI_SIZE_T"_Score", tRank);
	ThrowFirstOfRemain(fin, szTemp, strTemp);
	sscanf(strTemp.c_str(), "%lf", &lfTemp);
	mp.m_vlfScores.push_back(lfTemp);

	sprintf(szTemp, "NO"PRI_SIZE_T"_EValue", tRank);
	ThrowFirstOfRemain(fin, szTemp, strTemp);
	sscanf(strTemp.c_str(), "%lf", &lfTemp);
	//cout << lfTemp << " ";
	lfTemp = -10 * log10(lfTemp);
	//cout << lfTemp << endl;
	mp.m_vlfScores.push_back(lfTemp);

	///将EValue存到第一个里。
	lfTemp = mp.m_vlfScores[0];
	mp.m_vlfScores[0] = mp.m_vlfScores[1];
	mp.m_vlfScores[1] = lfTemp;

	GetLine(fin, strTemp);
	//cout << "tmp=" << strTemp << endl;
	string strResult;
	GetBeforeEqual(strTemp, strResult);
	//cout << "temp=" << strTemp << endl;
	//cout << "rest=" << strResult << endl;
	sscanf(strTemp.c_str(), "%lf", &mp.m_lfCalc_MH);
	bool old = false;
	if (strResult.find("MH") == string::npos)
		old = true;
	//cout << "d=" << mp.m_lfCalc_MH << endl;

	//mp.m_eEngineType = ST_PFIND;
	//mp.m_lfCalc_MH -= 1.00727647012;
	mp.m_lfDelta = mr.m_lfMH - mp.m_lfCalc_MH;
	mp.m_lfPPM = CBioMethods::TransDeltaMH2PPM(mp.m_lfDelta, mr.m_lfMH);
	//cout << mp.m_lfPPM << endl;
	//mp.m_lfCalc_MH -= protonH;
	//mp.m_lfCalc_M = mp.m_lfCalc_MH - protonH;
	sprintf(szTemp, "NO"PRI_SIZE_T"_SQ", tRank);
	ThrowFirstOfRemain(fin, szTemp, mp.m_strSQ);

	sprintf(szTemp, "NO"PRI_SIZE_T"_Proteins", tRank);
	ThrowFirstOfRemain(fin, szTemp, strTemp);
	//cout << strTemp << endl;
	GetNSplitStringByComma(strTemp, mp.m_vProteinAC);
//	printf("before\n");
//	for(size_t i = 0;i < mp.m_vProteinAC.size();++i){
//		printf("\t%s\n",mp.m_vProteinAC[i].c_str());
//	}
	DeleteRedundancy(mp.m_vProteinAC);
//	printf("after\n");
//	for(size_t i = 0;i < mp.m_vProteinAC.size();++i){
//		printf("\t%s\n",mp.m_vProteinAC[i].c_str());
//	}
//	system("pause");
	//	for (size_t t = 0; t < mp.m_vProteinAC.size(); t++)
	//	{
	//		GetBeforeDot(mp.m_vProteinAC[t]);
	//	}
	sprintf(szTemp, "NO"PRI_SIZE_T"_Modify_Pos", tRank);
	ThrowFirstOfRemain(fin, szTemp, strTemp);
	int nModifyNum = 0;
	sscanf(strTemp.c_str(), "%d", &nModifyNum);

	string strModPos;
//	string strTemp2 = strTemp.replace(strTemp.find("("), strTemp.find("-") - strTemp.find("(") + 1, "");
	size_t npos = strTemp.find_first_of(',');
	strTemp += ',';
	CModificationSiteInfo info;
	for (int j = 1; j <= nModifyNum; ++j)
	{
		size_t temp = npos + 1;
		npos = strTemp.find_first_of(',', temp);
		if (string::npos == npos)
			break;
		strModPos = strTemp.substr(temp, npos - temp);
		if (old == true)
			info.m_tPos = (atoi(strModPos.c_str()) + 1);
		else
			info.m_tPos = (atoi(strModPos.c_str()));
		mp.m_vMod.push_back(info);
	}

	sprintf(szTemp, "NO"PRI_SIZE_T"_Modify_Name", tRank);
	ThrowFirstOfRemain(fin, szTemp, strTemp);
	strTemp += ",";

	string strModName;
	npos = strTemp.find_first_of(',');
	for (int j = 1; j <= nModifyNum; ++j)
	{
		size_t temp = npos + 1;
		npos = strTemp.find_first_of(',', temp);
		if (string::npos == npos)
			break;
		strModName = strTemp.substr(temp, npos - temp);
		mp.m_vMod[j - 1].m_strModName = strModName;
	}
	mp.m_tRank = tRank;
	mp.m_cNext = mp.m_cPrev = '-';
}

void CpFindParser::_ParseMatchSpectraInfo(ifstream & fin, CMatchSpectraInfo & mr)
{
	string strTemp;
	CMatchPeptideInfo mp;
	RemainOfString(fin, "[Spectrum", strTemp);
	strTemp.erase(strTemp.end() - 1);
	sscanf(strTemp.c_str(), "%d", &mr.m_nSpectraNum);
	//cout << "spectranum=" << mr.m_nSpectraNum << endl;

	ThrowFirstOfRemain(fin, "Input", mr.m_strFileName);

	ThrowFirstOfRemain(fin, "Charge", strTemp);
	sscanf(strTemp.c_str(), "%d", &mr.m_nCharge);

	ThrowFirstOfRemain(fin, "MH", strTemp);
	sscanf(strTemp.c_str(), "%lf", &mr.m_lfMH);

	ThrowFirstOfRemain(fin, "ValidCandidate", strTemp);
	int nValid = 0;
	sscanf(strTemp.c_str(), "%d", &nValid);
	//cout << nValid << endl;'
	//	if (nValid > 0)
	//		nValid = 1;//只读第一个肽段
	for (int i = 1; i <= nValid; ++i)
	{
		_ParsePeptideInfor(fin, mr, mp, i);
		sort(mp.m_vMod.begin(), mp.m_vMod.end(), Mod_Sort);
		mr.m_vPeptides.push_back(mp);
		mp.m_vMod.clear();
		mp.m_vProteinAC.clear();
		mp.m_vlfScores.clear();
		mp.clear();
	}
	//cout << "return " << endl;
	//return true;
}

void CpFindParser::_ParseTXTFile(const string & strFilePath, CParseredDataSetInfo & pbres)
{
	ifstream fin(strFilePath.c_str());
	if (!fin.good())
	{
		CErrInfo info("CpFindParser.cpp", "ParseTXTFile", "Cannot open the file: " + strFilePath);
		throw runtime_error(info.Get());
	}
	//string strBuf;
	size_t nTotal;
	_ParseConditonInfo(fin, pbres, nTotal);//Parse condition of the file
	//	cout << nTotal << endl;
	//flog << nTotal << endl;
	//ostringstream oss;
	osspBuildLog << nTotal << endl;
#ifdef DEBUG
	cout << nTotal << endl;
#endif
	//pBuildLog(oss.str());
	CMatchSpectraInfo mr;
	for (size_t curr = 0; curr < nTotal; ++curr)
	{
#ifdef DEBUG
		cout << "cur = " << curr + 1 << " " << nTotal << endl;
#endif
		//cout << "cur = " << curr + 1 << " " << nTotal << endl;
		//if (_ParseMachSpectraInfo(fp, mr))没有候选太也都加进来
		_ParseMatchSpectraInfo(fin, mr);
		//cout << mr.m_strFileName << endl;
		//cout << "fsdjkl" << endl;
		pbres.m_mapResultData[mr.m_nCharge].push_back(mr);
		//cout << "yeah" << endl;
		mr.m_vPeptides.clear();
		//cout << "kfsl" << endl;
		//cout << "over" << endl;
		//cout << nTotal << " " << curr + 1 << endl;
	}
	//cout << "over" << endl;
	fin.close();
}

void CpFindParser::Parse(const string & strFilePath, CParseredDataSetInfo & pbres)
{
#ifdef DEBUG
	cout << "Start to Parse pFind..." << endl;
#endif
	if (!CheckFileValid(strFilePath))
	{
		CErrInfo info("pFindParser", "Parse", "Cannot access the file: " + strFilePath);
		//cout << "pbres.m_eEngineType" << pbres.m_eEngineType << endl;
		throw runtime_error(info.Get());
	}
	pbres.m_eEngineType = ST_PFIND;

	_ParseTXTFile(strFilePath, pbres);
	//cout << "pFind Parser Over" << endl;
}

