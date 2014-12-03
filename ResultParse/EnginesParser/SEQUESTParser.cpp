#include "SEQUESTParser.h"

using namespace bio_analysis;
//using namespace proteomics_sdk;

extern ostringstream osspBuildLog;

CSEQUESTParser::CSEQUESTParser()
{
}

CSEQUESTParser::~CSEQUESTParser()
{
}

void CSEQUESTParser::_GetSequestMod(string & strVal, vector<CModificationSiteInfo> & vModTemp)
{
	//cout << strVal << endl;
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
		else
		{
			ModTemp.m_strModName.clear();
			ModTemp.m_strModName = strVal[t];
			//ModTemp.m_tPos = k - 1; todo
			ModTemp.m_tPos = k;
			vModTemp.push_back(ModTemp);
		}
	}
	strVal.clear();
	strVal = strTemp;
}

void CSEQUESTParser::_GetVarMod(string & strTemp, CConditionInfo & cond, size_t & pos)
{
	//cout << strTemp << endl;
	string strVal;
	int f = 0;
	for (size_t i = 0; i < strTemp.size(); i++)
	{
		if (strTemp[i] == ')')
		{
			strVal += strTemp[i];
			cond.m_vVarMod.push_back(strVal);
			strVal.clear();
			f = 0;
		}
		else if (strTemp[i] == '(')
			f = 1;
		if (f == 1)
		{
			if (strTemp[i] == 32)
				strVal += '_';
			else
				strVal += strTemp[i];
		}
		if (f == 0 && strTemp[i] != 32)
		{
			pos = i;
			break;
		}
	}
	//cout << "over" << endl;
}

void CSEQUESTParser::_GetFixMod(string & strTemp, CConditionInfo & cond, size_t & pos)
{
	//cout << "fsa" << endl;
	string strVal;
	strVal = strTemp.substr(pos + 1);
	//cout << strVal << endl;
	istringstream stream(strVal);
	while (stream >> strTemp)
	{
		//cout << "strTemp=" << strTemp << endl;
		size_t pos = strTemp.find_first_of('=');
		//cout << pos << endl;
		if (pos == string::npos)
		{
			//if (strTemp.substr(0, 7) == "Enzyme:")
			//cond.m_Enzyme = strTemp.substr(8);
			if (strTemp.substr(0, 7) == "Enzyme:")
			{
				string strChange;
				while (stream >> strTemp)
				{
					strChange += strTemp;
				}
				//cout << "change=" << strChange << endl;
				size_t post = strChange.find_first_of('(');
				cond.m_Enzyme = strChange.substr(0, post);
				sscanf(strChange.substr(post).c_str(), "(%d)", &cond.m_nMaxMiss);

				//cout << cond.m_Enzyme << endl;
				//cout << cond.m_nMaxMiss << endl;
				//system("pause");
			}
			//
			//			stream >> strTemp;
			//			sscanf(strTemp.c_str(), "(%d)", &cond.m_nMaxMiss);
		}
		else
		{
			strTemp[pos] = '_';
			cond.m_vFixMod.push_back(strTemp);
		}
	}
	//cout << "fsd" << endl;
}

void CSEQUESTParser::_ParserFront(ifstream &fin, CMatchSpectraInfo & SpectraInfo,
		CConditionInfo & cond)
{
	size_t pos;
	string strTemp, strVal;
	while (1)
	{
		GetLine(fin, strTemp);
		//cout << strTemp << endl;
		DeleteFrontBlank(strTemp);
		if (strTemp.substr(0, 13) == "(M+H)+ mass =")
		{
			strVal = strTemp.substr(13);
			sscanf(strVal.c_str(), " %lf ~ %lf (%d), fragment tol = %lf", &SpectraInfo.m_lfMH,
					&cond.m_lfPepTol, &SpectraInfo.m_nCharge, &cond.m_lfFragTol);
			size_t pos = strVal.find_first_of(',');
			strTemp = strVal.substr(pos + 1);
			pos = strTemp.find_first_of(',');
			strTemp = strTemp.substr(pos + 1);
			FetchLetter(strTemp);
			if (strTemp.substr(0, 3) == "AVG")
				cond.m_bPepMono = false;
			else
				cond.m_bPepMono = true;
			if (strTemp.substr(3, 3) == "AVG")
				cond.m_bFragMono = false;
			else
				cond.m_bFragMono = true;
		}
		if (1 == SpectraInfo.m_nSpectraNum && strTemp.find(".fasta") != string::npos)
		{
			size_t p = strTemp.find_first_of(',');
			strTemp = strTemp.substr(p + 1);
			p = strTemp.find_first_of(',');
			strTemp = strTemp.substr(p + 1);
			DeleteFrontBlank(strTemp);
			pos = strTemp.find(".fasta");
			strTemp = strTemp.substr(0, pos + 6);
			cond.m_vDatabase.push_back(strTemp);
		}
		if (strTemp.substr(0, 10) == "ion series")
			continue;
		if (strTemp.substr(0, 11) == "display top")
		{
			GetLine(fin, strTemp);
			//cout << strTemp << endl;
			DeleteFrontBlank(strTemp);
			if (1 == SpectraInfo.m_nSpectraNum)
			{
				//cout << "over1" << endl;
				_GetVarMod(strTemp, cond, pos);
				_GetFixMod(strTemp, cond, pos);
			}
			//cout << "over" << endl;
			return;
		}
	}
}

void CSEQUESTParser::_ParserPeptide(ifstream & fin, CMatchSpectraInfo & SpectraInfo,
		string & InputPath)
{
	string strTemp, strVal, strPro;
	int ntmp, tag = 0;
	char protmp[1000], SQtmp[1000];
	double lftmp, score, E_Value;
	CMatchPeptideInfo mp;
	bool eof;
	GetLine(fin, strTemp);
	DeleteFrontBlank(strTemp);
	size_t pos = strTemp.find("Id#");
	if (pos == string::npos)
		tag = 0;
	else
		tag = 1;
	bool bPro = false;
	while (1)
	{
		if (bPro == false)
		{
			eof = GetLine(fin, strTemp);
			if (eof == false)
				return;
		}
		if (strTemp.substr(strTemp.length() - 4, 4) == ".out")
		{
			InputPath = strTemp;
			return;
		}
		DeleteFrontBlank(strTemp);
		string strValTemp;
		istringstream streamTemp(strTemp);
		streamTemp >> strValTemp;
		streamTemp >> strValTemp;
		streamTemp >> strValTemp;
		size_t pos = strValTemp.find_first_of('/');
		if (pos == string::npos)
		{
			bPro = false;
			continue;
		}
		pos = strTemp.find_first_of('+');
		int tsize = 0;
		if (1 == tag)
		{
			if (pos == string::npos)
				sscanf(strTemp.c_str(), "%d. %d / %d %d %lf %lf %lf %lf %d / %d %s %s",
						&mp.m_tRank, &ntmp, &ntmp, &ntmp, &mp.m_lfCalc_MH, &E_Value, &score,
						&lftmp, &ntmp, &ntmp, protmp, SQtmp);
			else
				sscanf(strTemp.c_str(), "%d. %d / %d %d %lf %lf %lf %lf %d / %d %s +%d %s",
						&mp.m_tRank, &ntmp, &ntmp, &ntmp, &mp.m_lfCalc_MH, &E_Value, &score,
						&lftmp, &ntmp, &ntmp, protmp, &tsize, SQtmp);
		}
		else
		{
			if (pos == string::npos)
				sscanf(strTemp.c_str(), "%d. %d / %d %lf %lf %lf %lf %d / %d %s %s", &mp.m_tRank,
						&ntmp, &ntmp, &mp.m_lfCalc_MH, &E_Value, &score, &lftmp, &ntmp, &ntmp,
						protmp, SQtmp);
			else
				sscanf(strTemp.c_str(), "%d. %d / %d %lf %lf %lf %lf %d / %d %s +%d %s",
						&mp.m_tRank, &ntmp, &ntmp, &mp.m_lfCalc_MH, &E_Value, &score, &lftmp,
						&ntmp, &ntmp, protmp, &tsize, SQtmp);
		}
		strVal = protmp;
		mp.m_vProteinAC.push_back(strVal);
		mp.m_vlfScores.push_back(score);
		mp.m_vlfScores.push_back(E_Value);
		//if (mp.m_tRank == 2) todo 这里注意
		if (SpectraInfo.m_vPeptides.size() == 1)
			SpectraInfo.m_vPeptides[0].m_vlfScores[1] = E_Value;
		mp.m_lfDelta = SpectraInfo.m_lfMH - mp.m_lfCalc_MH;//这里先没有算。
		//mp.m_lfDelta = mp.m_lfCalc_MH - SpectraInfo.m_lfMH;//这里先没有算
		mp.m_lfPPM = CBioMethods::TransDeltaMH2PPM(mp.m_lfDelta, SpectraInfo.m_lfMH);
		//mp.m_lfCalc_M = mp.m_lfCalc_MH - protonH;
		mp.m_cNext = mp.m_cPrev = '-';
		strVal = SQtmp;
		size_t posNP = strVal.find_first_of('.');
		if (posNP != string::npos)
			mp.m_cPrev = strVal[posNP - 1];
		strVal = strVal.substr(posNP + 1);
		//posNP = strVal.find_first_of('.');
		posNP = strVal.find_last_of('.');
		if (posNP != strVal.length() - 1 && posNP != string::npos)
			mp.m_cNext = strVal[posNP + 1];
		strVal = strVal.substr(0, posNP);
		vector<CModificationSiteInfo> vModTemp;
		_GetSequestMod(strVal, vModTemp);
		mp.m_vMod = vModTemp;
		mp.m_strSQ = strVal;
		while (1)
		{
			GetLine(fin, strTemp);
			istringstream stream(strTemp);
			strVal.clear();
			strPro.clear();
			stream >> strVal >> strPro;
			if (strVal.find_first_of('.') != string::npos)
			{
				bPro = true;
				break;
			}
			mp.m_vProteinAC.push_back(strPro);
		}
		DeleteRedundancy(mp.m_vProteinAC);
		sort(mp.m_vMod.begin(), mp.m_vMod.end(), Mod_Sort);
		SpectraInfo.m_vPeptides.push_back(mp);
		mp.m_vProteinAC.clear();
		mp.m_vMod.clear();
		mp.m_vlfScores.clear();
		mp.clear();
	}
}

void CSEQUESTParser::_SEQUSETParser(ifstream & fin, CMatchSpectraInfo & SpectraInfo,
		CConditionInfo & cond, string & InputPath)
{
	//read from the second line
	DeleteFrontBlank(InputPath);
	//InputPath.replace(InputPath.length() - 3, 3, "out");
	InputPath.replace(InputPath.length() - 3, 3, "dta");
	SpectraInfo.m_strFileName = InputPath;
	cond.m_lfFragBase = 0;
	cond.m_lfPepBase = 0;
	//cond.m_nMaxMiss = 0;
	_ParserFront(fin, SpectraInfo, cond);
	//cout << "_ParserFront(fin, SpectraInfo, cond);" << endl;
	_ParserPeptide(fin, SpectraInfo, InputPath);
}

void CSEQUESTParser::_SEQUESTParserFile(const string & strFilePath, CParseredDataSetInfo & pbres)
{
	ifstream fin(strFilePath.c_str());
	if (!fin.good())
	{
		CErrInfo info("SEQUESTParser.cpp", "_SEQUESTParserFile", "Cannot open the file: "
				+ strFilePath);
		throw runtime_error(info.Get());
	}
	//string strBuf;
	string strTemp, InputPath;
	ThrowFirstOfRemain(fin, "OUTFILE_COUNT", strTemp);
	CConditionInfo cond;
	sscanf(strTemp.c_str(), "%d", &cond.m_nTotalSpec);
	int nTotal = cond.m_nTotalSpec;
	//	cout << nTotal << endl;
	//flog << nTotal << endl;
	//ostringstream oss;
	osspBuildLog << nTotal << endl;
#ifdef DEBUG
	cout << nTotal << endl;
#endif
	//pBuildLog(oss.str());
	CMatchSpectraInfo SpectraInfo;
	while (1)
	{
		GetLine(fin, strTemp);
		DeleteFrontBlank(strTemp);
		if (strTemp.substr(strTemp.length() - 4, 4) == ".out")
		{
			InputPath = strTemp;
			break;
		}
	}
	for (int i = 1; i <= nTotal; i++)
	{
		SpectraInfo.m_nSpectraNum = i;
		_SEQUSETParser(fin, SpectraInfo, cond, InputPath);
		//if (SpectraInfo.m_vPeptides.size() > 0)
		pbres.m_mapResultData[SpectraInfo.m_nCharge].push_back(SpectraInfo);
		SpectraInfo.m_vPeptides.clear();
	}
	pbres.m_vConditions.push_back(cond);
	//pbres.m_vCondition = cond;
	cond.m_vDatabase.clear();
	cond.m_vFixMod.clear();
	cond.m_vVarMod.clear();
	fin.close();
}

void CSEQUESTParser::_SEQUESTParserFolder(const string & strFilePath, CParseredDataSetInfo & pbres)
{
	//readdir_sequence filenames(strFilePath.c_str(), readdir_sequence::files);
	DIR * dir = opendir(strFilePath.c_str());
	if (NULL == dir)
	{
		CErrInfo info("SEQUESTParser", "_SEQUESTParserFolder", "Cannot open the file: "
				+ strFilePath);
		throw runtime_error(info.Get());
	}
	struct dirent * de;
	string strTemp, InputPath, FileName;
	CMatchSpectraInfo SpectraInfo;
	CConditionInfo cond;
	int SpectraNum = 0;
	for (; NULL != (de = readdir(dir));)
	{
		if (de->d_name[0] == '.')
			continue;

		FileName = strFilePath;
		FileName += de->d_name;
		if (FileName.substr(FileName.length() - 4, 4) != ".out")
			continue;
		//cout << FileName << endl;
		SpectraNum++;
		//cout << FileName << endl;
		ifstream fin(FileName.c_str());
		if (!fin.good())
		{
			CErrInfo info("SEQUESTOutParser.cpp", "_SEQUESTParserFolder", "Cannot open the file: "
					+ FileName);
			throw runtime_error(info.Get());
		}
		//string strBuf;
		while (1)
		{
			GetLine(fin, strTemp);
			//cout << strTemp << endl;
			DeleteFrontBlank(strTemp);
			if (strTemp.substr(strTemp.length() - 4, 4) == ".out")
			{
				InputPath = strTemp;
				break;
			}
		}
		SpectraInfo.m_nSpectraNum = SpectraNum;
		_SEQUSETParser(fin, SpectraInfo, cond, InputPath);
		//if (SpectraInfo.m_vPeptides.size() > 0)
		pbres.m_mapResultData[SpectraInfo.m_nCharge].push_back(SpectraInfo);
		SpectraInfo.m_vPeptides.clear();
		fin.close();
	}
	cond.m_nTotalSpec = SpectraNum;
	//cout << SpectraNum << endl;
	//flog << SpectraNum << endl;
	//ostringstream oss;
	osspBuildLog << SpectraNum << endl;
#ifdef DEBUG
	cout << SpectraNum << endl;
#endif
	//pBuildLog(oss.str());
	pbres.m_vConditions.push_back(cond);
	//pbres.m_vCondition = cond;
	closedir(dir);
}

void CSEQUESTParser::Parse(const string & strFilePath, CParseredDataSetInfo & pbres)
{
	//cout << "d" << endl;
	pbres.m_eEngineType = ST_SEQUEST;
	if ('/' == strFilePath[strFilePath.length() - 1])
		_SEQUESTParserFolder(strFilePath, pbres);
	else
		_SEQUESTParserFile(strFilePath, pbres);
}
