#include "CommonProcess.h"

extern ostringstream osspBuildLog;

using namespace bio_analysis;

bool IsProteinReverse(const string & protein, const CConf & conf)
{
	for (size_t j = 0; j < conf.m_vDecoyTags.size(); j++)
	{
		if (protein.find(conf.m_vDecoyTags[j]) != string::npos)
			return true;
	}
	return false;
}

bool IsPepReverseNeedAll(const vector<string> & vPro, const CConf & conf)
{
	for (size_t t = 0; t < vPro.size(); t++)
	{
		if (!IsProteinReverse(vPro[t], conf))
			return false;
	}
	return true;
}

bool IsPepReverseNeedOne(const vector<string> & vPro, const CConf & conf)
{
	for (size_t t = 0; t < vPro.size(); t++)
	{
		if (IsProteinReverse(vPro[t], conf))
			return true;
	}
	return false;
}

bool IsPepFromDecoy(const vector<string> & vPro, const CConf & conf)
{
	if (conf.m_bCrossLink == 1 || conf.m_bCrossLink == 2)
	{
		return IsPepReverseNeedAll(vPro, conf);
	}
	else
	{
		return IsPepReverseNeedOne(vPro, conf);
	}
}

int IsReverseCrossLink(const CMatchSpectraInfo & mr, const CConf & conf)
{
	//1表示正库，2表示反库，3表示crosslink里一个来自正，一个来自反
	int nCrossType = GetCrossType(mr.m_vPeptides[0].m_strSQ);
	if (nCrossType == 4)
	{
		vector<string> Pro1;
		vector<string> Pro2;
		for (size_t t = 0; t < mr.m_vPeptides[0].m_vProteinAC.size(); t++)
		{
			const string & strPro = mr.m_vPeptides[0].m_vProteinAC[t];
			size_t pos = strPro.find_first_of('-');
			Pro1.push_back(strPro.substr(0, pos));
			Pro2.push_back(strPro.substr(pos + 1));
		}
		bool IsProteinRev1 = IsPepFromDecoy(Pro1, conf);
		bool IsProteinRev2 = IsPepFromDecoy(Pro2, conf);

		if (IsProteinRev1 == false && IsProteinRev2 == false) //两个来自正库
			return 1;
		else if (IsProteinRev1 == true && IsProteinRev2 == true) //两个来自反库
			return 2;
		else
			return 3;
	}
	else
	{
		bool nIsPepFromDecoy = IsPepFromDecoy(mr.m_vPeptides[0].m_vProteinAC, conf);
		if (nIsPepFromDecoy == true)
			return 2;
		else
			return 1;
	}
}


int IsReverseTriLink(const CMatchSpectraInfo & mr, const CConf & conf) //add by fan 2013.6.28
{
	//1表示正库，2表示全反库，3表示2个反1个正，4表示1个反两个正

		vector<string> Pro1;
		vector<string> Pro2;
		vector<string> Pro3;
		for (size_t t = 0; t < mr.m_vPeptides[0].m_vProteinAC.size(); t++)
		{
			const string & strPro = mr.m_vPeptides[0].m_vProteinAC[t];
			size_t pos = strPro.find_first_of('-');
			size_t pos2 = strPro.find_first_of("-", pos+1);

			Pro1.push_back(strPro.substr(0, pos));
			Pro2.push_back(strPro.substr(pos + 1, pos2));
			Pro3.push_back(strPro.substr(pos2+1));
		}
		bool IsProteinRev1 = IsPepFromDecoy(Pro1, conf);
		bool IsProteinRev2 = IsPepFromDecoy(Pro2, conf);
		bool IsProteinRev3 = IsPepFromDecoy(Pro3, conf);

		if (IsProteinRev1 == false && IsProteinRev2 == false && IsProteinRev3 == false) //三个来自正库
			return 1;
		else if (IsProteinRev1 == true && IsProteinRev2 == true && IsProteinRev3 == true) //三个来自反库
			return 2;
		else if ( (IsProteinRev1 == true && IsProteinRev2 == false && IsProteinRev3 == true) ||
					 (IsProteinRev1 == false && IsProteinRev2 == true && IsProteinRev3 == true) ||
					 (IsProteinRev1 == true && IsProteinRev2 == true && IsProteinRev3 == false ) ) //两反一正
			return 3;
		else if ( (IsProteinRev1 == false && IsProteinRev2 == false && IsProteinRev3 == true) ||
			 (IsProteinRev1 == false && IsProteinRev2 == true && IsProteinRev3 == false) ||
			 (IsProteinRev1 == true && IsProteinRev2 == false && IsProteinRev3 == false ) ) //一反两正
			return 4;
		else
		{
			cout<<"Impossible is nothing."<<endl;
			return 1;
		}
}


int IsReverse(const CMatchSpectraInfo & mr, const CConf & conf)
{
	if (conf.m_bCrossLink == 1)
		return IsReverseCrossLink(mr, conf);
	else if (conf.m_bCrossLink == 2)//add by fan 2013.6.28
		return IsReverseTriLink(mr, conf);
	bool nIsPepFromDecoy = IsPepFromDecoy(mr.m_vPeptides[0].m_vProteinAC, conf);
	if (nIsPepFromDecoy == true)
		return 2;
	else
		return 1;
}

bool CMatchSpectraInfo_SortScore(const CMatchSpectraInfo & a, const CMatchSpectraInfo & b)
{
	return a.m_vPeptides[0].m_vlfScores[0] > b.m_vPeptides[0].m_vlfScores[0];
}

bool CMatchSpectraInfo_SortSCANPATH(const CMatchSpectraInfo & a, const CMatchSpectraInfo & b)
{
	return a.m_strFileName < b.m_strFileName;
}

bool SPECTRAINFO_SortScore(const SPECTRAINFO & a, const SPECTRAINFO & b)
{
	return a.second.m_lfScore > b.second.m_lfScore;
}

bool SPECTRAINFO_SortInPutPath(const SPECTRAINFO & a, const SPECTRAINFO & b)
{
	return a.second.m_strInPutPath < b.second.m_strInPutPath;
}

bool Mod_Sort(const CModificationSiteInfo & s1, const CModificationSiteInfo & s2)
{
	if (s1.m_tPos == s2.m_tPos)
		return s1.m_strModName < s2.m_strModName;

	return s1.m_tPos < s2.m_tPos;
}

bool Peptide_SortScore(const CMatchPeptideInfo & mp1, const CMatchPeptideInfo & mp2)
{
	return mp1.m_vlfScores[0] > mp2.m_vlfScores[0];
}

bool ModificationSiteEqual(const CModificationSiteInfo & s1, const CModificationSiteInfo & s2)
{
	return (s1.m_strModName == s2.m_strModName) && (s1.m_tPos == s2.m_tPos);
}

bool VecotrModEqual(const vector<CModificationSiteInfo> & m1,
		const vector<CModificationSiteInfo> & m2)
{
	if (m1.size() != m2.size())
		return false;

	for (size_t t = 0; t < m1.size(); t++)
	{
		if (!ModificationSiteEqual(m1[t], m2[t]))
			return false;

	}//todo Masoct AND pFind and SEQUEST 固定修饰不一样
	return true;
}

bool PeptideEqual(const PEPTIDEINFO & PepA, const PEPTIDEINFO & PepB)
{
	if (PepA.first == PepB.first && VecotrModEqual(PepA.second, PepB.second))
		return true;

	return false;
}

bool PEP_SORT_SQ(const PEPTIDEINFO & pairPep1, const PEPTIDEINFO & pairPep2)
{
	if (pairPep1.first == pairPep2.first)
	{
		return pairPep1.second.size() < pairPep2.second.size();
	}

	return pairPep1.first < pairPep2.first;
}

void P123Count(const PEPTIDEINFO & PairPepTemp, int MODNUM[], const CConf & conf)
{
	size_t cnt = 0;

	for (size_t k = 0; k < PairPepTemp.second.size(); k++)
	{
		char chr = PairPepTemp.first[PairPepTemp.second[k].m_tPos];

		if (conf.m_Filter.m_strModSites.find(chr) != string::npos)
			cnt++;
	}

	if (cnt <= conf.m_Filter.m_strModSites.size())
		MODNUM[cnt - 1]++; //todo 这里这么弄不是很合理
}

int GetDifferNumOfPep(const vector<PEPTIDEINFO> & vPairPep, int P123[], const CConf & conf)
{
	for (size_t t = 0; t < conf.m_Filter.m_strModSites.size(); t++)
		P123[t] = 0;

	if (vPairPep.size() == 0)
		return 0;

	int num = 1;
	P123Count(vPairPep[0], P123, conf);

	for (size_t t = 1; t < vPairPep.size(); t++)
	{
		if (!PeptideEqual(vPairPep[t], vPairPep[t - 1]))
		{
			num++;
			P123Count(vPairPep[t], P123, conf);
		}
	}
	return num;
}

bool IsSameMod(const pair<string, CModificationSiteInfo> & m1, const pair<string,
		CModificationSiteInfo> & m2)
{
	if (m1.second.m_strModName != m2.second.m_strModName)
		return false;

	if (m1.first[m1.second.m_tPos] != m2.first[m2.second.m_tPos])
		return false;

	for (size_t s1 = m1.second.m_tPos, s2 = m2.second.m_tPos; s1 < m1.first.length() && s2
			< m2.first.length(); s1++, s2++)
	{
		if (m1.first[s1] != m2.first[s2])
			return false;
	}

	for (int s1 = m1.second.m_tPos - 1, s2 = m2.second.m_tPos - 1; s1 > 0 && s2 > 0; s1--, s2--)
	{
		if (m1.first[s1] != m2.first[s2])
			return false;
	}
	return true;
}

void STYCount(const pair<string, CModificationSiteInfo> & PairPepTemp, int STY[],
		const CConf & conf)
{
	char chr = PairPepTemp.first[PairPepTemp.second.m_tPos];

	size_t pos = conf.m_Filter.m_strModSites.find(chr);

	if (pos != string::npos)
		STY[pos]++;
}

int GetNumOfMod(FILE * fout, const vector<PEPTIDEINFO> & vPairPep, int STY[], const CConf & conf)
{
	for (size_t t = 0; t < conf.m_Filter.m_strModSites.size(); t++)
		STY[t] = 0;

	if (vPairPep.size() == 0)
		return 0;

	vector<pair<string, CModificationSiteInfo> > vModTemp;

	for (size_t t = 0; t < vPairPep.size(); t++)
	{
		pair<string, CModificationSiteInfo> pairTemp;
		pairTemp.first = vPairPep[t].first;

		for (size_t k = 0; k < vPairPep[t].second.size(); k++)
		{
			char chr = vPairPep[t].first[vPairPep[t].second[k].m_tPos];

			if (conf.m_Filter.m_strModSites.find(chr) != string::npos)
			{
				pairTemp.second = vPairPep[t].second[k];
				vModTemp.push_back(pairTemp);
			}
		}
	}
	int * tag = new int[vModTemp.size()];

	for (size_t t = 0; t < vModTemp.size(); t++)
		tag[t] = 0;

	for (size_t t = 0; t < vModTemp.size(); t++)
	{
		if (tag[t] == 1)
			continue;
		for (size_t k = t + 1; k < vModTemp.size(); k++)
		{
			if (vModTemp[t].first == vModTemp[k].first && ModificationSiteEqual(vModTemp[t].second,
					vModTemp[k].second))
			{
				tag[k] = 1;
			}
			else if (IsSameMod(vModTemp[t], vModTemp[k]))
			{
				tag[k] = 1;
			}
		}
	}
	int num = 0;
	for (size_t t = 0; t < vModTemp.size(); t++)
	{
		if (tag[t] == 0)
		{
			num++;
			STYCount(vModTemp[t], STY, conf);
		}
	}
	delete[] tag;
	return num;
}

void GetIndex(map<int, vector<SPECTRAINFO> > & DataSetSpectra,
		const CMatchSpectraInfo & SpectraTemp, const SPECTRA_PATH & SpecPathTemp)
{
	SPECTRAINFO SpectraInfoTemp;
	SpectraInfoTemp.first = SpecPathTemp;
	if (SpectraTemp.m_vPeptides.size() > 0)
	{
		SpectraInfoTemp.second.m_lfScore = SpectraTemp.m_vPeptides[0].m_vlfScores[0];
		SpectraInfoTemp.second.m_strFirstPeptide = SpectraTemp.m_vPeptides[0].m_strSQ;
		//		SpectraInfoTemp.second.m_vProtein
		//				= SpectraTemp.m_vPeptides[0].m_vProteinAC;
	}
	else
	{
		SpectraInfoTemp.second.m_lfScore = 0;
		SpectraInfoTemp.second.m_strFirstPeptide = "";
		//		SpectraInfoTemp.second.m_vProtein.clear();
	}
	SpectraInfoTemp.second.m_strInPutPath = SpectraTemp.m_strFileName;
	DataSetSpectra[SpectraTemp.m_nCharge].push_back(SpectraInfoTemp);
}

double EValueToLager(const double & EValue)
{
	return -10 * log10(EValue);
}

double EValueToSmall(const double & EValue)
{
	return pow(10, -0.1 * EValue);
}

SearchEngineType GetEnginType(const string & strTemp)
{
	SearchEngineType EnginType;

	if (strTemp == "pFind")
		EnginType = ST_PFIND;

	else if (strTemp == "Mascot")
		EnginType = ST_MASCOT;

	else if (strTemp == "SEQUEST")
		EnginType = ST_SEQUEST;

	else if (strTemp == "SQT")
		EnginType = ST_SQT;

	else
		EnginType = ST_NEXT_ENGINE;
	return EnginType;
}

bool InChargeState(const int & ChargeState, const CConf & conf)
{
	for (size_t t = 0; t < conf.m_Filter.m_vChargeState.size(); t++)
		if (ChargeState == conf.m_Filter.m_vChargeState[t])
			return true;

	return false;
}

string GetPeptideWithMod(const CMatchPeptideInfo & mp)
{
	string strSQMod = mp.m_strSQ;
	for (size_t t = 0; t < mp.m_vMod.size(); t++)
	{
		strSQMod += mp.m_vMod[t].m_strModName;
		char chr[PATH_MAX];
		sprintf(chr, PRI_SIZE_T, mp.m_vMod[t].m_tPos);
		strSQMod += chr;
	}
	return strSQMod;
}

bool ISinVector(const int & a, vector<int> m_vINT)
{
	for (size_t t = 0; t < m_vINT.size(); t++)
		if (a == m_vINT[t])
			return true;

	return false;
}

int getEXPCode(const int & a, const int & b)
{
	int res = a;
	for (int i = 1; i < b; i++)
	{
		res *= a;
	}
	return res;
}

int GetReagentType(const string & strSQ)
{
	size_t pos1 = strSQ.find_first_of(':');
	string strReagentType = strSQ.substr(pos1 + 1);
	int nReagentType = atoi(strReagentType.c_str());
	return nReagentType;
}

int GetCrossType(const string & strSQ)
{
	int nCrossType = 0;
	size_t pos2 = strSQ.find_first_of('-');
	if (pos2 != string::npos)
	{
		nCrossType = 4;
		return nCrossType;
	}

	size_t pos3 = strSQ.find_first_of('(');
	if (pos3 == string::npos)
	{
		nCrossType = 1;
		return nCrossType;
	}
	else
	{
		string strTemp = strSQ.substr(pos3 + 1);
		size_t pos4 = strTemp.find_first_of('(');
		if (pos4 == string::npos)
		{
			nCrossType = 2;
			return nCrossType;
		}
		else
		{
			nCrossType = 3;
			return nCrossType;
		}
	}
}

int GetCrossLinkType(const int & charge, const int & nCrossType, const int & nReagentType)
{
	return getEXPCode(2, charge) * getEXPCode(3, nCrossType);
}

int DeCrossLinkType(const int & charge)
{
	int res = 0;
	int tmp = charge;
	while (tmp % 2 == 0)
	{
		res++;
		tmp /= 2;
	}
	return res;
}

