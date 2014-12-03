#include "FilterUtility.h"

//#define SUMO

extern ostringstream osspBuildLog;

CFilterUtility::CFilterUtility()
{
	// TODO Auto-generated constructor stub

}

CFilterUtility::~CFilterUtility()
{
	// TODO Auto-generated destructor stub
}

double CFilterUtility::_Get_FDR(const int & F, const int & T, const int & U, const CConf & conf)
{
	double lfFDR;
	if (conf.m_bCrossLink == 1)
	{
		//U等于0时，计算结果相等
		if (T == 0)
			return 1.0;

#ifdef SUMO
		lfFDR = F / (double) T;
		if (lfFDR > 1.0)
			return 1.0;
		else
			return lfFDR;
#endif

		if (U - 2 * F > 0)
		{
			lfFDR = (U - F) / (double) T;
		}
		else
		{
			lfFDR = F / (double) T;
		}
		if (lfFDR > 1.0)
			return 1.0;
		else
			return lfFDR;
	}
	if (conf.m_Filter.m_FDRFormula == 1)
	{
		if (T == 0)
			return 1.0;
		lfFDR = (double) F / T;
		if (lfFDR > 1.0)
			return 1.0;
		else
			return lfFDR;
	}
	else
	{
		lfFDR = (double) 2 * F / (T + F);
		if (lfFDR > 1.0)
			return 1.0;
		else
			return lfFDR;
	}
}

double CFilterUtility::_Get_FDR(const int & T,  const int & F3, const int & F2, const int & F1, const CConf & conf)
{
	double lfFDR;

	//U等于0时，计算结果相等
	if (T == 0)
		return 1.0;

	if (F1-F2+F3 > 0)
	{
		lfFDR = (F1-F2+F3) / (double) T;
	}
	else
	{
		lfFDR = F3 / (double) T;
	}
	if (lfFDR > 1.0)
		return 1.0;
	else
		return lfFDR;
}

bool CFilterUtility::_IsExist(const CMatchSpectraInfo mr, set<string> & setStr)
{
	size_t tsize = setStr.size();
	string strSQMod = GetPeptideWithMod(mr.m_vPeptides[0]);
	setStr.insert(strSQMod);
	size_t tx = setStr.size();
	if (tx == tsize + 1)
		return false;
	return true;
}

void CFilterUtility::_FixedDeltCn(const vector<SPECTRAINFO> & vNewSpectraInfo,
		vector<SPECTRAINFO> & vpectraInfoTemp, double DeltCn)
{
	CMatchSpectraInfo SpectraTemp;
	for (size_t t = 0; t < vNewSpectraInfo.size(); t++)
	{
		ReadSingeSpectraFisrtPep(vNewSpectraInfo[t].first, SpectraTemp);
		if (SpectraTemp.m_vPeptides[0].m_vlfScores[1] >= DeltCn)
			vpectraInfoTemp.push_back(vNewSpectraInfo[t]);
	}
}

void CFilterUtility::_OUTPUT_FDR_Curve(vector<double> vlfFDR, const CConf & conf, const pair<int,
		int> & pairSimpleCharge)
{
	for (int t = vlfFDR.size() - 2; t >= 0; t--)
	{
		if (vlfFDR[t] > vlfFDR[t + 1])
			vlfFDR[t] = vlfFDR[t + 1];
	}
	char OutPathTemp[PATH_MAX] =
	{ 0 };
	sprintf(OutPathTemp, "Sample%d_Charge%d.txt", pairSimpleCharge.first, pairSimpleCharge.second);
	string ChargeOutPath = conf.m_outPutForder_Charge + OutPathTemp;

	ofstream fout(ChargeOutPath.c_str());
	for (size_t t = 0; t < vlfFDR.size(); t++)
	{
		fout << t + 1 << " " << vlfFDR[t] << endl;
	}
	fout.close();
}

size_t CFilterUtility::_Filter_FDR(vector<SPECTRAINFO> & vNewSpectraInfo, const CConf & conf,
		const pair<int, int> & pairSimpleCharge)
{
	double lfTemp;
	set<string> setStr;
	vector<double> vlfFDR;

	CFiltration m_Filter = conf.m_Filter;
	CMatchSpectraInfo SpectraTemp;

	int F = 0, T = 0, U = 0;
	int F2 = 0, F1 = 0;

	set<int> ndel;
	ndel.clear();
	for (size_t t = 0; t < vNewSpectraInfo.size(); t++)
	{
		ReadSingeSpectraFisrtPep(vNewSpectraInfo[t].first, SpectraTemp);
		if (m_Filter.m_bRedundant == 2 && _IsExist(SpectraTemp, setStr))
		{
			ndel.insert(t);
			continue;
		}
		int IsRev = IsReverse(SpectraTemp, conf);
		//osspBuildLog << SpectraTemp.m_strFileName << endl;
		//osspBuildLog << T << " " << " " << F << " " << U << endl;
		if (IsRev == 1)
			T++;
		else if (IsRev == 2)
			F++;
		else if (IsRev == 3) //表示一正一反 fan: 或一正两反
		{
			F2++;
			U++;
		}
		else if (IsRev == 4) // 表示两正一反
			F1++;

		if (! (conf.m_bCrossLink == 2))
			lfTemp = _Get_FDR(F, T, U, conf);
		else
			lfTemp = _Get_FDR(T, F, F2, F1, conf);

		vlfFDR.push_back(lfTemp);
	}

	if (ndel.size() > 0)
	{
		vector<SPECTRAINFO> vNewSpectraTemp;
		for (size_t t = 0; t < vNewSpectraInfo.size(); t++)
		{
			if (ndel.find(t) != ndel.end())
				continue;
			vNewSpectraTemp.push_back(vNewSpectraInfo[t]);
		}
		vNewSpectraInfo.clear();
		vNewSpectraInfo = vNewSpectraTemp;
		vNewSpectraTemp.clear();
	}
	ndel.clear();

	if (vlfFDR.size() == 0)
		return 0;
	if (vlfFDR.size() > 100)
	{
		_OUTPUT_FDR_Curve(vlfFDR, conf, pairSimpleCharge);
	}
	if (conf.m_bCrossLink == 1 || conf.m_bCrossLink == 2)
	{
		for (size_t t = 0; t < vlfFDR.size(); t++)
		{
			if (vlfFDR[t] > m_Filter.m_lfFDR)
				return t;
		}
		return vlfFDR.size();
	}
	size_t t = vlfFDR.size() - 1;
	for (; t >= 0; t--) //注意size_t 用--(减减操作符)时注意t==0时要写break，因为t不会小于0

	{
		if (vlfFDR[t] <= m_Filter.m_lfFDR)
		{
			return t + 1;
		}
		if (t == 0)
		{
			return 0;
		}
	}
}

void CFilterUtility::_FilterMLD(const vector<SPECTRAINFO> & vSpectraInfoTemp,
		vector<SPECTRAINFO> & vNewSpectraInfo, const CConf & conf)
{
	size_t lSQ;
	double massTemp;
	CMatchSpectraInfo SpectraTemp;
	for (size_t t = 0; t < vSpectraInfoTemp.size(); t++)
	{
		ReadSingeSpectraFisrtPep(vSpectraInfoTemp[t].first, SpectraTemp);

		if (SpectraTemp.m_vPeptides.size() == 0)
			continue;

		massTemp = SpectraTemp.m_vPeptides[0].m_lfCalc_MH - protonH;
		lSQ = SpectraTemp.m_vPeptides[0].m_strSQ.size();

		if (massTemp < conf.m_Filter.m_lfPepMassLowerBound || massTemp
				> conf.m_Filter.m_lfPepMassUpperBound || lSQ < conf.m_Filter.m_tLengthLowerBound
				|| lSQ > conf.m_Filter.m_tLengthUpperBound)
			continue;

		vector<double> DeltaTemp;
		if ("Da" == conf.m_Filter.m_strPepTolType)
		{
			for (size_t i = 0; i < conf.m_Filter.m_lfPepTolBase.size(); i++)
			{
				DeltaTemp.push_back(SpectraTemp.m_lfMH - (SpectraTemp.m_vPeptides[0].m_lfCalc_MH
						+ conf.m_Filter.m_lfPepTolBase[i]));
			}
		}

		else
		{
			for (size_t i = 0; i < conf.m_Filter.m_lfPepTolBase.size(); i++)
			{
				DeltaTemp.push_back(CBioMethods::TransDeltaMH2PPM(
						SpectraTemp.m_lfMH - (SpectraTemp.m_vPeptides[0].m_lfCalc_MH
								+ conf.m_Filter.m_lfPepTolBase[i]), SpectraTemp.m_lfMH));
			}
		}

		for (size_t i = 0; i < DeltaTemp.size(); i++)
		{
			if (DeltaTemp[i] >= conf.m_Filter.m_lfPepTolLowerBound[i] && DeltaTemp[i]
					<= conf.m_Filter.m_lfPepTolUpperBound[i])
			{
				vNewSpectraInfo.push_back(vSpectraInfoTemp[t]);
				break;
			}
		}
	}
}

void CFilterUtility::_FilterMLD_Separate(const vector<SPECTRAINFO> & vSpectraInfoTemp, vector<
		vector<SPECTRAINFO> > & vNewSpectraInfo, const CConf & conf)
{
	size_t lSQ;
	double massTemp;
	CMatchSpectraInfo SpectraTemp;

	vNewSpectraInfo.resize(conf.m_Filter.m_lfPepTolBase.size());

	for (size_t t = 0; t < vSpectraInfoTemp.size(); t++)
	{
		ReadSingeSpectraFisrtPep(vSpectraInfoTemp[t].first, SpectraTemp);
		if (SpectraTemp.m_vPeptides.size() == 0)
			continue;

		massTemp = SpectraTemp.m_vPeptides[0].m_lfCalc_MH - protonH;
		lSQ = SpectraTemp.m_vPeptides[0].m_strSQ.size();

		if (massTemp < conf.m_Filter.m_lfPepMassLowerBound || massTemp
				> conf.m_Filter.m_lfPepMassUpperBound || lSQ < conf.m_Filter.m_tLengthLowerBound
				|| lSQ > conf.m_Filter.m_tLengthUpperBound)
			continue;

		vector<double> DeltaTemp;
		if ("Da" == conf.m_Filter.m_strPepTolType)
		{
			for (size_t i = 0; i < conf.m_Filter.m_lfPepTolBase.size(); i++)
			{
				DeltaTemp.push_back(SpectraTemp.m_lfMH - (SpectraTemp.m_vPeptides[0].m_lfCalc_MH
						+ conf.m_Filter.m_lfPepTolBase[i]));
			}
		}
		else
		{
			for (size_t i = 0; i < conf.m_Filter.m_lfPepTolBase.size(); i++)
			{
				DeltaTemp.push_back(CBioMethods::TransDeltaMH2PPM(
						SpectraTemp.m_lfMH - (SpectraTemp.m_vPeptides[0].m_lfCalc_MH
								+ conf.m_Filter.m_lfPepTolBase[i]), SpectraTemp.m_lfMH));
			}
		}
		for (size_t i = 0; i < DeltaTemp.size(); i++)
		{
			if (DeltaTemp[i] >= conf.m_Filter.m_lfPepTolLowerBound[i] && DeltaTemp[i]
					<= conf.m_Filter.m_lfPepTolUpperBound[i])
			{
				vNewSpectraInfo[i].push_back(vSpectraInfoTemp[t]);
				break;
			}
		}
	}
}

void CFilterUtility::_OutTo(vector<SPECTRAINFO> & vSpectraOut, const vector<SPECTRAINFO> & vNEW,
		const size_t resSize, const CConf & conf, const pair<int, int> & pairSimpleCharge)
{
	CMatchSpectraInfo SpectraTemp;
#ifdef DEBUG
	///////////////////////////////////////////////////////////
	//This is test.
	//static int nx = 1;
	//ostringstream oss;
	//oss << conf.m_outPutForder_Index << "test" << nx++ << ".txt";
	char OutPathTemp[100];
	sprintf(OutPathTemp, "TEST_Sample%d_Charge%d.txt", pairSimpleCharge.first,
			pairSimpleCharge.second);
	string ChargeOutPath = conf.m_outPutForder_Charge + OutPathTemp;
	FILE * fout = fopen(ChargeOutPath.c_str(), "w");
	for (size_t t = 0; t < vNEW.size(); t++)
	{
		ReadSingeSpectraFisrtPep(vNEW[t].first, SpectraTemp);
		fprintf(fout, ""PRI_SIZE_T" %s   %10.10e   %lf %lf ", t + 1, SpectraTemp.m_strFileName.c_str(),
				EValueToSmall(SpectraTemp.m_vPeptides[0].m_vlfScores[0]),
				SpectraTemp.m_vPeptides[0].m_lfDelta, SpectraTemp.m_vPeptides[0].m_lfPPM);
		for (size_t i = 0; i < SpectraTemp.m_vPeptides[0].m_vProteinAC.size(); i++)
		{
			fprintf(fout, "%s; ", SpectraTemp.m_vPeptides[0].m_vProteinAC[i].c_str());
		}
		fprintf(fout, "\n");
	}
	fclose(fout);
#endif

	vSpectraOut.clear();

	for (size_t t = 0; t < resSize; t++)
	{
		ReadSingeSpectraFisrtPep(vNEW[t].first, SpectraTemp);
		int nIsRev = IsReverse(SpectraTemp, conf);
#ifdef DEBUG
		for (size_t j = 0; j < SpectraTemp.m_vPeptides[0].m_vProteinAC.size(); j++)
		{
			osspBuildLog << SpectraTemp.m_vPeptides[0].m_vProteinAC[j] << endl;
		}
		osspBuildLog << "nIsRev = " << nIsRev << endl;
#endif
		if (conf.m_bReserveDecoy != 1 && (nIsRev == 2 || nIsRev == 3 || nIsRev == 4))
			continue;
		vSpectraOut.push_back(vNEW[t]);
	}
}

void CFilterUtility::_FDR(vector<SPECTRAINFO> & vNewSpectraInfo, vector<SPECTRAINFO> & vSpectraOut,
		const CConf & conf, const SearchEngineType & m_EngineType,
		const pair<int, int> & pairSimpleCharge)
{
	size_t resSize = 0;

	if (m_EngineType == ST_PFIND || m_EngineType == ST_MASCOT)
	{
		sort(vNewSpectraInfo.begin(), vNewSpectraInfo.end(), SPECTRAINFO_SortScore);
#ifdef DEBUG
		cout << "vNewSpectraInfo = " << vNewSpectraInfo.size() << endl;
#endif
		resSize = _Filter_FDR(vNewSpectraInfo, conf, pairSimpleCharge);
#ifdef DEBUG
		cout << "resSize= " << resSize << endl;
#endif

		_OutTo(vSpectraOut, vNewSpectraInfo, resSize, conf, pairSimpleCharge);
	}
	else if (m_EngineType == ST_SEQUEST || m_EngineType == ST_SQT)
	{
		vector<SPECTRAINFO> vSpectraInfoTemp;
		if (conf.m_Filter.m_bFixedDeltCn == true)
		{
			_FixedDeltCn(vNewSpectraInfo, vSpectraInfoTemp, conf.m_Filter.m_lfFixedDeltCn);
			sort(vSpectraInfoTemp.begin(), vSpectraInfoTemp.end(), SPECTRAINFO_SortScore);
			resSize = _Filter_FDR(vSpectraInfoTemp, conf, pairSimpleCharge);

			_OutTo(vSpectraOut, vSpectraInfoTemp, resSize, conf, pairSimpleCharge);
		}
		else
		{
			size_t MaxSize = 0;
			double MaxDeltCn = 0.0;
			for (double DeltCn = 0; DeltCn <= 0.5; DeltCn += 0.01)
			{
				vSpectraInfoTemp.clear();
				_FixedDeltCn(vNewSpectraInfo, vSpectraInfoTemp, conf.m_Filter.m_lfFixedDeltCn);
				sort(vSpectraInfoTemp.begin(), vSpectraInfoTemp.end(), SPECTRAINFO_SortScore);
				resSize = _Filter_FDR(vSpectraInfoTemp, conf, pairSimpleCharge);
				if (resSize > MaxSize)
				{
					MaxSize = resSize;
					MaxDeltCn = DeltCn;
				}
			}
			vSpectraInfoTemp.clear();
			_FixedDeltCn(vNewSpectraInfo, vSpectraInfoTemp, conf.m_Filter.m_lfFixedDeltCn);
			sort(vSpectraInfoTemp.begin(), vSpectraInfoTemp.end(), SPECTRAINFO_SortScore);
			resSize = _Filter_FDR(vSpectraInfoTemp, conf, pairSimpleCharge);
			_OutTo(vSpectraOut, vSpectraInfoTemp, resSize, conf, pairSimpleCharge);
		}
	}
}

void CFilterUtility::_SepatrateOrNot(const vector<SPECTRAINFO> & vSpectraInfoIn,
		vector<SPECTRAINFO> & vSpectraOut, const CConf & conf,
		const SearchEngineType & m_EngineType, const pair<int, int> & pairSimpleCharge)
{
	if (conf.m_Filter.m_Sepatrate == false)
	{
		vector<SPECTRAINFO> vNewSpectraInfo;
		_FilterMLD(vSpectraInfoIn, vNewSpectraInfo, conf);

		_FDR(vNewSpectraInfo, vSpectraOut, conf, m_EngineType, pairSimpleCharge);
	}
	else
	{
		vector<vector<SPECTRAINFO> > vNewSpectraInfo;
		_FilterMLD_Separate(vSpectraInfoIn, vNewSpectraInfo, conf);
		vector<SPECTRAINFO> vSpectraOutTemp;
		for (size_t t = 0; t < vNewSpectraInfo.size(); t++)
		{
			_FDR(vNewSpectraInfo[t], vSpectraOutTemp, conf, m_EngineType, pairSimpleCharge);

			for (size_t k = 0; k < vSpectraOutTemp.size(); k++)
			{
				vSpectraOut.push_back(vSpectraOutTemp[k]);
			}
		}
	}
}

void CFilterUtility::_Score(const vector<SPECTRAINFO> & vSpectraInfoIn,
		vector<SPECTRAINFO> & vSpectraOut, const CConf & conf,
		const SearchEngineType & m_EngineType)
{
	vector<SPECTRAINFO> vNewSpectraInfo;
	_FilterMLD(vSpectraInfoIn, vNewSpectraInfo, conf);
	CMatchSpectraInfo SpectraTemp;

	for (size_t t = 0; t < vNewSpectraInfo.size(); t++)
	{
		ReadSingeSpectraFisrtPep(vNewSpectraInfo[t].first, SpectraTemp);
		if (SpectraTemp.m_vPeptides.size() == 0)
			continue;

		double Score = SpectraTemp.m_vPeptides[0].m_vlfScores[0];
		if (m_EngineType == ST_PFIND)
			Score = EValueToSmall(Score);

		if (Score >= conf.m_Filter.m_ScoreMin && Score <= conf.m_Filter.m_ScoreMax)
		{
			if ((m_EngineType == ST_SEQUEST || m_EngineType == ST_SQT)
					&& SpectraTemp.m_vPeptides[0].m_vlfScores[1] < conf.m_Filter.m_lfFixedDeltCn)
				continue;
			int IsRev = IsReverse(SpectraTemp, conf);
			if (IsRev == 2 || IsRev == 3 || IsRev == 4)
				continue;//todo pFind的分数都是0～1之间的
			vSpectraOut.push_back(vNewSpectraInfo[t]);
		}
	}
}

void CFilterUtility::Filter(const vector<SPECTRAINFO> & vSpectraInfoIn,
		vector<SPECTRAINFO> & vSpectraOut, const CConf & conf,
		const SearchEngineType & m_EngineType, const pair<int, int> & pairSimpleCharge)
{
	if (conf.m_Filter.m_bUseFDR == true)
	{
		_SepatrateOrNot(vSpectraInfoIn, vSpectraOut, conf, m_EngineType, pairSimpleCharge);
	}
	else
	{
		_Score(vSpectraInfoIn, vSpectraOut, conf, m_EngineType);
	}
}
