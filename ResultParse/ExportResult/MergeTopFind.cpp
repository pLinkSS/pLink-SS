/*
 * CMergeTopFind.cpp
 *
 *  Created on: 2010-10-12
 *      Author: Hailer
 */

#include "MergeTopFind.h"

using namespace bio_analysis;

extern ostringstream osspBuildLog;

CMergeTopFind::CMergeTopFind()
{
	// TODO Auto-generated constructor stub

}

CMergeTopFind::~CMergeTopFind()
{
	// TODO Auto-generated destructor stub
}

void CMergeTopFind::_OutPutSpectraASpFind(FILE * fout, const CMatchSpectraInfo & SpectraTemp,
		const size_t & SpectraNum)
{
	double lftmp = 0.0;
	int ntmp = 0;
	fprintf(fout, "[Spectrum"PRI_SIZE_T"]\n", SpectraNum);
	fprintf(fout, "Input=%s\n", SpectraTemp.m_strFileName.c_str());
	fprintf(fout, "Charge=%d\n", SpectraTemp.m_nCharge);
	fprintf(fout, "Intensity=%lf\n", lftmp);
	fprintf(fout, "MH=%lf\n", SpectraTemp.m_lfMH);
	fprintf(fout, "MZ=%lf\n", lftmp);
	fprintf(fout, "Candidate_Total=%d\n", ntmp);
	fprintf(fout, "ValidCandidate="PRI_SIZE_T"\n", SpectraTemp.m_vPeptides.size());
	for (size_t k = 0; k < SpectraTemp.m_vPeptides.size(); k++)
	{
		fprintf(fout, "NO"PRI_SIZE_T"_Score=%lf\n", k + 1, SpectraTemp.m_vPeptides[k].m_vlfScores[1]);//only
		fprintf(fout, "NO"PRI_SIZE_T"_EValue=%10.10e\n", k + 1, pow(10, -0.1
				* SpectraTemp.m_vPeptides[k].m_vlfScores[0]));//only pFind

		fprintf(fout, "NO"PRI_SIZE_T"_MH=%lf\n", k + 1, SpectraTemp.m_vPeptides[k].m_lfCalc_MH);
		fprintf(fout, "NO"PRI_SIZE_T"_SQ=%s\n", k + 1, SpectraTemp.m_vPeptides[k].m_strSQ.c_str());
		fprintf(fout, "NO"PRI_SIZE_T"_Proteins="PRI_SIZE_T"", k + 1, SpectraTemp.m_vPeptides[k].m_vProteinAC.size());
		for (size_t j = 0; j < SpectraTemp.m_vPeptides[k].m_vProteinAC.size(); j++)
		{
			fprintf(fout, ",%s", SpectraTemp.m_vPeptides[k].m_vProteinAC[j].c_str());
		}
		fprintf(fout, "\n");
		fprintf(fout, "NO"PRI_SIZE_T"_ProteinIDs="PRI_SIZE_T"", k + 1, SpectraTemp.m_vPeptides[k].m_vMod.size());
		for (size_t j = 0; j < SpectraTemp.m_vPeptides[k].m_vMod.size(); j++)
			fprintf(fout, ",%d", ntmp);
		fprintf(fout, "\n");
		fprintf(fout, "NO"PRI_SIZE_T"_Modify_Pos="PRI_SIZE_T"", k + 1, SpectraTemp.m_vPeptides[k].m_vMod.size());
		for (size_t j = 0; j < SpectraTemp.m_vPeptides[k].m_vMod.size(); j++)
		{
			fprintf(fout, ","PRI_SIZE_T"", SpectraTemp.m_vPeptides[k].m_vMod[j].m_tPos);
		}
		fprintf(fout, "\n");
		fprintf(fout, "NO"PRI_SIZE_T"_Modify_Name="PRI_SIZE_T"", k + 1, SpectraTemp.m_vPeptides[k].m_vMod.size());
		for (size_t j = 0; j < SpectraTemp.m_vPeptides[k].m_vMod.size(); j++)
		{
			fprintf(fout, ",%s", SpectraTemp.m_vPeptides[k].m_vMod[j].m_strModName.c_str());
		}
		fprintf(fout, "\n");
	}
}

void CMergeTopFind::_PrintfSpectraOut(FILE * fout, CMatchSpectraInfo & SpectraOut, bool & tag,
		size_t & SpectraNum)
{
	sort(SpectraOut.m_vPeptides.begin(), SpectraOut.m_vPeptides.end(), Peptide_SortScore);
	SpectraNum++;
	_OutPutSpectraASpFind(fout, SpectraOut, SpectraNum);
	SpectraOut.m_vPeptides.clear();
	SpectraOut.clear();
	tag = false;
}

void CMergeTopFind::_PeptideJoin(vector<CMatchPeptideInfo> & vPep,
		const CMatchPeptideInfo & pepTemp)
{
	for (size_t t = 0; t < vPep.size(); t++)
	{
		if (pepTemp.m_strSQ == vPep[t].m_strSQ && VecotrModEqual(pepTemp.m_vMod, vPep[t].m_vMod))
		{
			if (pepTemp.m_vlfScores[0] > vPep[t].m_vlfScores[0])
				vPep[t] = pepTemp;
			return;
		}
	}
	vPep.push_back(pepTemp);
	return;
}

void CMergeTopFind::_MergePeptide(CMatchSpectraInfo & SpectraOut,
		const CMatchSpectraInfo & SpectraTemp)
{
	for (size_t t = 0; t < SpectraTemp.m_vPeptides.size(); t++)
	{
		_PeptideJoin(SpectraOut.m_vPeptides, SpectraTemp.m_vPeptides[t]);
	}
}

void CMergeTopFind::_Merge(const OneDATASET & OneResult, const CConf & conf, int DataSetID)
{
	char szID[FILE_NAME_SIZE] = { 0 };
	sprintf(szID, "_DataSet%d_pfind.txt", DataSetID);
	string ExportpFind = conf.m_OutPutName;
	ExportpFind += szID;
	FILE * fout = fopen(ExportpFind.c_str(), "w");
	string strOutFile;
	set<string> strPath;
	vector<SPECTRAINFO> vSpectraTemp;
	vSpectraTemp.clear();
	set<string> SpectraPath;
	for (map<int, vector<SPECTRAINFO> >::const_iterator it = OneResult.second.begin(); it
			!= OneResult.second.end(); it++)
	{
		for (size_t t = 0; t < it->second.size(); t++)
		{
			vSpectraTemp.push_back(it->second[t]);
			SpectraPath.insert(it->second[t].second.m_strInPutPath);
		}
	}
	sort(vSpectraTemp.begin(), vSpectraTemp.end(), SPECTRAINFO_SortInPutPath);
	OutPutCondition(fout, OneResult.first, SpectraPath.size());
	CMatchSpectraInfo SpectraTemp, SpectraOut;
	bool tag = false;
	size_t SpectraNum = 0;
	for (size_t t = 0; t < vSpectraTemp.size(); t++)
	{
		ReadSingeSpectra(vSpectraTemp[t].first, SpectraTemp);
		if (false == tag)
		{
			SpectraOut = SpectraTemp;
			tag = true;
		}
		else
		{
			if (SpectraOut.m_strFileName == SpectraTemp.m_strFileName)
				_MergePeptide(SpectraOut, SpectraTemp);
			else
			{
				_PrintfSpectraOut(fout, SpectraOut, tag, SpectraNum);
				SpectraOut = SpectraTemp;
				tag = true;
			}
		}
	}
	_PrintfSpectraOut(fout, SpectraOut, tag, SpectraNum);
	fclose(fout);
}

void CMergeTopFind::pFindMerge_Export(const vector<OneDATASET> & ResultDataSet, const CConf & conf)
{
	//cout << "Merge to pFind..." << endl;
	osspBuildLog << "Merge to pFind..." << endl;
#ifdef DEBUG
	cout << "Merge to pFind..." << endl;
#endif
	for (size_t t = 0; t < ResultDataSet.size(); t++)
	{
		_Merge(ResultDataSet[t], conf, t + 1);
	}

}
