/*
 * CExportToIntermediateFile.cpp
 *
 *  Created on: 2010-10-25
 *      Author: chenhaifeng
 */

#include "ExportToIntermediateFile.h"

using namespace bio_analysis;

extern char SearchEngineName[4][20];

extern ostringstream osspBuildLog;
//#define DEBUG
CExportToIntermediateFile::CExportToIntermediateFile()
{
	// TODO Auto-generated constructor stub

}

CExportToIntermediateFile::~CExportToIntermediateFile()
{
	// TODO Auto-generated destructor stub
}

int structUniPepofProtein_cmp(const void * a, const void * b)
{
	structUniPepofProtein * c = (structUniPepofProtein *) a;
	structUniPepofProtein * d = (structUniPepofProtein *) b;

	if (c->UniPep > d->UniPep)
		return -1;
	return 1;
}

void CExportToIntermediateFile::_WriteSinglePeptide(FILE * fout, const CMatchPeptideInfo & PepTemp,
		const CMatchSpectraInfo & SpectraTemp)
{
	if (SpectraTemp.EnginType == ST_PFIND)
		fprintf(fout, "%2.2e\t", EValueToSmall(PepTemp.m_vlfScores[0]));
	else
		fprintf(fout, "%.2lf\t", PepTemp.m_vlfScores[0]);

	fprintf(fout, "%lf\t", PepTemp.m_lfCalc_MH - protonH);
	fprintf(fout, "%lf\t", PepTemp.m_lfDelta);
	fprintf(fout, "%lf\t", PepTemp.m_lfPPM);

	int nShiftForCxL = -1;
	if (PepTemp.m_strSQ.find("-") != -1 && PepTemp.m_strSQ.find("(") != -1)
		nShiftForCxL = PepTemp.m_strSQ.find("-") - PepTemp.m_strSQ.find("(") ;

	for (size_t t = 0; t < PepTemp.m_vMod.size(); t++)
	{
		if (PepTemp.m_vMod[t].m_tPos == 0)
		{
			fprintf(fout, ""PRI_SIZE_T",%c(%s);", PepTemp.m_vMod[t].m_tPos,
					PepTemp.m_strSQ[PepTemp.m_vMod[t].m_tPos],
					PepTemp.m_vMod[t].m_strModName.c_str());
		}
		else if (PepTemp.m_vMod[t].m_tPos >= PepTemp.m_strSQ.size() + 1)
		{
			fprintf(fout, ""PRI_SIZE_T",%c(%s);", PepTemp.m_vMod[t].m_tPos,
					PepTemp.m_strSQ[PepTemp.m_strSQ.size() - 1],
					PepTemp.m_vMod[t].m_strModName.c_str());
		}
		else
		{
			int nModSite = PepTemp.m_vMod[t].m_tPos - 1  + nShiftForCxL + 1; //加入（3）-这样的内部符号偏移后的修饰位置

			if (nModSite > PepTemp.m_strSQ.length())
				throw runtime_error("In _WriteSinglePeptide out of range");

			if (PepTemp.m_vMod[t].m_tPos <= PepTemp.m_strSQ.find("("))
				nModSite-=(nShiftForCxL+1);

			fprintf(fout, ""PRI_SIZE_T",%c(%s);", PepTemp.m_vMod[t].m_tPos,
					PepTemp.m_strSQ[nModSite],
					PepTemp.m_vMod[t].m_strModName.c_str());
		}
	}

	if (PepTemp.m_vMod.size() == 0)
		fprintf(fout, "null");

	fprintf(fout, "\t");
	fprintf(fout, "%d.%s\t", SpectraTemp.m_nDataSetID, SearchEngineName[SpectraTemp.EnginType]);
	fprintf(fout, "%d.%s\t", SpectraTemp.m_nFileID, SearchEngineName[SpectraTemp.EnginType]);

	fprintf(fout, "0\t");
	fprintf(fout, "-1\t");
	fprintf(fout, "1\t");

	if (PepTemp.m_vProteinAC.size()>0) //added at 2014.8.29
		fprintf(fout, "%s", PepTemp.m_vProteinAC[0].c_str());

	for (size_t p = 1; p < PepTemp.m_vProteinAC.size(); p++)
		fprintf(fout, "/%s", PepTemp.m_vProteinAC[p].c_str());
	fprintf(fout, "\n");
}

void CExportToIntermediateFile::_Spectra_IntermediateFile(const CConf & conf,
		const IntegratedInfo & HaveToInterData, FILE * fSS)
{
#ifdef DEBUG
	cout << "_Spectra_IntermediateFile..." << endl;
#endif
	string ExportSpectra = conf.m_outPutForder_Index + "IntermediateFile.spectra";
	FILE * fout = fopen(ExportSpectra.c_str(), "w");
	if (fout == NULL)
	{
		cout<<"Error: The file name: "<<ExportSpectra<<" is too long to generate."<<endl;
		return;
	}

	map<int, int> Sample_Spectra;
	Sample_Spectra.clear();

	fprintf(fout, "---Order\tTitle\tPeptideNum\tUniquePepNum\tSamples\tScore\tCondition\n");
	fprintf(
			fout,
			"---\tOrder\tSequence\tScore\tCalc_M\tDelta\tppm\tModification\tSampleID\tEngine\tMatchedIons\tMissCleaveNum\tRank\tProteins\n");

	int order = 1;
	CMatchSpectraInfo SpectraTemp;

	for (map<string, vector<SPECTRAINFO> >::const_iterator it =
			HaveToInterData.spec_peptide.begin(); it != HaveToInterData.spec_peptide.end(); it++)
	{
		if (HaveToInterData.m_deleteSpectra.find(it->first)
				!= HaveToInterData.m_deleteSpectra.end())
			continue;

		vector<CMatchSpectraInfo> vSpectra;
		set<string> setUniPep;
		set<int> setIDs;

		for (size_t t = 0; t < it->second.size(); t++)
		{
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);
			setUniPep.insert(SpectraTemp.m_vPeptides[0].m_strSQ);
			setIDs.insert(SpectraTemp.m_nDataSetID);
			vSpectra.push_back(SpectraTemp);
		}
		sort(vSpectra.begin(), vSpectra.end(), CMatchSpectraInfo_SortScore);
		fprintf(fout, "%d\t%s\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t", order++, it->first.c_str(), it->second.size(),
				setUniPep.size(), setIDs.size());

		if (vSpectra[0].EnginType == ST_PFIND)
			fprintf(fout, "%2.2e\t", EValueToSmall(vSpectra[0].m_vPeptides[0].m_vlfScores[0]));
		else
			fprintf(fout, "%.2lf\t", vSpectra[0].m_vPeptides[0].m_vlfScores[0]);

		if (setIDs.size() == 1)
			fprintf(fout, "Distinct\n");

		else if (setIDs.size() > 1 && setUniPep.size() == 1)
			fprintf(fout, "Consistent\n");

		else if (setIDs.size() > 1 && setUniPep.size() > 1)
			fprintf(fout, "Inconsistent\n");

		bitset<SIMPLE_MAX> m_bitSet;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			fprintf(fout, "*\t");
			fprintf(fout, ""PRI_SIZE_T",%d\t%c.%s.%c\t", t + 1, order - 1,
					vSpectra[t].m_vPeptides[0].m_cPrev, vSpectra[t].m_vPeptides[0].m_strSQ.c_str(),
					vSpectra[t].m_vPeptides[0].m_cNext);
			m_bitSet.set(vSpectra[t].m_nDataSetID - 1);
			_WriteSinglePeptide(fout, vSpectra[t].m_vPeptides[0], vSpectra[t]);
		}
		Sample_Spectra[m_bitSet.to_ulong()]++;
	}

	fclose(fout);
	nSpectra = order - 1;

	fprintf(fSS, "[Spectra]\n");
	for (map<int, int>::iterator it = Sample_Spectra.begin(); it != Sample_Spectra.end(); it++)
	{
		fprintf(fSS, "%d %d\n", it->first, it->second);
	}
	fprintf(fSS, "[END]\n");

}

void CExportToIntermediateFile::_Peptide_IntermediateFile(const CConf & conf,
		const IntegratedInfo & HaveToInterData, FILE * fSS)
{
#ifdef DEBUG
	cout << "_Peptide_IntermediateFile..." << endl;
#endif
	string ExportPeptide = conf.m_outPutForder_Index + "IntermediateFile.peptide";
	FILE * fout = fopen(ExportPeptide.c_str(), "w");

	fprintf(fout, "---#\tSQ\tSpectra\tSamples\tScore\tProtein\n");
	fprintf(
			fout,
			"---\tOrder\tSequence\tScore\tCalc_M\tDelta\tppm\tModification\tSampleID\tEngine\tMatchedIons\tMissCleaveNum\tRank\tProteins\n");

	map<int, int> Sample_Peptide;
	Sample_Peptide.clear();

	CMatchSpectraInfo SpectraTemp;
	int order = 1;

	for (map<string, vector<SPECTRAINFO> >::const_iterator it =
			HaveToInterData.peptide_spec.begin(); it != HaveToInterData.peptide_spec.end(); it++)
	{
		if (HaveToInterData.m_deletePeptide.find(it->first)
				!= HaveToInterData.m_deletePeptide.end())
			continue;

		vector<CMatchSpectraInfo> vSpectra;
		set<int> setIDs;
		set<string> UniProtein;

		for (size_t t = 0; t < it->second.size(); t++)
		{
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);
			setIDs.insert(SpectraTemp.m_nDataSetID);
			vSpectra.push_back(SpectraTemp);
			for (size_t k = 0; k < SpectraTemp.m_vPeptides[0].m_vProteinAC.size(); k++)
			{
				UniProtein.insert(SpectraTemp.m_vPeptides[0].m_vProteinAC[k]);
			}
		}
		sort(vSpectra.begin(), vSpectra.end(), CMatchSpectraInfo_SortScore);
		fprintf(fout, "%d\t%s\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t", order++, it->first.c_str(), it->second.size(),
				setIDs.size());

		if (vSpectra[0].EnginType == ST_PFIND)
			fprintf(fout, "%2.2e\t", EValueToSmall(vSpectra[0].m_vPeptides[0].m_vlfScores[0]));
		else
			fprintf(fout, "%.2lf\t", vSpectra[0].m_vPeptides[0].m_vlfScores[0]);

		set<string>::iterator k = UniProtein.begin();

		if (UniProtein.size() > 0)
		{
			fprintf(fout, "%s", (*k).c_str());

			for (k++; k != UniProtein.end(); k++)
			{
				fprintf(fout, "/%s", (*k).c_str());
			}
			fprintf(fout, "\n");
		}

		bitset<SIMPLE_MAX> m_bitSet;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			fprintf(fout, "*\t"PRI_SIZE_T",%d\t", t + 1, order - 1);
			fprintf(fout, "%s\t", vSpectra[t].m_strFileName.c_str());
			_WriteSinglePeptide(fout, vSpectra[t].m_vPeptides[0], vSpectra[t]);
			m_bitSet.set(vSpectra[t].m_nDataSetID - 1);
		}
		Sample_Peptide[m_bitSet.to_ulong()]++;
	}
	fclose(fout);
	nPeptide = order - 1;

	fprintf(fSS, "[Peptide]\n");
	for (map<int, int>::iterator it = Sample_Peptide.begin(); it != Sample_Peptide.end(); it++)
	{
		fprintf(fSS, "%d %d\n", it->first, it->second);
	}
	fprintf(fSS, "[END]\n");

}

void CExportToIntermediateFile::_Protein_IntermediateFile(const CConf & conf,
		const IntegratedInfo & HaveToInterData, FILE * fSS, map<string, CProteinInfo> & ProteinInfo)
{
#ifdef DEBUG
	cout << "_Protein_IntermediateFile..." << endl;
#endif
	////////////////////////////////////////////////////////////////////////////
	//	cout << "pQuant" << endl;
	string ExportpQuant = conf.m_outPutForder_Index + "IntermediateFile.pquant";
	FILE * fout_pQuant = fopen(ExportpQuant.c_str(), "w");

	map<int, int> Sample_Protein;
	Sample_Protein.clear();

	fprintf(fout_pQuant, "Order\tProteinAC\tCoverage\tPeptideNum\tUniquePepNum\n");
	fprintf(
			fout_pQuant,
			"\tOrder\tUniqueProNum\tSpectrum\tCharge\tSequence\tScore\tCalc_M\tDelta\tModification\tSampleID\tEngine\tRank\tProteins\n");

	map<string, int> peptide_protein_num;
	for (map<string, set<string> >::const_iterator it = HaveToInterData.peptide_protein_set.begin(); it
			!= HaveToInterData.peptide_protein_set.end(); it++)
	{
		peptide_protein_num[it->first] = it->second.size();
	}
	int order_pQuant = 1;
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//	cout << "Protein" << endl;
	string ExportProtein = conf.m_outPutForder_Index + "IntermediateFile.protein";
	FILE * fout = fopen(ExportProtein.c_str(), "w");

	fprintf(
			fout,
			"---Order\tProteinAC\tMW\tpI\tCoverage\tUniquePepNum\tSpecNum\tNonModifiedSpecNum\tModifiedSpecNum\tUniqueModifiedPepNum\tSamples\tDescription\n");
	fprintf(
			fout,
			"---\tOrder\tSpectrum\tSequence\tScore\tCalc_M\tDelta_M\tModification\tSampleID\tEngine\tMatchedIons\tMissCleaveNum\tRank\tProteins\n");

	CMatchSpectraInfo SpectraTemp;
	int order = 1;
	map<int, int, greater<int> > UniPepCount;
	set<string> setpQuantProtein;
	map<string, int> UniPepofProtein;
	map<string, set<string> > m_mapSameSet = HaveToInterData.m_mapSameSet;
	map<string, set<string> > m_mapSubSet = HaveToInterData.m_mapSubSet;
	map<string, vector<SPECTRAINFO> > protein_peptide = HaveToInterData.protein_peptide;

	int nDecoy = 0;
	int nTarget = 0;

	for (map<string, vector<SPECTRAINFO> >::const_iterator it =
			HaveToInterData.protein_peptide.begin(); it != HaveToInterData.protein_peptide.end(); it++)
	{
		if (HaveToInterData.m_setInvalidProteins.find(it->first)
				!= HaveToInterData.m_setInvalidProteins.end())
			continue;

		if (HaveToInterData.m_deleteProtein.find(it->first)
				!= HaveToInterData.m_deleteProtein.end())
			continue;

		set<string> setUniPep;
		set<int> setIDs;
		set<string> setUniPepWithMod;
		vector<CMatchSpectraInfo> vSpectra;
		size_t nModifiedSpecNum = 0;

		for (size_t t = 0; t < it->second.size(); t++)
		{
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);
			setUniPep.insert(SpectraTemp.m_vPeptides[0].m_strSQ);
			setIDs.insert(SpectraTemp.m_nDataSetID);

			if (SpectraTemp.m_vPeptides[0].m_vMod.size() > 0)
			{
				nModifiedSpecNum++;
				string strSQMod = GetPeptideWithMod(SpectraTemp.m_vPeptides[0]);
				setUniPepWithMod.insert(strSQMod);
			}
			vSpectra.push_back(SpectraTemp);
		}
		sort(vSpectra.begin(), vSpectra.end(), CMatchSpectraInfo_SortScore);

		string Description = "null";

		fprintf(fout, "[1Title]\n");
		if (ProteinInfo.find(it->first) == ProteinInfo.end())
		{
			fprintf(fout, "%d\t%s\t0\t0\t0\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\tnull\n", order++,
					it->first.c_str(), setUniPep.size(), it->second.size(), it->second.size()
							- nModifiedSpecNum, nModifiedSpecNum, setUniPepWithMod.size(),
					setIDs.size());
		}
		else
		{
			fprintf(fout, "%d\t%s\t%lf\t%lf\t%.2lf\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t%s\n", order++,
					it->first.c_str(), ProteinInfo[it->first].m_MW, ProteinInfo[it->first].m_pI,
					ProteinInfo[it->first].m_Coverage * 100, setUniPep.size(), it->second.size(),
					it->second.size() - nModifiedSpecNum, nModifiedSpecNum,
					setUniPepWithMod.size(), setIDs.size(), ProteinInfo[it->first].m_strDE.c_str());
		}

		UniPepCount[setUniPep.size()]++;
		UniPepofProtein[it->first] = setUniPep.size();

		fprintf(fout, "[2Title]\n");
		//////////////////////////////////////////////////
		bool bDecoy = false;
		if (IsProteinReverse(it->first, conf))
			bDecoy = true;
		for (set<string>::const_iterator t = m_mapSameSet[it->first].begin(); t
				!= m_mapSameSet[it->first].end(); t++)
		{
			if (IsProteinReverse(*t, conf))
				bDecoy = true;
		}
		//		for (set<string>::const_iterator t = m_mapSubSet[it->first].begin(); t
		//				!= m_mapSubSet[it->first].end(); t++)
		//		{
		//			if (IsProteinReverse(*t, conf))
		//				bDecoy = true;
		//		}
		if (bDecoy == true)
			nDecoy++;
		else
			nTarget++;
		//////////////////////////////////////////////////
		for (size_t t = 0; t < it->second.size(); t++)
		{
			fprintf(fout, "*\t"PRI_SIZE_T",%d\t", t + 1, order - 1);
			fprintf(fout, "%s\t", vSpectra[t].m_strFileName.c_str());
			fprintf(fout, "%c.%s.%c\t", vSpectra[t].m_vPeptides[0].m_cPrev,
					vSpectra[t].m_vPeptides[0].m_strSQ.c_str(), vSpectra[t].m_vPeptides[0].m_cNext);
			_WriteSinglePeptide(fout, vSpectra[t].m_vPeptides[0], vSpectra[t]);
		}
		//cout
		///////////////////////////////////////////////////////////////
		if (setpQuantProtein.find(it->first) == setpQuantProtein.end())
		{
			if (ProteinInfo.find(it->first) == ProteinInfo.end())
			{
				fprintf(fout_pQuant, "P%d\t%s\t0\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t\n", order_pQuant++, it->first.c_str(),
						it->second.size(), setUniPep.size());
			}
			else
			{
				fprintf(fout_pQuant, "P%d\t%s\t%lf%%\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t\n", order_pQuant++,
						it->first.c_str(), ProteinInfo[it->first].m_Coverage * 100,
						it->second.size(), setUniPep.size());
			}

			bitset<SIMPLE_MAX> m_bitSet;
			for (size_t t = 0; t < it->second.size(); t++)
			{
				fprintf(fout_pQuant, "\tS"PRI_SIZE_T"\t%d\t", t + 1,
						peptide_protein_num[vSpectra[t].m_vPeptides[0].m_strSQ]);

				fprintf(fout_pQuant, "%s\t%d\t", vSpectra[t].m_strFileName.c_str(),
						vSpectra[t].m_nCharge);

				fprintf(fout_pQuant, "%c.%s.%c\t", vSpectra[t].m_vPeptides[0].m_cPrev,
						vSpectra[t].m_vPeptides[0].m_strSQ.c_str(),
						vSpectra[t].m_vPeptides[0].m_cNext);

				m_bitSet.set(vSpectra[t].m_nDataSetID - 1);

				CMatchPeptideInfo PepTemp = vSpectra[t].m_vPeptides[0];
				WriteSingePeptide(fout_pQuant, PepTemp, vSpectra[t]);
			}
			//cout << m_bitSet.to_ulong() << endl;
			Sample_Protein[m_bitSet.to_ulong()]++;
			setpQuantProtein.insert(it->first);
		}
		/////////////////////////////////////////////////////////////////
		if (m_mapSameSet[it->first].size() > 0 || m_mapSubSet[it->first].size() > 0)
			fprintf(fout, "[3Title]\n");

		for (set<string>::const_iterator t = m_mapSameSet[it->first].begin(); t
				!= m_mapSameSet[it->first].end(); t++)
		{
			fprintf(fout, "SameSet");

			setUniPep.clear();
			setIDs.clear();
			setUniPepWithMod.clear();
			vSpectra.clear();
			nModifiedSpecNum = 0;

			for (size_t i = 0; i < protein_peptide[*t].size(); i++)
			{
				ReadSingeSpectraFisrtPep(protein_peptide[*t][i].first, SpectraTemp);
				setUniPep.insert(SpectraTemp.m_vPeptides[0].m_strSQ);
				setIDs.insert(SpectraTemp.m_nDataSetID);

				if (SpectraTemp.m_vPeptides[0].m_vMod.size() > 0)
				{
					nModifiedSpecNum++;
					string strSQMod = GetPeptideWithMod(SpectraTemp.m_vPeptides[0]);
					setUniPepWithMod.insert(strSQMod);
				}

				vSpectra.push_back(SpectraTemp);
			}
			sort(vSpectra.begin(), vSpectra.end(), CMatchSpectraInfo_SortScore);
			UniPepofProtein[*t] = setUniPep.size();
			if (ProteinInfo.find((*t)) == ProteinInfo.end())
			{
				fprintf(fout, "\t%s\t0\t0\t0\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\tnull\n", (*t).c_str(),
						setUniPep.size(), protein_peptide[*t].size(), protein_peptide[*t].size()
								- nModifiedSpecNum, nModifiedSpecNum, setUniPepWithMod.size(),
						setIDs.size());
			}
			else
			{
				fprintf(fout, "\t%s\t%lf\t%lf\t%.2lf\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t%s\n", (*t).c_str(),
						ProteinInfo[(*t)].m_MW, ProteinInfo[(*t)].m_pI,
						ProteinInfo[(*t)].m_Coverage * 100, setUniPep.size(),
						protein_peptide[*t].size(), protein_peptide[*t].size() - nModifiedSpecNum,
						nModifiedSpecNum, setUniPepWithMod.size(), setIDs.size(),
						ProteinInfo[(*t)].m_strDE.c_str());
			}

			for (size_t k = 0; k < vSpectra.size(); k++)
			{
				fprintf(fout, "*\t"PRI_SIZE_T",%d\t", k + 1, order - 1);
				fprintf(fout, "%s\t", vSpectra[k].m_strFileName.c_str());

				fprintf(fout, "%c.%s.%c\t", vSpectra[k].m_vPeptides[0].m_cPrev,
						vSpectra[k].m_vPeptides[0].m_strSQ.c_str(),
						vSpectra[k].m_vPeptides[0].m_cNext);

				_WriteSinglePeptide(fout, vSpectra[k].m_vPeptides[0], vSpectra[k]);
			}
			/////////////////////////////////////////////////////////////////////////////
			if (setpQuantProtein.find((*t)) == setpQuantProtein.end())
			{
				if (ProteinInfo.find((*t)) == ProteinInfo.end())
				{
					fprintf(fout_pQuant, "P%d\t%s\t0\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t\n", order_pQuant++, (*t).c_str(),
							vSpectra.size(), setUniPep.size());
				}
				else
				{
					fprintf(fout_pQuant, "P%d\t%s\t%lf%%\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t\n", order_pQuant++,
							(*t).c_str(), ProteinInfo[(*t)].m_Coverage * 100, vSpectra.size(),
							setUniPep.size());
				}

				bitset<SIMPLE_MAX> m_bitSet;

				for (size_t k = 0; k < vSpectra.size(); k++)
				{
					fprintf(fout_pQuant, "\tS"PRI_SIZE_T"\t%d\t", k + 1,
							peptide_protein_num[vSpectra[k].m_vPeptides[0].m_strSQ]);

					fprintf(fout_pQuant, "%s\t%d\t", vSpectra[k].m_strFileName.c_str(),
							vSpectra[k].m_nCharge);

					fprintf(fout_pQuant, "%c.%s.%c\t", vSpectra[k].m_vPeptides[0].m_cPrev,
							vSpectra[k].m_vPeptides[0].m_strSQ.c_str(),
							vSpectra[k].m_vPeptides[0].m_cNext);

					m_bitSet.set(vSpectra[k].m_nDataSetID - 1);

					CMatchPeptideInfo PepTemp = vSpectra[k].m_vPeptides[0];
					WriteSingePeptide(fout_pQuant, PepTemp, vSpectra[k]);
				}
				Sample_Protein[m_bitSet.to_ulong()]++;
				setpQuantProtein.insert((*t));
			}
			////////////////////////////////////////////////////////////////
		}
		for (set<string>::const_iterator t = m_mapSubSet[it->first].begin(); t
				!= m_mapSubSet[it->first].end(); t++)
		{
			fprintf(fout, "SubSet");

			setUniPep.clear();
			setIDs.clear();
			setUniPepWithMod.clear();
			vSpectra.clear();
			nModifiedSpecNum = 0;

			for (size_t i = 0; i < protein_peptide[*t].size(); i++)
			{
				ReadSingeSpectraFisrtPep(protein_peptide[*t][i].first, SpectraTemp);
				setUniPep.insert(SpectraTemp.m_vPeptides[0].m_strSQ);
				setIDs.insert(SpectraTemp.m_nDataSetID);

				if (SpectraTemp.m_vPeptides[0].m_vMod.size() > 0)
				{
					nModifiedSpecNum++;
					string strSQMod = GetPeptideWithMod(SpectraTemp.m_vPeptides[0]);
					setUniPepWithMod.insert(strSQMod);
				}

				vSpectra.push_back(SpectraTemp);
			}
			sort(vSpectra.begin(), vSpectra.end(), CMatchSpectraInfo_SortScore);

			UniPepofProtein[*t] = setUniPep.size();

			if (ProteinInfo.find((*t)) == ProteinInfo.end())
			{
				fprintf(fout, "\t%s\t0\t0\t0\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\tnull\n", (*t).c_str(),
						setUniPep.size(), protein_peptide[*t].size(), protein_peptide[*t].size()
								- nModifiedSpecNum, nModifiedSpecNum, setUniPepWithMod.size(),
						setIDs.size());
			}
			else
			{
				fprintf(fout, "\t%s\t%lf\t%lf\t%.2lf\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t%s\n", (*t).c_str(),
						ProteinInfo[(*t)].m_MW, ProteinInfo[(*t)].m_pI,
						ProteinInfo[(*t)].m_Coverage * 100, setUniPep.size(),
						protein_peptide[*t].size(), protein_peptide[*t].size() - nModifiedSpecNum,
						nModifiedSpecNum, setUniPepWithMod.size(), setIDs.size(),
						ProteinInfo[(*t)].m_strDE.c_str());
			}

			for (size_t k = 0; k < vSpectra.size(); k++)
			{
				fprintf(fout, "*\t"PRI_SIZE_T",%d\t", k + 1, order - 1);

				fprintf(fout, "%s\t", vSpectra[k].m_strFileName.c_str());

				fprintf(fout, "%c.%s.%c\t", vSpectra[k].m_vPeptides[0].m_cPrev,
						vSpectra[k].m_vPeptides[0].m_strSQ.c_str(),
						vSpectra[k].m_vPeptides[0].m_cNext);

				_WriteSinglePeptide(fout, vSpectra[k].m_vPeptides[0], vSpectra[k]);
			}
			/////////////////////////////////////////////////////////////////////////////
			if (setpQuantProtein.find((*t)) == setpQuantProtein.end())
			{
				if (ProteinInfo.find((*t)) == ProteinInfo.end())
				{

					fprintf(fout_pQuant, "P%d\t%s\t0\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t\n", order_pQuant++, (*t).c_str(),
							vSpectra.size(), setUniPep.size());
				}
				else
				{
					fprintf(fout_pQuant, "P%d\t%s\t%lf%%\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t\n", order_pQuant++,
							(*t).c_str(), ProteinInfo[(*t)].m_Coverage * 100, vSpectra.size(),
							setUniPep.size());
				}

				bitset<SIMPLE_MAX> m_bitSet;

				for (size_t k = 0; k < vSpectra.size(); k++)
				{
					fprintf(fout_pQuant, "\tS"PRI_SIZE_T"\t%d\t", k + 1,
							peptide_protein_num[vSpectra[k].m_vPeptides[0].m_strSQ]);

					fprintf(fout_pQuant, "%s\t%d\t", vSpectra[k].m_strFileName.c_str(),
							vSpectra[k].m_nCharge);

					fprintf(fout_pQuant, "%c.%s.%c\t", vSpectra[k].m_vPeptides[0].m_cPrev,
							vSpectra[k].m_vPeptides[0].m_strSQ.c_str(),
							vSpectra[k].m_vPeptides[0].m_cNext);

					m_bitSet.set(vSpectra[k].m_nDataSetID - 1);

					CMatchPeptideInfo PepTemp = vSpectra[k].m_vPeptides[0];

					WriteSingePeptide(fout_pQuant, PepTemp, vSpectra[k]);

				}
				Sample_Protein[m_bitSet.to_ulong()]++;
				setpQuantProtein.insert((*t));
			}
			////////////////////////////////////////////////////////////////
		}
	}

	fclose(fout_pQuant);

	fprintf(fout, "[END]\n");
	////////////////////////////////////

	ofstream TDout(string(conf.m_OutPutForder + "TargetDecoyNum.txt").c_str());
	TDout << "Target: " << nTarget << endl;
	TDout << "Decoy: " << nDecoy << endl;
	TDout << "Decoy/Target: " << (double) nDecoy / nTarget << endl;
	TDout << "2*Decoy/(Target+Decoy): " << 2.0 * nDecoy / (nTarget + nDecoy) << endl;
	TDout.close();
	////////////////////////////////////

	//todo 按从大到小排序//
	nProGroup = order - 1;
	fprintf(fout, "\n\n---summary---\n");
	//Protein group:
	fprintf(fout, "Protein groups: %d\n", nProGroup);
	fprintf(fout, "UniquePepNum\t");

	for (map<int, int>::iterator it = UniPepCount.begin(); it != UniPepCount.end(); it++)
	{
		fprintf(fout, "%d\t", it->first);
	}

	fprintf(fout, "\nNum of groups\t");
	for (map<int, int>::iterator it = UniPepCount.begin(); it != UniPepCount.end(); it++)
	{
		fprintf(fout, "%d\t", it->second);
	}

	fprintf(fout, "\n  %% in groups\t");
	for (map<int, int>::iterator it = UniPepCount.begin(); it != UniPepCount.end(); it++)
	{
		fprintf(fout, "%lf%%\t", 100.00 * it->second / (double) nProGroup);
	}
	//////////////////////////////////////////////////////////////////////
	//Total_proteins:
	fprintf(fout, "\n\nTotal_proteins: "PRI_SIZE_T"\n", UniPepofProtein.size());
	map<int, int, greater<int> > UniPepCountofProtein2;
	for (map<string, int>::iterator it = UniPepofProtein.begin(); it != UniPepofProtein.end(); it++)
	{
		UniPepCountofProtein2[it->second]++;
	}

	fprintf(fout, "UniquePepNum\t");
	for (map<int, int>::iterator it = UniPepCountofProtein2.begin(); it
			!= UniPepCountofProtein2.end(); it++)
	{
		fprintf(fout, "%d\t", it->first);
	}

	fprintf(fout, "\nNum of proteins\t");
	for (map<int, int>::iterator it = UniPepCountofProtein2.begin(); it
			!= UniPepCountofProtein2.end(); it++)
	{
		fprintf(fout, "%d\t", it->second);
	}

	fprintf(fout, "\n  %% in proteins\t");
	for (map<int, int>::iterator it = UniPepCountofProtein2.begin(); it
			!= UniPepCountofProtein2.end(); it++)
	{
		fprintf(fout, "%lf%%\t", 100.00 * it->second / (double) UniPepofProtein.size());
	}

	fprintf(fout, "\n");
	fclose(fout);

	string ExportProtein_simple = conf.m_outPutForder_Index + conf.m_OutPutFile
			+ ".proteins_sample";
	fout = fopen(ExportProtein_simple.c_str(), "w");

	fprintf(fout, "Order\tDescription\tUniquePepNum\tCoverage\tProteinAC\n");
	order = 0;
	structUniPepofProtein * m_structUniPepofProtein =
			new structUniPepofProtein[UniPepofProtein.size() + 1];
	for (map<string, int>::iterator it = UniPepofProtein.begin(); it != UniPepofProtein.end(); it++)
	{
		m_structUniPepofProtein[order].strPro = it->first;
		m_structUniPepofProtein[order].UniPep = it->second;
		order++;
	}
	nProtein = order;
	//sort(m_structUniPepofProtein, m_structUniPepofProtein + order);
	qsort(m_structUniPepofProtein, order, sizeof(m_structUniPepofProtein[0]),
			structUniPepofProtein_cmp);

	for (int i = 0; i < order; i++)
	{
		if (ProteinInfo.find(m_structUniPepofProtein[i].strPro) == ProteinInfo.end())
		{
			fprintf(fout, "%d\tnull\t%d\t0\t%s\n", i + 1, m_structUniPepofProtein[i].UniPep,
					m_structUniPepofProtein[i].strPro.c_str());
		}
		else
		{
			fprintf(fout, "%d\t%s\t%d\t%lf%%\t%s\n", i + 1,
					ProteinInfo[m_structUniPepofProtein[i].strPro].m_strDE.c_str(),
					m_structUniPepofProtein[i].UniPep,
					ProteinInfo[m_structUniPepofProtein[i].strPro].m_Coverage * 100,
					m_structUniPepofProtein[i].strPro.c_str());
		}
	}

	delete[] m_structUniPepofProtein;
	fclose(fout);

	fprintf(fSS, "[Protein]\n");
	for (map<int, int>::iterator it = Sample_Protein.begin(); it != Sample_Protein.end(); it++)
	{
		fprintf(fSS, "%d %d\n", it->first, it->second);
	}
	fprintf(fSS, "[END]\n");

}

void CExportToIntermediateFile::_plabel_IntermediateFile(const CConf & conf,
		const IntegratedInfo & HaveToInterData)
{
	/*如果是以前的索引已经存在，则这里不能再输出了，因为这里的m_vEngineInputInfo[t].m_vConditions[k]
	 * 信息已经在第一次运行完成之后没有保存，丢失了结果，所以，如果只读索引的话，这里不能再输出了
	 * */
#ifdef DEBUG
	cout << "_plabel_IntermediateFile..." << endl;
#endif
#ifdef WIN32
		mkdir(conf.m_outPutForder_pLabel.c_str());
#else			
		mkdir(conf.m_outPutForder_pLabel.c_str(),0775);
#endif

	set<string> vStrPath;

	for (size_t t = 0; t < conf.m_vEngineInputInfo.size(); t++)
	{
		for (size_t k = 0; k < conf.m_vEngineInputInfo[t].m_vConditions.size(); k++)
		{
			vStrPath.insert(conf.m_vEngineInputInfo[t].m_vConditions[k].m_strInputPath);
		}
	}

	/*第二次再过滤的时候，condition的信息全部丢失，不能再写*/
	//cout << "vStrPath.size()" << vStrPath.size() << endl;
	if (vStrPath.size() <= 0)
		return;

	string ExportpLabel_Name = conf.m_outPutForder_pLabel + "plabel_Name.txt";
	FILE * fout = fopen(ExportpLabel_Name.c_str(), "w");
	FILE ** foutlabel = new FILE*[vStrPath.size() + 1];

	int *nTotal = new int[vStrPath.size() + 1];
	int *move = new int[vStrPath.size() + 1];
	int cnt = 0;
	map<string, int> mapMod;

	for (set<string>::const_iterator it = vStrPath.begin(); it != vStrPath.end(); it++, cnt++)
	{
		char chrTmp[PATH_MAX];
		sprintf(chrTmp, "InputFILE%d.plabel", cnt + 1);

		string ExportpLabel = conf.m_outPutForder_pLabel + chrTmp;
		fprintf(fout, "%s\n", ExportpLabel.c_str());
		foutlabel[cnt] = fopen(ExportpLabel.c_str(), "w");

		fprintf(foutlabel[cnt], "[FilePath]\n");
		fprintf(foutlabel[cnt], "File_Path=%s\n", (*it).c_str());//todo 这里的路径没有添加
		fprintf(foutlabel[cnt], "[Modification]\n");

		int order = 1;
		for (set<string>::iterator it = HaveToInterData.m_setAllModName.begin(); it
				!= HaveToInterData.m_setAllModName.end(); ++it)
		{
			fprintf(foutlabel[cnt], "%d=%s\n", order, it->c_str());
			mapMod[*it] = order;
			order++;
		}

		fprintf(foutlabel[cnt], "[xlink]\nxlink=NULL\n");
		fprintf(foutlabel[cnt], "[Total]\ntotal=" /*HaveToInterData.spec_peptide.size()*/);//todo 这个是假的

		move[cnt] = ftell(foutlabel[cnt]);
		fprintf(foutlabel[cnt], "                             \n");
		nTotal[cnt] = 0;
	}
	fclose(fout);

	for (map<string, vector<SPECTRAINFO> >::const_iterator it =
			HaveToInterData.spec_peptide.begin(); it != HaveToInterData.spec_peptide.end(); it++)
	{
		//	now2++;
		////////////////////////////////
		//todo 2011/03/10添加的，郴哥发现怒nubmer of peptide >=2时过滤的结果是772，而pLabel确实1354
		if (HaveToInterData.m_deleteSpectra.find(it->first)
				!= HaveToInterData.m_deleteSpectra.end())
			continue;
		/////////////////////////////
		///////////////////////////////////////////////
		for (size_t t = 0; t < it->second.size(); t++)
		{
			//now++;
			CMatchSpectraInfo SpectraTemp;
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);
			int
					pos =
							setFindInt(vStrPath, conf.m_vEngineInputInfo[SpectraTemp.m_nDataSetID
									- 1].m_vConditions[SpectraTemp.m_nFileID - 1].m_strInputPath);

			nTotal[pos]++;

			fprintf(foutlabel[pos], "[Spectrum%d]\nname=%s\n", nTotal[pos],
					stringToUpper(it->first).c_str());

			fprintf(foutlabel[pos], "pep%d=0 %s", 1, SpectraTemp.m_vPeptides[0].m_strSQ.c_str());

			if (SpectraTemp.EnginType == ST_PFIND)
				fprintf(foutlabel[pos], " %10.10e", EValueToSmall(
						SpectraTemp.m_vPeptides[0].m_vlfScores[0]));
			else
				fprintf(foutlabel[pos], " %.6lf", SpectraTemp.m_vPeptides[0].m_vlfScores[0]);

			for (size_t j = 0; j < SpectraTemp.m_vPeptides[0].m_vMod.size(); j++)
			{
				fprintf(foutlabel[pos], " "PRI_SIZE_T",%d", SpectraTemp.m_vPeptides[0].m_vMod[j].m_tPos,
						mapMod[SpectraTemp.m_vPeptides[0].m_vMod[j].m_strModName]);
			}
			fprintf(foutlabel[pos], "\n");
		}
	}
	//cout << "now = " << now << endl;
	cnt = 0;
	//cout << "now2" << now2 << endl;
	for (set<string>::const_iterator it = vStrPath.begin(); it != vStrPath.end(); it++, cnt++)
	{
		fseek(foutlabel[cnt], move[cnt], SEEK_SET);
		fprintf(foutlabel[cnt], "%d", nTotal[cnt]);
		fclose(foutlabel[cnt]);
	}

	delete[] foutlabel;
	delete[] move;
	delete[] nTotal;
}

void CExportToIntermediateFile::_Statistics(const CConf & conf)
{
#ifdef DEBUG
	cout << "_Statistics..." << endl;
#endif
	string ExportStatistics = conf.m_OutPutName + ".statistics.txt";
	FILE * fout = fopen(ExportStatistics.c_str(), "w");

	fprintf(fout, "Spectra_Num = %d\n", nSpectra);
	fprintf(fout, "Peptide_Num = %d\n", nPeptide);
	fprintf(fout, "Protein_Num = %d\n", nProtein);
	fprintf(fout, "Protein_Group_Num = %d\n", nProGroup);

	osspBuildLog << "   Spectra_Num = " << nSpectra << endl;
	osspBuildLog << "   Peptide_Num = " << nPeptide << endl;
	osspBuildLog << "   Protein_Num = " << nProtein << endl;
	osspBuildLog << "   Protein_Group_Num = " << nProGroup << endl;
#ifdef DEBUG
	cout << "   Spectra_Num = " << nSpectra << endl;
	cout << "   Peptide_Num = " << nPeptide << endl;
	cout << "   Protein_Num = " << nProtein << endl;
	cout << "   Protein_Group_Num = " << nProGroup << endl;
#endif

	fclose(fout);
}

void CExportToIntermediateFile::_ProteinInferenceProblem(const CConf & conf,
		const IntegratedInfo & HaveToInterData, map<string, CProteinInfo> & ProteinInfo)
{
#ifdef DEBUG
	cout << "_ProteinInferenceProblem..." << endl;
#endif
	map<string, double> ProteinInferenceScore;
	for (map<string, vector<SPECTRAINFO> >::const_iterator it =
			HaveToInterData.protein_peptide.begin(); it != HaveToInterData.protein_peptide.end(); it++)
	{
		double ProteinScoreTemp = 0.0;
		ProteinScoreTemp = ProteinInfo[it->first].m_Coverage;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			ProteinScoreTemp += it->second[t].second.m_lfScore;
		}
	}
}

void CExportToIntermediateFile::ExportToIntermediateFile(const IntegratedInfo & HaveToInterData,
		const CConf & conf)
{
#ifdef DEBUG
	cout << "ExportToIntermediateFile..." << endl;
#endif
	map<string, CProteinInfo> ProteinInfo;
	ProteinInfo.clear();

	if (conf.m_vProDBPath.size() > 0)
	{
		ReadProteinANDCalculation(conf, ProteinInfo, HaveToInterData.protein_peptide);
#ifdef DEBUG
		cout << "ReadProteinANDCalculation Over" << endl;
#endif
	}

	else
	{
#ifdef DEBUG
		cout << "No Database" << endl;
#endif
		//开始这里忘记了，2010年8月6日在家
		//cout << "conf.m_outPutForder_Index" << endl;
		string strFASTA = conf.m_outPutForder_Index + conf.m_OutPutFile + ".fasta";
		FILE * fFASTA = fopen(strFASTA.c_str(), "w");
		for (map<string, vector<SPECTRAINFO> >::const_iterator it =
				HaveToInterData.protein_peptide.begin(); it
				!= HaveToInterData.protein_peptide.end(); it++)
		{
			//cout << it->first << endl;
			fprintf(fFASTA, ">%s\n", it->first.c_str());
		}
		fclose(fFASTA);
	}

	//cout << "Export To Intermediate File..." << endl;
	osspBuildLog << "Export To Intermediate File..." << endl;

	string strfSS = conf.m_outPutForder_Java + "Sample_VennDiagram.txt";
	FILE * fSS = fopen(strfSS.c_str(), "w");

	for (size_t t = 0; t < conf.m_vEngineInputInfo.size(); t++)
	{
		fprintf(fSS, ""PRI_SIZE_T".%s\n", t + 1, SearchEngineName[conf.m_vEngineInputInfo[t].m_Type]);
	}
	fprintf(fSS, "[END]\n");

	_Spectra_IntermediateFile(conf, HaveToInterData, fSS);
	_Peptide_IntermediateFile(conf, HaveToInterData, fSS);
	_Protein_IntermediateFile(conf, HaveToInterData, fSS, ProteinInfo);
	_plabel_IntermediateFile(conf, HaveToInterData);

	_ProteinInferenceProblem(conf, HaveToInterData, ProteinInfo);

	_Statistics(conf);

	fclose(fSS);
}
