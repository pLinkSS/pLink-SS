#include "DataIntegrate.h"

extern ostringstream osspBuildLog;
using namespace bio_analysis;

CDataIntegrate::CDataIntegrate()
{
	// TODO Auto-generated constructor stub

}

CDataIntegrate::~CDataIntegrate()
{
	// TODO Auto-generated destructor stub
}

int CDataIntegrate::_IsSameSetOrSubSet(const set<string> & a, const set<string> & b)
{
	for (set<string>::const_iterator it = b.begin(); it != b.end(); it++)
	{
		if (a.find(*it) == a.end())
			return 0;
	}

	if (a.size() == b.size())
		return 1;//same set

	else
		return 2;// sub set
}

void CDataIntegrate::_ProteinInfer(IntegratedInfo & HaveToInterData)
{
	for (map<string, vector<SPECTRAINFO> >::iterator it = HaveToInterData.protein_peptide.begin(); it
			!= HaveToInterData.protein_peptide.end(); it++)
	{
		string strTemp2 = it->first;

		if (HaveToInterData.m_setInvalidProteins.find(strTemp2)
				!= HaveToInterData.m_setInvalidProteins.end())
			continue;

		//subset,sameset的定义， 谱or肽段
		set<string> pep1;
		for (size_t i = 0; i < it->second.size(); i++)
			//pep1.insert(it->second[i].second.m_strFirstPeptide);
			pep1.insert(it->second[i].second.m_strInPutPath);

		set<string> proteins;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			CMatchSpectraInfo SpectraTemp;
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);

			const vector<string> vstrProtein = SpectraTemp.m_vPeptides[0].m_vProteinAC;

			for (size_t k = 0; k < vstrProtein.size(); k++)
			{
				string strTemp1 = vstrProtein[k];
				if (strTemp1 == strTemp2)
					continue;

				if (proteins.find(strTemp1) != proteins.end())
					continue;

				if (it->second.size() < HaveToInterData.protein_peptide[vstrProtein[k]].size())
					continue;

				proteins.insert(strTemp1);
				set<string> pep2;

				for (size_t i = 0; i < HaveToInterData.protein_peptide[vstrProtein[k]].size(); i++)
					pep2.insert(
							HaveToInterData.protein_peptide[vstrProtein[k]][i].second.m_strInPutPath);//todo
				//pep2.insert(
				//	HaveToInterData.protein_peptide[vstrProtein[k]][i].second.m_strFirstPeptide);

				int nRes = _IsSameSetOrSubSet(pep1, pep2);
				if (1 == nRes)
				{
					HaveToInterData.m_mapSameSet[strTemp2].insert(strTemp1);
					HaveToInterData.m_setInvalidProteins.insert(strTemp1);
				}

				else if (2 == nRes)
				{
					HaveToInterData.m_mapSubSet[strTemp2].insert(strTemp1);
					HaveToInterData.m_setInvalidProteins.insert(strTemp1);
				}
			}
		}
	}
}

void CDataIntegrate::_DistinctPeptide(FILE * fout, pair<string, vector<SPECTRAINFO> > pairProTemp,
		const CConf & conf, IntegratedInfo & HaveToInterData)
{
	vector<PEPTIDEINFO> vPairPep;

	for (size_t t = 0; t < pairProTemp.second.size(); t++)
	{
		PEPTIDEINFO pairPepTemp;
		CMatchSpectraInfo SpectraTemp;

		ReadSingeSpectraFisrtPep(pairProTemp.second[t].first, SpectraTemp);
		pairPepTemp.first = SpectraTemp.m_vPeptides[0].m_strSQ;
		pairPepTemp.second = SpectraTemp.m_vPeptides[0].m_vMod;

		vPairPep.push_back(pairPepTemp);
	}

	sort(vPairPep.begin(), vPairPep.end(), PEP_SORT_SQ);
	fprintf(fout, "%s\n-------\n", pairProTemp.first.c_str());
	//OutSTY(fout, vPairPep, conf, HaveToInterData);
}

void CDataIntegrate::_CExportIndex(const vector<OneDATASET> & ResultDataSet, const CConf & conf,
		IntegratedInfo & HaveToInterData)
{
	for (size_t Rt = 0; Rt < ResultDataSet.size(); Rt++)
	{
		for (map<int, vector<SPECTRAINFO> >::const_iterator it = ResultDataSet[Rt].second.begin(); it
				!= ResultDataSet[Rt].second.end(); it++)
		{
			if (!InChargeState(it->first, conf))
				continue;

			CMatchSpectraInfo SpectraTemp;
			for (size_t t = 0; t < it->second.size(); t++)
			{
				ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);
				HaveToInterData.peptide_spec[SpectraTemp.m_vPeptides[0].m_strSQ].push_back(
						it->second[t]);

				HaveToInterData.spec_peptide[SpectraTemp.m_strFileName].push_back(it->second[t]);
				for (size_t k = 0; k < SpectraTemp.m_vPeptides[0].m_vProteinAC.size(); k++)
				{
					HaveToInterData.protein_peptide[SpectraTemp.m_vPeptides[0].m_vProteinAC[k]].push_back(
							it->second[t]);

					HaveToInterData.peptide_protein_set[SpectraTemp.m_vPeptides[0].m_strSQ].insert(
							SpectraTemp.m_vPeptides[0].m_vProteinAC[k]);

					for (size_t t = 0; t < SpectraTemp.m_vPeptides[0].m_vMod.size(); t++)
					{
						HaveToInterData.m_setAllModName.insert(
								SpectraTemp.m_vPeptides[0].m_vMod[t].m_strModName);
					}
				}
			}
		}
	}
#ifdef DEBUG
	cout << HaveToInterData.spec_peptide.size() << endl;
	cout << HaveToInterData.peptide_spec.size() << endl;
	cout << HaveToInterData.protein_peptide.size() << endl;
#endif
}

bool CDataIntegrate::_IsDeleteSpec(const vector<SPECTRAINFO> & vPeptide,
		IntegratedInfo & HaveToInterData)
{
	for (size_t t = 0; t < vPeptide.size(); t++)
	{
		CMatchSpectraInfo SpectraTemp;
		ReadSingeSpectraFisrtPep(vPeptide[t].first, SpectraTemp);

		for (size_t k = 0; k < SpectraTemp.m_vPeptides[0].m_vProteinAC.size(); k++)
		{
			if (HaveToInterData.m_deleteProtein.find(SpectraTemp.m_vPeptides[0].m_vProteinAC[k])
					== HaveToInterData.m_deleteProtein.end())
			{
				return false;
			}
		}
	}
	return true;
}

bool CDataIntegrate::_IsDeletePeptide(const vector<SPECTRAINFO> & vSpectra,
		IntegratedInfo & HaveToInterData)
{
	for (size_t t = 0; t < vSpectra.size(); t++)
	{
		if (HaveToInterData.m_deleteSpectra.find(vSpectra[t].second.m_strInPutPath)
				== HaveToInterData.m_deleteSpectra.end())
			return false;
	}
	return true;
}

void CDataIntegrate::_DeleteSet(const CConf & conf, IntegratedInfo & HaveToInterData)
{
	for (map<string, vector<SPECTRAINFO> >::const_iterator it =
			HaveToInterData.protein_peptide.begin(); it != HaveToInterData.protein_peptide.end(); it++)
	{
		set<string> setStrTemp;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			setStrTemp.insert(it->second[t].second.m_strFirstPeptide);
		}
		if (setStrTemp.size() < conf.m_Filter.m_nDistinctPepLimit)
			HaveToInterData.m_deleteProtein.insert(it->first);
	}
	for (map<string, vector<SPECTRAINFO> >::iterator it = HaveToInterData.spec_peptide.begin(); it
			!= HaveToInterData.spec_peptide.end(); it++)
	{
		if (_IsDeleteSpec(it->second, HaveToInterData))
		{
			HaveToInterData.m_deleteSpectra.insert(it->first);
		}
	}
	for (map<string, vector<SPECTRAINFO> >::iterator it = HaveToInterData.peptide_spec.begin(); it
			!= HaveToInterData.peptide_spec.end(); it++)
	{
		if (_IsDeletePeptide(it->second, HaveToInterData))
		{
			HaveToInterData.m_deletePeptide.insert(it->first);
		}
	}
}

void CDataIntegrate::_Initialize(IntegratedInfo & HaveToInterData)
{
	HaveToInterData.spec_peptide.clear();
	HaveToInterData.peptide_spec.clear();
	HaveToInterData.protein_peptide.clear();
	HaveToInterData.peptide_protein_set.clear();

	HaveToInterData.m_setAllModName.clear();
	HaveToInterData.m_setInvalidProteins.clear();
	HaveToInterData.m_mapSameSet.clear();
	HaveToInterData.m_mapSubSet.clear();

	HaveToInterData.m_deleteProtein.clear();
	HaveToInterData.m_deletePeptide.clear();
	HaveToInterData.m_deleteSpectra.clear();
}

void CDataIntegrate::DataTOIntegrate(const vector<OneDATASET> & ResultDataSet, const CConf & conf,
		IntegratedInfo & HaveToInterData)
{
#ifdef DEBUG
	cout << "Data Integrate..." << endl;
#endif
	osspBuildLog << "Data Integrate..." << endl;

	_Initialize(HaveToInterData);

	_CExportIndex(ResultDataSet, conf, HaveToInterData);

	if (conf.m_Filter.m_nDistinctPepLimit > 1)
		_DeleteSet(conf, HaveToInterData);

	_ProteinInfer(HaveToInterData);

	//TODO 注意这里DeleteSet和ProteinInferd的顺序，对结果没有影响，但是最终显示需要删除的没有在map中删除，最近
	//导致某些显示应该删除的蛋白可能出现在其他蛋白的sameset,subset中
}
