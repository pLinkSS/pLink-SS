/*
 * ProteinInference.cpp
 *
 *  Created on: 2010-10-5
 *      Author: hfchen
 */

#include "ProteinInference.h"

CProteinInference::CProteinInference()
{
	// TODO Auto-generated constructor stub

}

CProteinInference::~CProteinInference()
{
	// TODO Auto-generated destructor stub
}

void CProteinInference::_Score_Normalization(
		const map<string, vector<SPECTRAINFO> > & HaveToInterData_peptide_spec)
{
	map<string, vector<SPECTRAINFO> > peptide_spec = HaveToInterData_peptide_spec;
	/*同一个肽段，可能来自不同的引擎，可能来自同一个引擎
	 *
	 */
	double Max[20] =
	{ -1.0e100, -1.0e100, -1.0e100, -1.0e100, -1.0e100, -1.0e100, -1.0e100, -1.0e100, -1.0e100,
			-1.0e100, -1.0e100, -1.0e100, -1.0e100, -1.0e100, -1.0e100 };
	double Min[20] =
	{ 1.0e100, 1.0e100, 1.0e100, 1.0e100, 1.0e100, 1.0e100, 1.0e100, 1.0e100, 1.0e100, 1.0e100,
			1.0e100, 1.0e100, 1.0e100, 1.0e100, 1.0e100 };

	map<string, vector<pair<double, int> > > Peptide_Score_Engine;//记录某个肽段在各个搜索引擎下的分数

	for (map<string, vector<SPECTRAINFO> >::iterator it = peptide_spec.begin(); it
			!= peptide_spec.end(); it++)
	{
		vector<pair<double, int> > vectorTemp;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			CMatchSpectraInfo SpectraTemp;
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);
			pair<double, int> pairTemp;
			pairTemp.first = SpectraTemp.m_vPeptides[0].m_vlfScores[0];
			pairTemp.second = SpectraTemp.m_nDataSetID;
			if (pairTemp.first > Max[pairTemp.second])
				Max[pairTemp.second] = pairTemp.first;
			if (pairTemp.first < Min[pairTemp.second])
				Min[pairTemp.second] = pairTemp.first;
			vectorTemp.push_back(pairTemp);
		}
		Peptide_Score_Engine[it->first] = vectorTemp;
	}

	double sum[20] =
	{ 0 };
	for (map<string, vector<pair<double, int> > >::iterator it = Peptide_Score_Engine.begin(); it
			!= Peptide_Score_Engine.end(); it++)
	{
		for (size_t t = 0; t < it->second.size(); t++)
		{
			double tmp = (it->second[t].first - Min[it->second[t].second])
					/ (Max[it->second[t].second] - Min[it->second[t].second]);
			sum[it->second[t].second] += tmp;
			it->second[t].first = tmp;
		}
	}

	for (map<string, vector<pair<double, int> > >::iterator it = Peptide_Score_Engine.begin(); it
			!= Peptide_Score_Engine.end(); it++)
	{
		for (size_t t = 0; t < it->second.size(); t++)
		{
			it->second[t].first = it->second[t].first / sum[it->second[t].second];
		}
	}

	for (map<string, vector<pair<double, int> > >::iterator it = Peptide_Score_Engine.begin(); it
			!= Peptide_Score_Engine.end(); it++)
	{
		double score = 1.0;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			score *= it->second[t].first;
		}
		score = pow(score, 1.0 / it->second.size());

		peptide_score[it->first] = score;

		/*这里的公式还得再讨论，暂时参考
		 *  Improving sensitivity in proteome studies by analysis of false discovery rates for multiple search engines
		 */
	}

	double max = -1.0e100, min = 1.0e100;
	for (map<string, double>::iterator it = peptide_score.begin(); it != peptide_score.end(); it++)
	{
		if (it->second > max)
			max = it->second;
		if (it->second < min)
			min = it->second;
	}

	ofstream fout1("peptide_score1.txt");
	double sum_2 = 0;
	for (map<string, double>::iterator it = peptide_score.begin(); it != peptide_score.end(); it++)
	{
		fout1 << it->first << " " << it->second << endl;
		if (it->second > 0.1)
			fout1 << endl;
		sum_2 += it->second;
	}
	fout1 << sum_2 << " " << sum_2 / peptide_score.size() << endl;
	fout1.close();

	double sum_score = 0.0;
	for (map<string, double>::iterator it = peptide_score.begin(); it != peptide_score.end(); it++)
	{
		it->second = (it->second - min) / (max - min);
		sum_score += it->second;
	}

	for (map<string, double>::iterator it = peptide_score.begin(); it != peptide_score.end(); it++)
	{
		it->second /= sum_score;
	}
	ofstream fout("peptide_score.txt");
	double sum_1 = 0;
	for (map<string, double>::iterator it = peptide_score.begin(); it != peptide_score.end(); it++)
	{
		fout << it->first << " " << it->second << endl;
		if (it->second > 0.1)
			fout << endl;
		sum_1 += it->second;
	}
	fout << sum_1 << " " << sum_1 / peptide_score.size() << endl;
	fout.close();
}

int CProteinInference::_CountPeptideNum(const string & ProteinAC, const string & StrPeptide)
{
	int num = 0;
	size_t pos = -1;
	while ((pos = ProteinAC.find(StrPeptide, pos + 1)) != string::npos)
		num++;
	return num;
}

void CProteinInference::_Calculate_pij(const CConf & conf,
		const map<string, vector<SPECTRAINFO> > & HaveToInterData_protein_peptide)
{
	map<string, vector<SPECTRAINFO> > protein_peptide = HaveToInterData_protein_peptide;

	for (size_t t = 0; t < conf.m_vProDBPath.size(); t++)
	{
		ifstream fin(conf.m_vProDBPath[t].c_str());
		if (!fin.good())
		{
			CErrInfo info("ProteinInference.cpp", "_prob_TFIDF_cosine", "Cannot open the file: "
					+ conf.m_vProDBPath[t]);
			throw runtime_error(info.Get());
		}

		string ProteinSQ = "";
		string ProteinAC = "";
		string ProteinDE = "";
		string strRead = "";

		GetLine(fin, strRead);

		size_t pos = strRead.find_first_of(' ');
		ProteinAC = strRead.substr(1, pos - 1);
		ProteinDE = strRead.substr(pos + 1);

		peptide_count.clear();
		protein_peptide_sum = 0;
		while (GetLine(fin, strRead))
		{
			if (strRead[0] == '>')
			{
				vector<pair<string, int> > vectorTemp;
				for (size_t pt = 0; pt < protein_peptide[ProteinAC].size(); pt++)
				{
					pair<string, int> pairTemp;
					pairTemp.second = _CountPeptideNum(ProteinSQ,
							protein_peptide[ProteinAC][pt].second.m_strFirstPeptide);
					pairTemp.first = protein_peptide[ProteinAC][pt].second.m_strFirstPeptide;
					peptide_count[pairTemp.first] += pairTemp.second;
					protein_peptide_sum += pairTemp.second;
					vectorTemp.push_back(pairTemp);
				}
				protein_peptide_count[ProteinAC] = vectorTemp;
				////////////////////////////
				ProteinSQ = "";
				pos = strRead.find_first_of(' ');
				ProteinAC = strRead.substr(1, pos - 1);
				ProteinDE = strRead.substr(pos + 1);
			}
			else
			{
				//ProteinSQ += strRead.substr(0, strRead.size());
				//cout << strRead << endl;
				for (size_t strReadt = 0; strReadt < strRead.size(); strReadt++)
				{
					if (strRead[strReadt] >= 'A' && strRead[strReadt] <= 'Z')
						ProteinSQ += strRead[strReadt];
				}
				//cout << ProteinSQ << endl;
			}
		}
		//The last protein
		vector<pair<string, int> > vectorTemp;
		for (size_t pt = 0; pt < protein_peptide[ProteinAC].size(); pt++)
		{
			pair<string, int> pairTemp;
			pairTemp.second = _CountPeptideNum(ProteinSQ,
					protein_peptide[ProteinAC][pt].second.m_strFirstPeptide);
			pairTemp .first = protein_peptide[ProteinAC][pt].second.m_strFirstPeptide;
			peptide_count[pairTemp.first] += pairTemp.second;
			protein_peptide_sum += pairTemp.second;
			vectorTemp.push_back(pairTemp);
		}
		protein_peptide_count[ProteinAC] = vectorTemp;

		fin.close();
	}
}

void CProteinInference::_prob_AND()
{
	vector<double> pi;
	for (map<string, int>::iterator it = peptide_count.begin(); it != peptide_count.end(); it++)
	{
		pi.push_back(peptide_count[it->first] / protein_peptide_sum);
	}

	for (map<string, vector<pair<string, int> > >::iterator it = protein_peptide_count.begin(); it
			!= protein_peptide_count.end(); it++)
	{
		double score = 1.0;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			double pij = (protein_peptide_count[it->first][t].second + u * pi[t])
					/ protein_peptide_sum;
			score *= pow(pij, peptide_score[it->second[t].first]);
		}
		protein_score_prob_AND[it->first] = score; //注意这里没有取负号
	}
}

void CProteinInference::_prob_OR()
{
	for (map<string, vector<pair<string, int> > >::iterator it = protein_peptide_count.begin(); it
			!= protein_peptide_count.end(); it++)
	{
		double score = 1.0;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			score *= (1 - peptide_score[it->second[t].first]);
		}
		protein_score_prob_OR[it->first] = 1 - score;
	}
}

void CProteinInference::_prob_TFIDF_cosine()
{
	for (map<string, vector<pair<string, int> > >::iterator it = protein_peptide_count.begin(); it
			!= protein_peptide_count.end(); it++)
	{
		double sum = 0.0;
		vector<double> dij;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			double dijTemp = (log(protein_peptide_count[it->first][t].second) + 1.0) * log(
					protein_peptide_sum / peptide_count[protein_peptide_count[it->first][t].first]);
			dij.push_back(dijTemp);
			sum += (dijTemp * dijTemp);
		}
		for (size_t t = 0; t < it->second.size(); t++)
		{
			dij[t] /= sum;
		}
		double sum1 = 0.0;
		double sum2 = 0.0;
		double sum3 = 0.0;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			sum1 += dij[t] * peptide_score[protein_peptide_count[it->first][t].first];
			sum2 += dij[t] * dij[t];
			sum3 += peptide_score[protein_peptide_count[it->first][t].first]
					* peptide_score[protein_peptide_count[it->first][t].first];
		}
		protein_score_prob_TFIDF_cosine[it->first] = sum1 / sqrt(sum2 * sum3);
	}
}

void CProteinInference::_OutPut(const CConf & conf)
{
	string ExportPathAndName = conf.m_OutPutForder + conf.m_OutPutFile;
	string proteininference_xls = ExportPathAndName + "_ProteinInference.xls";
	ofstream fout(proteininference_xls.c_str());
	fout << "ProteinAC\t";
	fout << "prob_AND\t";
	fout << "prob_OR\t";
	fout << "TFIDF_cosine" << endl;

	for (map<string, double>::iterator it = protein_score_prob_AND.begin(); it
			!= protein_score_prob_AND.end(); it++)
	{
		fout << it->first << "\t";
		fout << it->second << "\t";
		fout << protein_score_prob_OR[it->first] << "\t";
		fout << protein_score_prob_TFIDF_cosine[it->first] << endl;
	}
	fout.close();
}

void CProteinInference::Run(const IntegratedInfo & HaveToInterData, const CConf & conf)
{
	//cout << "CProteinInference Run" << endl;
	if (conf.m_vProDBPath.size() == 0)
		return;
	//cout << "CProteinInference Run2" << endl;
	u = 5000.0;
	_Score_Normalization(HaveToInterData.peptide_spec);
	_Calculate_pij(conf, HaveToInterData.protein_peptide);

	_prob_OR();
	_prob_AND();
	_prob_TFIDF_cosine();
	_OutPut(conf);
}
