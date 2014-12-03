/*
 * DataAnalysis.cpp
 *
 *  Created on: 2011-1-10
 *      Author: chenhaifeng
 */

#include "DataAnalysis.h"

using namespace bio_analysis;

extern ostringstream osspBuildLog;
//extern vector<CMatchSpectraInfo> SpectraOnlyforhfchen; //todo only for hfchen这里

//set<int> already;


void CDataAnalysis::_OutPutOneSpectra(FILE * fout, const CMatchSpectraInfo & SpectraTemp,
		const set<string> & SpectraScan)
{
	//	string type_hfchen[10] =
	//	{ "HCD-FTMS-FROM50", "HCD-FTMS-FROM100", "ETD-FTMS", "ETD-ITM", "CID-ITMS" };
	string type_hfchen[10] =
	{ "HCD-Merged", "ETD-FTMS" };

	fprintf(fout, "%s\t", SpectraTemp.m_strFileName.c_str());

	//	if (SpectraTemp.m_vPeptides.size() > 0)
	//	{
	//		fprintf(fout, "%s\t", SpectraTemp.m_vPeptides[0].m_strSQ.c_str());
	//		if (SpectraTemp.EnginType == ST_PFIND)
	//			fprintf(fout, "%2.2e\t", EValueToSmall(SpectraTemp.m_vPeptides[0].m_vlfScores[0]));
	//		else
	//			fprintf(fout, "%.2lf\t", SpectraTemp.m_vPeptides[0].m_vlfScores[0]);
	//	}
	//	else
	//		fprintf(fout, "****\t0.000\t");


	if (SpectraTemp.m_vPeptides.size() > 0)
	{
		//		for (size_t t = 0; t < SpectraTemp.m_vPeptides.size(); t++)
		//		{
		//			fprintf(fout, "%s\t", SpectraTemp.m_vPeptides[t].m_strSQ.c_str());
		//			if (SpectraTemp.EnginType == ST_PFIND)
		//				fprintf(fout, "%2.2e\t", EValueToSmall(SpectraTemp.m_vPeptides[t].m_vlfScores[0]));
		//			else
		//				fprintf(fout, "%.2lf\t", SpectraTemp.m_vPeptides[t].m_vlfScores[0]);
		//		}

		fprintf(fout, "%s\t", SpectraTemp.m_vPeptides[0].m_strSQ.c_str());
		if (SpectraTemp.EnginType == ST_PFIND)
			fprintf(fout, "%2.2e\t", EValueToSmall(SpectraTemp.m_vPeptides[0].m_vlfScores[0]));
		else
			fprintf(fout, "%.2lf\t", SpectraTemp.m_vPeptides[0].m_vlfScores[0]);

	}
	else
		fprintf(fout, "****\t0.000\t");

	if (SpectraScan.find(SpectraTemp.m_strFileName) != SpectraScan.end())
	{
		fprintf(fout, "YES\t");
	}
	else
	{
		fprintf(fout, "NO\t");
		//		for (size_t t = 0; t < SpectraTemp.m_vPeptides.size(); t++)
		//		{
		//			if (SpectraTemp.m_vPeptides[t].m_strSQ == )
		//				fprintf(fout, "%s\t", SpectraTemp.m_vPeptides[t].m_strSQ.c_str());
		//			if (SpectraTemp.EnginType == ST_PFIND)
		//				fprintf(fout, "%2.2e\t", EValueToSmall(SpectraTemp.m_vPeptides[t].m_vlfScores[0]));
		//			else
		//				fprintf(fout, "%.2lf\t", SpectraTemp.m_vPeptides[t].m_vlfScores[0]);
		//		}
	}

	fprintf(fout, "%lf\t", SpectraTemp.m_lfMH);

	if (SpectraTemp.m_vPeptides.size() == 0)
		fprintf(fout, "0.000\t0.000\t0.000\tNULL\t");
	else if (SpectraTemp.m_vPeptides.size() > 0)
	{
		fprintf(fout, "%lf\t", SpectraTemp.m_vPeptides[0].m_lfCalc_MH - protonH);
		fprintf(fout, "%lf\t", SpectraTemp.m_vPeptides[0].m_lfDelta);
		fprintf(fout, "%lf\t", SpectraTemp.m_vPeptides[0].m_lfPPM);

		if (SpectraTemp.m_vPeptides[0].m_vMod.size() == 0)
			fprintf(fout, "NULL");
		CMatchPeptideInfo PepTemp = SpectraTemp.m_vPeptides[0];
		for (size_t t = 0; t < SpectraTemp.m_vPeptides[0].m_vMod.size(); t++)
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
				fprintf(fout, ""PRI_SIZE_T",%c(%s);", PepTemp.m_vMod[t].m_tPos,
						PepTemp.m_strSQ[PepTemp.m_vMod[t].m_tPos - 1],
						PepTemp.m_vMod[t].m_strModName.c_str());
			}
		}
		fprintf(fout, "\t");
	}

	fprintf(fout, "trypsin_%s\t", type_hfchen[SpectraTemp.m_nDataSetID - 1].c_str());
	fprintf(fout, "NULL\n");
}

void CDataAnalysis::OutPutEachResult(IntegratedInfo & HaveToInterData, const CConf & conf)
{

	int fragment = 2;
	string ExportResult = strDataAnalysis + "ExportResult.txt";
	cout << ExportResult << endl;
	FILE * fout = fopen(ExportResult.c_str(), "w");

	if (fout == NULL)
	{
		printf("NO\n");
	}
	//fprintf(fout_hfchen, "SCAN号\t序列\t分数\t1%%FDR是否鉴定到\tMH\tCalc_M\tDelta_M\tppm\t修饰\t仪器\n");
	set<string> SpectraScan;
	for (map<string, vector<SPECTRAINFO> >::iterator it = HaveToInterData.spec_peptide.begin(); it
			!= HaveToInterData.spec_peptide.end(); it++)

	{
		SpectraScan.insert(it->first);
	}

	sort(CompleteSpectraInfo.begin(), CompleteSpectraInfo.end(), CMatchSpectraInfo_SortSCANPATH);

	int cnt = 1;
	for (size_t t = 0; t < CompleteSpectraInfo.size(); t++)
	{
		fprintf(fout, ""PRI_SIZE_T"\t", (t + fragment) / fragment);
		_OutPutOneSpectra(fout, CompleteSpectraInfo[t], SpectraScan);
		cnt++;
	}
	fclose(fout);
}

void CDataAnalysis::_ReadProteinANDGetSQ(const CConf & conf, map<string, string> & ProteinACSQ)
{
	osspBuildLog << "Read Protein Database..." << endl;

	for (size_t t = 0; t < conf.m_vProDBPath.size(); t++)
	{
		ifstream fin(conf.m_vProDBPath[t].c_str());
		if (!fin.good())
		{
			CErrInfo info("ReadWrite.cpp", "ReadProteinANDGetSQ",
					"Cannot open the file: " + conf.m_vProDBPath[t]);
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

		while (GetLine(fin, strRead))
		{
			if (strRead[0] == '>')
			{
				ProteinACSQ[ProteinAC] = ProteinSQ;

				ProteinSQ = "";
				pos = strRead.find_first_of(' ');
				ProteinAC = strRead.substr(1, pos - 1);
				ProteinDE = strRead.substr(pos + 1);
			}
			else
			{
				for (size_t strReadt = 0; strReadt < strRead.size(); strReadt++)
				{
					if (strRead[strReadt] >= 'A' && strRead[strReadt] <= 'Z')
						ProteinSQ += strRead[strReadt];
				}
			}
		}
		ProteinACSQ[ProteinAC] = ProteinSQ;
		fin.close();
	}
}

bool CDataAnalysis::_IsInDatabase(const string & sequence, map<string, string> & ProteinACSQ)
{
	for (map<string, string>::iterator it = ProteinACSQ.begin(); it != ProteinACSQ.end(); it++)
	{
		if (it->second.find(sequence) != string::npos)
			return true;
	}
	return false;
}

void CDataAnalysis::_Out_three_result(ofstream & fout3,
		const vector<CSpectrum_forhfchen> & spectra1, const vector<CSpectrum_forhfchen> & spectra2,
		map<string, vector<pair<string, double> > > devovo_res, const int & i)
{
	//already.insert(spectra1[i].order);

	spectra1[i].Output(fout3);
	spectra2[i].Output(fout3);

	for (size_t t = 0; t < devovo_res[spectra1[i].scan].size(); t++)
	{
		fout3 << spectra1[i].scan << " " << devovo_res[spectra1[i].scan][t].first << " "
				<< devovo_res[spectra1[i].scan][t].second << endl;
	}
	fout3 << endl;
}

void CDataAnalysis::_get_another_ans1(const CConf & conf,
		const vector<CSpectrum_forhfchen> & spectra1, const vector<CSpectrum_forhfchen> & spectra2,
		map<string, vector<pair<string, double> > > devovo_res)
{
	string str_fout3 = strDataAnalysis + "get_another_result1.txt";
	cout << str_fout3 << endl;
	ofstream fout3(str_fout3.c_str());

	cout << "get_another_result.txt" << endl;
	for (size_t i = 0; i < spectra1.size(); i++)
	{
		//if (already.find(spectra1[i].order) != already.end())
		//	continue;
		if (spectra1[i].YESorNOT == "YES" && (spectra1[i].score < 0.00001 && spectra1[i].score > 0))
			continue;

		int flag = 0;
		if (spectra2[i].YESorNOT == "YES" || (spectra2[i].score < 0.0001 && spectra2[i].score > 0))
			flag = 1;

		if (devovo_res[spectra1[i].scan].size() == 0)
			continue;
		if (devovo_res[spectra1[i].scan].size() >= 1 && devovo_res[spectra1[i].scan][0].second
				> 0.6)
			flag = 1;
		if (flag == 0)
			continue;

		_Out_three_result(fout3, spectra1, spectra2, devovo_res, i);
	}

	fout3.close();
}

void CDataAnalysis::_get_true_ans2(const CConf & conf,
		const vector<CSpectrum_forhfchen> & spectra1, const vector<CSpectrum_forhfchen> & spectra2,
		map<string, vector<pair<string, double> > > devovo_res)
{
	/*三个有一个正确*/
	string str_fout3 = strDataAnalysis + "get_true_ans2.txt";
	ofstream fout3(str_fout3.c_str());

	cout << "get_true_ans2.txt" << endl;
	for (size_t i = 0; i < spectra1.size(); i++)
	{
		//if (already.find(spectra1[i].order) != already.end())
		//continue;
		int flag = 0;
		if (spectra1[i].YESorNOT == "YES" && (spectra1[i].score < 0.00001 && spectra1[i].score > 0))
			flag = 1;

		if (spectra2[i].YESorNOT == "YES" && (spectra2[i].score < 0.00001 && spectra2[i].score > 0))
			flag = 1;

		if (devovo_res[spectra1[i].scan].size() >= 1 && devovo_res[spectra1[i].scan][0].second
				> 0.6)
			flag = 1;

		if (flag == 0)
			continue;

		_Out_three_result(fout3, spectra1, spectra2, devovo_res, i);

	}
	fout3.close();
}

void CDataAnalysis::_get_true_ans(const CConf & conf, const vector<CSpectrum_forhfchen> & spectra1,
		const vector<CSpectrum_forhfchen> & spectra2,
		map<string, vector<pair<string, double> > > devovo_res)
{
	/*三个都要正确*/
	string str_fout3 = strDataAnalysis + "get_true_result.txt";
	ofstream fout3(str_fout3.c_str());

	cout << "get_true_result.txt" << endl;
	for (size_t i = 0; i < spectra1.size(); i++)
	{
		//if (already.find(spectra1[i].order) != already.end())
		//continue;
		if (spectra1[i].TRUEORNOTTRUE == "TRUE" || spectra1[i].TRUEORNOTTRUE == "other")
			continue;
		if (spectra1[i].YESorNOT == "NO" || spectra1[i].score > 0.0001)
			continue;

		if (spectra2[i].YESorNOT == "NO" || spectra2[i].score > 0.0001)
			continue;

		if (devovo_res[spectra1[i].scan].size() == 0)
			continue;
		if (devovo_res[spectra1[i].scan][0].second < 0.5)
			continue;

		_Out_three_result(fout3, spectra1, spectra2, devovo_res, i);
	}
	fout3.close();
	cout << "hahlkj" << endl;
}

void CDataAnalysis::getwhy(IntegratedInfo & HaveToInterData, const CConf & conf, const int & xx)
{
	string type_hfchen[10] =
	{ "trypsin_HCD-FTMS-FROM50", "trypsin_HCD-FTMS-FROM100", "trypsin_ETD-FTMS", "trypsin_ETD-ITM",
			"trypsin_CID-ITMS" };
	set<string> SpectraScan;
	for (map<string, vector<SPECTRAINFO> >::iterator it = HaveToInterData.spec_peptide.begin(); it
			!= HaveToInterData.spec_peptide.end(); it++)

	{
		SpectraScan.insert(it->first);
	}

	sort(CompleteSpectraInfo.begin(), CompleteSpectraInfo.end(), CMatchSpectraInfo_SortSCANPATH);

	int numberoneyestrue = 0;
	int numberoneyes = 0;
	int numberoneno = 0;
	int number33 = 0;
	int number33not = 0;
	int number3 = 0;
	int notin = 0;
	int cnt = 0;
	int cnt2 = 0;
	ofstream see("see.txt");
	for (size_t i = 0; i < spectra1.size(); i++)
	{
		//cout << spectra1[i].TRUEORNOTTRUE << endl;
		if (spectra1[i].TRUEORNOTTRUE == "FALSE")
			continue;
		cnt2++;
		cout << spectra1[i].instrument << endl;
		if (spectra1[i].instrument != type_hfchen[xx])
		{
			//see << spectra1[i].scan << endl;
			//see << spectra1[i].instrument << endl;
			continue;
		}
		cnt++;
		for (size_t t = 0; t < CompleteSpectraInfo.size(); t++)
		{
			//cout << "t = " << t << endl;
			if (CompleteSpectraInfo[t].m_strFileName != spectra1[i].scan)
				continue;

			//cout << "NOW" << endl;
			const CMatchSpectraInfo * SpectraTemp = &CompleteSpectraInfo[t];
			if (SpectraScan.find(SpectraTemp->m_strFileName) != SpectraScan.end())
			{

				if (SpectraTemp->m_vPeptides[0].m_strSQ == spectra1[i].truesequence)
					numberoneyestrue++;
				else
				{
					numberoneyes++;
					int exist = false;
					for (size_t t = 0; t < SpectraTemp->m_vPeptides.size(); t++)
					{
						if (SpectraTemp->m_vPeptides[t].m_strSQ == spectra1[i].truesequence)
						{
							number33++;
							exist = true;
							break;
						}
					}
					if (exist == false)
						number33not++;
				}
				break;
			}
			else
			{
				if (SpectraTemp->m_vPeptides.size() == 0)
				{
					notin++;
					break;
				}
				if (SpectraTemp->m_vPeptides[0].m_strSQ == spectra1[i].truesequence)
				{
					numberoneno++;
					break;
				}
				else
				{
					bool exist = false;
					for (size_t t = 0; t < SpectraTemp->m_vPeptides.size(); t++)
					{

						if (SpectraTemp->m_vPeptides[t].m_strSQ == spectra1[i].truesequence)
						{
							number3++;
							exist = true;
							break;
						}
					}
					if (exist == false)
					{
						notin++;
					}
					break;
				}
				break;
			}
		}
	}
	char chr[100];
	sprintf(chr, "%s_OStat_%d.txt", strDataAnalysis.c_str(), xx);
	ofstream OStat(chr);
	OStat << numberoneyestrue << endl;
	OStat << numberoneyes << endl;
	OStat << numberoneno << endl;
	OStat << number3 << endl;
	OStat << notin << endl;
	OStat << cnt << endl;
	OStat << cnt2 << endl;
	OStat << number33 << endl;
	OStat << number33not << endl;
	OStat.close();
	see.close();
}

void CDataAnalysis::_getinter(const vector<int> & index, int inter[], int intertrue[])
{
	for (int i = 31; i >= 1; i--)
	{
		bitset<5> bit(i);
		//cout << bit << endl;
		set<string> zero;
		set<string> one;

		for (size_t j = 0; j < index.size(); j++)
		{
			if (bit[j] == 1)
				one.insert(spectra1[index[j]].sequence);
			else
				zero.insert(spectra1[index[j]].sequence);
		}

		if (one.size() == 1)
		{
			set<string>::iterator it = one.begin();
			if (*it == "****")
				continue;
			bool tag = true;
			for (set<string>::iterator it2 = zero.begin(); it2 != zero.end(); it2++)
			{
				if (*it == *it2)
					tag = false;
			}
			if (tag == true)
			{
				if (i == 31)
				{
					cout << *it << " " << spectra1[index[0]].truesequence << endl;
					if (*it != spectra1[index[0]].truesequence)
					{
						cout << " klsjdflkfjkl  = " << *it << " "
								<< spectra1[index[0]].truesequence << endl;
						cout << spectra1[index[0]].order << endl;
					}
				}

				inter[i]++;

				//cout << *it << " " << spectra1[index[0]].truesequence << endl;
				if (*it == spectra1[index[0]].truesequence)
					intertrue[i]++;
			}
		}
	}
	//system("pause");
}

//void CDataAnalysis::_getinter(const vector<int> & index, int inter[], int intertrue[])
//{
//	if(spectra1[index[0]] == spectra1[index])
//}

struct stat1
{
	string str;
	int a;
	int b;
};

int cmp(const struct stat1 & a, const struct stat1 & b)
{
	if (a.str > b.str)
		return 1;
	return -1;

}
void CDataAnalysis::getVenn()
{
	if (spectra1.size() == 0)
		return;

	int inter[64] =
	{ 0 };
	int intertrue[64] =
	{ 0 };
	int order = spectra1[0].order;
	vector<int> index;
	index.push_back(0);
	int cnt = 0;
	for (size_t t = 0; t < spectra1.size(); t++)
	{
		//cout << "order=" << order << endl;
		//cout << spectra1[t].order << endl;
		//system("pause");
		//		if (spectra1[t].sequence == "****")
		//			continue;
		if (spectra1[t].order == order)
			index.push_back(t);
		else
		{
			if (spectra1[index[0]].TRUEORNOTTRUE == "FALSE")
			{
				index.clear();
				index.push_back(t);
				order = spectra1[t].order;
				continue;
			}
			_getinter(index, inter, intertrue);
			cout << "cnt = " << ++cnt << endl;
			index.clear();
			index.push_back(t);
			order = spectra1[t].order;
		}
	}
	if (spectra1[index[0]].TRUEORNOTTRUE != "FALSE")
	{
		_getinter(index, inter, intertrue);
	}
	ofstream fout(string(strDataAnalysis + "venn.xls").c_str());

	//int sum1 = 0, sum2 = 0;
	vector<stat1> m_stat;
	for (int i = 31; i >= 1; i--)
	{
		bitset<5> bit(i);
		stat1 tmp;
		string strtmp;
		for (int j = 0; j < 5; j++)
		{
			if (bit[j] == 1)
			{
				char chr = 'A' + j;
				strtmp += chr;
			}
		}
		tmp.str = strtmp;
		tmp.a = inter[i];
		tmp.b = intertrue[i];
		m_stat.push_back(tmp);
	}
	sort(m_stat.begin(), m_stat.end(), cmp);
	for (size_t t = 0; t < m_stat.size(); t++)
	{
		fout << m_stat[t].str << "\t" << m_stat[t].a << "\t" << m_stat[t].b << endl;
	}
	//fout << "\t";
	//fout << inter[i] << "\t" << intertrue[i] << endl;
	//sum1 += inter[i];
	////sum2 += intertrue[i];
	//fout << sum1 << endl;
	//fout << sum2 << endl;
	fout.close();
}

void CDataAnalysis::_GetSameAsMunal()
{
	int nsame = 0;
	for (size_t t = 0; t < spectra1.size(); t++)
	{
		for (size_t tj = 0; tj < spectra2.size(); tj++)
		{
			if (spectra1[t].scan == spectra2[tj].scan)
			{
				cout << "spectra1[t] = " << spectra1[t].scan << endl;
				cout << spectra1[t].truesequence << endl;
				cout << spectra1[t].sequence << endl;
				cout << spectra1[t].TRUEORNOTTRUE << endl;
				cout << spectra1[t].instrument << endl;

				cout << spectra2[tj].sequence << endl;
				if (spectra1[t].truesequence == spectra2[tj].sequence)
				{
					nsame++;
				}
			}
		}
	}
	cout << "nsame=" << nsame << endl;
}
void CDataAnalysis::_getHCDETDsame()
{
	string SameResult = strDataAnalysis + "same.txt";
	ofstream fout(SameResult.c_str());
	for (size_t t = 0; t < spectra2.size(); t++)
	{
		if (spectra2[t].sequence != "****" && spectra2[t].sequence == spectra2[t + 1].sequence)
		{
			fout << spectra2[t].order << " " << spectra2[t].scan << " " << spectra2[t].sequence
					<< " " << spectra2[t].YESorNOT << " " << spectra2[t].score << " "
					<< spectra2[t].instrument << endl;
			fout << spectra2[t + 1].order << " " << spectra2[t + 1].scan << " "
					<< spectra2[t + 1].sequence << " " << spectra2[t + 1].YESorNOT << " "
					<< spectra2[t + 1].score << " " << spectra2[t + 1].instrument << endl;
			fout << endl;
			t++;
		}
	}
	fout.close();
	cout << "over" << endl;

}
void CDataAnalysis::DataAnalysis(const CConf & conf)
{
	osspBuildLog << "_DataAnalysis_forhfchen" << endl;
	cout << "_DataAnalysis_forhfchen" << endl;
	map<string, string> ProteinACSQ;

	vector<string> dev_res_path;
	//dev_res_path.push_back("param_cid_res_1.txt");
	//dev_res_path.push_back("param_hcd_res_1.txt");
	//dev_res_path.push_back("param_hcd_res_2.txt");

	//	ifstream fin1("std8_trypsin_denovo_20090729_result.txt");
	//
	//string ExportResult = strDataAnalysis + "ExportResult.txt";
	//ifstream fin2(ExportResult .c_str());
	//ifstream fin3("param_res_1_chihao.txt");

	//ofstream fout1("out1.txt");
	//ofstream fout2("out2.txt");

	//	vector<CSpectrum_forhfchen> spectra1;
	//	vector<CSpectrum_forhfchen> spectra2;
	map<string, vector<pair<string, double> > > devovo_res;
	map<string, bool> InDatabase;
	char ss[10000];
//	while (fin1.getline(ss, 1000, '\n'))
//	{
//		istringstream iss(ss, istringstream::in);
//		CSpectrum_forhfchen spectraTemp;
//		spectraTemp.Input(iss);
//		//spectraTemp.Output(fout1);
//		spectra1.push_back(spectraTemp);
//	}
	//_getwhy();
	return;
//	while (fin2.getline(ss, 1000, '\n'))
//	{
//		istringstream iss(ss, istringstream::in);
//		CSpectrum_forhfchen spectraTemp;
//		spectraTemp.Input(iss);
//		//spectraTemp.Output(fout2);
//		spectra2.push_back(spectraTemp);
//	}

	//_ReadProteinANDGetSQ(conf, ProteinACSQ);

	_GetSameAsMunal();
	_getHCDETDsame();
	return;
	for (size_t dev_res_path_t = 0; dev_res_path_t < dev_res_path.size(); dev_res_path_t++)
	{
		ifstream fin3(dev_res_path[dev_res_path_t].c_str());
		while (fin3.getline(ss, 1000, '\n'))
		{
			istringstream iss(ss, istringstream::in);
			string scan;
			double MH, spscore;
			int total, charge;
			vector<pair<string, double> > seq;
			iss >> scan >> MH >> charge >> spscore >> total;
			for (int t = 0; t < total; t++)
			{
				int order;
				string sequence;
				double score;
				fin3.getline(ss, 1000, '\n');
				istringstream iss2(ss, istringstream::in);
				iss2 >> order >> sequence >> score;
				pair<string, double> pairtemp;
				pairtemp.first = sequence;
				//InDatabase[sequence] = _IsInDatabase(sequence, ProteinACSQ);
				//cout << sequence << endl;
				//if (InDatabase[sequence] == true)
				//{
				//	fout3 << scan << " " << sequence << endl;
				//}
				pairtemp.second = score;
				seq.push_back(pairtemp);
			}
			devovo_res[scan] = seq;
		}
		fin3.close();
	}

	_get_true_ans(conf, spectra1, spectra2, devovo_res);
	_get_true_ans2(conf, spectra1, spectra2, devovo_res);
	_get_another_ans1(conf, spectra1, spectra2, devovo_res);

	//fin1.close();
	//fin2.close();
	//fout1.close();
	//fout2.close();
	//system("pause");
}
