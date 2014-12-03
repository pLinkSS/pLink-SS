///*
// * DataAnalysis.h
// *
// *  Created on: 2011-1-10
// *      Author: chenhaifeng
// */

#ifndef DATAANALYSIS_H_
#define DATAANALYSIS_H_

#include <map>
#include <set>
#include <bitset>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <string.h>
#include <fstream>
#include <sstream>

#include "ReadWrite.h"
#include "CommonProcess.h"
#include "StringUtility.h"

#include "../bio_analysis_sdk.h"

using namespace bio_analysis;

class CSpectrum_forhfchen
{
public:
	int order;
	double score, MH, Cal_H, Delt, ppm;
	string scan, sequence, YESorNOT, modification, instrument, truesequence;
	string TRUEORNOTTRUE;

	void Input(istringstream & iss)
	{
		iss >> order >> scan >> sequence >> score >> YESorNOT >> MH >> Cal_H >> Delt >> ppm
				>> modification >> instrument >> truesequence >> TRUEORNOTTRUE;
	}
	void Output(ofstream & fout) const
	{
		fout << order << " " << scan << " " << sequence << " " << score << " " << YESorNOT << " "
				<< MH << " " << Cal_H << " " << Delt << " " << ppm << " " << modification << " "
				<< instrument << " " << truesequence << endl;
		//		cout << order << " " << scan << " " << sequence << " " << score << " " << YESorNOT << " "
		//				<< MH << " " << Cal_H << " " << Delt << " " << ppm << " " << modification << " "
		//				<< instrument << " " << truesequence << endl;
	}
	bool operator==(const CSpectrum_forhfchen & b)
	{
		if (order != b.order || scan != b.scan || sequence != b.sequence
				|| /* score != b.score ||*/YESorNOT != b.YESorNOT || MH != b.MH || Cal_H != b.Cal_H
				|| Delt != b.Delt || ppm != b.ppm || modification != b.modification || instrument
				!= b.instrument)
			return false;
		return true;
	}
};

class CDataAnalysis
{
public:
	CDataAnalysis(const CConf & conf)
	{
		//strDataAnalysis = conf.m_OutPutName + "_DataAnalysis/";
		strDataAnalysis = conf.m_outPutForder_Java + "DataAnalysis_Result/";
		//cout << strDataAnalysis << endl;

#ifdef WIN32
		mkdir(strDataAnalysis.c_str());
#else			
		mkdir(strDataAnalysis.c_str(),0775);
#endif

		ifstream fin1("std8_trypsin_denovo_20090729_result.txt");

		string ExportResult = strDataAnalysis + "ExportResult.txt";
		char ss[10000];
		while (fin1.getline(ss, 1000, '\n'))
		{
			istringstream iss(ss, istringstream::in);
			CSpectrum_forhfchen spectraTemp;
			spectraTemp.Input(iss);
			//spectraTemp.Output(fout1);
			spectra1.push_back(spectraTemp);
			//cout << "order=" << spectraTemp.order << endl;
			//system("pause");
		}
		fin1.close();
		CompleteSpectraInfo.clear();
	}
	virtual ~CDataAnalysis()
	{
		CompleteSpectraInfo.clear();
	}

	vector<CMatchSpectraInfo> CompleteSpectraInfo;
	string strDataAnalysis;
	vector<CSpectrum_forhfchen> spectra1;
	vector<CSpectrum_forhfchen> spectra2;
	void _GetSameAsMunal();
	void DataAnalysis(const CConf & conf);
	void OutPutEachResult(IntegratedInfo & HaveToInterData, const CConf & conf);
	void getwhy(IntegratedInfo & HaveToInterData, const CConf & conf,const int & xx);
	void getVenn();
	void _getHCDETDsame();
private:
	void _OutPutOneSpectra(FILE * fout, const CMatchSpectraInfo & SpectraTemp,
			const set<string> & SpectraScan);

	void _ReadProteinANDGetSQ(const CConf & conf, map<string, string> & ProteinACSQ);
	bool _IsInDatabase(const string & sequence, map<string, string> & ProteinACSQ);

	void _get_true_ans(const CConf & conf, const vector<CSpectrum_forhfchen> & spectra1,
			const vector<CSpectrum_forhfchen> & spectra2,
			map<string, vector<pair<string, double> > > devovo_res);
	void _get_true_ans2(const CConf & conf, const vector<CSpectrum_forhfchen> & spectra1,
			const vector<CSpectrum_forhfchen> & spectra2,
			map<string, vector<pair<string, double> > > devovo_res);
	void _get_another_ans1(const CConf & conf, const vector<CSpectrum_forhfchen> & spectra1,
			const vector<CSpectrum_forhfchen> & spectra2,
			map<string, vector<pair<string, double> > > devovo_res);
	void _Out_three_result(ofstream & fout3, const vector<CSpectrum_forhfchen> & spectra1,
			const vector<CSpectrum_forhfchen> & spectra2,
			map<string, vector<pair<string, double> > > devovo_res, const int & i);
	void _getinter(const vector<int> & index, int inter[], int intertrue[]);

};

#endif /* DATAANALYSIS_H_ */

