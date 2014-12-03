/*
 * ProteinInference.h
 *
 *  Created on: 2010-10-5
 *      Author: hfchen
 */

#ifndef PROTEININFERENCE_H_
#define PROTEININFERENCE_H_

#include "../BasicFunction/ReadWrite.h"
#include "../BasicFunction/CommonProcess.h"
#include "../BasicFunction/StringUtility.h"
#include "../bio_analysis_sdk.h"

using namespace bio_analysis;

class CProteinInference
{
public:
	CProteinInference();
	virtual ~CProteinInference();

	void Run(const IntegratedInfo & HaveToInterData, const CConf & conf);

private:
	void _Score_Normalization(
			const map<string, vector<SPECTRAINFO> > & HaveToInterData_peptide_spec);
	int _CountPeptideNum(const string & ProteinAC, const string & StrPeptide);
	void _Calculate_pij(const CConf & conf,
			const map<string, vector<SPECTRAINFO> > & HaveToInterData_protein_peptide);
	void _OutPut(const CConf & conf);
	void _prob_AND();
	void _prob_OR();
	void _prob_TFIDF_cosine();

private:
	double u;
	map<string, double> peptide_score;
	map<string, double> protein_score_prob_AND;
	map<string, double> protein_score_prob_OR;
	map<string, double> protein_score_prob_TFIDF_cosine;//
	map<string, vector<pair<string, int> > > protein_peptide_count;//在某种蛋白中，各肽段的数目
	map<string, int> peptide_count;//每种肽段在所有蛋白中的总数
	int protein_peptide_sum;//鉴定出肽段在蛋白中的总个数

};

#endif /* PROTEININFERENCE_H_ */
