#ifndef COMMOMPROCESS_H_
#define COMMOMPROCESS_H_

#include "ReadWrite.h"

#include "../bio_analysis_sdk.h"

using namespace bio_analysis;

int IsReverse(const CMatchSpectraInfo & mr, const CConf & conf);
bool IsProteinReverse(const string & protein, const CConf & conf);

bool
Mod_Sort(const CModificationSiteInfo & s1, const CModificationSiteInfo &s2);
bool SPECTRAINFO_SortScore(const SPECTRAINFO & a, const SPECTRAINFO & b);
bool SPECTRAINFO_SortInPutPath(const SPECTRAINFO & a, const SPECTRAINFO & b);
bool Peptide_SortScore(const CMatchPeptideInfo & mp1, const CMatchPeptideInfo & mp2);
bool ModificationSiteEqual(const CModificationSiteInfo & s1, const CModificationSiteInfo & s2);

bool CMatchSpectraInfo_SortScore(const CMatchSpectraInfo & a, const CMatchSpectraInfo & b);
bool CMatchSpectraInfo_SortSCANPATH(const CMatchSpectraInfo & a, const CMatchSpectraInfo & b);

bool VecotrModEqual(const vector<CModificationSiteInfo> & m1,
		const vector<CModificationSiteInfo> & m2);
bool PeptideEqual(const PEPTIDEINFO & PepA, const PEPTIDEINFO & PepB);
bool PEP_SORT_SQ(const PEPTIDEINFO & pairPep1, const PEPTIDEINFO & pairPep2);

bool IsSameMod(const pair<string, CModificationSiteInfo> & m1, const pair<string,
		CModificationSiteInfo> & m2);

void STYCount(const pair<string, CModificationSiteInfo> & PairPepTemp, int STY[],
		const CConf & conf);
int GetNumOfMod(FILE * fout, const vector<PEPTIDEINFO> & vPairPep, int STY[], const CConf & conf);
void P123Count(const PEPTIDEINFO & PairPepTemp, int MODNUM[], const CConf & conf);
int GetDifferNumOfPep(const vector<PEPTIDEINFO> & vPairPep, int MODNUM[], const CConf & conf);

void GetIndex(map<int, vector<SPECTRAINFO> > & DataSetSpectra,
		const CMatchSpectraInfo & SpectraTemp, const SPECTRA_PATH & SpecPathTemp);

SearchEngineType GetEnginType(const string & strTemp);

double EValueToLager(const double & EValue);
double EValueToSmall(const double & EValue);

bool InChargeState(const int & ChargeState, const CConf & conf);

string GetPeptideWithMod(const CMatchPeptideInfo & mp);

bool ISinVector(const int & a, vector<int> m_vINT);
int GetCrossType(const string & strSQ);
int GetReagentType(const string & strSQ);
int getEXPCode(const int & a, const int & b);
int GetCrossLinkType(const int & charge, const int & nCrossType, const int & nReagentType);
int DeCrossLinkType(const int & charge);

void Test(const vector<OneDATASET> & ResultDataSet);

#endif /* SORTFUNCTION_H_ */
