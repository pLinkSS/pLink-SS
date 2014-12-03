#ifndef READWRITE_H_
#define READWRITE_H_

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/stat.h>
#include "CommonProcess.h"
#include "StringUtility.h"
#include "CalcBioVariable.h"
#include "../bio_analysis_sdk.h"

using namespace bio_analysis;

void pBuildLog(const string & strOut);
string GetTimePre(const CConf & conf);
void ReadCConf(CConf & conf, const string & m_strParamFilePath);
void WriteCConf(const CConf & conf, const string & strFilePath);

void RemoveFolderOrFile(string strPath);

void ReadCFiltration(CFiltration & m_Filter, FILE * fp);
void WriteCFiltration(const CFiltration & m_Filter, FILE * fp);

void SetDefault(CConf & conf);
void ReadSingeSpectra(const SPECTRA_PATH & pairPath, CMatchSpectraInfo & SpectraTemp);
void ReadSingeSpectraFisrtPep(const SPECTRA_PATH & pairPath, CMatchSpectraInfo & SpectraTemp);

int ReadSpectra(ifstream & fin, string & strBuf, CMatchSpectraInfo & SpectraTemp);

void WriteCConfHTML(FILE * fp, const CConf & conf);
void WriteFiltrationHTML(const CConf & conf, FILE * fp);

int ReadSpectra(ifstream & fin, CMatchSpectraInfo & SpectraTemp);

//void OutPutSpectra(FILE * fout, const CMatchSpectraInfo & SpectraTemp);
void OutPutSpectra(FILE * fout, const CMatchSpectraInfo & SpectraTemp, FILE * ferrorfigure_Da,
		FILE * ferrorfigure_ppm, const CConf & conf);
void OutPutSpectraFirstPep(FILE * fout, const CMatchSpectraInfo & SpectraTemp);
void
OutPutCondition(FILE * fp, const CConditionInfo & cond, const int & TotalSpectra);

void WriteIndexFile(const string & strIndexFile, const vector<OneDATASET> & AllDataSetParser);
void ReadIndexFile(const string & strIndexFile, vector<OneDATASET> & AllDataSetParser);
void WriteParserSimple(const string & strSimpleFile, const CConf & conf);

void CheckIndexFile(const string & strIndexFile);
bool CheckFileExist(const string & strFilePath);
bool CheckFolderExist(const string & strForderPath);

void WriteSingePeptide(FILE * fout, const CMatchPeptideInfo & PepTemp,
		const CMatchSpectraInfo & SpectraTemp);
void WritePeptide_Spec(FILE * fout, const CMatchSpectraInfo & SpectraTemp, const size_t & t);

void GetProtein_Peptide(FILE * fout, const CMatchSpectraInfo & SpectraTemp, const int & i);
//void WriteSingePeptide(FILE * fout, const CMatchPeptideInfo & PepTemp,
//		const CMatchSpectraInfo & SpectraTemp);

void GetProteinHead(FILE * fout, const string & strPro, const vector<SPECTRAINFO> & SpecPath);
//void WriteSingePeptideforCSV(FILE * fout, const CMatchPeptideInfo & PepTemp,
//		const CMatchSpectraInfo & SpectraTemp);

//void OutPutProteinDatabase(PROTEINIDATABASEINFO & ProDatabaseInfo, const CConf & conf);
//CMatchProteinInfo ReadOneProteinFromDatabase(const string & strPro,
//		PROTEINIDATABASEINFO & ProDatabaseInfo, const CConf & conf);
//void OutPutProteinDatabase(PROTEINIDATABASEINFO & ProDatabaseInfo, const CConf & conf);
//void ReadProteinDATABASE(const CConf & conf, PROTEINIDATABASEINFO & ProDatabaseInfo);
double CalcCoverage(const CMatchProteinInfo & ProteinTemp,
		map<string, vector<SPECTRAINFO> > & protein_peptide);

void ReadProteinANDCalculation(const CConf & conf, map<string, CProteinInfo> & ProteinInfo,
		const map<string, vector<SPECTRAINFO> > & protein_peptide);
bool CheckSimpleSame(const CConf & conf, const string & strParserSimple);

void PrintProteinDatabase(FILE * fout, const string & strPro);

void ReadFromIntermediateFromWriteSinglePeptide(istringstream & iss, string & Score,
		string & Calc_M, string & Delta_M, string & ppm, string & Mod_Sites, string & SampleID,
		string & Engine, string & MatchedIons, string & MissCleaveNum, string & Rank,
		string & Proteins);

void ReadFromIntermediateFileSpectra1(istringstream & iss, string & order, string & Spectrum,
		string & PeptideNum, string & UniquePepNum, string & Samples, string & Score,
		string & Condition);
void ReadFromIntermediateFileSpectra2(istringstream & iss, string & strOrder, string & Peptide,
		string & Score, string & Calc_M, string & Delta_M, string & ppm, string & Mod_Sites,
		string & SampleID, string & Engine, string & MatchedIons, string & MissCleaveNum,
		string & Rank, string & Proteins);

void ReadFromIntermediateFilePeptide1(istringstream & iss, string & order, string & SQ,
		string & Spectra, string & Samples, string & Score, string & Proteins);
void ReadFromIntermediateFilePeptide2(istringstream & iss, string & strOrder, string & Spectra,
		string & Score, string & Calc_M, string & Delta_M, string & ppm, string & Mod_Sites,
		string & SampleID, string & Engine, string & MatchedIons, string & MissCleaveNum,
		string & Rank, string & Proteins);

void ReadFromIntermediateFileProtein1(istringstream & iss, string & order, string & ProteinAC,
		string & MW, string & pI, string & Coverage, string & UniquePepNum, string & SpecNum,
		string & NonModifiedSpecNum, string & ModifiedSpecNum, string & UniqueModifiedPepNum,
		string& Samples, string & Description);
void ReadFromIntermediateFileProtein2(istringstream & iss, string & strOrder, string & Spectrum,
		string & Sequence, string & Score, string & Calc_M, string & Delta_M, string & ppm,
		string & Mod_Sites, string & SampleID, string & Engine, string & MatchedIons,
		string & MissCleaveNum, string & Rank, string & Proteins);

void OutPutSpectra_onlyforhfchen(FILE * fout, const CMatchSpectraInfo & SpectraTemp,
		const set<string> & SpectraScan);

#endif /* READPARAM_H_ */
