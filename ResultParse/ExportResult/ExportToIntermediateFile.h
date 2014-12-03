/*
 * ExportToIntermediateFile.h
 *
 *  Created on: 2010-10-25
 *      Author: chenhaifeng
 */

#ifndef EXPORTTOINTERMEDIATEFILE_H_
#define EXPORTTOINTERMEDIATEFILE_H_

#include "../BasicFunction/ReadWrite.h"
#include "../BasicFunction/CommonProcess.h"
#include "../BasicFunction/StringUtility.h"
#include "../bio_analysis_sdk.h"
#include <bitset>

using namespace bio_analysis;
namespace bio_analysis
{

class CExportToIntermediateFile
{
public:
	CExportToIntermediateFile();
	virtual ~CExportToIntermediateFile();

	void ExportToIntermediateFile(const IntegratedInfo & HaveToInterData, const CConf & conf);

private:
	void _WriteSinglePeptide(FILE * fout, const CMatchPeptideInfo & PepTemp,
			const CMatchSpectraInfo & SpectraTemp);
	void _Spectra_IntermediateFile(const CConf & conf, const IntegratedInfo & HaveToInterData,
			FILE * fSS);
	void _Peptide_IntermediateFile(const CConf & conf, const IntegratedInfo & HaveToInterData,
			FILE * fSS);
	void _Protein_IntermediateFile(const CConf & conf, const IntegratedInfo & HaveToInterData,
			FILE * fSS, map<string, CProteinInfo> & ProteinInfo);
	void _plabel_IntermediateFile(const CConf & conf, const IntegratedInfo & HaveToInterData);
	void _Statistics(const CConf & conf);
	void _ProteinInferenceProblem(const CConf & conf, const IntegratedInfo & HaveToInterData, map<
			string, CProteinInfo> & ProteinInfo);
private:
	int nSpectra;
	int nPeptide;
	int nProtein;
	int nProGroup;

};

}

struct structUniPepofProtein
{
	string strPro;
	int UniPep;
};

#endif /* EXPORTTOINTERMEDIATEFILE_H_ */
