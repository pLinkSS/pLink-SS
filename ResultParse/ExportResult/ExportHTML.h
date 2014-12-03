#ifndef EXPORTHTML_H_
#define EXPORTHTML_H_

#include "ResultExport.h"

#include "../bio_analysis_sdk.h"
#include "../BasicFunction/ReadWrite.h"
#include "../BasicFunction/CommonProcess.h"
#include "../BasicFunction/StringUtility.h"

namespace bio_analysis
{

class CExportHTML: public CResultExport
{
public:
	CExportHTML();
	virtual ~CExportHTML();

	virtual void Export(const IntegratedInfo & HaveToInterData, const CConf & conf);
private:
	void _GetLink(FILE * fout, const CConf & conf);
	void _GetProteinHead(FILE * fout, const string & strPro, const vector<SPECTRAINFO> & SpecPath);
	void _Protein(const CConf & conf, const IntegratedInfo & HaveToInterData);
	void _Peptide(const CConf & conf, const IntegratedInfo & HaveToInterData);
	void _Spectra(const CConf & conf, const IntegratedInfo & HaveToInterData);
	void _Summary(const CConf & conf, const IntegratedInfo & HaveToInterData);
private:
	size_t nProtein_Group_Num;
};

}

#endif /* EXPORTHTML_H_ */
