#ifndef EXPORTFORINTERFACE_H_
#define EXPORTFORINTERFACE_H_

#include "ResultExport.h"

#include "../BasicFunction/ReadWrite.h"
#include "../BasicFunction/CommonProcess.h"
#include "../BasicFunction/StringUtility.h"
#include <sstream>

#include <set>

namespace bio_analysis
{

class CExportForInterface
{
public:
	CExportForInterface();
	virtual ~CExportForInterface();
	void Export(const CConf & conf);
private:
	void
	_Statistics(const CConf & conf, const IntegratedInfo & HaveToInterData);
	void _Protein_java(const CConf & conf);
	void _Spectra_java(const CConf & conf);
	void _Peptide_java(const CConf & conf);
	void _GetProtein_Peptide(FILE * fout, const vector<CMatchSpectraInfo> & vSpectra,
			const int & order);

private:
	size_t nProtein_Group_Num;
};

}

#endif /* EXPORTFORINTERFACE_H_ */
