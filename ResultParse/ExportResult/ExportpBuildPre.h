#ifndef EXPORTPBUILDPRE_H_
#define EXPORTPBUILDPRE_H_

#include "ResultExport.h"

#include "../BasicFunction/ReadWrite.h"
#include "../BasicFunction/CommonProcess.h"
#include "../BasicFunction/StringUtility.h"

#include <set>

namespace bio_analysis
{

class CExportpBuildPre
{
public:
	CExportpBuildPre();
	virtual ~CExportpBuildPre();
	void Export(const CConf & conf, const string & ExportPath, const string & ExportFileType);

private:
	void _Protein_Export_NoUniPepSort(const CConf & conf, const string & ExportPathAndName);
	void _Spectra_Export(const CConf & conf, const string & ExportPathAndName);
	void _Peptide_Export(const CConf & conf, const string & ExportPathAndName);
	void _Protein_Export(const CConf & conf, const string & ExportPathAndName);
	void _pQuant_Export(const CConf & conf, const string & ExportPathAndName);
	void _plabel_Export(const CConf & conf, const string & ExportPathAndName);
	void _FASTA_Export(const CConf & conf, const string & ExportPathAndName);

};

}

#endif /* EXPORTPBUILDPRE_H_ */
