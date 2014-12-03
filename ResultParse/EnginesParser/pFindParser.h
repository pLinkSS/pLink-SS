#ifndef PFINDPARSER_H_
#define PFINDPARSER_H_

#include "../bio_analysis_sdk.h"
#include "ResultParser.h"
#include "../../include/sdk.h"

#include "../BasicFunction/ReadWrite.h"
#include "../BasicFunction/StringUtility.h"
#include "../BasicFunction/CommonProcess.h"

namespace bio_analysis
{
class CpFindParser: public CResultParser
{
public:
	CpFindParser();
	virtual ~CpFindParser();

	virtual void
	Parse(const string & strFilePath, CParseredDataSetInfo & pbres);
private:
	void _ParseMatchSpectraInfo(ifstream & fin, CMatchSpectraInfo & mr);
	void _ParsePeptideInfor(ifstream & fin, CMatchSpectraInfo & mr, CMatchPeptideInfo & mp,
			size_t tRank);
	void _ParseConditonInfo(ifstream & fin, CParseredDataSetInfo & pbres, size_t & nTotal);
	void
	_ParseTXTFile(const string & strFilePath, CParseredDataSetInfo & pbres);
};
}
;
#endif /*PFINDPARSER_H_*/
