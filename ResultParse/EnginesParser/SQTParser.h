#ifndef SQTPARSER_H_
#define SQTPARSER_H_

#include "../bio_analysis_sdk.h"
#include "ResultParser.h"

#include "../../include/sdk.h"

#include "../BasicFunction/StringUtility.h"
#include "../BasicFunction/CommonProcess.h"

namespace bio_analysis
{

class CSQTParser: public CResultParser
{
public:
	CSQTParser();
	virtual ~CSQTParser();

	virtual void
	Parse(const string & strFilePath, CParseredDataSetInfo & pbres);
private:
	void _GetSQTMod(string & strVal, vector<CModificationSiteInfo> & vModTemp);
	void _ParserFront(ifstream & fin, CParseredDataSetInfo & pbres);
	void _ParseSpectra(ifstream & fin, CParseredDataSetInfo & pbres);
	string m_strFilePath;
};

}

#endif /* SQTPARSER_H_ */
