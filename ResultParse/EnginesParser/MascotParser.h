#ifndef MASCOTPARSER_H_
#define MASCOTPARSER_H_

#include "ResultParser.h"

#include "../../include/sdk.h"
#include "../BasicFunction/ReadWrite.h"
#include "../BasicFunction/StringUtility.h"
#include "../BasicFunction/CommonProcess.h"

#include <fstream>
#include <sstream>

namespace bio_analysis
{

class CMascotParser: public CResultParser
{
	string m_strBoundary;
public:
	CMascotParser();
	virtual ~CMascotParser();

	virtual void
	Parse(const string & strFilePath, CParseredDataSetInfo & pbres);
private:
	void _ParseParameters(ifstream & fin, string & strBuf, CConditionInfo & cond);
	void _ParseHeader(ifstream & fin, string & strBuf, CConditionInfo & cond, int & nTotal);
	void _ParseMasses(ifstream & fin, string & strBuf, CConditionInfo & cond);
	void _ParseSummary(ifstream & fin, string & strBuf, int *naCharge);
	void _ParsePeptides(ifstream & fin, string & strBuf, vector<CMatchSpectraInfo> & vSpectraInfo,
			CConditionInfo & cond, const int * anCharge);
	void _ParseQuery(ifstream & fin, string & strBuf, vector<CMatchSpectraInfo> & vSpectraInfo,
			const int & nTotal);
};

}

#endif /* MASCOTPARSER_H_ */
