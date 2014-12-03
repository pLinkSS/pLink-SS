#ifndef SEQUESTPARSER_H_
#define SEQUESTPARSER_H_

#include "../bio_analysis_sdk.h"
#include "ResultParser.h"
#include <dirent.h>
#include "../BasicFunction/ReadWrite.h"
#include "../BasicFunction/StringUtility.h"
#include "../BasicFunction/CommonProcess.h"

#include <fstream>
#include <sstream>

namespace bio_analysis
{

class CSEQUESTParser: public CResultParser
{
public:
	CSEQUESTParser();
	virtual ~CSEQUESTParser();
	virtual void
	Parse(const string & strFilePath, CParseredDataSetInfo & pbres);
private:
	void _GetVarMod(string & strTemp, CConditionInfo & cond, size_t & pos);
	void _GetFixMod(string & strTemp, CConditionInfo & cond, size_t & pos);
	void _ParserFront(ifstream & fin, CMatchSpectraInfo & SpectraInfo, CConditionInfo & cond);
	void _ParserPeptide(ifstream & fin, CMatchSpectraInfo & SpectraInfo, string & InputPath);
	void _SEQUSETParser(ifstream & fin, CMatchSpectraInfo & SpectraInfo, CConditionInfo & cond,
			string & InputPath);
	void _SEQUESTParserFile(const string & strFilePath, CParseredDataSetInfo & pbres);
	void _SEQUESTParserFolder(const string & strFilePath, CParseredDataSetInfo & pbres);
	void _GetSequestMod(string & strVal, vector<CModificationSiteInfo> & vModTemp);
};

}

#endif /* SEQUESTPARSER_H_ */
