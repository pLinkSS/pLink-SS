#ifndef RESULTPARSER_H_
#define RESULTPARSER_H_

#include "../bio_analysis_sdk.h"

namespace bio_analysis
{
class CResultParser
{
public:
	virtual ~CResultParser()
	{
	}
	;
	virtual void
	Parse(const string & strFilePath, CParseredDataSetInfo & pbres) = 0;
};
}

#endif /* RESULTPARSER_H_ */
