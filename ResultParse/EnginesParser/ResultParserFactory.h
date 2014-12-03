/*
 * ResultParserFactory.h
 *
 *  Created on: 2009-3-23
 *      Author: Haifeng Chen & Hao Chi
 */

#ifndef RESULTPARSERFACTORY_H_
#define RESULTPARSERFACTORY_H_

//#include "bio_analysis_sdk.h"

#include "pFindParser.h"
#include "MascotParser.h"
#include "SEQUESTParser.h"
#include "SQTParser.h"

using namespace bio_analysis;

class CResultParserFactory
{
public:
	CResultParser * GetParser(SearchEngineType eType)
	{
		switch (eType)
		{
		case ST_PFIND:
			return new CpFindParser;
		case ST_MASCOT:
			return new CMascotParser;
		case ST_SEQUEST:
			return new CSEQUESTParser;
		case ST_SQT:
			return new CSQTParser;

		default:
		{
			CErrInfo info("ResultParserFactory", "GetParser",
					"Con't find the SearchEngineType you want!");
			throw runtime_error(info.Get());
		}
		}
	}
};

#endif /* RESULTPARSERFACTORY_H_ */
