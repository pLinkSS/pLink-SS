/*
 * MS2InputFactory.cpp
 *
 *  Created on: 2009-7-7
 *      Author: fan
 */


#include "DTAInput.h"
#include "DTASInput.h"
#include "MGFInput.h"
#include "PKLInput.h"
#include "MS2TypeInput.h"
#include "RAWInput.h"
#include "MS2InputFactory.h"
#include "DTASingleInput.h"


CMS2InputFactory::CMS2InputFactory() {
	// TODO Auto-generated constructor stub
}

CMS2InputFactory::~CMS2InputFactory() {
	// TODO Auto-generated destructor stub
}


CMS2Input * CMS2InputFactory::GetImporter(MS2FormatType eType)
{

	switch(eType)
		{
		case(PFF_DTA):
			return new CDTAInput;
		case(PFF_DTAS):
			return new CDTASInput;
		case(PFF_MGF):
			return new CMGFInput;
		case(PFF_PKL):
			return new CPKLInput;
		case(PFF_MS2):
			return new CMS2TypeInput;
		case(PFF_RAW):
			return new CRAWInput;
		case(PFF_SDTA):
			return new CDTASingleInput;
		default:
			return NULL;
		}
}


