/*
 * MS2InputFactory.cpp
 *
 *  Created on: 2009-7-7
 *      Author: fan
 */


//#include "DTAInput.h"
//#include "DTASInput.h"
//#include "MGFInput.h"
//#include "PKLInput.h"
//#include "MS2TypeInput.h"
//#include "RAWInput.h"
//#include "MS2InputFactory.h"

#include "MS2OutputFactory.h"
#include "MGFOutput.h"
#include "MS2TypeOutput.h"
#include "DTAOutput.h"

CMS2OutputFactory::CMS2OutputFactory() {
	// TODO Auto-generated constructor stub

}

CMS2OutputFactory::~CMS2OutputFactory() {
	// TODO Auto-generated destructor stub
}

CMS2Output * CMS2OutputFactory::GetImporter(MS2FormatType eType)
{

	switch(eType)
		{
		case(PFF_DTA):
			return new CDTAOutput;
//		case(PFF_DTAS):
//			return new CDTASOutput;
		case(PFF_MGF):
			return new CMGFOutput;
//		case(PFF_PKL):
//			return new CPKLOutput;
		case(PFF_MS2):
			return new CMS2TypeOutput;
//		case(PFF_RAW):
//			return new CRAWOutput;
		default:
			return NULL;
		}
}
