#include "Mass2PepIndexCreatorFactory.h"

#include "SPSingleScanIndexCreator.h"
//#include "SPMulitiScanIndexCreator.cpp.h"
//#include "MPMulitiScanIndexCreator.cpp.h"

namespace Mass2PepIndex{
CMass2PepIndexCreatorFactory::CMass2PepIndexCreatorFactory()
{
}

CMass2PepIndexCreatorFactory::~CMass2PepIndexCreatorFactory()
{
}

CMass2PepIndexCreator * CMass2PepIndexCreatorFactory::GetCreator(IndexCreatorType eType)
{
	switch(eType)
	{
		case ICT_Default: 
		case ICT_SPMulitiScan:
		case ICT_MPMulitiScan: 	
		case ICT_SPSingleScan:
			return new CSPSingleScanIndexCreator();
		case ICT_MPSingleScan:
		default:
			return new CSPSingleScanIndexCreator();
	}
}

}
