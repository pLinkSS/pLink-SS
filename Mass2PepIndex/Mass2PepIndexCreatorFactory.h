#ifndef MASS2PEPINDEXCREATORFACTORY_H_
#define MASS2PEPINDEXCREATORFACTORY_H_

#include "Mass2PepIndexCreator.h"

namespace Mass2PepIndex{

enum IndexCreatorType
{
	ICT_Default = 0,
	ICT_SPMulitiScan = 1,
	ICT_MPMulitiScan = 2,	
	ICT_SPSingleScan = 3,
	ICT_MPSingleScan = 4
};

class CMass2PepIndexCreatorFactory
{
public:
	CMass2PepIndexCreatorFactory();
	virtual ~CMass2PepIndexCreatorFactory();
	
	CMass2PepIndexCreator * GetCreator(IndexCreatorType eType);
};

}
#endif /*MASS2PEPINDEXCREATORFACTORY_H_*/
