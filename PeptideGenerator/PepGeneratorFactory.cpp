#include <string>
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/predefine.h"
#include "common.h"
#include "SimpleGenerator.h"
#include "PepGeneratorFactory.h"

using namespace std;
using namespace proteomics_sdk;

CPeptideGenerator * CPepGeneratorFactory::GetGenerator(PepGeneType eGeneType) const
{

	switch (eGeneType)
	{
	case PEPGENE_DEFAULT:
		return new CSimpleGenerator();
		
	default:
		return new CSimpleGenerator();
		
	}
	
};

