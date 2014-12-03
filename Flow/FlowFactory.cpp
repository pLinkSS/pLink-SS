#include <set>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "../Mass2PepIndex/Mass2PepIndex.h"
#include "../Mass2PepIndex/PeptideReaderFactory.h"
#include "SpectraSearcher.h"
#include "XLinkPepFlow.h"
#include "XLinkOpenFlow.h"
#include "FlowFactory.h"

namespace proteomics_search
{


CFlowFactory::CFlowFactory()
{
}

CFlowFactory::~CFlowFactory()
{
}

CFlow * CFlowFactory::GetFlow(IndexContent eType, vector<CSpectrum> & vSpec)
{			
	if(PEPTIDE_PAIR == eType)
	{
		// pep flow
		// a pep pair list is built firstly
		return new CXLinkPepSalvoFlow(vSpec);
	}
	else if (PEPTIDE_PAIR_OPEN == eType)
	{
		// open flow for xlink to deal with large database
		// no pep pair list is built
		// first : use non-stringent scoring to filter out un-relevant peptides
		//			high sensitivity is required
		// second : use stringent scoring to refine the search results
		//			high precision is required

		return new CXLinkOpenFlow(vSpec);
	}
	else
	{
		return new CXLinkPepSalvoFlow(vSpec);
	}

}

}
