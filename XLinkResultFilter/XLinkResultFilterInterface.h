#ifndef XLINKRESULTFILTERINTERFACE_H_
#define XLINKRESULTFILTERINTERFACE_H_

#include <iostream>
#include <string>

#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "../XLinkPepResultFilter/XLinkPepResultFilter.h"
#include "FilterConf.h"
using namespace std ;
using namespace proteomics_sdk;

class CXLinkResultFilterInterface
{
public:
	CXLinkResultFilterInterface(){};
	virtual ~CXLinkResultFilterInterface(){};
	virtual void Init(string strConf) = 0;
	virtual void Run() = 0;
	virtual void Close() = 0;
};


#endif /*XLINKRESULTFILTERINTERFACE_H_*/
