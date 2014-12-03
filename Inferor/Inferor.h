#ifndef INFEROR_H_
#define INFEROR_H_

#include <string>
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
using namespace std;

enum INFEROR_TYPE
{
	COMMON_INFEROR,
	XLINK_INFEROR,
	TRI_INFEROR,
};

class CInferor
{
public:
	CInferor(){};
	virtual ~CInferor(){};
	
	virtual bool Init(string strOption, time_t tmStartTime = 0) = 0;
	virtual void Run(string strParams, string strOutputPath, int nFileTotal, int nSpectraTotal, string strIdentifier) = 0;
	
};

#endif /*INFEROR_H_*/
