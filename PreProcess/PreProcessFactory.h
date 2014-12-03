#ifndef _PRE_PROCESS_FACTORY_H
#define _PRE_PROCESS_FACTORY_H

using namespace std;

namespace proteomics_sdk
{
class CPreProcess;

class CPreProcessFactory
{
public:
	CPreProcessFactory(void);
	~CPreProcessFactory(void);

	CPreProcess * GetPreProcessor(PreProcType ePP) const;
};


}

#endif//_PRE_PROCESS_FACTORY_H
