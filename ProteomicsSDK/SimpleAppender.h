#ifndef SIMPLEAPPENDER_H_
#define SIMPLEAPPENDER_H_

namespace proteomics_sdk
{

class CSimpleAppender:public CAppender
{
public:
	
	CSimpleAppender();
	virtual ~CSimpleAppender();
	
	virtual void Init(const CCondition &condition) ;
	virtual void Out(const string &strOunt) ;
	virtual void Close(void) ;
};

}

#endif /*SIMPLEAPPENDER_H_*/
