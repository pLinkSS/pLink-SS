#ifndef XLINKRESULTFILTERFACTORY_H_
#define XLINKRESULTFILTERFACTORY_H_

class CXLinkResultFilterFactory
{
public:
	CXLinkResultFilterFactory();
	virtual ~CXLinkResultFilterFactory();
	CXLinkResultFilterInterface * GetFilter(int nType);
};

#endif /*XLINKRESULTFILTERFACTORY_H_*/
