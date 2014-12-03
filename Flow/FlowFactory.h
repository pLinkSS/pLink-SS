#ifndef FLOW_FACTORY_H_
#define FLOW_FACTORY_H_

namespace proteomics_search
{
class CSearchThread;
class CFlowFactory
{
public:
	CFlowFactory();
	virtual ~CFlowFactory();
	CFlow * GetFlow(IndexContent eType, vector<CSpectrum> & vSpec);
};

}

#endif /*FLOW_FACTORY_H_*/
