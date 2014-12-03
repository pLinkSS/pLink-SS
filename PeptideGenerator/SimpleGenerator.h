#ifndef PEPTIDEGENERATOR_H_
#define PEPTIDEGENERATOR_H_
#include "../include/interface.h"
using namespace proteomics_sdk;
namespace proteomics_sdk
{
class CAAConf;

class CSimpleGenerator:
	public CPeptideGenerator{
		
public:
	virtual void Init(const CAAConf & aaconf, bool bMono);
	virtual void GetPeptide(double lfMass, string & strOut);
	virtual void Batch(double lfMass, size_t tNum, vector<string> & vstrOut);
    virtual void Close(void);

protected:
	double m_lfMass[26];
	char m_cAA[26];
};

}

#endif /*PEPTIDEGENERATOR_H_*/
