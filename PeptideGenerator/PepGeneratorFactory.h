#ifndef _PEP_GENERATOR__FACTORY_H_INCLUDED_
#define _PEP_GENERATOR__FACTORY_H_INCLUDED_

using namespace proteomics_sdk;
using namespace std;

namespace proteomics_sdk
{
class CPeptideGenerator;

class CPepGeneratorFactory
{
public:
	CPepGeneratorFactory(void){;};
	~CPepGeneratorFactory(void){;};

	CPeptideGenerator * GetGenerator(PepGeneType eGeneType) const;
};


}


#endif /*_PEP_GENERATOR__FACTORY_H_INCLUDED_*/
