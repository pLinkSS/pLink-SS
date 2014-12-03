#ifndef RESULTOUT_H_
#define RESULTOUT_H_

#include "../bio_analysis_sdk.h"

namespace bio_analysis
{

class CResultExport
{
public:
	virtual ~CResultExport()
	{
	}
	virtual void
	Export(const IntegratedInfo & HaveToInterData, const CConf & conf) = 0;
};
}

#endif /* RESULTOUT_H_ */
