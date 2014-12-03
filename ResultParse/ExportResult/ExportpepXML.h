#ifndef EXPORTPEPXML_H_
#define EXPORTPEPXML_H_

#include "ResultExport.h"

namespace bio_analysis
{

class CExportpepXML: public CResultExport
{
public:
	CExportpepXML();
	virtual ~CExportpepXML();

	virtual void Export(const IntegratedInfo & HaveToInterData,
			const CConf & conf);
};

}

#endif /* EXPORTPEPXML_H_ */
