#include <iostream>
#include <string>
#include <time.h>
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"

#include "XLinkResultReport.h"
#include "XLinkpLabelReport.h"
#include "XLinkpBuildReport.h"
#include "XLinkpXBuildReport.h"
#include "XLinkProteinReport.h"
#include "XLinkPepReport.h"
#include "XLinkResultReportFactory.h"

using namespace std;
using namespace proteomics_sdk;

CXLinkResultReportFactory::CXLinkResultReportFactory()
{
	m_pTrace = CTrace::GetInstance();
	m_pTrace->Debug("Factory initialization.");
}

CXLinkResultReportFactory::~CXLinkResultReportFactory()
{
}

CXLinkResultReport * CXLinkResultReportFactory::GetReport(XLINK_RESULT_REPORT_TYPE reportType)
{
	m_pTrace->Debug("in the factory");
	switch(reportType)
	{
	case Xlink_STANDARD:
	{
		m_pTrace->Debug("in the xlink standard");
		return new CXLinkProteinReport();
	}
	case XLINK_PBUILD:
		return new CXLinkpBuildReport();
	case XLINK_PLABEL:
		return new CXLinkpLabelReport();
	case XLINK_PXBUILD:
		return new CXLinkpXBuildReport();
	case XLINK_PEPTIDE:
		return new CXLinkPepReport();
	case XLINK_PQUANT:
	
		
	default:
		return new CXLinkProteinReport(); 
	}
}
