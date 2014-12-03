#ifndef XLINKRESULTREPORTFACTORY_H_
#define XLINKRESULTREPORTFACTORY_H_


class CXLinkResultReportFactory
{
public:
	CXLinkResultReportFactory();
	virtual ~CXLinkResultReportFactory();
	
	CXLinkResultReport * GetReport(XLINK_RESULT_REPORT_TYPE reportType);
private:
	CTrace *m_pTrace;
};

#endif /*XLINKRESULTREPORTFACTORY_H_*/
