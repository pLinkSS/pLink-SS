#ifndef TRACE_H_
#define TRACE_H_


namespace proteomics_sdk
{
class CAppender;

class CTrace
{
protected:	
	CTrace(CAppender * pAppender);
		
	static CTrace * m_pInstance;

public:
	
	static CTrace * GetInstance();
	
	static CTrace * GetInstance(CCondition & condition);
	
	static void DeleteInstance();
	void Set(CCondition & condition);
		
	void Alert(string strAlert, ModuleType eType=MODULE_ALL);
	
	void Info(string strInfo, ModuleType eType=MODULE_ALL);

	void Debug(string strDebug, ModuleType eType=MODULE_ALL);
	
	void Break();

protected:		
	CAppender * m_pAppender;
	ModuleType m_eModule;
	LogRankType m_eLogRank;	
};


}
#endif /*TRACE_H_*/
