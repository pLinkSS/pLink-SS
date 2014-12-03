#ifdef WIN32
#include "windows.h"	
#else
#include "stdio.h"
#endif

#include <string>
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "AppenderFactory.h"
#include "Trace.h"

using namespace std;

namespace proteomics_sdk
{
CTrace * CTrace::m_pInstance(NULL);

CTrace::CTrace(CAppender * pAppender)
		:m_pAppender(pAppender), m_eModule(MODULE_ALL), m_eLogRank(LOG_RANK_INFO)
{
	;
};
		
CTrace * CTrace::GetInstance()
{
	if (NULL == m_pInstance)
	{
		CAppenderFactory factory;
		CAppender * pAppender = factory.GetAppender(LOG_APPENDER_CMD);
		m_pInstance = new CTrace(pAppender);
	}
	return m_pInstance;
}
	
CTrace * CTrace::GetInstance(CCondition & condition)
{
	if (NULL == m_pInstance)
	{
		CAppenderFactory factory;
		CAppender * pAppender = factory.GetAppender(condition.m_eLogAppender);
		m_pInstance = new CTrace(pAppender);
		m_pInstance->Set(condition);
	}
	else
	{
		// modify by emily 
		m_pInstance->Set(condition);
	}
	return m_pInstance;
}

void CTrace::DeleteInstance()
{
	// do not use this function
	if (NULL != m_pInstance)
	{
		delete m_pInstance;
	}
	m_pInstance = NULL;
}
void CTrace::Set(CCondition & condition)
{
	m_eModule=condition.m_eModule;
	m_eLogRank=condition.m_eLogRank;
};
		
void CTrace::Alert(string strAlert, ModuleType eType)
{
	if((MODULE_ALL == m_eModule || MODULE_ALL == eType || eType == m_eModule) &&
			LOG_RANK_ALERT <= m_eLogRank)
	{
#ifdef WIN32
	    HANDLE consol = GetStdHandle(STD_OUTPUT_HANDLE);
	    SetConsoleTextAttribute(consol,FOREGROUND_RED|FOREGROUND_INTENSITY);
	    m_pAppender->Out("[alert] "+strAlert);
	    SetConsoleTextAttribute(consol,FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE); 
#else
        /* set foreground color red */
	    m_pAppender->Out("\033[1;31m[alert] "+strAlert+"\033[1;37m");
#endif

	}
};
	
void CTrace::Info(string strInfo, ModuleType eType)
{
	if((MODULE_ALL == m_eModule || MODULE_ALL == eType || eType == m_eModule) &&
			LOG_RANK_INFO <= m_eLogRank)
	{		
		m_pAppender->Out("[info] "+strInfo);
	}	

};

void CTrace::Debug(string strDebug, ModuleType eType)
{
	if((MODULE_ALL == m_eModule || MODULE_ALL == eType || eType == m_eModule) &&
			LOG_RANK_DEBUG <= m_eLogRank)
	{		
#ifdef WIN32
		HANDLE consol = GetStdHandle(STD_OUTPUT_HANDLE);
		SetConsoleTextAttribute(consol,FOREGROUND_GREEN|FOREGROUND_INTENSITY);
		m_pAppender->Out("[debug] "+strDebug);
		SetConsoleTextAttribute(consol,FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE);
#else
		/* set foreground color light green */
		m_pAppender->Out("\033[1;32m[debug] "+strDebug+"\033[1;37m");
#endif
	}
};

void CTrace::Break()
{
//	if(LOG_RANK_DEBUG <= m_eLogRank)
//	{
//		cout << "[stop] Press any key to continue ..." << endl;
//		getchar();
//	}
}


}

