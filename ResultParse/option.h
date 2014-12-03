#ifndef _OPTION_H_INCLUDED_
#define _OPTION_H_INCLUDED_

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <stdexcept>
#include <limits.h>

#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <sstream>
#include<set>
#include<fstream>
#include<iostream>
#include<iomanip>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <dirent.h>
#include <algorithm>
using namespace std;


const int BUFFER_SIZE = 81960;

using namespace std;

class CErrInfo
{
	string m_strClass;
	string m_strMethod;
	string m_strDetail;
	string m_strInfo;
	string m_strException;
public:
	CErrInfo(const string &strClass, const string &strMethod, const string &strDetail = "")
	{
		m_strClass = strClass;
		m_strMethod = strMethod;
		m_strDetail = strDetail;
	}

	CErrInfo(const string &strClass, const string &strMethod, const string &strDetail,
			const exception & e)
	{
		m_strClass = strClass;
		m_strMethod = strMethod;
		m_strDetail = strDetail;
		SetException(e);
	}
	void Append(const string &strInfo)
	{
		if (strInfo.empty())
			return;
		else
			m_strInfo += "\t\t  " + strInfo + "\n";
	}

	string Get() const
	{
		string strError = m_strException;
		strError += "\t  at " + m_strClass + "::" + m_strMethod + "() " + m_strDetail + "\n";
		strError += m_strInfo;
		return strError;
	}

	string Get(const exception& e)
	{
		SetException(e);
		return Get();
	}

	inline void SetException(const exception & e)
	{
		m_strException = e.what();
	}

	friend ofstream& operator<<(ofstream& os, const CErrInfo& info)
	{
		os << endl << "==========================" << endl;
		time_t current_time;
		time(&current_time);
		os << ctime(&current_time) << endl;
		os << info.Get() << endl;
		return os;
	}
	friend ostream& operator<<(ostream& os, const CErrInfo& info)
	{
		os << endl << "==========================" << endl;
		time_t current_time;
		time(&current_time);
		os << ctime(&current_time) << endl;
		os << info.Get() << endl;
		return os;
	}

};
class COptionTool
{
protected:
	void _GetPrivateProfileString(const string & strApp, const char * szKey,
			const char * szDefault, char szValue[], int nBufSize, const string & strFile)
	{
		szValue[0] = 0;
		FILE * fp = fopen(strFile.c_str(), "r");
		if (!fp)
		{
			strcpy(szValue, szDefault);
			CErrInfo info("COptionTool", "_GetPrivateProfileString", "cannot find the file.");
			throw runtime_error(info.Get().c_str());
		}
		size_t len = strlen(szKey);
		char szBuf[BUFFER_SIZE] =
		{ 0 };
		bool bRange = false;
		string str;
		str = "[" + strApp;
		str += "]";
		size_t lenApp = str.length();
		while (1)
		{
			if (0 == fgets(szBuf, BUFFER_SIZE - 1, fp))
			{
				break;
			}

			szBuf[lenApp] = 0;

			if (strcmp(str.c_str(), szBuf) == 0)
			{
				bRange = true;
				break;
			}
		}
		if (bRange)
		{
			while (1)
			{
				if (0 == fgets(szBuf, BUFFER_SIZE - 1, fp))
				{
					break;
				}
				size_t tCurrLen = strlen(szBuf) - 1;
				if (szBuf[len] != '=')
				{
					continue;
				}
				szBuf[len] = 0;

				if (0 == strcmp(szBuf, szKey))
				{
					while (tCurrLen >= 0 && (szBuf[tCurrLen] == 0xa || szBuf[tCurrLen] == 0xd))
					{
						szBuf[tCurrLen--] = 0;
					}
					strcpy(szValue, szBuf + len + 1);
					fclose(fp);
					return;
				}
			}
			strcpy(szValue, szDefault);
			fclose(fp);
			return;
		}
		else
		{
			strcpy(szValue, szDefault);
			fclose(fp);
			return;
		}

	}
	int _GetPrivateProfileInt(const string & strApp, const char * szKey, const int & nDefault,
			string & strFile)
	{
		char szValue[BUFFER_SIZE] =
		{ 0 };
		_GetPrivateProfileString(strApp, szKey, "", szValue, BUFFER_SIZE, strFile);
		if (0 == strlen(szValue))
		{
			return nDefault;
		}
		else
		{
			return atoi(szValue);
		}
	}

public:
	COptionTool(const char * szApp, const char * szFile = "option.ini") :
		m_strFile(szFile), m_strApp(szApp)
	{
	}
	;

	string GetString(const char * szKey, const char * szDefault = "")
	{
		char szValue[PATH_MAX] =
		{ 0 };
		_GetPrivateProfileString(m_strApp, szKey, szDefault, szValue, PATH_MAX, m_strFile);
		return string(szValue);
	}
	;

	int GetInteger(const char * szKey, int nDefault = 0)
	{
		return _GetPrivateProfileInt(m_strApp, szKey, nDefault, m_strFile);
	}
	;

protected:
	string m_strFile;
	string m_strApp;
};

#endif
