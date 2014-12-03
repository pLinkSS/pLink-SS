#ifndef _OPTION_H_INCLUDED_
#define _OPTION_H_INCLUDED_

#include <string>
#include <iostream>
#include <vector>
#include <dirent.h>
#include <string.h>

using namespace std;

namespace proteomics_sdk {
/* 81920 is too much */
const size_t BUFFER_SIZE = 1000;
class COptionTool {
protected:
	void _GetPrivateProfileString(const string & strApp, const char * szKey,
			const char * szDefault, char szValue[], int nBufSize,
			const string & strFile) {
		szValue[0] = '\0';
		FILE * fp = fopen(strFile.c_str(), "r");
		if (!fp) {
			strcpy(szValue, szDefault);
			CErrInfo info("COptionTool", "_GetPrivateProfileString",
					"cannot find the file.");
			throw runtime_error(info.Get().c_str());
		}
		size_t len = strlen(szKey);
		char szBuf[BUFFER_SIZE] = { 0 };
		bool bRange = false;
		string str;
		str = "[" + strApp;
		str += "]";
		size_t lenApp = str.length();
		while (1) {
			if (0 == fgets(szBuf, BUFFER_SIZE - 1, fp)) {
				break;
			}

			szBuf[lenApp] = 0;

			if (strcmp(str.c_str(), szBuf) == 0) {
				bRange = true;
				break;
			}
		}
		if (bRange) {
			while (1) {
				if (0 == fgets(szBuf, BUFFER_SIZE - 1, fp)) {
					break;
				}
				size_t tCurrLen = strlen(szBuf) - 1;
				if (szBuf[len] != '=') {
					continue;
				}
				szBuf[len] = 0;

				if (0 == strcmp(szBuf, szKey)) {
					while (tCurrLen >= 0
							&& (szBuf[tCurrLen] == 0xa || szBuf[tCurrLen] == 0xd)) {
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
		} else {
			strcpy(szValue, szDefault);
			fclose(fp);
			return;
		}

	}
	;
	int _GetPrivateProfileInt(const string & strApp, const char * szKey,
			const int & nDefault, string & strFile) {
		char szValue[BUFFER_SIZE] = { '\0' };
		_GetPrivateProfileString(strApp, szKey, "", szValue, BUFFER_SIZE,
				strFile);
		if (0 == strlen(szValue)) {
			return nDefault;
		} else {
			return atoi(szValue);
		}
	}
	;

public:
	COptionTool(const char * szApp, const char * szFile = "option.ini") :
			m_strFile(szFile), m_strApp(szApp) {
	}

	string GetString(const char * szKey, const char * szDefault = "") {
		char szValue[PATH_MAX] = { 0 };
		_GetPrivateProfileString(m_strApp, szKey, szDefault, szValue, PATH_MAX,
				m_strFile);
		return string(szValue);
	}

	string GetString(const char * szApp, const char * szKey,
			const char * szDefault) {
		char szValue[PATH_MAX] = { 0 };
		_GetPrivateProfileString(szApp, szKey, szDefault, szValue, PATH_MAX,
				m_strFile);
		return string(szValue);
	}

	int GetInteger(const char * szKey, int nDefault = 0) {
		return _GetPrivateProfileInt(m_strApp, szKey, nDefault, m_strFile);
	}

	int GetInteger(const char * szApp, const char * szKey, int nDefault) {
		return _GetPrivateProfileInt(szApp, szKey, nDefault, m_strFile);
	}

	size_t GetSizeT(const char * strSection, const char * strKey,
			const size_t tDefault) {
		return (size_t) GetInteger(strSection, strKey, (int) tDefault);
	}

	bool GetBool(const char * strSection, const char * strKey,
			const bool bDefault) {
		return (bool) GetInteger(strSection, strKey, (int) bDefault);
	}

protected:
	string m_strFile;
	string m_strApp;
};

class CDBConf {
public:
	CDBConf(void) {
	}
	;
	CDBConf(const char * szDB) {
		Load(szDB);
	}
	;

	bool Load(void) {
		FILE * pDBConf = NULL;
		pDBConf = fopen(m_strDBConfName.c_str(), "r");

		if (NULL == pDBConf) {
			CErrInfo info("CDBConf", "Load", "fopen returns NULL");
			info.Append("m_strDBConfName=" + m_strDBConfName);
			throw runtime_error(info.Get().c_str());
		}

		fclose(pDBConf);

		COptionTool * pOption = new COptionTool("total",
				m_strDBConfName.c_str());

		m_strEnzymeListPath = pOption->GetString("enzyme_list", "enzyme.ini");
		m_strModifyListPath = pOption->GetString("modify_list", "modify.ini");

		int iIndexNum = (int) (pOption->GetInteger("index_num", 1));

		delete pOption;

		char str[20];

		for (int i = 1; i <= iIndexNum; ++i) {
			sprintf(str, "%d", i);
			string strIndex = "index" + string(str);

			COptionTool * pOption = new COptionTool(strIndex.c_str(),
					m_strDBConfName.c_str());

			string strDBName = pOption->GetString("db_name", "Null");
			string strEnzyme = pOption->GetString("enzyme", "Null");
			string strPath = pOption->GetString("path", "Null");
			string strMetaName = pOption->GetString("meta_name", "Null");

			m_vstrDBName.push_back(strDBName);
			m_vstrEZName.push_back(strEnzyme);
			m_vstrPath.push_back(strPath);
			m_vstrMetaName.push_back(strMetaName);

			delete pOption;
		}
		return true;
	}
	;

	bool Load(const char * szDB) {
		m_strDBConfName = szDB;

		return Load();
	}
	;

	vector<string> GetDBNameList(void) {
		vector<string> vDB;
		for (size_t i = 0; i < m_vstrDBName.size(); ++i) {
			bool bFind(false);
			for (size_t j = 0; j < vDB.size(); j++)
				if (vDB[j] == m_vstrDBName[i])
					bFind = true;

			if (!bFind)
				vDB.push_back(m_vstrDBName[i]);
		}

		return vDB;
	}
	;

	vector<string> GetEnzymeNameList(string strDB) {
		vector<string> vEnzyme;
		for (size_t i = 0; i < m_vstrDBName.size(); ++i) {
			if (strDB == m_vstrDBName[i])
				vEnzyme.push_back(m_vstrEZName[i]);
		}

		return vEnzyme;
	}
	;

	const char * GetPath(const char * szDBName, const char * szEnzymeName) {
		for (int i = 0; i < (int) (m_vstrDBName.size()); ++i) {
			if (m_vstrDBName[i] == szDBName && m_vstrEZName[i] == szEnzymeName)
				return m_vstrPath[i].c_str();
		}

		return "NULL";
	}
	;

	const char * GetMetaName(const char * szDBName, const char * szEnzymeName,
			bool bMono) {
		string strMono = bMono ? "mono.meta" : "avrg.meta";
		int i = 0;
		for (i = 0; i < (int) (m_vstrDBName.size()); ++i) {
			if (m_vstrDBName[i] == szDBName
					&& m_vstrEZName[i] == szEnzymeName) {
				if (m_vstrMetaName[i].find(strMono) != string::npos)
					return m_vstrMetaName[i].c_str();
			}
		}

		return "NULL";
	}
	;

	bool Insert(const char * szDBName, const char * szEnzymeName,
			const char * szPath, const char * szMetaName) {
		int i = 0;
		for (i = 0; i < (int) (m_vstrDBName.size()); ++i) {
			//if the index is already exist,then update it
			if (m_vstrDBName[i] == szDBName
					&& m_vstrEZName[i] == szEnzymeName) {
				m_vstrPath[i] = szPath;
				m_vstrMetaName[i] = szMetaName;
				break;
			}
		}

		//new index
		if (i == (int) (m_vstrDBName.size())) {

			m_vstrDBName.push_back(szDBName);
			m_vstrEZName.push_back(szEnzymeName);
			m_vstrPath.push_back(szPath);
			m_vstrMetaName.push_back(szMetaName);

		}
		return true;
	}
	;

	bool SetEnzymeListPath(const char * szEnzymeListPath) {
		m_strEnzymeListPath = szEnzymeListPath;
		return true;
	}
	;

	const char * GetEnzymeListPath(void) {
		return m_strEnzymeListPath.c_str();
	}
	;

	bool SetModifyListPath(const char * szModifyListPath) {
		m_strModifyListPath = szModifyListPath;
		return true;
	}
	;

	const char * GetModifyListPath(void) {
		return m_strModifyListPath.c_str();
	}
	;

	bool Save(void) {
		FILE * pDBConf = fopen(m_strDBConfName.c_str(), "w");

		fprintf(pDBConf, "[total]\n");
		fprintf(pDBConf, "enzyme_list = %s\n", m_strEnzymeListPath.c_str());
		fprintf(pDBConf, "modify_list = %s\n", m_strModifyListPath.c_str());

		int iIndexNum = (int) (m_vstrDBName.size());
		fprintf(pDBConf, "index_num=%d\n", iIndexNum);

		char str[20];

		for (int i = 0; i < iIndexNum; ++i) {
			sprintf(str, "%d", i + 1);
			string strIndex = "[index" + string(str) + "]";
			fprintf(pDBConf, "%s\n", strIndex.c_str());

			fprintf(pDBConf, "db_name=%s\n", m_vstrDBName[i].c_str());
			fprintf(pDBConf, "enzyme=%s\n", m_vstrEZName[i].c_str());
			fprintf(pDBConf, "path=%s\n", m_vstrPath[i].c_str());
			fprintf(pDBConf, "meta_name=%s\n", m_vstrMetaName[i].c_str());
		}

		fclose(pDBConf);
		return true;
	}
	;

private:

	string m_strDBConfName;
	string m_strEnzymeListPath;
	string m_strModifyListPath;

	vector<string> m_vstrDBName;
	vector<string> m_vstrEZName;
	vector<string> m_vstrPath;
	vector<string> m_vstrMetaName;
};

class CEnzyme;

class CEnzymeConf {
public:
	CEnzymeConf(string strEnzymeListPath) :
			m_strEnzymeListPath(strEnzymeListPath) {
	}
	;

	bool GetEnzyme(string strName, CEnzyme & enzyme) const {
		COptionTool option("enzyme", m_strEnzymeListPath.c_str());
		string strEnzymeValue = option.GetString(strName.c_str(), "Null");

		if ("Null" == strEnzymeValue)
			return false;

		enzyme.m_strName = strName;
		enzyme.InitByString(strEnzymeValue);

		return true;
	}
	;

	vector<CEnzyme> GetEnzymeList(void) const {
		vector<CEnzyme> vEnzyme;

		COptionTool option("enzyme", m_strEnzymeListPath.c_str());
		char szValue[BUFFER_SIZE] = { 0 };
		int nEnzyme = option.GetInteger("total", 0);
		for (int i = 0; i < nEnzyme; ++i) {
			sprintf(szValue, "name%d", i + 1);
			string strKey(szValue);

			string strEnzymeName = option.GetString(strKey.c_str(), "Null");
			if ("Null" == strEnzymeName)
				continue;

			string strEnzymeValue = option.GetString(strEnzymeName.c_str(),
					"Null");
			if ("Null" == strEnzymeValue)
				continue;

			CEnzyme enzyme;
			enzyme.m_strName = strEnzymeName;
			enzyme.InitByString(strEnzymeValue);
			vEnzyme.push_back(enzyme);
		}

		return vEnzyme;
	}
	;

	string m_strEnzymeListPath;
};

class CXLinker;

// pfind-xlink
class CXLinkerConf {
public:

	string m_strXLinkerListPath;

	CXLinkerConf(string strXLinkerListPath) :
			m_strXLinkerListPath(strXLinkerListPath) {
	}
	;

	bool GetXLinker(string strName, CXLinker & linker) const {
		COptionTool option("xlink", m_strXLinkerListPath.c_str());

		string strXLinkerValue = option.GetString(strName.c_str(), "Null");

		if ("Null" == strXLinkerValue)
			return false;

		linker.m_strName = strName;
		linker.InitByString(strXLinkerValue);
		return true;
	}
	;

	vector<CXLinker> GetXLinkerList(void) const {
		vector<CXLinker> vLinker;
		COptionTool option("xlink", m_strXLinkerListPath.c_str());
		char szValue[BUFFER_SIZE] = { 0 };
		int nLinker = option.GetInteger("total", 0);

		for (int i = 0; i < nLinker; ++i) {
			sprintf(szValue, "name%d", i + 1);
			string strKey(szValue);

			string strLinkerName = option.GetString(strKey.c_str(), "Null");
			if ("Null" == strLinkerName)
				continue;

			string strLinkerValue = option.GetString(strLinkerName.c_str(),
					"Null");
			if ("Null" == strLinkerValue)
				continue;

			CXLinker linker;
			linker.m_strName = strLinkerName;
			linker.InitByString(strLinkerValue);
			vLinker.push_back(linker);
		}
		return vLinker;
	}
	;
};

class CModification;

class CModifyConf {
public:
	CModifyConf(string strModifyListPath) :
			m_strModifyListPath(strModifyListPath) {
	}
	;

	bool GetModify(string strName, CModification & modify) const {
		COptionTool option("modify", m_strModifyListPath.c_str());

		string strModifyValue = option.GetString(strName.c_str(), "Null");

		if ("Null" == strModifyValue)
			return false;

		modify.m_strName = strName;
		modify.InitByString(strModifyValue);

		return true;
	}
	;

	vector<CModification> GetModifyList(void) const {
		vector<CModification> vMod;

		COptionTool option("modify", m_strModifyListPath.c_str());
		char szValue[BUFFER_SIZE] = { 0 };
		int nModify = option.GetInteger("total", 0);
		for (int i = 0; i < nModify; ++i) {
			sprintf(szValue, "name%d", i + 1);
			string strKey(szValue);

			string strModifyName = option.GetString(strKey.c_str(), "Null");
			if ("Null" == strModifyName)
				continue;

			string strModifyValue = option.GetString(strModifyName.c_str(),
					"Null");
			if ("Null" == strModifyValue)
				continue;

			CModification mod;
			mod.m_strName = strModifyName;
			mod.InitByString(strModifyValue);
			vMod.push_back(mod);
		}

		return vMod;
	}
	;

	string m_strModifyListPath;
};

class CMapAAMass;

class CAAConf {
public:
	CAAConf(string strAAListPath) :
			m_strAAListPath(strAAListPath) {
	}
	;

	bool GetAA(string strName, CAA & aa) const {
		COptionTool option("aa", m_strAAListPath.c_str());

		string strAAValue = option.GetString(strName.c_str(), "Null");

		if ("Null" == strAAValue)
			return false;

		aa.m_cAA = strName[0];
		aa.InitByString(strAAValue);

		return true;
	}
	;

	CMapAAMass GetMapAAMass(void) const {
		CMapAAMass map;
		COptionTool option("aa", m_strAAListPath.c_str());

		char szValue[BUFFER_SIZE] = { 0 };
		int nAA = option.GetInteger("total", 0);
		for (int i = 0; i < nAA; ++i) {
			sprintf(szValue, "name%d", i + 1);
			string strKey(szValue);

			string strAAName = option.GetString(strKey.c_str(), "Null");
			if ("Null" == strAAName)
				continue;

			string strAAValue = option.GetString(strAAName.c_str(), "Null");
			if ("Null" == strAAValue)
				continue;

			CAA aa;

			aa.m_cAA = strAAName[0];
			aa.InitByString(strAAValue);

			map.SetAAMap(aa.m_cAA, aa);
		}

		return map;
	}
	;

	string m_strAAListPath;
};

}

#endif
