#ifndef DBCONF_H_
#define DBCONF_H_

#include <vector>
#include <string>
using namespace std;

typedef struct _DBCONF_HEAD
{
	string  enzymeList;	    
	int	  index_num;                  // The num of index existing
}DBCONF_HEAD;

typedef struct _DBCONF_ITEM
{
	string strDBName;   	    // The name of database, type string?
	string strEnzymeName;	// The name of enzyme, type string?
	string strPath;			// The path of index
	string strMetaName;		// The meta name of this index
	string strRemark;
	//time
}DBCONF_ITEM;


class CDBConfManage
{
public:
	CDBConfManage(void);
	virtual ~CDBConfManage(void);

	bool Load();
	bool Load(string strDBConfName);

	bool Open(string strOpenMode);
	bool Close();
	bool ReadDBConf();

	bool AddItem(DBCONF_ITEM sDBConfItem, bool bCheck = true, int iAddPos = -1);
	bool DelItem(DBCONF_ITEM sDBConfItem);

	bool CheckValidity();
	int CheckItemValidity(DBCONF_ITEM sDBConfItem);
	int CheckItemValidity(string strDBName, string strEnzyme,string strPath, string strMetaName);

	DBCONF_HEAD GetDBConfHead();
	bool GetDBConfItem(DBCONF_ITEM & sDBConfItem, string strDBName, string strEnzyme,string strPath, string strMetaName);
	vector<DBCONF_ITEM> GetDBConfAllItems();

	string GetDBConfName();
	bool Save();

private:
	DBCONF_HEAD m_sDBConfHead;
	vector<DBCONF_ITEM> m_vsDBConfItem;
	//string m_strDBConfPath;
	string m_strDBConfName;

	FILE * m_pDBConf;
};


#endif /*DBCONF_H_*/
