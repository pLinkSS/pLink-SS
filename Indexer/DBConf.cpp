//#include <stdlib.h>
#include <string>
#include <vector>
#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/option.h"
#include "ACEConfigTool.h"

using namespace std;
using namespace proteomics_sdk;

#include "DBConf.h"

CDBConfManage::CDBConfManage(void):
m_strDBConfName(".//dbconf.ini")
{
;
}

CDBConfManage::~CDBConfManage(void)
{
}

bool CDBConfManage::Load()
{
	if (!Open("r"))//check if the dbconf.ini is exist
	{
		FILE * f = fopen("dbconf.ini","w");
		fclose(f);
		return false;
	}

	Close();

	return ReadDBConf();
}
bool CDBConfManage::Load(string strDBConfName)
{
		m_strDBConfName = strDBConfName;

		return Load();
}

bool CDBConfManage::Open(string strOpenMode)
{
	m_pDBConf = fopen(m_strDBConfName.c_str(), strOpenMode.c_str());

	if( 0 == m_pDBConf)
	{
		return false;
	}

	return true;
}

bool CDBConfManage::Close()
{
	if (0 != m_pDBConf)
	{
		fclose(m_pDBConf);
		return true;
	}
	
	return false;
}

bool CDBConfManage::ReadDBConf()
{
	m_vsDBConfItem.clear();

	CACEConfigTool config;
	config.Open(m_strDBConfName);
	
	const string strSection("total");
	m_sDBConfHead.enzymeList = config.GetString(strSection, string("enzyme_list"));
	int iIndexNum = config.GetInteger(strSection, string("index_num"));

	char str[20];

	for(int i=1; i<=iIndexNum; i++)
	{
//		itoa( i, str, 10);
		sprintf(str,"%d",i);
		string strIndex = "index" + string(str);

		DBCONF_ITEM sDBConfItem;

		sDBConfItem.strDBName = config.GetString(strIndex, string("db_name"));
		sDBConfItem.strEnzymeName = config.GetString(strIndex, string("enzyme"));
		sDBConfItem.strPath = config.GetString(strIndex, string("path"));
		sDBConfItem.strMetaName = config.GetString(strIndex, string("meta_name"));
		sDBConfItem.strRemark = config.GetString(strIndex, string("remark"));

		m_vsDBConfItem.push_back(sDBConfItem);	

	}

	return true;
}

bool CDBConfManage::AddItem(DBCONF_ITEM sDBConfItem, bool bCheck, int iAddPos)
{
	if (bCheck)
	{
		iAddPos = CheckItemValidity(sDBConfItem);
		
		if (-1 == iAddPos)//new index	
		{	
			m_vsDBConfItem.push_back(sDBConfItem);
			return true;
		}

		else 
		{	
			return false;
		}
	}

	else
	{
		m_vsDBConfItem[iAddPos] = sDBConfItem;
		return true;
	}
}

bool CDBConfManage::DelItem(DBCONF_ITEM sDBConfItem)
{
	int iDelPos = CheckItemValidity(sDBConfItem);

	if (-1 == iDelPos)
	{
		return false;
	}

	else
	{
		m_vsDBConfItem.erase(m_vsDBConfItem.begin() + iDelPos);

		return true;
	}
	
}

bool CDBConfManage::CheckValidity()
{
	return true;
}

int CDBConfManage::CheckItemValidity(DBCONF_ITEM sDBConfItem)
{	
	for(int i=0; i<(int)(m_vsDBConfItem.size()); i++)
	{
		//if(m_vsDBConfItem[i].strDBName == sDBConfItem.strDBName && m_vsDBConfItem[i].strEnzymeName == sDBConfItem.strEnzymeName
				//&&m_vsDBConfItem[i].strPath == sDBConfItem.strPath && m_vsDBConfItem[i].strMetaName == sDBConfItem.strMetaName)
//		if(m_vsDBConfItem[i].strDBName == sDBConfItem.strDBName && m_vsDBConfItem[i].strEnzymeName == sDBConfItem.strEnzymeName)
//		{
//			return i;
//		}
		if(m_vsDBConfItem[i].strMetaName == sDBConfItem.strMetaName)
		{
			return i;
		}
	}

	return -1;
}

int CDBConfManage::CheckItemValidity(string strDBName, string strEnzyme,string strPath, string strMetaName)
{	
	for(int i=0; i<(int)(m_vsDBConfItem.size()); i++)
	{
		//if(m_vsDBConfItem[i].strDBName == strDBName && m_vsDBConfItem[i].strEnzymeName == strEnzyme&&m_vsDBConfItem[i].strPath == strPath && m_vsDBConfItem[i].strMetaName == strMetaName)
		if(m_vsDBConfItem[i].strDBName == strDBName && m_vsDBConfItem[i].strEnzymeName == strEnzyme)
		{
			return i;
		}
	}

	return -1;
}

DBCONF_HEAD CDBConfManage::GetDBConfHead()
{
	return m_sDBConfHead;
}

 bool CDBConfManage::GetDBConfItem(DBCONF_ITEM & sDBConfItem, string strDBName, string strEnzyme,string strPath, string strMetaName)
{
	int iPos = CheckItemValidity(strDBName, strEnzyme, strPath, strMetaName);

	if (-1 != iPos)
	{
		sDBConfItem = m_vsDBConfItem[iPos];
		return true;
	}

	return false;
}

vector<DBCONF_ITEM> CDBConfManage::GetDBConfAllItems()
{
	return m_vsDBConfItem;
}

string CDBConfManage::GetDBConfName()
{
	return m_strDBConfName;
}

bool CDBConfManage::Save()
{
	//Close();
	//Open("w");
	m_pDBConf = fopen(m_strDBConfName.c_str(), "w");

	fprintf(m_pDBConf, "[total]\n");
	//fprintf(m_pDBConf, "enzyme_list = %s\n", m_sDBConfHead.enzymeList.c_str());

	int iIndexNum = (int)(m_vsDBConfItem.size());
	fprintf(m_pDBConf,"index_num=%d\n", iIndexNum);

	char str[20];

	for(int i=0; i<iIndexNum; i++)
	{
//		itoa( i+1, str, 10);
		sprintf(str,"%d",i+1);
		string strIndex = "[index" + string(str) + "]";
		fprintf(m_pDBConf,"%s\n", strIndex.c_str());

		fprintf(m_pDBConf,"db_name=%s\n", m_vsDBConfItem[i].strDBName.c_str());
		fprintf(m_pDBConf,"enzyme=%s\n", m_vsDBConfItem[i].strEnzymeName.c_str());
		fprintf(m_pDBConf,"path=%s\n", m_vsDBConfItem[i].strPath.c_str());
		fprintf(m_pDBConf,"meta_name=%s\n", m_vsDBConfItem[i].strMetaName.c_str());
		fprintf(m_pDBConf,"remark=%s\n", m_vsDBConfItem[i].strRemark.c_str());
	}

	Close();
	
	return true;
}

