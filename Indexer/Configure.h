#ifndef CONFIGURE_H_
#define CONFIGURE_H_

#include "../include/sdk.h"
using namespace Mass2PepIndex;


class CConfigure
{
public:
	CConfigure();
	virtual ~CConfigure();
	
	void Init(std::string strPIndexName);
	void ReadHead();
	void ReadItems();
	
	string m_strPIndexName;
	PINDEX_HEAD m_PIndexHead;
	vector<PINDEX_ITEM> m_PIndexItems;
};

#endif /*CONFIGURE_H_*/
