#ifndef METEREADER_H_
#define METEREADER_H_

//#include <string>
//#include <vector>
#include "../include/sdk.h"
#include "Mass2PepIndex.h"

using namespace std;
class sting;

namespace Mass2PepIndex
{

class CMetaReader
{
public:
	CMetaReader();
	virtual ~CMetaReader();
	
	void	Load( const string strWorkPath, const string strMetaName);
	
	void	OpenFile(void);
	void	CloseFile(void);

	inline const string & GetMetaName(void) const 
	{
		return m_strMetaName;
	}
	
	inline const META_HEAD & GetMetaHead(void) const 
	{
		return m_stMetaHead;
	}
	
	inline const vector<META_ITEM> & GetMetaItems(void) const 
	{
		return m_vstMetaItems;
	}
	
	inline string GetWorkDir(void) const 
	{
		return m_strWorkDir;
	}

	void	SetMetaName(string strMetaName);
	void	SetMetaHead();
	void	SetMetaItems();	
	void 	SetWorkDir(const char * szDestDir);

protected:
	
	size_t	_GetMetaHeadSize(void)const;
	size_t	_GetMetaItemSize(void)const;
	
	META_HEAD m_stMetaHead;
	vector<META_ITEM > m_vstMetaItems;

private:
	string m_strWorkDir;
	string m_strMetaName;
	FILE * m_fMeta;
};

}

#endif /*METEREADER_H_*/
