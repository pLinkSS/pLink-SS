#ifndef METACREATER_H_
#define METACREATER_H_

#include<stdlib.h>
#include<string>
#include "../include/sdk.h"
#include "Mass2PepIndex.h"

using namespace std;
//class std::string;

namespace Mass2PepIndex
{

class CMetaCreater
{
public:
	CMetaCreater();
	virtual ~CMetaCreater();	
	
	void	Load( const std::string strWorkPath, const std::string strDBName, const std::string strEnzymeName);
	
	void	SetMetaHead(META_HEAD stMetaHead);
	void	SetMetaItems(vector<META_ITEM > vstMetaItems);
	
	void	SetDBName(std::string strDBName);
	void	SetEnzymeName(std::string strEnzymeName);
	void	SetIdxNum(unsigned char tIdxNum);
	void	SetMissCleaves(unsigned char tMissCleaveSites);
	void	SetMultiplier(size_t	tMultiplier);
	void	SetRPepSQNum(size_t	tPepSQNum);
	void	SetRPepSQNum(size_t	tPepSQNum, size_t tCurrent);
	void	SetPepSQNum(size_t	tUniquePepSQNum);
	void	SetPepSQNum(size_t	tUniquePepSQNum, size_t tCurrent);
	void	SetPepMassNum(size_t	tUniquePepMassNum);
	void	SetPepMassNum(size_t	tUniquePepMassNum, size_t tCurrent);	
	void	SetMinPepMass(size_t tMinMass);
	void	SetMinPepMass(size_t tMinMass, size_t tCurrent);
	void	SetMaxPepMass(size_t tMaxMass);
	void	SetMaxPepMass(size_t tMaxMass, size_t tCurrent);
	void 	SetWorkDir(const char * szDestDir);
	void    WriteMetaItemMass(const int nIndex, const size_t tMassNum);
	void    WriteTotalMass(const size_t tMassNum);
	
	inline void SetIsMono(bool bMono)
	{
		m_stMetaHead.bMono = bMono;
	}
	
	inline void SetIsPep2Pro(bool bPep2Pro)
	{
		m_stMetaHead.bPep2Pro = bPep2Pro;
	}
	
	void	WriteMeta();
	
	inline	std::string	GetDBName() const
	{
		return std::string(m_stMetaHead.cDBName);
	}
	
	inline	std::string	GetEnzymeName()const
	{
		return std::string(m_stMetaHead.cEnzyme);
	}
	
	inline bool	GetIsMono()const
	{
		return m_stMetaHead.bMono;
	}

	inline std::string	GetWorkDir()const
	{
		return m_strWorkDir;
	}
	
	META_HEAD m_stMetaHead;
	vector<META_ITEM > m_vstMetaItems;
	
protected:
	
	void	_WriteMetaHead(FILE * fMeta);
	void	_WriteMetaItems(FILE * fMeta);	
	
	size_t	_GetMetaHeadSize(void)const;
	size_t	_GetMetaItemSize(void)const;
	
private:
	
	std::string m_strWorkDir;
	std::string m_strMetaName;
	std::string m_strMetaNameTxt;

	FILE * m_fMeta;
	FILE * m_fMetaTxt;
};

}

#endif /*METACREATER_H_*/
