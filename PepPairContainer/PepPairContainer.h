#ifndef PEPPAIRCONTAINER_H_
#define PEPPAIRCONTAINER_H_

#include "common.h"
#include "../include/sdk.h"
using namespace proteomics_search;
using namespace proteomics_sdk;

namespace proteomics_search{

class CPepPairContainer
{	
	friend class CPepPairContainerReader;

public:
	CPepPairContainer();
	virtual ~CPepPairContainer();
	
	/*
	 * construction of the list
	 */
	
	// initialize the list
	void Init();
	// attach one pair at the end of the list
	void Attach(bool bPair,int tPepId1,int tPepId2,double lfMass);
	// finish attaching 
	void EndAttach();
	// destroy the list
	void Clear();
	// merge pList into this list , pList will vanish
	// todo: bug still remains
	void Merge(CPepPairContainer * pList);
	// sort the list
	void Sort();
	 
	void Print();
	
	size_t GetLength()
	{
		return m_tLength;
	}
protected :
	void _WriteToFile();
	void _RemoveErrorFiles();
protected :
	
	size_t m_tLength;
	size_t m_tBufLength;
	struct PEP_PAIR_ITEM * m_buffer;
	string m_strFileName;
	FILE * m_fp;

	CTrace *m_pTrace;
};

class CPepPairContainerReader
{
public:
	CPepPairContainerReader();
	~CPepPairContainerReader();
	
	void Init(CPepPairContainer * pList);
	void Begin();
	bool GetNext(struct PEP_PAIR_ITEM & stPepPair);
	
protected:
	
	void _ReadFromFile();
	
protected:
	struct PEP_PAIR_ITEM * m_buffer;
	size_t m_tCurId;
	size_t m_tCurSize;
	size_t m_tCurRead;
	size_t m_tTotalNum;
	FILE * m_fp;
};

}

#endif /*PEPPAIRCONTAINER_H_*/
