#ifndef PEPPAIRLIST_H_
#define PEPPAIRLIST_H_

namespace proteomics_search{
struct PEP_PAIR
{
	bool bPair;
	int tPepId1;
	int tPepId2;
	double	lfMass;
	struct PEP_PAIR * pNextPepPair;
	
	PEP_PAIR()
	{
		bPair = true;
		tPepId1 = 0;
		tPepId2 = 0;
		lfMass = 0.0;
		pNextPepPair = NULL;
	}
};

class CPepPairList
{
	friend class CPepPairListReader;
	
public:
	CPepPairList();
	virtual ~CPepPairList();
	
	/*
	 * construction of the list
	 */
	
	// initialize the list
	void Init();
	// attach one pair at the end of the list
	void Attach(bool bPair,int tPepId1,int tPepId2,double lfMass);
	// destroy the list
	void Clear();
	// merge pList into this list , pList will vanish 
	void Merge(CPepPairList * pList);
	 
	void Print();
protected :
	
	size_t m_tLength;
	struct PEP_PAIR * m_pHead;
	struct PEP_PAIR * m_pTail;

};

/*
 * list reader 
 */

class CPepPairListReader
{
public:
	CPepPairListReader();
	~CPepPairListReader();
	
	void Init(CPepPairList * pList);
	void Begin();
	bool GetNext(struct PEP_PAIR & stPepPair);
	
protected:
	struct PEP_PAIR * m_pHead;
	struct PEP_PAIR * m_pCurrentPair;
	
};
}

#endif /*PEPPAIRLIST_H_*/
