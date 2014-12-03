#include <stdio.h>
#include "../include/sdk.h"
#include "PepPairList.h"
#include <iostream>
#include <string>
using namespace std;
using namespace proteomics_search;
using namespace proteomics_sdk;

CPepPairList::CPepPairList()
{
	m_pHead = NULL;
	m_pTail = NULL;
	m_tLength = 0;
}

CPepPairList::~CPepPairList()
{
	Clear();
}


// initialize the list
void CPepPairList::Init()
{
	Clear();
}

// attach one pair at the end of the list
void CPepPairList::Attach(bool bPair,int tPepId1,int tPepId2,double lfMass)
{
	struct PEP_PAIR * pPepPair;
	pPepPair = new struct PEP_PAIR();
	
	if(pPepPair == NULL)
	{
		CErrInfo info("CPepPairList", "Attach", "fail");
		throw runtime_error(info.Get().c_str());
	}
	
	pPepPair->bPair = bPair;
	pPepPair->tPepId1 = tPepId1;
	pPepPair->tPepId2 = tPepId2;
	pPepPair->lfMass = lfMass;
	pPepPair->pNextPepPair = NULL;
	
	m_tLength ++ ;
	if(m_pTail)
	{
		m_pTail->pNextPepPair = pPepPair;
		m_pTail = pPepPair;
	}
	else
	{
		m_pHead = m_pTail = pPepPair;
	}
}

// destroy the list
void CPepPairList::Clear()
{
	struct PEP_PAIR * pPepPair;
	while(m_pHead)
	{
		pPepPair = m_pHead;
		m_pHead = m_pHead->pNextPepPair;
		delete pPepPair;
	}
	m_pHead = NULL;
	m_pTail = NULL;
	m_tLength = 0;
}

void CPepPairList::Print()
{
	string f = "peppair.txt";
	ofstream os(f.c_str(),ios::app);

	struct PEP_PAIR * pCurPair;
	pCurPair = m_pHead;
	while(pCurPair)
	{
		os << pCurPair->lfMass ;
		if(pCurPair->bPair)
		{
			os << "	" << pCurPair->tPepId1 << "	" << pCurPair->tPepId2 << endl;
		}
		else
		{
			os << "	" << pCurPair->tPepId1 << endl;
		}
		pCurPair = pCurPair->pNextPepPair;
	}
	os << "length = " << m_tLength;
	
	os.close();
}

// merge pList into this list , pList will vanish 
void CPepPairList::Merge(CPepPairList * pList)
{
	if(pList == NULL)
		return ;
	
	m_tLength += pList->m_tLength;
	
	struct PEP_PAIR * pLeftPair , * pRightPair, * pPreLeftPair, * pNextRightPair;
	pLeftPair = m_pHead;
	pRightPair = pList->m_pHead;
	pPreLeftPair = NULL;
	pNextRightPair = NULL;
	
	while(pLeftPair && pRightPair)
	{
		if(pLeftPair->lfMass < pRightPair->lfMass)
		{
			pPreLeftPair = pLeftPair;
			pLeftPair = pLeftPair->pNextPepPair;
		}
		else
		{
			pNextRightPair = pRightPair->pNextPepPair;
			if(pPreLeftPair)
			{
				pRightPair->pNextPepPair = pPreLeftPair->pNextPepPair;
				pPreLeftPair->pNextPepPair = pRightPair;
				
				
				 
			}
			else
			{
				pRightPair->pNextPepPair = m_pHead;
				m_pHead = pRightPair; 
			}
			
			pPreLeftPair = pRightPair;
			pRightPair = pNextRightPair;
		}
	}
	
	if(pRightPair)
	{
		if(pPreLeftPair)
		{
			pPreLeftPair->pNextPepPair = pRightPair;
		}
		else
		{
			m_pHead = pRightPair; 
		}
		
		m_pTail = pList->m_pTail;
	}
	
	pList->m_pHead = NULL;
	pList->Clear();
}
 
/*
 * list reader 
 */

CPepPairListReader::CPepPairListReader()
{
	m_pHead = NULL;
	m_pCurrentPair = NULL;
}

CPepPairListReader::~CPepPairListReader()
{
	m_pHead = NULL;
	m_pCurrentPair = NULL;
}

void CPepPairListReader::Init(CPepPairList * pList)
{
	m_pHead = m_pCurrentPair = pList->m_pHead;
}

void CPepPairListReader::Begin()
{
	m_pCurrentPair = m_pHead;
}

bool CPepPairListReader::GetNext(struct PEP_PAIR & stPepPair)
{
	if(m_pCurrentPair)
	{
		stPepPair = *m_pCurrentPair;
		m_pCurrentPair = m_pCurrentPair->pNextPepPair;
		return true;
	}
	else
	{
		return false;
	}
}

