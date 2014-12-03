#ifndef MYSTACK_H_
#define MYSTACK_H_

#include <stdio.h>
//using namespace proteomics_search;

namespace proteomics_search{

template <class ELEMENT_TYPE>
struct CMyStack_Node
{
	CMyStack_Node(){
		m_pNext=NULL;
	}
	CMyStack_Node(ELEMENT_TYPE data){
		m_data=data;
		m_pNext=NULL;
	}
	ELEMENT_TYPE m_data;
	struct CMyStack_Node<ELEMENT_TYPE> * m_pNext;
	
};

template <class ELEMENT_TYPE>
class CMyStack
{
public:
	CMyStack(void);
	bool Push(ELEMENT_TYPE data);
	bool Pop(ELEMENT_TYPE & data);
	bool Init();
	bool Destroy();
	int GetDepth(){
		return m_depth;
	}
	~CMyStack(void);

protected:
	CMyStack_Node<ELEMENT_TYPE> * m_pTopEle;
	int m_depth;

};

template <class ELEMENT_TYPE>
CMyStack<ELEMENT_TYPE>::CMyStack(void)
{
	m_pTopEle=NULL;
	m_depth=0;
}

template <class ELEMENT_TYPE>
CMyStack<ELEMENT_TYPE>::~CMyStack(void)
{
	Destroy();
}

template <class ELEMENT_TYPE>
bool CMyStack<ELEMENT_TYPE>::Push(ELEMENT_TYPE data)
{
	struct CMyStack_Node<ELEMENT_TYPE> * pNewNode = new struct CMyStack_Node<ELEMENT_TYPE>(data);
	if(pNewNode==NULL)
		return false;
	pNewNode->m_pNext = this->m_pTopEle;
	m_pTopEle = pNewNode;
	m_depth++;
	return true;
}

template <class ELEMENT_TYPE>
bool CMyStack<ELEMENT_TYPE>::Pop(ELEMENT_TYPE & data)
{
	if(m_pTopEle==NULL)
		return false;
	struct CMyStack_Node<ELEMENT_TYPE> * pDelNode;
	pDelNode = m_pTopEle;
	m_pTopEle = m_pTopEle->m_pNext;
    data = pDelNode->m_data;
	delete pDelNode;
	m_depth--;
	return true;
}

template <class ELEMENT_TYPE>
bool CMyStack<ELEMENT_TYPE>::Init()
{
	Destroy();
	return true;
}

template <class ELEMENT_TYPE>
bool CMyStack<ELEMENT_TYPE>::Destroy()
{
	struct CMyStack_Node<ELEMENT_TYPE> * pTmpNode;
	if(this->m_pTopEle){
		pTmpNode=m_pTopEle;
		m_pTopEle=m_pTopEle->m_pNext;
		delete pTmpNode;
	}
	m_pTopEle=NULL;
	m_depth=0;
	return true;
}

}



#endif /*MYSTACK_H_*/
