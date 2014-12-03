#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "../include/sdk.h"
#include "../ProteomicsSDK/ProteomicsSDK.h"
#include "common.h"
#include "PepPairContainer.h"
#include "PepPairSorter.h" 
using namespace std;
using namespace proteomics_search;
using namespace proteomics_sdk;

CPepPairContainer::CPepPairContainer()
{
	m_tLength = 0;
	m_tBufLength = 0;
	m_strFileName = "";
	m_fp = NULL;
	m_pTrace = CTrace::GetInstance();

	// srand(time(NULL));
	ostringstream oss;
	oss<<(unsigned long)(this)<<".list";
	m_strFileName = oss.str();

	string strDebug = m_strFileName;
	strDebug += ": ";
	strDebug += "in CPepPairContainer Constructor.";
	m_pTrace->Debug(strDebug.c_str());

	m_buffer = new PEP_PAIR_ITEM[WRITE_BUFFER_SIZE];
}

CPepPairContainer::~CPepPairContainer()
{
	m_tLength = 0;
	m_tBufLength = 0;
	if(m_fp) {
		fclose(m_fp);
		m_fp = NULL;
	}
	/*
#ifdef WIN32
	if(0 == access(m_strFileName.c_str(), 0)) {
		remove(m_strFileName.c_str());
	}
#else
	m_strFileName = "./" + m_strFileName;
	if(0 == access(m_strFileName.c_str(), 0)) {
		remove(m_strFileName.c_str());
	}
#endif
*/
	remove(m_strFileName.c_str());
	if(m_buffer)
		delete [] m_buffer;
}

void CPepPairContainer::Init()
{
	if ( (m_fp = fopen(m_strFileName.c_str(),"wb")) == NULL)
	{
		CErrInfo info("CPepPairContainer", "Init()", "Failed to open.");
		throw runtime_error(info.Get().c_str());
	}
	m_tLength = 0;
	m_tBufLength = 0;
	
}

void CPepPairContainer::Attach(bool bPair,int tPepId1,int tPepId2,double lfMass)
{
	if(m_tBufLength == WRITE_BUFFER_SIZE)
	{
		_WriteToFile();
	}
	m_buffer[m_tBufLength].bPair = bPair;
	m_buffer[m_tBufLength].tPepId1 = tPepId1;
	m_buffer[m_tBufLength].tPepId2 = tPepId2;
	m_buffer[m_tBufLength].lfMass = lfMass;
	++m_tBufLength;
	++m_tLength;
}

void CPepPairContainer::EndAttach()
{
	_WriteToFile();
	if(m_fp) {
		fclose(m_fp);
		m_fp = NULL;
	}
}

void CPepPairContainer::Clear()
{
	m_tLength = 0;
	m_tBufLength = 0;
	
	if(m_fp) {
		fclose(m_fp);
		m_fp = NULL;
	}
	
	m_fp = fopen(m_strFileName.c_str(),"wb");
	
	if(m_fp) {
		fclose(m_fp);
		m_fp = NULL;
	}
		
}
 
void CPepPairContainer::Merge(CPepPairContainer * pList)
{
	char tmpch;
	
	CPepPairContainer tempContainer;
	
	tempContainer.Init();
	
	CPepPairContainerReader reader1,reader2;
	
	reader1.Init(this);
	
	reader2.Init(pList);
	
	reader1.Begin();
	
	reader2.Begin();
	
	PEP_PAIR_ITEM	item1,item2;
	bool b1,b2;
	
	b1 = reader1.GetNext(item1);
	b2 = reader2.GetNext(item2);
	while(1)
	{
		if(!b1 || !b2)
			break;
		if(item1.lfMass < item2.lfMass)
		{
			tempContainer.Attach(item1.bPair,item1.tPepId1,item1.tPepId2,item1.lfMass);
			b1 = reader1.GetNext(item1);
		}
		else
		{
			tempContainer.Attach(item2.bPair,item2.tPepId1,item2.tPepId2,item2.lfMass);
			b2 = reader2.GetNext(item2);
		}
	}
	
	
	
	while(b1)
	{
		tempContainer.Attach(item1.bPair,item1.tPepId1,item1.tPepId2,item1.lfMass);
		b1 = reader1.GetNext(item1);
	}
	
	while(b2)
	{
		tempContainer.Attach(item2.bPair,item2.tPepId1,item2.tPepId2,item2.lfMass);
		b2 = reader2.GetNext(item2);
	}
	
	tempContainer.EndAttach();
	Clear();
	pList->Clear();
	
	CPepPairContainerReader reader3;
	reader3.Init(&tempContainer);
	reader3.Begin();
	Init();
	while(reader3.GetNext(item1))
	{
		Attach(item1.bPair,item1.tPepId1,item1.tPepId2,item1.lfMass);
	}
	EndAttach();
}
	 

void CPepPairContainer::Print()
{
	CPepPairContainerReader reader;
	reader.Init(this);
	PEP_PAIR_ITEM item;
	while(reader.GetNext(item))
	{
		if(item.bPair)
			cout << item.lfMass << "	" << item.tPepId1 << "	" << item.tPepId2 << endl;
		else
			cout << item.lfMass << "	" << item.tPepId1 << endl;
			
	}
	
}

void CPepPairContainer::_WriteToFile()
{
	// todo: exception
	char tmpch;
	if(m_fp)
	{
		if(m_tBufLength!=fwrite(m_buffer,sizeof(PEP_PAIR_ITEM),m_tBufLength,m_fp))
		{
			CErrInfo info("CPepPairContainer", "_WriteToFile()", "Failed to write.");
			throw runtime_error(info.Get().c_str());
		}
	}
	m_tBufLength = 0;
}

void CPepPairContainer::Sort()
{
	CPepPairTwoWaySorter sorter(m_strFileName, m_tLength);
//	CPepPairMoreWaySorter sorter(m_strFileName, m_tLength);
	if(false == sorter.BatchDivision()) {
		_RemoveErrorFiles();
		CErrInfo info("CPepPairContainer", "Sort", "Batch division failed!");
		throw new runtime_error(info.Get().c_str());
	}
	if(false == sorter.BatchMergeSort()) {
		_RemoveErrorFiles();
		CErrInfo info("CPepPairContainer", "Sort", "Batch merge sort failed!");
		throw new runtime_error(info.Get().c_str());
	}
}

void CPepPairContainer::_RemoveErrorFiles()
{
	string strCommand = "";
#ifdef WIN32
	strCommand = "del *.list";
#else
	strCommand = "rm *.list";
#endif
	system(strCommand.c_str());
}

CPepPairContainerReader::CPepPairContainerReader()
{
	m_tCurId = 0;
	m_tCurSize = 0;
	m_tCurRead = 0;
	m_tTotalNum = 0;
	m_fp = NULL;
	m_buffer = new PEP_PAIR_ITEM[READ_BUFFER_SIZE];
}

CPepPairContainerReader::~CPepPairContainerReader()
{
	m_tCurId = 0;
	m_tCurSize = 0;
	m_tCurRead = 0;
	m_tTotalNum = 0;
	if(m_fp) {
		fclose(m_fp);
		m_fp = NULL;
	}
	if(m_buffer)
		delete [] m_buffer;
}
	
void CPepPairContainerReader::Init(CPepPairContainer * pList)
{
	m_tCurId = 0;
	m_tCurSize = 0;
	m_tCurRead = 0;
	m_tTotalNum = pList->m_tLength;
	if( (m_fp = fopen(pList->m_strFileName.c_str(),"rb")) == NULL )
	{
		CErrInfo info("CPepPairContainerReader", "Init()", "Failed to open.");
		throw runtime_error(info.Get().c_str());
	}
}

void CPepPairContainerReader::Begin()
{
	m_tCurId = 0;
	m_tCurSize = 0;
	m_tCurRead = 0;
	// todo : exception
	if(m_fp)
		rewind(m_fp);
}

bool CPepPairContainerReader::GetNext(struct PEP_PAIR_ITEM & stPepPair)
{
	if(m_tCurId == m_tCurSize)
	{
		_ReadFromFile();
	}
	
	if(!m_tCurSize)
		return false;
	
	stPepPair = m_buffer[m_tCurId];
	++m_tCurId;
	++m_tCurRead;
	
	return true;
}

void CPepPairContainerReader::_ReadFromFile()
{
	size_t tReadNum;
	tReadNum =  m_tTotalNum - m_tCurRead;
	if(tReadNum > READ_BUFFER_SIZE)
	{
		tReadNum = READ_BUFFER_SIZE;
	}
	
	m_tCurId = 0;
	if(tReadNum <= 0)
	{
		m_tCurSize = 0;
	}
	else
	{
		m_tCurSize = tReadNum;
		if(m_fp)
		{
			if(m_tCurSize != fread(m_buffer,sizeof(PEP_PAIR_ITEM),m_tCurSize,m_fp))
			{
				CErrInfo info("CPepPairContainerReader", "_ReadFromFile()", "Failed to read.");
				throw runtime_error(info.Get().c_str());
			}
		}
	}
}
	

