#ifdef WIN32
#include <direct.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <utility>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include "common.h"
#include "MyStack.h"
#include "PepPairSorter.h"
#include "../include/sdk.h"

using namespace std;
using namespace proteomics_sdk;
using namespace proteomics_search;

CPepPairSorter::CPepPairSorter()
{
	m_strFileName = "";
	m_tDataNum = 0;
}

CPepPairSorter::CPepPairSorter(string strFile, size_t tDataNum)
{
	m_strFileName = strFile;
	m_tDataNum = tDataNum;
}

void CPepPairSorter::SetDataInfo(string strFile, size_t tDataNum)
{
	m_strFileName = strFile;
	m_tDataNum = tDataNum;
}

CPepPairMoreWaySorter::CPepPairMoreWaySorter()
{
	m_pTrace = CTrace::GetInstance();
}

CPepPairMoreWaySorter::CPepPairMoreWaySorter(string strFile, size_t tDataNum)
: CPepPairSorter(strFile, tDataNum)
{
	m_pTrace = CTrace::GetInstance();
}

CPepPairMoreWaySorter::~CPepPairMoreWaySorter()
{
	if(m_arrData)
		delete []m_arrData;
}

bool CPepPairMoreWaySorter::BatchDivision()
{
	m_pTrace->Info("BatchDivision enter.");
	m_pTrace->Debug(m_strFileName);

	if(m_tDataNum == 0)
		return false;

	/* initialization */
	m_vTmpFileName.clear();
	m_tTmpFileNum = int(m_tDataNum / MAX_SORT_UNIT);
	if(m_tDataNum % MAX_SORT_UNIT)
		m_tTmpFileNum++;
	m_arrData = new PEP_PAIR_ITEM[MAX_SORT_UNIT];

	FILE *fpIn = NULL, *fpOut = NULL;
	if(!(fpIn = fopen(m_strFileName.c_str(), "rb"))) {
		CErrInfo info("CPepPairSorter", "BatchDivision", "Failed to open m_strFileName!");
		throw runtime_error(info.Get());
	}

	size_t tLeftDataNum = m_tDataNum;
	size_t tCurDataNum;
	srand(time(0));
	string strRand;
	ostringstream oss;
	oss<<(unsigned long)this;
	strRand = oss.str();
	for(unsigned int i = 0; i < m_tTmpFileNum; i++) {
		/* randomly generate temporary files names */
		string strFileName = "";
		ostringstream oss;
		oss<<i<<"_"<<strRand.c_str()<<".tmp";
		strFileName = oss.str();
		m_pTrace->Debug(strFileName);

		/* store files names */
		m_vTmpFileName.push_back(strFileName);

		tCurDataNum = (tLeftDataNum < MAX_SORT_UNIT) ? tLeftDataNum : MAX_SORT_UNIT;
		if(tCurDataNum != fread(m_arrData, sizeof(PEP_PAIR_ITEM), tCurDataNum, fpIn)) {
			CErrInfo info("CPepPairSorter", "BatchDivision", "Failed to read m_strFileName!");
			throw runtime_error(info.Get());
		}

		/* sorting */
		sort(m_arrData, m_arrData+tCurDataNum);

		/* open file */
		fpOut = fopen(strFileName.c_str(), "wb");

		/* write sorted data */
		if(tCurDataNum != fwrite(m_arrData, sizeof(PEP_PAIR_ITEM), tCurDataNum, fpOut)) {
			CErrInfo info("CPepPairSorter", "BatchDivision", "Failed to write sorted data!");
			throw runtime_error(info.Get());
		}

		/* close write file */
		fclose(fpOut);
		fpOut = NULL;
		tLeftDataNum -= tCurDataNum;
	}
	fclose(fpIn);
	fpIn = NULL;
	m_pTrace->Info("BatchDivision exit.");
	return true;
}

bool CPepPairMoreWaySorter::BatchMergeSort()
{
	m_pTrace->Info("BatchMergeSort enter.");
	m_pTrace->Debug(m_strFileName);
	LoserTree loserTree(m_tTmpFileNum);
	loserTree.ConnectToChunks(m_vTmpFileName, m_strFileName);
	m_pTrace->Debug("Initialize LoserTree...");
	loserTree.CreateLoserTree();
	m_pTrace->Debug("External sorting, merging...");
	loserTree.Merge();

	m_pTrace->Info("Remove temporary files in BatchMergeSort().");
	for(int i = 0; i < m_tTmpFileNum; i++)
		remove(m_vTmpFileName[i].c_str());
	m_pTrace->Info("BatchMergeSort exit.");
	return true;
}

CPepPairMoreWaySorter::LoserTree::LoserTree(const size_t tWay)
{
	m_tChunkNum = tWay;
	if(m_tChunkNum > 128) { /* maximum limits of file descriptor for one time */
		CErrInfo info("CPepPairMoreWaySorter::LoserTree", "LoserTree", "File descriptor error!");
		throw runtime_error(info.Get());
	}
	m_fpIn = new FILE*[m_tChunkNum];
	for(unsigned int i = 0; i < m_tChunkNum; i++)
		m_fpIn[i] = NULL;
	m_fpOut = NULL;

	m_arrBuf = new PEP_PAIR_ITEM[m_tChunkNum+1];
	m_nTree = new int[m_tChunkNum];
}

CPepPairMoreWaySorter::LoserTree::~LoserTree()
{
	if(m_fpIn) {
		for(unsigned int i = 0; i < m_tChunkNum; i++)
			if(m_fpIn[i]) {
				fclose(m_fpIn[i]);
				m_fpIn[i] = NULL;
			}
		delete []m_fpIn;
		m_fpIn = NULL;
	}
	if(m_fpOut) {
		fclose(m_fpOut);
		m_fpOut = NULL;
	}
	delete []m_arrBuf;
	delete []m_nTree;
}

void CPepPairMoreWaySorter::LoserTree::ConnectToChunks(const vector<string> &vInFileName,
		const string &strOutFileName)
{
	if(m_tChunkNum != vInFileName.size()) {
		CErrInfo info("CPepPairMoreWaySorter::LoserTree", "ConnectToChunks", "Not agree in file number!");
		throw runtime_error(info.Get());
	}
	for(unsigned int i = 0; i < m_tChunkNum; i++) {
		m_fpIn[i] = fopen(vInFileName[i].c_str(), "rb");
		if(!m_fpIn[i]) {
			CErrInfo info("CPepPairMoreWaySorter::LoserTree", "ConnectToChunks", "Failed to open file!");
			throw runtime_error(info.Get());
		}
		if(1 != fread(m_arrBuf+i, sizeof(PEP_PAIR_ITEM), 1, m_fpIn[i])) {
			CErrInfo info("CPepPairMoreWaySorter::LoserTree", "ConnectToChunks", "Failed to read item!");
			throw runtime_error(info.Get());
		}
	}
	m_fpOut = NULL;
	m_fpOut = fopen(strOutFileName.c_str(), "wb");
}

void CPepPairMoreWaySorter::LoserTree::CreateLoserTree()
{
	int i;

	m_arrBuf[m_tChunkNum].lfMass = MIN_MASS;
	for(i = 0; i < m_tChunkNum; i++)
		m_nTree[i] = m_tChunkNum;
	for(i = m_tChunkNum-1; i >= 0; i--)
		_Adjust(i);
}

void CPepPairMoreWaySorter::LoserTree::Merge()
{
	int nWay;

	do {
		nWay = m_nTree[0];
		if(1 != fwrite(m_arrBuf+nWay, sizeof(PEP_PAIR_ITEM), 1, m_fpOut)) {
			CErrInfo info("CPepPairMoreWaySorter::LoserTree", "Merge()", "Failed to write item!");
			throw runtime_error(info.Get());
		}

		if(NULL != m_fpIn[nWay]) {
			if(1 != fread(m_arrBuf+nWay, sizeof(PEP_PAIR_ITEM), 1, m_fpIn[nWay])) {
				m_arrBuf[nWay].lfMass = MAX_MASS;
				fclose(m_fpIn[nWay]);
				m_fpIn[nWay] = NULL;
			}
		}
		_Adjust(nWay);
	} while(m_arrBuf[m_nTree[0]].lfMass != MAX_MASS);

	fclose(m_fpOut);
	m_fpOut = NULL;
}

void CPepPairMoreWaySorter::LoserTree::_Adjust(int nItem)
{
	int nTmp;
	int nPre = ((int)m_tChunkNum + nItem) / 2;
	while(nPre > 0) {
		if(m_arrBuf[nItem].lfMass > m_arrBuf[m_nTree[nPre]].lfMass) {
			nTmp = nItem;
			nItem = m_nTree[nPre];
			m_nTree[nPre] = nTmp;
		}
		nPre /= 2;
	}
	m_nTree[0] = nItem;
}

/*
int main()
{
	time_t start, end;
	start = clock();
	CPepPairMoreWaySorter sorter("1859043368_1.list", 3752429);
	sorter.BatchDivision();
	sorter.BatchMergeSort();
	end = clock();
	cout<<"Time consumption: "<<(end-start)*1.0/CLOCKS_PER_SEC<<endl;

	start = clock();
	CPepPairTwoWaySorter sorter2("1859043368_2.list", 3752429);
	sorter2.BatchDivision();
	sorter2.BatchMergeSort();
	end = clock();
	cout<<"Time consumption: "<<(end-start)*1.0/CLOCKS_PER_SEC<<endl;
	return 0;
}
*/

CPepPairTwoWaySorter::CPepPairTwoWaySorter()
{
	m_Data = NULL;
	m_tCurDataNum = 0;
	m_nTmpFileNum = 0;
	m_strTmpFileName = NULL;
	m_pTrace = CTrace::GetInstance();
	srand(time(NULL));
}

CPepPairTwoWaySorter::CPepPairTwoWaySorter(string strFile, size_t tDataNum)
: CPepPairSorter(strFile, tDataNum)
{
	m_Data = NULL;
	m_tCurDataNum = 0;
	m_nTmpFileNum = 0;
	m_strTmpFileName = NULL;
	m_pTrace = CTrace::GetInstance();
	srand(time(NULL));
}

CPepPairTwoWaySorter::~CPepPairTwoWaySorter()
{
	if(m_Data)
		delete []m_Data;

	for(int i = 0;i < m_nTmpFileNum; i++){
		delete [] m_strTmpFileName[i];
		m_strTmpFileName[i] = NULL;
	}
	delete []m_strTmpFileName;
}

bool CPepPairTwoWaySorter::BatchDivision()
{
	if(false == BatchReadFile())
	{
		m_pTrace->Alert("batch read file fail");
		m_pTrace->Break();
		return false;
	}
	return true;
}
bool CPepPairTwoWaySorter::BatchMergeSort()
{
	if(false == BatchSort())
	{
		m_pTrace->Alert("batch sort file fail");
		m_pTrace->Break();
		return false;
	}
	if(false == BatchWriteFile())
	{
		m_pTrace->Alert("batch write file fail");
		m_pTrace->Break();
		return false;
	}
	return true;
}

void CPepPairTwoWaySorter::ClearTmpFile()
{
	for(int i=0;i<m_nTmpFileNum;i++)
		remove(m_strTmpFileName[i]);
}

bool CPepPairTwoWaySorter::BatchReadFile()
{
	FILE * fp;

	m_pTrace->Debug(m_strFileName.c_str());
	if(m_strFileName.empty() || (fp = fopen(m_strFileName.c_str(),"rb")) == NULL){
		// can't open input file
		string info = "Can't open input file: " + m_strFileName;
		m_pTrace->Info(info.c_str());
		return false;
	}

	m_nTmpFileNum = 0;
	//存储临时文件名的临时空间
	char szBuf[256];
	//存储当前时间
	time_t tTmptime;
	//临时文件指针
	FILE * fpTmpFile=NULL;
	//临时存放数据
	PEP_PAIR_ITEM * pTmpData = NULL;
	pTmpData = new PEP_PAIR_ITEM[MAX_SORT_UNIT];


	if(pTmpData == NULL)
	{
		// allocate space failure
		m_pTrace->Info("Allocate space failure");
		m_pTrace->Break();
		return false;
	}

	//归并段数目
	m_nTmpFileNum = int(m_tDataNum / MAX_SORT_UNIT);
	if(m_tDataNum % MAX_SORT_UNIT)
		m_nTmpFileNum ++ ;

	m_strTmpFileName = new char *[m_nTmpFileNum];
	if(!m_strTmpFileName)
	{
		delete [] pTmpData;
		fclose(fp);
		m_pTrace->Info("Too many temporary files");
		m_pTrace->Break();
		return false;
	}
	for(int i = 0;i < m_nTmpFileNum; ++i)
		m_strTmpFileName[i] = NULL;

	string strInfo = "";
	ostringstream oss;
	oss<<"Temporary file number = "<<m_nTmpFileNum;
	strInfo += oss.str();
	m_pTrace->Info(strInfo.c_str());

	size_t tLeftDataNum = m_tDataNum;
	size_t tCurDataNum;
	for(int i = 0;i < m_nTmpFileNum; ++i)
	{
		if(tLeftDataNum < MAX_SORT_UNIT)
		{
			tCurDataNum = tLeftDataNum;
		}
		else
		{
			tCurDataNum = MAX_SORT_UNIT;
		}

		if(tCurDataNum != fread(pTmpData,sizeof(PEP_PAIR_ITEM),tCurDataNum,fp))
		{
			delete [] pTmpData;
			fclose(fp);
			m_pTrace->Info("fread fail.");
			m_pTrace->Break();
			return false;
		}

		// time(&tTmptime);
		sprintf(szBuf,"%d_%d%d.tmp",i, rand(), rand());
		m_strTmpFileName[i]=new char[strlen(szBuf)+1];
		strcpy(m_strTmpFileName[i],szBuf);

		if(fpTmpFile){
			fclose(fpTmpFile);
			fpTmpFile = NULL;
		}

		if((fpTmpFile=fopen(m_strTmpFileName[i],"wb"))==NULL)
		{
			delete [] pTmpData;
			fclose(fp);

			string strInfo = "";
			ostringstream oss;
			oss<<"Open temporary file "<< m_strTmpFileName[i] << " failed.";
			strInfo += oss.str();
			m_pTrace->Info(strInfo);
			m_pTrace->Break();

			return false;
		}

		if(tCurDataNum!=fwrite(pTmpData,sizeof(PEP_PAIR_ITEM),tCurDataNum,fpTmpFile)){

			delete [] pTmpData;
			fclose(fp);
			if(fpTmpFile)
				fclose(fpTmpFile);

			string strInfo = "";
			ostringstream oss;
			oss<<"Write temporary file "<< m_strTmpFileName[i] << " failed.";
			strInfo += oss.str();
			m_pTrace->Info(strInfo);
			m_pTrace->Break();

			return false;

		}
		tLeftDataNum -= tCurDataNum;
	}

	delete []pTmpData;
	pTmpData = NULL;

	fclose(fp);
	fclose(fpTmpFile);

	return true;

}

//批量排序-排序部分
//排序次数：m_nTmpFileNum
//归并次数：m_nTmpFileNum-1
//流程：
//1. 读入第i个子文件
//2. 快排
//3. 与第i-1个子文件归并到第i个子文件中
//程序结束后，结果在最后一个子文件中

bool CPepPairTwoWaySorter::BatchSort()
{
	FILE *fp, *fp1;
	int nSubDataNum;
	for(int i = 0; i < m_nTmpFileNum; i++)
	{
		//读入第i个子文件
		if((fp=fopen(m_strTmpFileName[i],"rb")) == NULL){
			string strInfo = "Failed to open file: ";
			strInfo += m_strTmpFileName[i];
			m_pTrace->Info(strInfo);
			m_pTrace->Break();

			if(m_Data)
				delete [] m_Data;
			m_Data = NULL;
			return false;
		}

		//计算读入数据量
		if(fseek(fp,0,SEEK_END)){
			m_pTrace->Info("Reach end of file.");
			m_pTrace->Break();
			fclose(fp);

			if(m_Data)
				delete [] m_Data;
			m_Data = NULL;
			return false;
		}
		nSubDataNum = ftell(fp)/sizeof(PEP_PAIR_ITEM);

		//文件指针复位, 指向文件开始
		rewind(fp);

		//读入第一个文件时分配数组空间
		if(i == 0){
			m_Data = new PEP_PAIR_ITEM[nSubDataNum];
		}

		//读入数据
		if(nSubDataNum!=fread(m_Data,sizeof(PEP_PAIR_ITEM),nSubDataNum,fp))
		{
			m_pTrace->Info("Read data failed.");
			m_pTrace->Break();
			fclose(fp);

			if(m_Data)
				delete [] m_Data;
			m_Data = NULL;
			return false;
		}

		if(fp) {
			fclose(fp);
			fp = NULL;
		}

		//排序阶段, 快排
		m_tCurDataNum = nSubDataNum;
		QuickSort_NoneRecursive();

		//排序阶段, 归并
		if(i >= 1)
		{
			//打开第i-1个子文件
			if((fp=fopen(m_strTmpFileName[i-1],"rb"))==NULL){
				string strInfo = "Failed to open file(rb): ";
				strInfo += m_strTmpFileName[i-1];
				m_pTrace->Info(strInfo);
				m_pTrace->Break();

				if(m_Data)
					delete [] m_Data;
				return false;
			}

			// 以写的方式，打开第i个文件
			/* 这里的代码可以提速, 将前一个文件排好序的结果也一次读入内存, 在内存中归并 : by Jeremy*/
			if((fp1=fopen(m_strTmpFileName[i],"wb")) == NULL){
				string strInfo = "Failed to open file(wb): ";
				strInfo += m_strTmpFileName[i];
				m_pTrace->Info(strInfo);
				m_pTrace->Break();

				fclose(fp);
				if(m_Data)
					delete [] m_Data;
				m_Data = NULL;
				return false;
			}

			int j = 0;
			PEP_PAIR_ITEM tPepPair1,tPepPair2;

			int nReadRet = fread(&tPepPair1, sizeof(PEP_PAIR_ITEM), 1, fp);

			while(j < nSubDataNum && nReadRet == 1){
				if(m_Data[j].lfMass<=tPepPair1.lfMass){
					tPepPair2 = m_Data[j];
					j++;
				}
				else{
					tPepPair2 = tPepPair1;
					nReadRet = fread(&tPepPair1,sizeof(PEP_PAIR_ITEM),1,fp);
				}
				fwrite(&tPepPair2,sizeof(PEP_PAIR_ITEM),1,fp1);
			}

			if(j >= nSubDataNum)
			{
				//如果是数组中的数据输出完毕
				while( nReadRet == 1 )
				{
					nReadRet = fread(&tPepPair1,sizeof(PEP_PAIR_ITEM),1,fp);
					fwrite(&tPepPair1,sizeof(PEP_PAIR_ITEM),1,fp1);
				}
			}
			else
			{
				fwrite(&m_Data[j],sizeof(PEP_PAIR_ITEM),nSubDataNum-j,fp1);
			}


			fclose(fp);
			fclose(fp1);

			//删除上一个临时文件
			string strDel = "";
#ifdef WIN32
			strDel = m_strTmpFileName[i-1];
#else
			strDel += "./";
			strDel += m_strTmpFileName[i-1];
#endif
			remove(strDel.c_str());
		}
		else
		{
			//如果当前正在处理第0个子文件，则直接将内存中数组的排序结果输出到第0个子文件中

			//以写的方式，打开第i个文件
			if((fp1=fopen(m_strTmpFileName[0],"wb"))==NULL){
				//打开子文件失败
				delete [] m_Data;
				m_Data = NULL;
				return false;
			}

			fwrite(m_Data,sizeof(PEP_PAIR_ITEM),nSubDataNum,fp1);

			fclose(fp1);
		}
	}

	delete [] m_Data;
	m_Data = NULL;
	return true;
}


//批量排序-输出文件
bool CPepPairTwoWaySorter::BatchWriteFile(){

	FILE * fp, * fp1;
	if(m_nTmpFileNum<=0)
		return true;

	if((fp1 = fopen(m_strTmpFileName[m_nTmpFileNum-1],"rb"))==NULL ){
		return false;
	}

	if((fp=fopen(m_strFileName.c_str(),"wb"))==NULL){
		return false;
	}

	PEP_PAIR_ITEM tPepPair;
	for(int i = 0;i < m_tDataNum; i++)
	{
		if(1!=fread(&tPepPair,sizeof(PEP_PAIR_ITEM),1,fp1)){
			fclose(fp1);
			fclose(fp);
			return false;
		}

		if(1!=fwrite(&tPepPair,sizeof(PEP_PAIR_ITEM),1,fp)){
			fclose(fp1);
			fclose(fp);
			return false;
		}
	}

	fclose(fp);
	fclose(fp1);

	// delete the last tmp file
	string strDel = "";
#ifdef WIN32
	strDel = m_strTmpFileName[m_nTmpFileNum-1];
#else
	strDel += "./";
	strDel += m_strTmpFileName[m_nTmpFileNum-1];
#endif
	remove(strDel.c_str());

	return true;
}

void CPepPairTwoWaySorter::bubblesort(int low,int high)
{
	bool bChanged;
	PEP_PAIR_ITEM tPepPair;

	for(int i=0;i<=high-low;i++){
		bChanged=false;
		for(int j=low;j<=high-i-1;j++)
		{
			if(m_Data[j].lfMass>m_Data[j+1].lfMass)
			{
				bChanged=true;
				PEP_PAIR_ITEM tPepPair=m_Data[j+1];
				m_Data[j+1]=m_Data[j];
				m_Data[j]=tPepPair;
			}
		}
		if(bChanged==false){
			break;
		}
	}
}

void CPepPairTwoWaySorter::QuickSort_NoneRecursive()
{
	//非递归的快排函数
	//当待排序数组长度等于3时，用比较法直接排序
	//否则先partition再压栈
	CMyStack<pair<int,int> > mystk;
	int nLow,nHigh,nMid;
	nLow=0;
	nHigh=m_tCurDataNum-1;

	pair<int,int> limitPair;


	while(1){
		if(nHigh-nLow<=LOW_LIMIT_FOR_QSORT){
			bubblesort(nLow,nHigh);
			if(mystk.Pop(limitPair)==false)
				break;
			nLow = limitPair.first;
			nHigh = limitPair.second;
		}
		else{
			nMid=quicksort_part(nLow,nHigh);
			if(nHigh-nMid>nMid-nLow){
				limitPair.first=nMid+1;
				limitPair.second=nHigh;
				mystk.Push(limitPair);
				nHigh=nMid-1;
			}
			else{
				limitPair.first=nLow;
				limitPair.second=nMid-1;
				mystk.Push(limitPair);
				nLow=nMid+1;
			}
		}
	}
}

int CPepPairTwoWaySorter::quicksort_part(int low,int high)
{
	// srand((unsigned)time(NULL));
	int pos=rand()%(high-low+1);

	PEP_PAIR_ITEM pepPair=m_Data[low];
	m_Data[low]=m_Data[low+pos];
	m_Data[low+pos]=pepPair;

	pepPair=m_Data[low];
	while(low<high){
		while(low<high && m_Data[high].lfMass>=pepPair.lfMass)
			high--;
		m_Data[low]=m_Data[high];
		while(low<high && m_Data[low].lfMass<=pepPair.lfMass)
			low++;
		m_Data[high]=m_Data[low];
	}
	m_Data[low]=pepPair;
	return low;
}

void CPepPairTwoWaySorter::quicksort(int low,int high)
{
	if(low<high){
		int mid=quicksort_part(low,high);
		if(mid-low<high-mid){
			quicksort(low,mid-1);
			quicksort(mid+1,high);
		}
		else{
			quicksort(mid+1,high);
			quicksort(low,mid-1);
		}
	}
}
