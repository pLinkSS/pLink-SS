#ifndef PEPPAIRSORTER_H_
#define PEPPAIRSORTER_H_

#include <string>
#include <cfloat>
#include "../include/sdk.h"
using namespace std;
using namespace proteomics_sdk;

namespace proteomics_search {

class CPepPairSorter
{
public:
	CPepPairSorter();
	CPepPairSorter(string strFile, size_t tDataNum);
	virtual ~CPepPairSorter() {}

	void SetDataInfo(string strFile, size_t tDataNum);
	/* divide original file into several chunks, meanwhile sort them */
	virtual bool BatchDivision() = 0;
	/* merge all those chunks */
	virtual bool BatchMergeSort() = 0;
protected:
	/* file to read */
	string m_strFileName;
	/* how many elements */
	size_t m_tDataNum;
};

/* Written by Jeremy */
class CPepPairMoreWaySorter : public CPepPairSorter
{
public:
	CPepPairMoreWaySorter();
	CPepPairMoreWaySorter(string strFile, size_t tDataNum);
	virtual ~CPepPairMoreWaySorter();

	virtual bool BatchDivision();
	virtual bool BatchMergeSort();
private:
	vector<string> m_vTmpFileName;
	size_t m_tTmpFileNum;

	PEP_PAIR_ITEM *m_arrData;

	CTrace *m_pTrace;
#define MIN_MASS (-DBL_MAX)
#define MAX_MASS DBL_MAX
private:
	class LoserTree {
	public:
		LoserTree(const size_t tWay);
		~LoserTree();
		void ConnectToChunks(const vector<string> &vInFileName,
				const string &strOutFileName);
		void CreateLoserTree();
		void Merge();
	private:
		void _Adjust(int nItem);
	private:
		FILE **m_fpIn;
		FILE *m_fpOut;

		size_t m_tChunkNum; /* how many chunks to merge */
		PEP_PAIR_ITEM *m_arrBuf; /* m_tChunkNum+1 */
		int *m_nTree; /* the number of nodes is m_tChunkNum */
	};
};

/* Designed by Wu */
class CPepPairTwoWaySorter : public CPepPairSorter
{
public:
	CPepPairTwoWaySorter();
	CPepPairTwoWaySorter(string strFile, size_t tDataNum);
	virtual ~CPepPairTwoWaySorter();
	
	virtual bool BatchDivision();
	virtual bool BatchMergeSort();
private:
	void QuickSort_NoneRecursive();

	//批次排序程序
	//读取数据文件，按块分割若干子文件
	//子文件数目由m_nTmpFileNum记录
	//子文件名由m_strTmpFileName数组记录
	//每个子文件由8Bdouble型数据组成
	bool BatchReadFile();
	//子文件排序=单个文件快排+多个文件归并
	bool BatchSort();
	//最终将归并完成的文件导入输出文件
	bool BatchWriteFile();
	//临时文件清理：文件删除+回收存储文件名的空间
	void ClearTmpFile();

	//快速排序内调函数
	int quicksort_part(int low,int high);
	void quicksort(int low,int high);
	void bubblesort(int low,int high);

private:
	char **m_strTmpFileName; //临时文件名
	int m_nTmpFileNum; // 临时文件个数

	PEP_PAIR_ITEM * m_Data; //数据存储区
	size_t m_tCurDataNum;

	CTrace *m_pTrace;
};

}

#endif /*PEPPAIRSORTER_H_*/
