#ifndef ACALGORITHM_H_
#define ACALGORITHM_H_
#include<vector>
#include<string>
#include "../include/sdk.h"
using namespace std;
using namespace proteomics_sdk;
namespace AC
{
struct Node
{
	size_t m_tLeft;
	size_t m_tRight;
	char m_tLeftEdge;
	char m_tRightEdge;
	size_t m_tF;
	size_t inline GetChild(){return m_tLeft;};
	size_t inline GetBrother(){return m_tRight;};
	size_t inline GetF(){return m_tF;};
	void inline SetChild(size_t t, char c = 0){m_tLeft = t;m_tLeftEdge = c;};
	void inline SetBrother(size_t t, char c = 0){m_tRight = t;m_tRightEdge = c;};
	void inline SetF(size_t t) {m_tF = t;};
};
void Init(CCondition * pC);
size_t Find(size_t tCurr, char cVal);
void MergeOutput(size_t tDstState, size_t tSrcState);
void CreateFailFunction();
void Create(const vector<string> & vSrc);
bool Transfer(size_t & tCurr, char cVal);
int Output(const string &strSrc, size_t tPepEnd, size_t tCurr, set<size_t> & setIDs);
void Query(const string & strSrc, set<size_t> & setIDs);
bool Empty();

void Query(const string & strSrc, set< pair<size_t, size_t> > & setIDs);
int Output(const string &strSrc, size_t tPepEnd, size_t tCurr, set< pair<size_t, size_t> > & setIDs);

}
#endif /*ACALGORITHM_H_*/
