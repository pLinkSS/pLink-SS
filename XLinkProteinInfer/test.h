#ifndef ACALGORITHM_H_
#define ACALGORITHM_H_
#include<vector>
#include<string>
#include "../include/sdk.h"
using namespace std;
using namespace proteomics_sdk;

#define MAX_STATE 10000000

struct Node
{
	Node()
	{
		m_tLeft = 0;
		m_tRight = 0;
		m_tLeftEdge = 0;
		m_tRightEdge = 0;
		m_tF = 0;
		};
	~Node(){};
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

class AC
{
	public :
		AC(){pCond=NULL;};
		~AC(){};
		
		void Init_(CCondition * pC);
		size_t Find_(size_t tCurr, char cVal);
		void MergeOutput_(size_t tDstState, size_t tSrcState);
		void CreateFailFunction_();
		void Create_(const vector<string> & vSrc);
		bool Transfer_(size_t & tCurr, char cVal);
		int Output_(const string &strSrc, size_t tPepEnd, size_t tCurr, set<pair<size_t, size_t> > & setIDs);
		void Query_(const string & strSrc, set<pair<size_t,size_t> > & setIDs);
		bool Empty_();
private:
	//Node g_states[MAX_STATE];
	vector<Node> g_states; 
	map<size_t, vector<size_t> > mapOutput;
	vector<string> vTotalPep;
	size_t tCurrentState;
	size_t tTotQuery;
	CCondition * pCond;
};

#endif /*ACALGORITHM_H_*/
