#ifndef SEARCHER_H_
#define SEARCHER_H_

//#include <pthread.h>
namespace proteomics_search
{
class CSearchEngine;
class CSearcher
{
public:
	CSearcher();
	
	virtual ~CSearcher();
	
	void Init(int nArgc, char ** pArgv);
	
	void Run(void);
	
	void Close(void);

protected:
	bool _ParseCmdLine(int nArgc, char ** pArgv);
	
	string m_strParamFile;
	
	CSearchEngine m_SearchEngine;

};
}
#endif /*Searcher_H_*/
