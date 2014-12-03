#ifndef XLINKRESULTFILTER_H_
#define XLINKRESULTFILTER_H_

static bool EVLesser( const pair<double,char>& elem1, const pair<double,char>& elem2 )
{
	return elem1.first < elem2.first;
}


class CXLinkResultFilter : public CXLinkResultFilterInterface
{
public:
	CXLinkResultFilter();
	virtual ~CXLinkResultFilter();
	virtual void Init(string strConf);
	virtual void Run();
	virtual void Close();
	
protected :
	int _IsInTolWindow(double lfExpMH,double lfCalMH);
	char _GetReverseTag(const CXLinkPepResult & pep_res);
	double _GetEvaThreashold(vector<pair<double,char> > & vEvaHist);
	void _FilterByFDR(vector<CXLinkMatchResult> & vResults, vector<CSpectrum> & vSpectra,int nTolBase);
	void _LoadInclusionList();
	
	CFilterConf m_filerconf;
	CXLinkPepResultFilter m_pepResultFilter;
	vector<char > m_vBeOutput;
	map<string,bool> m_mapInclusionList;
	
};

#endif /*XLINKRESULTFILTER_H_*/
