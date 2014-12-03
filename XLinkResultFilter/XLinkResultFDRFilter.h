#ifndef XLINKRESULTFDRFILTER_H_
#define XLINKRESULTFDRFILTER_H_

#define MIN_WINDOW 50
#define SPECIFICITY 0.01


static bool _FeatureLesser( const pair<double,size_t>& elem1, const pair<double,size_t>& elem2 )
{
	return elem1.first < elem2.first;
}

static bool _FeatureGreater( const pair<double,size_t>& elem1, const pair<double,size_t>& elem2 )
{
	return elem1.first > elem2.first;
}

class CXLinkResultFDRFilter : public  CXLinkResultFilterInterface
{
public:
	CXLinkResultFDRFilter();
	virtual ~CXLinkResultFDRFilter();
	virtual void Init(string strConf);
	virtual void Run();
	virtual void Close();
	
protected :
	
	char _GetReverseTag(const CXLinkPepResult & pep_res);
	void _SetTUFTag();
	void _EstimateTargetNum();
	int _GetTagLength(const CXLinkPepResult & pep_res);
	int _CalCharge(bool bSet = false);
	int _CalDeltaMass(bool bSet = false);
	int _CalIntenRatio(bool bSet = false);
	int _CalTagLength(bool bSet = false);
	int _CalScore(bool bSet = false);
	int _CalFeature(bool bSet);
	int _CalEvalue(bool bSet = false);
	void _FilterByFDR();
	bool _IsInTolWindow(double lfExpMH,double lfCalMH);
	
	CFilterConf m_filerconf;
	CXLinkPepResultFilter m_pepResultFilter;
	vector<int > m_vBeOutput;
	vector<char > m_vTUFTag;
	vector<pair<double,size_t> > m_vFeatureList;
	
	int m_nTargetNum ;
	vector<CXLinkMatchResult> m_vResults;
	vector<CSpectrum> m_vSpectra;
	
};

#endif /*XLINKRESULTFDRFILTER_H_*/
