#ifndef SPECTRASEARCHER_H_
#define SPECTRASEARCHER_H_

namespace proteomics_search
{

class CSpectraSearcher
{
public:
	int * m_pData;
	int m_nStart;
	void Construct(const vector<SPEC_SORTER_INDEX_INFO> & vInfo);
	void Close();
	int FindLowerBound(double & lfMass);
	CSpectraSearcher();
	~CSpectraSearcher();
};

}

#endif /*SPECTRASEARCHER_H_*/
