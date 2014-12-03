#ifndef SPECTRAWINDOW_H_
#define SPECTRAWINDOW_H_

namespace proteomics_search
{

struct SPECTRA_WINDOW_INFO
{
	size_t tSpectraNum;
	int nStartId;
	double lfMin;
	double lfMax;
	
	SPECTRA_WINDOW_INFO()
	{
		tSpectraNum = 0;
		nStartId = 0;
		lfMin = 0.0;
		lfMax = 0.0;
	}
};

class CSpectraWindow
{
public:
	CSpectraWindow();
	virtual ~CSpectraWindow();
	bool Construct(const vector<SPEC_SORTER_INDEX_INFO> & vInfo);
	void Begin();
	bool GetCurrentWindowInfo(struct SPECTRA_WINDOW_INFO & stInfo);
	
	
private:
	vector<SPECTRA_WINDOW_INFO> vSpectraWindows;
	size_t m_tCurrentWindow;
	
	
};

}

#endif /*SPECTRAWINDOW_H_*/
