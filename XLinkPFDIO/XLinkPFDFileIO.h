#ifndef XLINK_PFDFILEIO_H_
#define XLINK_PFDFILEIO_H_

class CXLinkPFDFileIO
{
public:
	CXLinkPFDFileIO();
	virtual ~CXLinkPFDFileIO();
	bool WriteFile(string strOutputFile, const vector<CXLinkMatchResult> & vResults,
			const vector<CSpectrum> & vSpectra);
	bool LoadFile(string strOutputPath, vector<CXLinkMatchResult> & vResults,
			vector<CSpectrum> & vSpectra, size_t &tLoad);
	bool LoadAll(string strOutput, int nSpectraTotal, int nTotalFile, vector<CXLinkMatchResult> & vResults,
			vector<CSpectrum> & vSpectra, size_t &nTotalSize, const string & strIdentifier);
private:
	CTrace *m_pTrace;
};

#endif /*XLINK_PFDFILEIO_H_*/
