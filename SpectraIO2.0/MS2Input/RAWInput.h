/*
 * RAWInput.h
 *
 *  Created on: 2009-7-21
 *      Author: fan
 */

#ifndef RAWINPUT_H_
#define RAWINPUT_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"
#include <string>
#include <vector>
#include <fstream>
using namespace std;

namespace proteomics_sdk {

class CMS2Input;
class CSpectrum;

class CRAWInput:
	public CMS2Input{
public:
	CRAWInput();
	virtual ~CRAWInput();
	virtual void LoadAll(string strPath, vector<CSpectrum> & S);
	virtual void StartLoad(string strPath, int & nTotal);
	virtual void LoadNext(CSpectrum & spec, int & idx);
	virtual void LoadPrev(CSpectrum & spec, int &idx) ;
	virtual void LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx);
	virtual void LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo);
	virtual vector<string> GetAllSpecName(int begin, int step, int num);
	virtual void EndLoad(void);
	virtual MS2FormatType GetType(void) {return PFF_RAW;};
protected:
	void OpenInFile(string strFilePath);
	void CloseInFile(void);
	void ReadPFDHead(CSpectrum& spec);
	void ReadMZAndItensity(CSpectrum& spec, bool bAdmitBlank = false);
	bool GetLine(string& str, bool bAdmitBlank = false);

	void ReplaceSlash(string & command);

	bool m_bMono;
	ifstream m_ifIn;
	string m_strInPath;
	string m_strFileName;
	int m_nInCount;
	size_t m_tMaxPeaksNum;
	string strBFmark;
};

}

#endif /* RAWINPUT_H_ */
