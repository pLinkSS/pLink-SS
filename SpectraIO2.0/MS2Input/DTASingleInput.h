#pragma once
#ifndef SDTAINPUT_H_
#define SDTAINPUT_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"
#include <string>
#include <vector>
#include <fstream>
using namespace std;

namespace proteomics_sdk {

class CDTASingleInput:
	public CMS2Input
{
public:
	virtual void LoadAll(string strPath, vector<CSpectrum> & S);
	virtual void StartLoad(string strPath, int & nTotal);
	virtual void LoadNext(CSpectrum & spec, int & idx);
	virtual void LoadPrev(CSpectrum & spec, int &idx) ;
	virtual void LoadNextBatch(int nNumber, vector<CSpectrum> & S, int &idx);
	virtual void LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo);
	virtual vector<string> GetAllSpecName(int begin, int step, int num);
	virtual void EndLoad(void);
	virtual MS2FormatType GetType(void) {return PFF_SDTA;};
	CDTASingleInput(void);
	~CDTASingleInput(void);
protected:
	void OpenInFile(string strFilePath);
	void CloseInFile(void);
	void ReadMZAndItensity(CSpectrum& spec);
	bool GetLine(string& str, bool bAdmitBlank);
	void ReadMHAndCharge(CSpectrum& spec);

	ifstream m_ifIn;
	string m_strInPath;
};

}
#endif
