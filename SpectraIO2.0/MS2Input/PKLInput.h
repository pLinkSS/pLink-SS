/*
 * PKLInput.h
 *
 *  Created on: 2009-7-13
 *      Author: fan
 */

//--------------------------------------------------------------
//              PKL文件格式
//      若干单独的文件合并。文件头:MZ,INTENSITY,CHARGE
//      后面若干行与DTA文件类似
//      相邻文件以空行分隔。
//
//--------------------------------------------------------------

#ifndef PKLINPUT_H_
#define PKLINPUT_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include "SpecUtility.h"

using namespace std;


namespace proteomics_sdk {

class CSpectrum;
class CMS2Input;

class CPKLInput: public CMS2Input {
public:
	CPKLInput();
	virtual ~CPKLInput();

	virtual void LoadAll(string strPath, vector<CSpectrum> & S);

	virtual void StartLoad(string strPath, int & nTotal);
	virtual void LoadNext(CSpectrum & spec, int & idx);
	virtual void LoadPrev(CSpectrum & spec, int &idx) ;
	virtual void LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx);
	virtual void LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo);
	virtual vector<string> GetAllSpecName(int begin, int step, int num);
	virtual void EndLoad(void);

	virtual MS2FormatType GetType(void) {return PFF_PKL;};
protected:
	void OpenInFile(string strFilePath);
	void CloseInFile(void);
	void ReadPKLHead(CSpectrum& spec);
	void ReadMZAndItensity(CSpectrum& spec);
	bool GetLine(string& str, bool bAdmitBlank = false);
	bool IsPKLHead(string strLine);

	bool m_bMono;

	ifstream m_ifIn;
	string m_strInPath;
	int m_nInCount;
	size_t m_tMaxPeaksNum;
};

}

#endif /* PKLINPUT_H_ */
