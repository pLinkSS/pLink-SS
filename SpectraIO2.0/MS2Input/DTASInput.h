/*
 * DTASInput.h
 *
 *  Created on: 2009-7-13
 *      Author: fan
 */

//--------------------------------------------------------------
//              DTAS文件格式
//      DTAS文件头部可以有任意行的附加信息，后面是若干DTA文件。
//      每个DTA文件以文件名开头。后跟一个合法的DTA文件的内容。
//      相邻DTA文件以空行分隔。
//--------------------------------------------------------------

//------------------------------
//         class CDtasIO
//       对DTAS文件的读写类
//------------------------------

#ifndef DTASINPUT_H_
#define DTASINPUT_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include "SpecUtility.h"

using namespace std;

namespace proteomics_sdk {

class CMS2Input;
//class CMS2Info;

class CDTASInput:
	public CMS2Input {
public:
	CDTASInput();
	virtual ~CDTASInput();

	virtual MS2FormatType GetType(void);

	virtual void LoadAll(string strPath, vector<CSpectrum> & S);

	virtual void StartLoad(string strPath, int & nTotal);
	virtual void LoadNext(CSpectrum & spec, int & idx);
	virtual void LoadPrev(CSpectrum & spec, int &idx) ;
	virtual void LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx);
	virtual void LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo);
	virtual vector<string> GetAllSpecName(int begin, int step, int num);
	virtual void EndLoad(void);

protected:
	void OpenInFile(string strFilePath);
	void CloseInFile(void);
	void ReadDTAHead(CSpectrum& spec);
	void ReadMZAndItensity(CSpectrum& spec, bool bAdmitBlank = false);
	bool GetLine(string& str, bool bAdmitBlank = false);
	string GetStrWithStrSplit(string strParent , string::size_type& idx);

	bool m_bMono;
	bool GetDtaFileName(string & strFileName);
	bool IsDTAFileName(string strFileName);

	ifstream m_ifIn;

	string m_strInPath;

	int m_nInCount;

	size_t m_tMaxPeaksNum;
};

}

#endif /* CDTASInput_H_ */
