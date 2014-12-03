/*
 * MS2TypeOutput.h
 *
 *  Created on: 2010-3-11
 *      Author: fan
 */


#ifndef MS2TYPEOUTPUT_H_
#define MS2TYPEOUTPUT_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include "../MS2Input/SpecUtility.h"

using namespace std;

namespace proteomics_sdk {

class CMS2Output;
//class CMS2Info;

class CMS2TypeOutput :
	public CMS2Output{
public:
	CMS2TypeOutput();
	virtual ~CMS2TypeOutput();
//	virtual void WriteAll(string strPath, vector<CSpectrum> & S);
	virtual void WriteAll(string strPath, vector<CSpectrum> & S);
	virtual void StartWrite(string strPath, int & nTotal);
	virtual void WriteNext(CSpectrum & spec, int & idx);
	virtual void WriteNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx);
	virtual void EndWrite(void);
	virtual MS2FormatType GetType(void) {return PFF_MS2;};
	
protected:
	void OpenOutFile(string strFilePath);
	void CloseOutFile(void);
	void WriteMS2Head(CSpectrum& spec);
	void WriteMZAndItensity(CSpectrum& spec, bool bAdmitBlank = false);

	bool m_bMono;

	FILE * m_ofOut;
	ifstream m_ifIn;

	string m_strOutPath;

	int m_nOutCount;

	size_t m_tMaxPeaksNum;


	////END用于文件头有CHARGE信息

	void WriteFileHead();

};

}

#endif /* MGFINPUT_H_ */
