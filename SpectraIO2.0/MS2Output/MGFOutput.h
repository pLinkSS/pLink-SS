/*
 * MGFOutput.h
 *
 *  Created on: 2010-2-25
 *      Author: fan
 */

//--------------------------------------------------------------
//              MGF文件格式
//      若干单独的文件，起始于"BEGIN IONS",终止于"END IONS"
//      每个文件头:TITLE,CHARGE and PEPMASS，后面类似DTA。
//
//--------------------------------------------------------------

//------------------------------
//         class CMgfIO
//       对MGF文件的读写类
//------------------------------

#ifndef MGFOUTPUT_H_
#define MGFOUTPUT_H_
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

class CMGFOutput :
	public CMS2Output{
public:
	CMGFOutput();
	virtual ~CMGFOutput();
//	virtual void WriteAll(string strPath, vector<CSpectrum> & S);
	virtual void WriteAll(string strPath, vector<CSpectrum> & S);
	virtual void StartWrite(string strPath, int & nTotal);
	virtual void WriteNext(CSpectrum & spec, int & idx);
	virtual void WriteNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx);
	virtual void EndWrite(void);

	/*virtual void StartWrite(string strPath);
	virtual void WriteNext(CSpectrum & spec);
	virtual void EndWrite(void);
	*/
	virtual MS2FormatType GetType(void) {return PFF_MGF;};

protected:
	void OpenOutFile(string strFilePath);
	void CloseOutFile(void);
	void WriteMGFHead(CSpectrum& spec);
	void WriteMZAndItensity(CSpectrum& spec, bool bAdmitBlank = false);
//	bool GetLine(string& str, bool bAdmitBlank = false);
//	void OpenOutFile(string strFilePath);
//	void CloseOutFile(void);
	bool m_bMono;
//	void WriteSingleMGF(CSpectrum & spec);
	FILE * m_ofOut;
	ifstream m_ifIn;
//	ofstream m_ofOut;
	string m_strOutPath;
//	string m_strOutPath;
	int m_nOutCount;
//	int m_nOutCount;
	size_t m_tMaxPeaksNum;

	////用于文件头有CHARGE信息
	vector<int>  m_vCharge;//CHARGE的个数。以AND分开，如2+ and 3+
	vector<int> m_vCurrCharge;//当前的CHARGE
	bool m_bMissCharge;//看每个谱是否缺失CHARGE，缺失则以m_vCharge中的数目填补,将当前行记录到m_strMissCharge中。
	string m_buffString;//用于记录读到的第一对谱图数目项
	string m_strMissCharge;
	CSpectrum m_specTemp;
	////END用于文件头有CHARGE信息

	void WriteFileHead();
	string GetTitle(CSpectrum& spec);
//	void WriteFloat(double lfNum, int precision = PRECISION);
};

}

#endif /* MGFINPUT_H_ */
