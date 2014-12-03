/*
 * MGFInput.h
 *
 *  Created on: 2009-7-13
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

#ifndef MGFINPUT_H_
#define MGFINPUT_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include "SpecUtility.h"
#include <sys/stat.h>

using namespace std;

namespace proteomics_sdk {

class CMS2Input;
//class CMS2Info;

class CMGFInput :
	public CMS2Input{
public:
	CMGFInput();
	virtual ~CMGFInput();
//	virtual void WriteAll(string strPath, vector<CSpectrum> & S);
	virtual void LoadAll(string strPath, vector<CSpectrum> & S);
	virtual void StartLoad(string strPath, int & nTotal);
	virtual void LoadNext(CSpectrum & spec, int & idx);
	virtual void LoadPrev(CSpectrum & spec, int &idx) ;
	virtual void LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx);
	virtual void LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo);
	virtual vector<string> GetAllSpecName(int begin, int step, int num);
	virtual void EndLoad(void);

	/*virtual void StartWrite(string strPath);
	virtual void WriteNext(CSpectrum & spec);
	virtual void EndWrite(void);
	*/
	virtual MS2FormatType GetType(void) {return PFF_MGF;};

protected:
	void OpenInFile(string strFilePath);
	void CloseInFile(void);
	void ReadMGFHead(CSpectrum& spec);
	void ReadMZAndItensity(CSpectrum& spec, bool bAdmitBlank = false);
	bool GetLine(string& str, bool bAdmitBlank = false);
//	void OpenOutFile(string strFilePath);
//	void CloseOutFile(void);
	bool m_bMono;
//	void WriteSingleMGF(CSpectrum & spec);
	ifstream m_ifIn;
//	ofstream m_ofOut;
	string m_strInPath;
//	string m_strOutPath;
	int m_nInCount;
//	int m_nOutCount;
	size_t m_tMaxPeaksNum;
	int m_nTotal;

	////用于文件头有CHARGE信息
	vector<int>  m_vCharge;//CHARGE的个数。以AND分开，如2+ and 3+
	vector<int> m_vCurrCharge;//当前的CHARGE
	bool m_bMissCharge;//看每个谱是否缺失CHARGE，缺失则以m_vCharge中的数目填补,将当前行记录到m_strMissCharge中。
	string m_buffString;//用于记录读到的第一对谱图数目项
	string m_strMissCharge;
	CSpectrum m_specTemp;
	map<string, int> m_mapTitle2Loc;
	vector<int> m_vSpecLoc;
	////END用于文件头有CHARGE信息

	void ReadFileHead();
	
	double m_lfLineCount;//记录读到第几行
//	void WriteFloat(double lfNum, int precision = PRECISION);
};

}

#endif /* MGFINPUT_H_ */
