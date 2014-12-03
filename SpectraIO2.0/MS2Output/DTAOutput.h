/*
 * DTAOutput.h
 *
 *  Created on: 2010-3-22
 *      Author: fan
 */

//--------------------------------------------------------------
//              DTA file format
//      file head: MH+ CHARGE
//      several pairs followed: MZ INTENSITY
//--------------------------------------------------------------
///////
#ifndef DTAOUTPUT_H_
#define DTAOUTPUT_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"
#include <vector>
#include <fstream>
#include <string>
#include <dirent.h>
using namespace std;
namespace proteomics_sdk
{

class CMS2Output;
//class CMS2Info;

class CDTAOutput:
	public CMS2Output{
public:
	CDTAOutput();
	virtual ~CDTAOutput();

	virtual MS2FormatType GetType(void);
	virtual void WriteAll(string strPath, vector<CSpectrum> & S) ;

	virtual void StartWrite(string strPath, int & nTotal) ;
	virtual void WriteNext(CSpectrum & spec, int &idx) ;
	virtual void WriteNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx);
	virtual void EndWrite(void) ;

protected:

	//------------Write-----------//
	void WriteDTA(string strPath, CSpectrum & spec);
//	void WriteFloat(double lfNum);

	void WriteDTAHead(CSpectrum& spec);
	void WriteMZAndItensity(CSpectrum& spec, bool bAdmitBlank = false);
//	bool GetLine(string& str, bool bAdmitBlank = false);
	void OpenOutFile(string strFilePath);
	void CloseOutFile(void);

public:
	bool m_bMono;
protected:

	FILE * m_ofOut;
	string m_strOutPath;
	int m_nOutCount;
//	int m_nOutCount;
	CFileFind m_ffFinder;
	vector<string> m_vstrTotalFilePath;

	//store the current reserved space for peaks
	size_t m_tMaxPeaksNum;
};
}
#endif /* DTAOUTPUT_H_ */

