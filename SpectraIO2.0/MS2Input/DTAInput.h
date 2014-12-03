/*
 * DTAInput.h
 *
 *  Created on: 2009-7-7
 *      Author: fan
 */

//--------------------------------------------------------------
//              DTA file format
//      file head: MH+ CHARGE
//      several pairs followed: MZ INTENSITY
//--------------------------------------------------------------
///////
#ifndef DTAINPUT_H_
#define DTAINPUT_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"
#include <vector>
#include <fstream>
#include <string>
#include <dirent.h>
using namespace std;
namespace proteomics_sdk
{

class CMS2Input;
//class CMS2Info;

class CDTAInput:
	public CMS2Input{
public:
	CDTAInput();
	virtual ~CDTAInput();

	virtual MS2FormatType GetType(void);
	virtual void LoadAll(string strPath, vector<CSpectrum> & S) ;

	virtual void StartLoad(string strPath, int & nTotal) ;
	virtual void LoadNext(CSpectrum & spec, int &idx) ;
	virtual void LoadPrev(CSpectrum & spec, int &idx) ;
	virtual void LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx);
	virtual void LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo);
	virtual vector<string> GetAllSpecName(int begin, int step, int num);
	virtual void EndLoad(void) ;

protected:

	//------------LOAD-----------//
	void LoadDTA(const char* pszDTAFilePath, const char* pszFileOnly, vector<CSpectrum> & S);
	void ReadDTAHead(CSpectrum& spec);
	void ReadMZAndItensity(CSpectrum& spec, bool bAdmitBlank = false);
	bool GetLine(string& str, bool bAdmitBlank = false);
	void OpenInFile(string strFilePath);
	void CloseInFile(void);

public:
	bool m_bMono;
protected:
	ifstream m_ifIn;
//	ofstream m_ofOut;

	string m_strInPath;
//	string m_strOutPath;

	int m_nInCount;
	int m_nTotal;
//	int m_nOutCount;

	CFileFind m_ffFinder;
	vector<string> m_vstrTotalFilePath;
	vector<string> m_vstrPartFilePath;
	map<string, int> m_mapTitle2Loc;

	//store the current reserved space for peaks
	size_t m_tMaxPeaksNum;
};
}
#endif /* DTAINPUT_H_ */

