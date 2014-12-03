/*
 * MS2TypeInput.h
 *
 *  Created on: 2009-7-14
 *      Author: fan
 */

#ifndef MS2TYPEINPUT_H_
#define MS2TYPEINPUT_H_
#include "../../include/sdk.h"
#include "../../include/interface.h"
#include <vector>
#include <fstream>
#include <string>

using namespace std;

namespace proteomics_sdk {

class CMS2Input;
class CSpectrum;

class CMS2TypeInput
	:public CMS2Input{
public:
	CMS2TypeInput();
	virtual ~CMS2TypeInput();
	virtual MS2FormatType GetType(void){return PFF_MS2;}
	virtual void LoadAll(string strPath, vector<CSpectrum> & S) ;

	virtual void StartLoad(string strPath, int & nTotal) ;
	virtual void LoadNext(CSpectrum & spec, int &idx) ;
	virtual void LoadPrev(CSpectrum & spec, int &idx) ;
	virtual void LoadNextBatch(int nNumber ,vector<CSpectrum> & S, int &idx);
	virtual void LoadSpec(CSpectrum & spec, string strTitle, int& nSpecNo);
	virtual vector<string> GetAllSpecName(int begin, int step, int num);
	virtual void EndLoad(void) ;

public:
	bool m_bMono;
	ifstream m_ifIn;
	string m_strInPath;
	string m_strFileName;
	int m_nInCount;
	int m_nTotal;
	size_t m_tMaxPeaksNum;
	map<string, int> m_mapTitle2Loc;
	vector<int> m_vSpecLoc;

	// test
//	FILE * m_ifFile; 
	
//	string strBFmark;

protected:
	void OpenInFile(string strFilePath);
	void ReadFileHead();
	bool GetLine(string& str, bool bAdmitBlank);
	void ReadMS2Head(CSpectrum & spec);
	void ReadMZAndItensity(CSpectrum& spec, bool bAdmitBlank = false);
	void GetScanNo(char *p, int &, int &);
	int  GetCharge(char *p);
	void CloseInFile();

	string strBFmark;
//	bool b;

};

}

#endif /* MS2TYPEINPUT_H_ */
