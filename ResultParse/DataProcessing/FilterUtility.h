#ifndef CFILTERUTILITY_H_
#define CFILTERUTILITY_H_

#include "../BasicFunction/ReadWrite.h"
#include "FilterUtility.h"
#include "../BasicFunction/CommonProcess.h"

#include <set>

class CFilterUtility
{
public:
	CFilterUtility();
	virtual ~CFilterUtility();

	void Filter(const vector<SPECTRAINFO> & vSpectraInfoIn, vector<SPECTRAINFO> & vSpectraOut,
			const CConf & conf, const SearchEngineType & m_EngineType,
			const pair<int, int> & pairSimpleCharge);

private:
	void _Score(const vector<SPECTRAINFO> & vSpectraInfoIn, vector<SPECTRAINFO> & vSpectraOut,
			const CConf & conf, const SearchEngineType & m_EngineType);
	void _SepatrateOrNot(const vector<SPECTRAINFO> & vSpectraInfoIn,
			vector<SPECTRAINFO> & vSpectraOut, const CConf & conf,
			const SearchEngineType & m_EngineType, const pair<int, int> & pairSimpleCharge);
	void _FDR(vector<SPECTRAINFO> & vNewSpectraInfo, vector<SPECTRAINFO> & vSpectraOut,
			const CConf & conf, const SearchEngineType & m_EngineType,
			const pair<int, int> & pairSimpleCharge);
	void _OutTo(vector<SPECTRAINFO> & vSpectraOut, const vector<SPECTRAINFO> & vNEW,
			const size_t resSize, const CConf & conf, const pair<int, int> & pairSimpleCharge);
	void _FilterMLD_Separate(const vector<SPECTRAINFO> & vSpectraInfoTemp, vector<vector<
			SPECTRAINFO> > & vNewSpectraInfo, const CConf & conf);

	void _FilterMLD(const vector<SPECTRAINFO> & vSpectraInfoTemp,
			vector<SPECTRAINFO> & vNewSpectraInfo, const CConf & conf);
	size_t _Filter_FDR(vector<SPECTRAINFO> & vNewSpectraInfo, const CConf & conf, const pair<int,
			int> & pairSimpleCharge);
	void _FixedDeltCn(const vector<SPECTRAINFO> & vNewSpectraInfo,
			vector<SPECTRAINFO> & vpectraInfoTemp, double DeltCn);
	//bool _IsAllReverse(const CMatchSpectraInfo & mr, const CConf & conf);
	bool _IsExist(const CMatchSpectraInfo mr, set<string> & setStr);
	double _Get_FDR(const int & F, const int & T, const int & U, const CConf & conf);
	double _Get_FDR(const int & T,  const int & F3, const int & F2, const int & F1, const CConf & conf);
	void _OUTPUT_FDR_Curve(vector<double> vlfFDR, const CConf & conf,
			const pair<int, int> & pairSimpleCharge);

};

#endif /* CFILTERUTILITY_H_ */
