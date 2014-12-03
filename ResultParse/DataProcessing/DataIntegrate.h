#ifndef CDATAINTEGRATE_H_
#define CDATAINTEGRATE_H_

#include "DataIntegrate.h"
#include "../BasicFunction/ReadWrite.h"
#include "../BasicFunction/CommonProcess.h"

#include "../bio_analysis_sdk.h"

class CDataIntegrate
{
public:
	CDataIntegrate();
	virtual ~CDataIntegrate();

	void DataTOIntegrate(const vector<OneDATASET> & ResultDataSet, const CConf & conf,
			IntegratedInfo & HaveToInterData);

private:

	void _Initialize(IntegratedInfo & HaveToInterData);
	void _DeleteSet(const CConf & conf, IntegratedInfo & HaveToInterData);
	bool _IsDeletePeptide(const vector<SPECTRAINFO> & vSpectra, IntegratedInfo & HaveToInterData);
	bool _IsDeleteSpec(const vector<SPECTRAINFO> & vPeptide, IntegratedInfo & HaveToInterData);
	void _CExportIndex(const vector<OneDATASET> & ResultDataSet, const CConf & conf,
			IntegratedInfo & HaveToInterData);
	void _DistinctPeptide(FILE * fout, pair<string, vector<SPECTRAINFO> > pairProTemp,
			const CConf & conf, IntegratedInfo & HaveToInterData);
	void _ProteinInfer(IntegratedInfo & HaveToInterData);
	int _IsSameSetOrSubSet(const set<string> & a, const set<string> & b);
};

#endif /* CDATAINTEGRATE_H_ */
