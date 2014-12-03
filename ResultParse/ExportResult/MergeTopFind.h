/*
 * CMergeTopFind.h
 *
 *  Created on: 2010-10-12
 *      Author: Hailer
 */

#ifndef MERGETOPFIND_H_
#define MERGETOPFIND_H_

#include "../bio_analysis_sdk.h"

#include "../BasicFunction/CommonProcess.h"
#include "../BasicFunction/ReadWrite.h"

using namespace bio_analysis;

class CMergeTopFind
{
public:
	CMergeTopFind();
	virtual ~CMergeTopFind();

	void pFindMerge_Export(const vector<OneDATASET> & ResultDataSet, const CConf & conf);

private:
	void _OutPutSpectraASpFind(FILE * fout, const CMatchSpectraInfo & SpectraTemp,
			const size_t & SpectraNum);
	void _PrintfSpectraOut(FILE * fout, CMatchSpectraInfo & SpectraOut, bool & tag,
			size_t & SpectraNum);
	void _PeptideJoin(vector<CMatchPeptideInfo> & vPep, const CMatchPeptideInfo & pepTemp);
	void _MergePeptide(CMatchSpectraInfo & SpectraOut, const CMatchSpectraInfo & SpectraTemp);
	void _Merge(const OneDATASET & OneResult, const CConf & conf, int DataSetID);
};

#endif /* MERGETOPFIND_H_ */
