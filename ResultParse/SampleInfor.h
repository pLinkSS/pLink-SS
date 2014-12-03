#ifndef SAMPLEINFOR_H_
#define SAMPLEINFOR_H_

#include<string>
#include <ctime>
#include "bio_analysis_sdk.h"

#include <ctime>
#include <sstream>
#include "SampleInfor.h"
//#include "PSM\Denovo.h"
#include "ExportResult/ExportHTML.h"
#include "BasicFunction/ReadWrite.h"
#include "ExportResult/MergeTopFind.h"
#include "BasicFunction/DataAnalysis.h"
#include "BasicFunction/CommonProcess.h"
#include "BasicFunction/StringUtility.h"
#include "ExportResult/ExportpBuildPre.h"
#include "DataProcessing/FilterUtility.h"
#include "DataProcessing/DataIntegrate.h"
#include "ExportResult/ProteinInference.h"
#include "ExportResult/ExportForInterface.h"
#include "EnginesParser/ResultParserFactory.h"

#include "ExportResult/ExportToIntermediateFile.h"

using namespace bio_analysis;

class CSampleInfor
{

public:
	void Run(const string & strParamFilePath, const string & strIndexFile = "");

private:

	clock_t start, end;

//	CDeNovo m_CDeNovo;
	CFilterUtility m_FilterUtility;
	CDataIntegrate m_DataIntegrate;

	//CDataAnalysis m_DataAnalysis(conf);

	void _OutPutCConf(CConf & conf);
	void _Close(const CConf & conf);
	void _Makdir(CConf & conf);
	void _ReadparametersAndConditions(CConf & conf);

	void _Filter(vector<OneDATASET> & AllDataSetParser, vector<OneDATASET> & AllDataSetFilter,
			const CConf & conf);
	void _Parser(vector<OneDATASET> & AllDataSetParser, CConf & conf,
			CDataAnalysis & m_DataAnalysis);

	void _ExportResult(const CConf & conf);

	void _ModSiteFilter(OneDATASET & ResultDataSet, const CConf & conf);
	void _TerminalFilter(OneDATASET & ResultDataSet, const CConf & conf);
	void _CrossLinkFilter(vector<OneDATASET> & AllDataSetParserconst, const CConf & conf);
	void _DeCrossLinkFilter(vector<OneDATASET> & AllDataSetParser, const CConf & conf);

};

#endif /* SAMPLEINFOR_H_ */
