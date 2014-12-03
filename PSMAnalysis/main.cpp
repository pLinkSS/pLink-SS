
#include <iostream>
#include <stdio.h>
#include <string>
//#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "../XLinkResultReport/XLinkResultReport.h"
#include "../XLinkResultReport/XLinkResultReportFactory.h"
#include "../SpectraIO2.0/MS2Input/MS2InputFactory.h"
#include "../PreProcess/PreProcessFactory.h"
#include "string.h"
#include "PSMConf.h"
#include "PSMatch.h"
using namespace std;
using namespace proteomics_sdk;



// analyse the PSM result to output the unmatched peaks and find out the reason 
// input the PSM , ion types , mass scope to consider
// output the unmatched peaks , gaps between the missed peaks and its matched counterparts

/*
 *	format of *.psm:
 *  [pfd] 
 * 	pFindFile=
 *  ReportFile=
	OutputPath=
 *	[ion]
 * 	tolerance=0.5	// Da ; fragment tolerance
 * 	ion_type_total=4
	ion_type1=b 1 0 0 0
	ion_type2=b 2 0 0 0
	ion_type3=y 1 0 0 0
	ion_type4=y 2 0 0 0
 * 	[gap]
 * 	mass_scope=40 // Da ; mass gap to consider between unmatched and matched peaks
 *	
 */	
 
#ifdef WIN32
const char cSlash = '\\';
#else
const char cSlash = '/';
#endif
 
bool _LoadSpectraList(string strSpectraList, string strPrefix, map<string , bool> & mapSpectraList)
{
	FILE * fp;
	if ( (fp = fopen(strSpectraList.c_str(),"r")) == NULL )
	{
		return false;
	}
	
	char szbuf[1024];
	while(fgets(szbuf,1024,fp))
	{
		if(szbuf[strlen(szbuf)-1]==0x0a)
			szbuf[strlen(szbuf)-1]=0;
		if(szbuf[strlen(szbuf)-1]==0x0d)
			szbuf[strlen(szbuf)-1]=0;
		
		if(strPrefix.empty() || strncmp(szbuf,strPrefix.c_str(),strPrefix.length()) == 0)
		{
			mapSpectraList[&szbuf[strPrefix.length()]] = true;	
		}
	}
	fclose(fp);
	
	return true;
}

bool _LoadSpectra(string strpFindFile,vector<CSpectrum> & vSpectra)
{
	COptionTool * pOption = new COptionTool("spectrum", strpFindFile.c_str());
	string strSpectraPath = pOption->GetString("spec_path", "null");
    if("null" == strSpectraPath)
    {
        delete pOption;
        return false;
    }

    string strType = pOption->GetString("spec_type", "PKL");
    MS2FormatType eType=PFF_PKL;
    
    if("DTA" == strType)
    {
    	eType=PFF_DTA;
        if(strSpectraPath[strSpectraPath.length() - 1]
                          != cSlash)
        	strSpectraPath += cSlash;
    }
    else if("DTAS" == strType)
    	eType=PFF_DTAS;
    else if("MGF" == strType)
    {
    	eType=PFF_MGF;
    }
    else if("RAW" == strType)
    	eType=PFF_RAW;
    else if("MS2" == strType)
    	eType=PFF_MS2;
    
    delete pOption;
    pOption = NULL;
    
    
   
    CMS2InputFactory IOFactory;
    CMS2Input * pIO = IOFactory.GetImporter(eType);
    
    try
    {
		vSpectra.clear();
    	pIO->LoadAll(strSpectraPath,vSpectra);
    }
    catch(exception & e)
	{
		CErrInfo info("PSMAnalysis", "_LoadSpectra", "in function CMS2Input::LoadAll.");
		info.Append("strSpectraPath=" + strSpectraPath);
		throw runtime_error(info.Get(e).c_str());					
	}
	catch(...)
	{
		CErrInfo info("PSMAnalysis", "_LoadSpectra", "in function CMS2Input::LoadAll.");
		info.Append("strSpectraPath=" + strSpectraPath);
		throw runtime_error(info.Get().c_str());					
	}
	
	delete pIO;
}

string GetTitle(string strPath)
{
	char szbuf[1024];
	int nStart = 0,nEnd = strPath.length()-1;
	int i;
	for(i=strPath.length()-1;i>=0;--i)
	{
		if(strPath[i] == '.')
			break;
	}
	nEnd = i-1;
	for(;i>=0;--i)
	{
		if(strPath[i] == cSlash)
			break;
	}
	nStart = i+1;
	
	for(i=nStart;i<=nEnd;i++)
	{
		szbuf[i-nStart] = strPath[i];
	}
	
	szbuf[i-nStart] = 0;
	return string(szbuf);
}

struct MatchInfoMeasure
{
	double lfCountRatio;
	double lfIntenRatio;
	double lfTolRatio;
	double lfIsotopRatio;
	
	void Clear()
	{
		lfCountRatio = 0;
		lfIntenRatio = 0;
		lfTolRatio = 0;
		lfIsotopRatio = 0;
	}
};

bool ion_significance_lesser(const struct IonMatchInfo & elem1, const struct IonMatchInfo & elem2 ) 
{
	return elem1.lfSignificance > elem2.lfSignificance;
}

int main(int argc , char * argv[])
{
	
	if(argc != 2)
	{
		cout << "usage : PSMAnalysis.exe *.psm" << endl;
		return 1;
	}
	
	CTrace *pTrace = CTrace::GetInstance();
	CPSMConf psmconf;
	pTrace->Info("Loading PSM configuration information...");
	psmconf.Load(argv[1]);
	pTrace->Info("Successfully loaded PSM configuration information.");
	
	map<string , bool> mapSpectraList;
	bool bSpectraListAvail;
	bSpectraListAvail = _LoadSpectraList(psmconf.m_strSpectraList, psmconf.m_strSpectraListPrefix , mapSpectraList);
	
	vector<CXLinkMatchResult> vResults;
	vector<CSpectrum> vSpectra;

	// load results and spectra (with output peak info)
	CXLinkResultReportFactory reportFactory;
	
	CXLinkResultReport * report = NULL;
	if(psmconf.m_strInputType == "pbuild")
		report = reportFactory.GetReport(XLINK_PBUILD);
	else
		report = reportFactory.GetReport(Xlink_STANDARD);
	report->Init(psmconf.m_strpFindFile);
	
	vector<CSpectrum> vSpectraOfResult;
	report->LoadFile(psmconf.m_strReportFile,vResults,vSpectraOfResult);
	report->Close();
	// load spectra with peak info and construct index on peaks
	// load spectra  
	_LoadSpectra(psmconf.m_strpFindFile,vSpectra);
	
	map<string,size_t> mapSpectra;
	for(size_t i = 0 ; i < vSpectra.size() ; ++ i)
	{
		//cout << "read in " << vSpectra[i].m_strFilePath << endl;
		mapSpectra[vSpectra[i].m_strFilePath] = i;
	}
	/*
	if(vSpectra.size() != tSize)
	{
		cout << "spectra size doesn't math with result size .." << endl;
		return -1;
	}
	*/
	// construct peak indexes
	CPreProcessFactory PreProcessFactory;
	
	CPreProcess * pPreProc ;
	if(psmconf.m_bNoiseRemove)
		pPreProc = PreProcessFactory.GetPreProcessor(/*PRE_PROC_PSM PRE_PROC_NULL*/psmconf.m_Condition.m_ePreProcMethod);
	else
		pPreProc = PreProcessFactory.GetPreProcessor(/*PRE_PROC_PSM PRE_PROC_NULL*/psmconf.m_Condition.m_ePreProcMethod);

	pPreProc->Init(psmconf.m_Condition);
	for(size_t i=0; i< vSpectra.size(); ++i)
	{
		if(0 == vSpectra[i].m_tPeaksNum)
			continue;
		CSpectrum Output = vSpectra[i];
		try
		{
			if(!bSpectraListAvail || mapSpectraList.find(vSpectra[i].m_strFilePath) != mapSpectraList.end() )
				pPreProc->Run(vSpectra[i], Output);
		}
		catch(exception & e)
		{
			CErrInfo info("PSMAnalysis", "main", "in the function CPreProcess::Run.");
			throw runtime_error(info.Get(e).c_str());			
		}
		catch(...)
		{
			CErrInfo info("PSMAnalysis", "main", "caught an unknown exception in the function CPreProcess::Run.");
			throw runtime_error(info.Get().c_str());			
		}	
		vSpectra[i]=Output;
	}
	pPreProc->Close();
	
	CPSMatch match;
	FILE * fp;
	
	char szbuf[1024];
	
	sprintf(szbuf,"%s%c%s.psm.output",psmconf.m_strOutputPath.c_str(),cSlash,GetTitle(argv[1]).c_str());
	fp = fopen(szbuf,"w");
	if(fp == NULL)
	{
		cout << "in PSMAnalysis can't write to " << szbuf << endl;
		return 0;
	}
	match.Initialize(psmconf,fp);
	// cout << "init complete" << endl;
	
	// for xlink and common ion
	vector<MatchInfoMeasure> vXLinkIons;
	vector<MatchInfoMeasure> vCommonIons;
	
	vector<struct IonMatchInfo> vIonTypeMatchInfo;
	vector<struct IonMatchInfo> vCommonIonTypeMatchInfo;
	vector<struct IonMatchInfo> vXlinkIonTypeMatchInfo;
	
	vIonTypeMatchInfo.resize(psmconf.m_vIonTypes.size()+2);
	vCommonIonTypeMatchInfo.resize(psmconf.m_vIonTypes.size()+2);
	vXlinkIonTypeMatchInfo.resize(psmconf.m_vIonTypes.size()+2);
	
	for(size_t i = 0 ; i< psmconf.m_vIonTypes.size()+2; ++i)
	{
		vIonTypeMatchInfo[i].clear();
		vCommonIonTypeMatchInfo[i].clear();
		vXlinkIonTypeMatchInfo[i].clear();
	}

	double lfTop10Ratio = 0.0;
	double lf10PercentRatio = 0.0;
	double lfTotalRatio = 0.0;
	
	int nSpectraCount = 0;
	double lfNRErrorRatio = 0;
	for(int i=0;i<vResults.size();++i)
	{
		//cout << vSpectraOfResult[i].m_strFilePath << endl;
		if( mapSpectra.end()== mapSpectra.find(vSpectraOfResult[i].m_strFilePath))
			continue;
		//cout << "ok" << mapSpectra[vSpectraOfResult[i].m_strFilePath] << endl;
		match.SetSpectrum(vSpectra[mapSpectra[vSpectraOfResult[i].m_strFilePath]]);
		match.Preprocess();
		if(!bSpectraListAvail || mapSpectraList.find(vSpectra[i].m_strFilePath) != mapSpectraList.end() )
		{
			if(vResults[i].m_vPeptideResults.size()<=0)
			{
				continue;
			}
	
			match.SetPeptide(vResults[i].m_vPeptideResults[0]);
			//cout << "compute mz" << endl;
			match.ComputeMZ();
			//cout << "match" << endl;
			match.Match();
			//cout << "output " << endl;
			match.Output();
			
			lfNRErrorRatio += match.m_lfNRErrorRatio;
			
			double lfOdd1 = match.GetMatchOdd(0);
			double lfOdd2 = 0;
			//cout << vSpectra[i].m_strFilePath << endl;
			//cout << vResults[i].m_vPeptideResults[0].m_AlphaPeptide.m_szSequence << "	" << -lfOdd1 ;
			
			//cout << -lfOdd1 << "	" ;
			
			if(vResults[i].m_vPeptideResults[0].m_bPair)
			{
				lfOdd2 = match.GetMatchOdd(1);
				//cout << "	" << vResults[i].m_vPeptideResults[0].m_BetaPeptide.m_szSequence << "	" << -lfOdd2 << endl;
				//cout << -lfOdd2 << endl;
			}
			else
			{
				//cout << endl;
			}
			
			
			//if(lfOdd1 > -2 || lfOdd2 > -2)
				//cout << vSpectra[i].m_strFilePath << endl;
			//cout << "get psm match info" << endl;
			match.GetPSMMatchInfo();

			nSpectraCount++;

			for(size_t j = 0 ;j <= psmconf.m_vIonTypes.size() + 1; ++ j)
			{
				vIonTypeMatchInfo[j].nIonTypeOrder = j;
				vIonTypeMatchInfo[j].nMatchIonCount += match.m_vIonTypeMatchInfo[j].nMatchIonCount;
				vIonTypeMatchInfo[j].nTheoIonCount += match.m_vIonTypeMatchInfo[j].nTheoIonCount;
				vIonTypeMatchInfo[j].lfMatchCountRatio += match.m_vIonTypeMatchInfo[j].lfMatchCountRatio;
				vIonTypeMatchInfo[j].lfMatchGainRatio += match.m_vIonTypeMatchInfo[j].lfMatchGainRatio;
				vIonTypeMatchInfo[j].lfMatchAvrgItensity += match.m_vIonTypeMatchInfo[j].lfMatchAvrgItensity;
				vIonTypeMatchInfo[j].lfMatchAvrgTolerance += match.m_vIonTypeMatchInfo[j].lfMatchAvrgTolerance;
				vIonTypeMatchInfo[j].lfMatchContiRatio += match.m_vIonTypeMatchInfo[j].lfMatchContiRatio;
				vIonTypeMatchInfo[j].lfAvrgAAnum += match.m_vIonTypeMatchInfo[j].lfAvrgAAnum;
				
				if(match.m_vIonTypeMatchInfo[j].lfMatchAvrgItensity * match.m_vIonTypeMatchInfo[j].nMatchIonCount > vIonTypeMatchInfo[j].lfMaxItensity)
				{
					vIonTypeMatchInfo[j].lfMaxItensity = match.m_vIonTypeMatchInfo[j].lfMatchAvrgItensity * match.m_vIonTypeMatchInfo[j].nMatchIonCount;
					vIonTypeMatchInfo[j].strMaxItensityFileName = vSpectra[mapSpectra[vSpectraOfResult[i].m_strFilePath]].m_strFilePath;
					//vIonTypeMatchInfo[j].strMaxItensityFileName = vSpectra[i].m_strFilePath;
				}
				
				//if(match.m_vIonTypeMatchInfo[j].nTheoIonCount > 0)
				if(match.m_vIonTypeMatchInfo[j].nMatchIonCount > 0)
				{
					vIonTypeMatchInfo[j].nSpectraCnt ++ ;
				}
								
				if( j >= psmconf.m_vIonTypes.size())
					continue;
				
				vCommonIonTypeMatchInfo[j].nIonTypeOrder = j;
				vCommonIonTypeMatchInfo[j].nMatchIonCount += match.m_vCommonIonTypeMatchInfo[j].nMatchIonCount;
				vCommonIonTypeMatchInfo[j].nTheoIonCount += match.m_vCommonIonTypeMatchInfo[j].nTheoIonCount;
				vCommonIonTypeMatchInfo[j].lfMatchCountRatio += match.m_vCommonIonTypeMatchInfo[j].lfMatchCountRatio;
				vCommonIonTypeMatchInfo[j].lfMatchGainRatio += match.m_vCommonIonTypeMatchInfo[j].lfMatchGainRatio;
				vCommonIonTypeMatchInfo[j].lfMatchAvrgItensity += match.m_vCommonIonTypeMatchInfo[j].lfMatchAvrgItensity;
				vCommonIonTypeMatchInfo[j].lfMatchAvrgTolerance += match.m_vCommonIonTypeMatchInfo[j].lfMatchAvrgTolerance;
				vCommonIonTypeMatchInfo[j].lfMatchContiRatio += match.m_vCommonIonTypeMatchInfo[j].lfMatchContiRatio;
				vCommonIonTypeMatchInfo[j].lfAvrgAAnum += match.m_vCommonIonTypeMatchInfo[j].lfAvrgAAnum;
	
				vXlinkIonTypeMatchInfo[j].nIonTypeOrder = j;
				vXlinkIonTypeMatchInfo[j].nMatchIonCount += match.m_vXlinkIonTypeMatchInfo[j].nMatchIonCount;
				vXlinkIonTypeMatchInfo[j].nTheoIonCount += match.m_vXlinkIonTypeMatchInfo[j].nTheoIonCount;
				vXlinkIonTypeMatchInfo[j].lfMatchCountRatio += match.m_vXlinkIonTypeMatchInfo[j].lfMatchCountRatio;
				vXlinkIonTypeMatchInfo[j].lfMatchGainRatio += match.m_vXlinkIonTypeMatchInfo[j].lfMatchGainRatio;
				vXlinkIonTypeMatchInfo[j].lfMatchAvrgItensity += match.m_vXlinkIonTypeMatchInfo[j].lfMatchAvrgItensity;
				vXlinkIonTypeMatchInfo[j].lfMatchAvrgTolerance += match.m_vXlinkIonTypeMatchInfo[j].lfMatchAvrgTolerance;
				vXlinkIonTypeMatchInfo[j].lfMatchContiRatio += match.m_vXlinkIonTypeMatchInfo[j].lfMatchContiRatio;
				vXlinkIonTypeMatchInfo[j].lfAvrgAAnum += match.m_vXlinkIonTypeMatchInfo[j].lfAvrgAAnum;
			
						
				
				//if(match.m_vCommonIonTypeMatchInfo[j].nTheoIonCount > 0)
				if(match.m_vCommonIonTypeMatchInfo[j].nMatchIonCount > 0)
					vCommonIonTypeMatchInfo[j].nSpectraCnt ++ ;
				if(match.m_vXlinkIonTypeMatchInfo[j].nMatchIonCount > 0)
					vXlinkIonTypeMatchInfo[j].nSpectraCnt ++ ;
			}
			
			lfTop10Ratio += match.m_lfTop10Ratio;
			lf10PercentRatio += match.m_lf10PercentRatio;
			lfTotalRatio += match.m_lfTotalRatio;
		}
	}
	
	if(nSpectraCount > 0)
	{
		size_t i =0;
		for(i = 0 ;i <= psmconf.m_vIonTypes.size()+1; ++ i)
		{
			if(vIonTypeMatchInfo[i].nSpectraCnt > 0)
			{
				vIonTypeMatchInfo[i].nMatchIonCount /= vIonTypeMatchInfo[i].nSpectraCnt;
				vIonTypeMatchInfo[i].nTheoIonCount /= vIonTypeMatchInfo[i].nSpectraCnt;
				vIonTypeMatchInfo[i].lfMatchCountRatio /= vIonTypeMatchInfo[i].nSpectraCnt;
				vIonTypeMatchInfo[i].lfMatchGainRatio /= vIonTypeMatchInfo[i].nSpectraCnt;
				vIonTypeMatchInfo[i].lfMatchAvrgItensity /= vIonTypeMatchInfo[i].nSpectraCnt;
				vIonTypeMatchInfo[i].lfMatchAvrgTolerance /= vIonTypeMatchInfo[i].nSpectraCnt;
				vIonTypeMatchInfo[i].lfMatchContiRatio /= vIonTypeMatchInfo[i].nSpectraCnt;
				vIonTypeMatchInfo[i].lfAvrgAAnum /= vIonTypeMatchInfo[i].nSpectraCnt;
			}
			vIonTypeMatchInfo[i].lfSignificance = vIonTypeMatchInfo[i].lfMatchCountRatio*vIonTypeMatchInfo[i].lfMatchGainRatio*vIonTypeMatchInfo[i].lfMatchAvrgItensity;
			
			if( i >= psmconf.m_vIonTypes.size())
				continue;
			
			if(vCommonIonTypeMatchInfo[i].nSpectraCnt > 0)
			{
				vCommonIonTypeMatchInfo[i].nMatchIonCount /= vCommonIonTypeMatchInfo[i].nSpectraCnt;
				vCommonIonTypeMatchInfo[i].nTheoIonCount /= vCommonIonTypeMatchInfo[i].nSpectraCnt;
				vCommonIonTypeMatchInfo[i].lfMatchCountRatio /= vCommonIonTypeMatchInfo[i].nSpectraCnt;
				vCommonIonTypeMatchInfo[i].lfMatchGainRatio /= vCommonIonTypeMatchInfo[i].nSpectraCnt;
				vCommonIonTypeMatchInfo[i].lfMatchAvrgItensity /= vCommonIonTypeMatchInfo[i].nSpectraCnt;
				vCommonIonTypeMatchInfo[i].lfMatchAvrgTolerance /= vCommonIonTypeMatchInfo[i].nSpectraCnt;
				vCommonIonTypeMatchInfo[i].lfMatchContiRatio /= vCommonIonTypeMatchInfo[i].nSpectraCnt;
				vCommonIonTypeMatchInfo[i].lfAvrgAAnum /= vCommonIonTypeMatchInfo[i].nSpectraCnt;
			}
			vCommonIonTypeMatchInfo[i].lfSignificance = vCommonIonTypeMatchInfo[i].lfMatchCountRatio*vCommonIonTypeMatchInfo[i].lfMatchGainRatio*vCommonIonTypeMatchInfo[i].lfMatchAvrgItensity;

			if(vXlinkIonTypeMatchInfo[i].nSpectraCnt > 0)
			{
				vXlinkIonTypeMatchInfo[i].nMatchIonCount /= vXlinkIonTypeMatchInfo[i].nSpectraCnt;
				vXlinkIonTypeMatchInfo[i].nTheoIonCount /= vXlinkIonTypeMatchInfo[i].nSpectraCnt;
				vXlinkIonTypeMatchInfo[i].lfMatchCountRatio /= vXlinkIonTypeMatchInfo[i].nSpectraCnt;
				vXlinkIonTypeMatchInfo[i].lfMatchGainRatio /= vXlinkIonTypeMatchInfo[i].nSpectraCnt;
				vXlinkIonTypeMatchInfo[i].lfMatchAvrgItensity /= vXlinkIonTypeMatchInfo[i].nSpectraCnt;
				vXlinkIonTypeMatchInfo[i].lfMatchAvrgTolerance /= vXlinkIonTypeMatchInfo[i].nSpectraCnt;
				vXlinkIonTypeMatchInfo[i].lfMatchContiRatio /= vXlinkIonTypeMatchInfo[i].nSpectraCnt;
				vXlinkIonTypeMatchInfo[i].lfAvrgAAnum /= vXlinkIonTypeMatchInfo[i].nSpectraCnt;
			}
			vXlinkIonTypeMatchInfo[i].lfSignificance = vXlinkIonTypeMatchInfo[i].lfMatchCountRatio*vXlinkIonTypeMatchInfo[i].lfMatchGainRatio*vXlinkIonTypeMatchInfo[i].lfMatchAvrgItensity;

		}
		
		lfTop10Ratio /= nSpectraCount;
		lf10PercentRatio /= nSpectraCount;
		lfTotalRatio /= nSpectraCount;
		lfNRErrorRatio /= nSpectraCount;
	}

	sort(vIonTypeMatchInfo.begin(),vIonTypeMatchInfo.end(),ion_significance_lesser);
	sort(vCommonIonTypeMatchInfo.begin(),vCommonIonTypeMatchInfo.end(),ion_significance_lesser);
	sort(vXlinkIonTypeMatchInfo.begin(),vXlinkIonTypeMatchInfo.end(),ion_significance_lesser);
	
/*
	cout << "MatchCountRatio	" 
		<< "MatchGainRatio	"
		<< "MatchAverageIntensity	"
		<< "MatchAverageTolerance	"
		<< endl;
	for(size_t i = 0 ;i <= psmconf.m_vIonTypes.size()+1; ++ i)
	{
		cout << vIonTypeMatchInfo[i].nIonTypeOrder + 1<< "	"
			<< vIonTypeMatchInfo[i].lfSignificance << "	"
			<< vIonTypeMatchInfo[i].lfMatchCountRatio << "	"
			<< vIonTypeMatchInfo[i].lfMatchGainRatio << "	"
			<< vIonTypeMatchInfo[i].lfMatchAvrgItensity << "	"
			<< vIonTypeMatchInfo[i].lfMatchAvrgTolerance << "	"
			<< vIonTypeMatchInfo[i].lfMatchContiRatio << "	"
			<< vIonTypeMatchInfo[i].lfAvrgAAnum << "	"
			<< vIonTypeMatchInfo[i].strMaxItensityFileName << "	" 
			<< vIonTypeMatchInfo[i].nSpectraCnt 
			<< endl;
	}
	
	cout << "Summery : " << endl
		<< "10% : " << lf10PercentRatio << endl
		<< "Total : " << lfTotalRatio << endl
		<< "Preprocessing error ratio :" << lfNRErrorRatio << endl;
	
	cout << "FOR COMMON ION TYPES============================================" << endl;
	cout << "MatchCountRatio	" 
		<< "MatchGainRatio	"
		<< "MatchAverageIntensity	"
		<< "MatchAverageTolerance	"
		<< endl;
	for(size_t i = 0 ;i < psmconf.m_vIonTypes.size(); ++ i)
	{
		cout << "COMMON" << "	"
			<< vCommonIonTypeMatchInfo[i].nIonTypeOrder + 1<< "	"
			<< vCommonIonTypeMatchInfo[i].lfSignificance << "	"
			<< vCommonIonTypeMatchInfo[i].lfMatchCountRatio << "	"
			<< vCommonIonTypeMatchInfo[i].lfMatchGainRatio << "	"
			<< vCommonIonTypeMatchInfo[i].lfMatchAvrgItensity << "	"
			<< vCommonIonTypeMatchInfo[i].lfMatchAvrgTolerance << "	"
			<< vCommonIonTypeMatchInfo[i].lfMatchContiRatio << "	"
			<< vCommonIonTypeMatchInfo[i].lfAvrgAAnum << "	"
			//<< vCommonIonTypeMatchInfo[i].strMaxItensityFileName << "	" 
			<< vCommonIonTypeMatchInfo[i].nSpectraCnt 
			<< endl;
	}
	cout << "FOR XLINK ION TYPES============================================" << endl;
	cout << "MatchCountRatio	" 
		<< "MatchGainRatio	"
		<< "MatchAverageIntensity	"
		<< "MatchAverageTolerance	"
		<< endl;
	
	for(size_t i = 0 ;i < psmconf.m_vIonTypes.size(); ++ i)
	{
		cout << "XLINK" << "	"
			<< vXlinkIonTypeMatchInfo[i].nIonTypeOrder + 1<< "	"
			<< vXlinkIonTypeMatchInfo[i].lfSignificance << "	"
			<< vXlinkIonTypeMatchInfo[i].lfMatchCountRatio << "	"
			<< vXlinkIonTypeMatchInfo[i].lfMatchGainRatio << "	"
			<< vXlinkIonTypeMatchInfo[i].lfMatchAvrgItensity << "	"
			<< vXlinkIonTypeMatchInfo[i].lfMatchAvrgTolerance << "	"
			<< vXlinkIonTypeMatchInfo[i].lfMatchContiRatio << "	"
			<< vXlinkIonTypeMatchInfo[i].lfAvrgAAnum << "	"
			//<< vCommonIonTypeMatchInfo[i].strMaxItensityFileName << "	" 
			<< vXlinkIonTypeMatchInfo[i].nSpectraCnt 
			<< endl;
	}
	*/
	fclose(fp);
	
}
