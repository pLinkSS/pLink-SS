#include <iostream>
#include <string>
 
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"

#include "XLinkResultReport.h"
#include "XLinkProteinReport.h"

using namespace std;
using namespace proteomics_sdk;

#define MAXCHARINLINE 102400


CXLinkProteinReport::CXLinkProteinReport()
{
	char szBuf[MAXCHARINLINE] = {0};
	getcwd(szBuf, MAXCHARINLINE);
    m_strWorkDir = szBuf;
    if(m_strWorkDir[m_strWorkDir.length() - 1] != SLASH)
    {
    	m_strWorkDir += SLASH;
    }
    m_pTrace = CTrace::GetInstance();
}

CXLinkProteinReport::~CXLinkProteinReport()
{
}

void CXLinkProteinReport::Init(string strpFindFile,time_t tmStartTime)
{
	m_strpFindFile = strpFindFile;
	m_tmStartTime = tmStartTime;

	m_pTrace->Info("Protein report initializing...");
	CConditionReader reader(m_strpFindFile, m_strWorkDir);
	try
	{
		reader.Read();
	}
	catch(exception & e)
	{
		CErrInfo info("CXLinkpBuildReport", "Init", "in the function CConditionReader::Read.");
		info.Append("strOption=" + strpFindFile);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CXLinkpBuildReport", "Init", "caught an unknown exception in the function CConditionReader::Read.");
		info.Append("strOption=" + strpFindFile);
		info.Append("m_strWorkDir=" + m_strWorkDir);
		throw runtime_error(info.Get().c_str());
	}
	
	m_Condition = reader.m_Condition;
		
    for(size_t i = 0;i < m_Condition.m_vSelectedFixMod.size();++i)
    	m_Condition.m_vSelectedVarMod.push_back(m_Condition.m_vSelectedFixMod[i]);
    
    COptionTool * pOption = new COptionTool("spectrum", m_strpFindFile.c_str());
	m_strSpectraPath = pOption->GetString("spec_path", "null");
	delete pOption;
}

void CXLinkProteinReport::Close(void)
{
	m_Condition.clear();
}

double CXLinkProteinReport::_GetXLinkerMass(const CXLinkPepResult & pep_res)
{
	switch(pep_res.m_XLink.m_eXLinkType)
	{
	case NONE_LINK:
		return 0.0;
	case MONO_LINK:
		if(m_Condition.m_bPepMono)
		{
			return m_Condition.m_vSelectedXLinker[pep_res.m_XLink.m_nLinkerId].m_lfMLMonoMass_dif;
		}
		else
		{
			return m_Condition.m_vSelectedXLinker[pep_res.m_XLink.m_nLinkerId].m_lfMLAvrgMass_dif;
		}
		
	case LOOP_LINK:
		
	case X_LINK:
		if(m_Condition.m_bPepMono)
		{
			return m_Condition.m_vSelectedXLinker[pep_res.m_XLink.m_nLinkerId].m_lfMonoMass_dif;
		}
		else
		{
			return m_Condition.m_vSelectedXLinker[pep_res.m_XLink.m_nLinkerId].m_lfAvrgMass_dif;
		}
	default:
		return 0.0;
	}
}

int CXLinkProteinReport::_getModifyId(string strModName)
{
	for(int i=0;i<m_Condition.m_vSelectedVarMod.size();++i)
	{
		if(m_Condition.m_vSelectedVarMod[i].m_strName == strModName)
		{
			return i; 
		}
	}
	return -1;
}

void CXLinkProteinReport::_parseline(string strLine,char cSep,vector<string > & vStrs)
{
	vStrs.clear();
	int n = 0;
	char szbuf[MAXCHARINLINE];
	for(int k=0;k<strLine.length();++k)
	{
		if(strLine[k] == cSep)
		{
			szbuf[n] = 0;
			n = 0;
			vStrs.push_back(szbuf);
		}
		else
			szbuf[n++] = strLine[k];
	}
	szbuf[n] = 0;
	if(n)
	{
		vStrs.push_back(szbuf);
	}
}

bool CXLinkProteinReport::_getline(FILE * fp, string strTitle,string & strValue)
{
	strValue = "";
	if(fp == NULL)
		return false;
	
	char szbuf[MAXCHARINLINE];
	if(fgets(szbuf,MAXCHARINLINE,fp) == NULL)
	{
		return false;
	}
	
	if(szbuf[strlen(szbuf)-1] == 0x0a)
		szbuf[strlen(szbuf)-1] = 0;
	if(szbuf[strlen(szbuf)-1] == 0x0d)
		szbuf[strlen(szbuf)-1] = 0;

	if(strncmp(szbuf,strTitle.c_str(),strTitle.length())==0)
	{
		strValue = &szbuf[strTitle.length()];
		return true;
	}
	else
	{
		return false;
	}
}

bool CXLinkProteinReport::LoadFile(string strReportFile , vector<CXLinkMatchResult> & vResults,vector<CSpectrum> & vSpectra)
{
	cout<<"In the load file."<<endl;
	vResults.clear();
	vSpectra.clear();
	
	char szline[MAXCHARINLINE];
	char szbuf[MAXCHARINLINE];
	FILE * fp = fopen(strReportFile.c_str(),"r");
	if(fp == NULL)
	{
		cout << "can't open report file " << strReportFile << endl;
		return false;
	}

	bool bConti = false;
	while((bConti = fgets(szline,MAXCHARINLINE,fp)))
	{
		if(strncmp(szline,"Spectra=",strlen("Spectra="))==0)
			break;
	}
	
	if(szline[strlen(szline)-1] == 0x0a)
		szline[strlen(szline)-1] = 0;
	if(szline[strlen(szline)-1] == 0x0d)
		szline[strlen(szline)-1] = 0;
	
	size_t tSpectraCount = 0;
	if(bConti)
	{
		tSpectraCount = atoi(&szline[strlen("Spectra=")]);
		
	}
	
	size_t i;
	
	CXLinkMatchResult mr;
	CXLinkPepResult pep_res;
	pep_res.m_bEV = true;
	CSpectrum spec;

	for(i = 0 ; i < tSpectraCount ; ++i)
	{
		bool bGood = true;
		string strVal;
		vector<string > vStrs;

		sprintf(szbuf,"[Spectrum%d]",i+1);
		while((bConti = fgets(szline,MAXCHARINLINE,fp)))
		{
			if(strncmp(szline,szbuf,strlen(szbuf))==0)
				break;
		}
		
		if(!bConti)
			break;
		
		bGood=_getline(fp,"Input=",strVal);
		if(!bGood)
			break;
		spec.m_strFilePath = strVal;
		
		bGood=_getline(fp,"Charge=",strVal);
		if(!bGood)
			break;
		spec.m_nCharge = atoi(strVal.c_str());
		
		bGood=_getline(fp,"Intensity=",strVal);
		if(!bGood)
			break;
		spec.m_lfIntensity = atof(strVal.c_str());
		
		bGood=_getline(fp,"MH=",strVal);
		if(!bGood)
			break;
		spec.m_lfMH = atof(strVal.c_str());
		
		bGood=_getline(fp,"MZ=",strVal);
		if(!bGood)
			break;
		spec.m_lfMZ = atof(strVal.c_str());
		
		bGood=_getline(fp,"Candidate_Total=",strVal);
		if(!bGood)
			break;
		mr.m_tScore = atoi(strVal.c_str());
		
		bGood=_getline(fp,"ValidCandidate=",strVal);
		if(!bGood)
			break;
		mr.m_tRealCandidate = atoi(strVal.c_str());
		mr.m_vPeptideResults.clear();
			
		for(size_t j = 0;j<mr.m_tRealCandidate;++j)
		{
			sprintf(szbuf,"NO%d_Score=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_lfScore = atof(strVal.c_str());
			
			sprintf(szbuf,"NO%d_EValue=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_lfEvalue = atof(strVal.c_str());
			
			//added at 2014.9.10
			sprintf(szbuf,"NO%d_AlphaEValue=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_lfAlphaEvalue = atof(strVal.c_str());
			sprintf(szbuf,"NO%d_BetaEValue=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_lfBetaEvalue = atof(strVal.c_str());

			sprintf(szbuf,"NO%d_ContinousTag=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			_parseline(strVal,' ',vStrs);
			for(int k=0;k<vStrs.size();++k)
			{
				pep_res.m_stMatchInfo.aPepConf[k] = atoi(vStrs[k].c_str());
			}
			
			sprintf(szbuf,"NO%d_MatchedIntensity=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_stMatchInfo.lfMatchedSpecInt= atof(strVal.c_str());
			
			sprintf(szbuf,"NO%d_UnMatchedIntensity=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_stMatchInfo.lfUnMatchedSpecInt= atof(strVal.c_str());
			
			sprintf(szbuf,"NO%d_Mass=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_lfCalc_MH = atof(strVal.c_str());

			
			sprintf(szbuf,"NO%d_XLink_Type=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_XLink.m_eXLinkType = (XLinkType)atoi(strVal.c_str());
			
			if(pep_res.m_XLink.m_eXLinkType == 3)
				pep_res.m_bPair = true;
			else
				pep_res.m_bPair = false;
			
			sprintf(szbuf,"NO%d_Linker_Id=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_XLink.m_nLinkerId = atoi(strVal.c_str()); 

			sprintf(szbuf,"NO%d_XLink_Pos1=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_XLink.m_tAlphaSite = atoi(strVal.c_str());
			
			pep_res.m_XLink.m_tBetaSite = -1;
			if(pep_res.m_XLink.m_eXLinkType>=2)
			{
				sprintf(szbuf,"NO%d_XLink_Pos2=",j+1);
				
				bGood=_getline(fp,szbuf,strVal);
				if(!bGood)
					break;
				pep_res.m_XLink.m_tBetaSite = atoi(strVal.c_str());
			}
			
			sprintf(szbuf,"NO%d_Alpha_SQ=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_AlphaPeptide.SetPeptideInfor(strVal.c_str(),strVal.length(),0,0,0,false);
			
			
			sprintf(szbuf,"NO%d_Alpha_Proteins=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			_parseline(strVal,',',vStrs);
			for(int k=1;k<vStrs.size();++k)
				pep_res.m_vAlphaProteinAC.push_back(vStrs[k]);
			
			sprintf(szbuf,"NO%d_Alpha_ProteinIDs=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			_parseline(strVal,',',vStrs);
			for(int k=1;k<vStrs.size();++k)
				pep_res.m_vAlphaProteinID.push_back(atoi(vStrs[k].c_str()));
			
			sprintf(szbuf,"NO%d_Alpha_ProteinSites=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			_parseline(strVal,',',vStrs);
			for(int k=1;k<vStrs.size();++k)
				pep_res.m_vAlphaProteinSite.push_back(atoi(vStrs[k].c_str()));

			sprintf(szbuf,"NO%d_Alpha_Modify_Pos=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			_parseline(strVal,',',vStrs);
			pep_res.m_AlphaPeptide.m_tModCnt = atoi(vStrs[0].c_str());
			for(int k=1;k<vStrs.size();++k)
			{
				pep_res.m_AlphaPeptide.m_tModSites[k-1][0] = atoi(vStrs[k].c_str());
			}
			
			sprintf(szbuf,"NO%d_Alpha_Modify_Name=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			_parseline(strVal,',',vStrs);
			for(int k=1;k<vStrs.size();++k)
			{
				pep_res.m_AlphaPeptide.m_tModSites[k-1][1] = _getModifyId(vStrs[k]);
			}
			
			pep_res.m_BetaPeptide.SetPeptideInfor("",0,0,0,0,false);
			pep_res.m_BetaPeptide.m_tModCnt = 0;
			if(pep_res.m_bPair)
			{
				sprintf(szbuf,"NO%d_Beta_SQ=",j+1);
				bGood=_getline(fp,szbuf,strVal);
				if(!bGood)
					break;
				pep_res.m_BetaPeptide.SetPeptideInfor(strVal.c_str(),strVal.length(),0,0,0,false);
				
				vector<string > vStrs;
				sprintf(szbuf,"NO%d_Beta_Proteins=",j+1);
				bGood=_getline(fp,szbuf,strVal);
				if(!bGood)
					break;
				_parseline(strVal,',',vStrs);
				for(int k=1;k<vStrs.size();++k)
					pep_res.m_vBetaProteinAC.push_back(vStrs[k]);
				
				sprintf(szbuf,"NO%d_Beta_ProteinIDs=",j+1);
				bGood=_getline(fp,szbuf,strVal);
				if(!bGood)
					break;
				_parseline(strVal,',',vStrs);
				for(int k=1;k<vStrs.size();++k)
					pep_res.m_vBetaProteinID.push_back(atoi(vStrs[k].c_str()));
				
				sprintf(szbuf,"NO%d_Beta_ProteinSites=",j+1);
				bGood=_getline(fp,szbuf,strVal);
				if(!bGood)
					break;
				_parseline(strVal,',',vStrs);
				for(int k=1;k<vStrs.size();++k)
					pep_res.m_vBetaProteinSite.push_back(atoi(vStrs[k].c_str()));

				sprintf(szbuf,"NO%d_Beta_Modify_Pos=",j+1);
				bGood=_getline(fp,szbuf,strVal);
				if(!bGood)
					break;
				_parseline(strVal,',',vStrs);
				pep_res.m_BetaPeptide.m_tModCnt = atoi(vStrs[0].c_str());
				for(int k=1;k<vStrs.size();++k)
				{
					pep_res.m_BetaPeptide.m_tModSites[k-1][0] = atoi(vStrs[k].c_str());
				}
				
				sprintf(szbuf,"NO%d_Beta_Modify_Name=",j+1);
				bGood=_getline(fp,szbuf,strVal);
				if(!bGood)
					break;
				_parseline(strVal,',',vStrs);
				for(int k=1;k<vStrs.size();++k)
				{
					pep_res.m_BetaPeptide.m_tModSites[k-1][1] = _getModifyId(vStrs[k]);
				}
				
			}
			
			mr.m_vPeptideResults.push_back(pep_res);
			
			vector<string>().swap(pep_res.m_vAlphaProteinAC);
			vector<size_t>().swap(pep_res.m_vAlphaProteinID);
			vector<string>().swap(pep_res.m_vBetaProteinAC);
			vector<size_t>().swap(pep_res.m_vBetaProteinID);
			vector<size_t>().swap(pep_res.m_vAlphaProteinSite);
			vector<size_t>().swap(pep_res.m_vBetaProteinSite);
		}
		
		if(!bGood)
		{
			cout << i << endl;
			cout << szbuf << endl << "val = " << strVal << endl;
			break;
		}
		

		vResults.push_back(mr);
		vSpectra.push_back(spec);
	}
	

	cout<<"In the load file."<<endl;
	
	fclose(fp);
	return (i == tSpectraCount);

}


void CXLinkProteinReport::GetLines(const vector<CXLinkMatchResult> & vResults, const vector<CSpectrum> & vSpectra, string& strTXT)
{
	
	char szBuf[PATH_MAX * 100 ] = {'\0'};

	strTXT.clear();

	strTXT += "[Meta]\n";

	strTXT += "Enzyme_List=" + m_Condition.m_strEnzymeListPath + "\n";

	strTXT += "Modify_List=" + m_Condition.m_strModifyPath + "\n";

	strTXT += "AA_List=" + m_Condition.m_strAAListPath + "\n";

	strTXT += "Index_List=" + m_Condition.m_strDBConfPath + "\n";
	
	strTXT += "xlink_List=" + m_Condition.m_strXLinkerPath + "\n";

	strTXT += "[Search]\n";

	strTXT += "InputPath=" + m_strSpectraPath + "\n";

	sprintf(szBuf, "Database=%d", m_Condition.m_vSelectedDBName.size());
	strTXT += szBuf;
	for(size_t k = 0; k < m_Condition.m_vSelectedDBName.size(); k++)
		strTXT += "," + m_Condition.m_vSelectedDBName[k];
	strTXT += "\n";

	time_t ltime;
	time(&ltime); 
	sprintf(szBuf,"Time=%s",ctime( &ltime ));
	strTXT +=  szBuf;
	
	// output time spent
	time_t tmEndTime = time(NULL);
	double lfDur = difftime(tmEndTime,m_tmStartTime);
	sprintf(szBuf,"TimeCost=%f sec\n",lfDur);
	strTXT += szBuf;
	
	strTXT +=  "Enzyme=" + m_Condition.m_SelectedEnzyme.GetName() + "\n";

	sprintf(szBuf, "Fixed_modifications=%d", m_Condition.m_vSelectedFixMod.size());
	strTXT += szBuf;
	for(size_t k = 0; k < m_Condition.m_vSelectedFixMod.size(); k++)
		strTXT += "," + m_Condition.m_vSelectedFixMod[k].m_strName;
	strTXT += "\n";
	
	sprintf(szBuf, "Variable_modifications=%d", m_Condition.m_vSelectedVarMod.size() - m_Condition.m_vSelectedFixMod.size());
	strTXT += szBuf;
	
	size_t tNum = 0;
	if(m_Condition.m_vSelectedVarMod.size() < m_Condition.m_vSelectedFixMod.size())
		tNum = m_Condition.m_vSelectedVarMod.size();
	else
		tNum = m_Condition.m_vSelectedVarMod.size()- m_Condition.m_vSelectedFixMod.size();
		
	for(size_t k = 0; k < tNum ; k++)
		strTXT += "," + m_Condition.m_vSelectedVarMod[k].m_strName;
	strTXT += "\n";
	
	sprintf(szBuf,"XLinkerNum=%d\n",m_Condition.m_vSelectedXLinker.size());
	strTXT += szBuf;
	for(size_t k = 0 ; k< m_Condition.m_vSelectedXLinker.size() ; ++ k)
	{
		sprintf(szBuf,"XLinker%d=%s\n",k,m_Condition.m_vSelectedXLinker[k].m_strName.c_str());
		strTXT += szBuf;
	}

	if(m_Condition.m_bPepMono)
		strTXT +=  "Peptide_Mass=Monoisotopic\n";
	else
		strTXT +=  "Peptide_Mass=Average\n";

	// modify by emily
	sprintf(szBuf,"Peptide_Mass_Tolerance_Wnd_Num=%d\n",m_Condition.m_vPepTolWnds.size());
	for(size_t k = 0 ; k < m_Condition.m_vPepTolWnds.size() ; ++ k)
	{
		sprintf(szBuf, "Peptide_Mass_Tolerance_Base%d=%f\n", k,-m_Condition.m_vPepTolWnds[k].m_lfPepTolBase);
		strTXT += szBuf;
		sprintf(szBuf, "Peptide_Mass_Tolerance%d=%f\n", k,m_Condition.m_vPepTolWnds[k].m_lfPepTol);
		strTXT += szBuf;
		sprintf(szBuf,"Peptide_Mass_Tolerance_Base_Type%d=%s\n",k,m_Condition.m_vPepTolWnds[k].m_strPepTolBaseType.c_str());
		strTXT += szBuf;
		sprintf(szBuf,"Peptide_Mass_Tolerance_Type%d=%s\n",k,m_Condition.m_vPepTolWnds[k].m_strPepTolType.c_str());
		strTXT += szBuf;
	}

	if(m_Condition.m_bFragmentMono)
		strTXT +=  "Fragment_Mass=Monoisotopic\n";
	else
		strTXT +=  "Fragment_Mass=Average\n";

	// modify by emily
	sprintf(szBuf, "Fragment_Mass_Tolerance_Base=%f\n", -m_Condition.m_lfFragmentTolBase);
	strTXT +=  szBuf;
	sprintf(szBuf, "Fragment_Mass_Tolerance=%f\n", m_Condition.m_lfFragmentTol);
	strTXT +=  szBuf;
	strTXT +=  "Fragment_Mass_Tolerance_Base_Type=" + m_Condition.m_strFragmentTolBaseType + "\n";
	strTXT +=  "Fragment_Mass_Tolerance_Type=" + m_Condition.m_strFragmentTolType + "\n";

	sprintf(szBuf, "Max_Missed_Cleavages=%d\n", m_Condition.m_nMaxMissCleaves);
	strTXT +=  szBuf;
		
	strTXT +=  "Instrument_type=" + m_Condition.GetInstrumentType() + "\n";

	strTXT += "[Total]\n";

	sprintf(szBuf, "%d", vResults.size());
	strTXT += "Spectra=" + string(szBuf) + "\n";

	strTXT += "Proteins=0\n";

	sprintf(szBuf, "FPR=%f\n", m_Condition.m_lfFPRThreshold);
	strTXT += szBuf;

	double lfThreshold = 1;
	if(!vResults.empty())
	{
		vector<CXLinkPepResult>::const_iterator itBegin = vResults[vResults.size()-1].m_vPeptideResults.begin();
		if(!vResults.empty() && !vResults[vResults.size()-1].m_vPeptideResults.empty())
			lfThreshold = itBegin->m_lfEvalue;
	}

	sprintf(szBuf, "Threshold=%E\n", lfThreshold);
	strTXT += szBuf;
	


	for(size_t i=0; i<vResults.size(); ++i)
	{
		
		sprintf(szBuf, "[Spectrum%d]\n", i+1);
		strTXT += szBuf;
		
		sprintf(szBuf, "Input=%s\n", vSpectra[i].m_strFilePath.c_str());
		strTXT += szBuf;
		
		sprintf(szBuf, "Charge=%d\n", vSpectra[i].m_nCharge);
		strTXT += szBuf;

		sprintf(szBuf, "Intensity=%.5f\n", vSpectra[i].m_lfIntensity);
		strTXT += szBuf;

		sprintf(szBuf, "MH=%.5f\n", vSpectra[i].m_lfMH);
		strTXT += szBuf;

		sprintf(szBuf, "MZ=%.5f\n", vSpectra[i].m_lfMZ);
		strTXT += szBuf;

		sprintf(szBuf, "Candidate_Total=%u\n", vResults[i].m_tScore);
		strTXT += szBuf;

		sprintf(szBuf, "ValidCandidate=%d\n", vResults[i].m_vPeptideResults.size());
		strTXT += szBuf;

#ifdef EMILY_DEBUG
		sprintf(szBuf, "RealCandidate=%d\n", vResults[i].m_tRealCandidate);
		strTXT += szBuf;
#endif

		if(vResults[i].m_vPeptideResults.empty())
			continue;
		
		for(size_t rank_i = 0;rank_i < vResults[i].m_vPeptideResults.size();++rank_i)
		//for(size_t rank_i = 0;rank_i < 1;++rank_i)
		{
			const CXLinkPepResult & pep_res = vResults[i].m_vPeptideResults[rank_i];
	
			// todo : emily debug
			sprintf(szBuf, "NO%d_Score=%.5f\n", rank_i + 1, pep_res.m_lfScore);
			strTXT += szBuf;
			
			sprintf(szBuf, "NO%d_EValue=%E\n", rank_i + 1, pep_res.m_lfEvalue);
			strTXT += szBuf;
			
			//added at 2014.9.10 Alpha/Beta e-value
			sprintf(szBuf, "NO%d_AlphaEValue=%E\n", rank_i + 1, pep_res.m_lfAlphaEvalue);
			strTXT += szBuf;
			sprintf(szBuf, "NO%d_BetaEValue=%E\n", rank_i + 1, pep_res.m_lfBetaEvalue);
			strTXT += szBuf;

			// matching info
			// calculate continous amino num (the minimum of pep1 and pep2)
			sprintf(szBuf,"NO%d_ContinousTag=",rank_i + 1);
			strTXT += szBuf;
			size_t tLength = pep_res.m_AlphaPeptide.m_tLength;
			if(pep_res.m_bPair)
				tLength += pep_res.m_BetaPeptide.m_tLength;
			
			for(size_t k = 0 ; k<tLength; ++ k)
			{
				sprintf(szBuf,"%d ",pep_res.m_stMatchInfo.aPepConf[k]);
				strTXT += szBuf;
			}
			strTXT += "\n";

			
			sprintf(szBuf,"NO%d_MatchedIntensity=%f\n",rank_i + 1 ,pep_res.m_stMatchInfo.lfMatchedSpecInt);
			strTXT += szBuf;
			
			sprintf(szBuf,"NO%d_UnMatchedIntensity=%f\n",rank_i + 1 ,pep_res.m_stMatchInfo.lfUnMatchedSpecInt);
			strTXT += szBuf;
		
			double lfMass1 = 0,lfMass2 = 0;
			
			if(pep_res.m_lfCalc_MH <= 0.000001)
			{
				lfMass1 = CXLinkMatchResult::Calc_Theoretical_MH(pep_res, m_Condition.m_bPepMono);
				lfMass1 += _GetXLinkerMass(pep_res);
			}
			else
				lfMass1 = pep_res.m_lfCalc_MH;
			//lfMass2 = _GetDeltaMass(lfMass1,vSpectra[i].m_lfMH);
			
			sprintf(szBuf, "NO%d_Mass=%.5f\n", rank_i + 1, lfMass1);
			strTXT += szBuf;
			
			sprintf(szBuf,"NO%d_XLink_Type=%d\n", rank_i + 1, pep_res.m_XLink.m_eXLinkType);
			strTXT += szBuf;

			sprintf(szBuf,"NO%d_Linker_Id=%d\n", rank_i + 1, pep_res.m_XLink.m_nLinkerId);
			strTXT += szBuf;
			
			sprintf(szBuf,"NO%d_XLink_Pos1=%d\n", rank_i + 1, pep_res.m_XLink.m_tAlphaSite);
			strTXT += szBuf;
			
			if(pep_res.m_XLink.m_eXLinkType >=2)
			{
				sprintf(szBuf,"NO%d_XLink_Pos2=%d\n", rank_i + 1, pep_res.m_XLink.m_tBetaSite);
				strTXT += szBuf;
			}
			
			sprintf(szBuf, "NO%d_Alpha_SQ=%s\n", rank_i + 1, pep_res.m_AlphaPeptide.m_szSequence);
			//sprintf(szBuf, "NO1_SQ=%c.%s.%c\n", pep_res.m_peptide.m_cPrev, pep_res.m_peptide.m_szSequence, pep_res.m_peptide.m_cNext);
			strTXT += szBuf;
			
			sprintf(szBuf, "NO%d_Alpha_Proteins=%d", rank_i + 1, pep_res.m_vAlphaProteinAC.size());
			strTXT += szBuf;
	
			for(size_t k = 0; k < pep_res.m_vAlphaProteinAC.size(); k++)
			{
				sprintf(szBuf, ",%s", pep_res.m_vAlphaProteinAC[k].c_str());
				strTXT += szBuf;			
			}
	
			strTXT += "\n";
	
			sprintf(szBuf, "NO%d_Alpha_ProteinIDs=%d", rank_i + 1, pep_res.m_vAlphaProteinID.size());
			strTXT += szBuf;
			for(size_t k = 0; k < pep_res.m_vAlphaProteinID.size(); k++)
			{
				sprintf(szBuf, ",%d",pep_res.m_vAlphaProteinID[k]);
				strTXT += szBuf;			
			}	
	
			strTXT += "\n";
	
			// site of the first AA of peptide in protein sequence
			sprintf(szBuf, "NO%d_Alpha_ProteinSites=%d", rank_i + 1, pep_res.m_vAlphaProteinSite.size());
			strTXT += szBuf;
			for(size_t k = 0; k < pep_res.m_vAlphaProteinSite.size(); k++)
			{
				sprintf(szBuf, ",%d",pep_res.m_vAlphaProteinSite[k]);
				strTXT += szBuf;			
			}	
			strTXT += "\n";

			sprintf(szBuf, "NO%d_Alpha_Modify_Pos=%d", rank_i + 1, pep_res.m_AlphaPeptide.m_tModCnt);
			strTXT += szBuf;
			for(size_t k = 0; k < pep_res.m_AlphaPeptide.m_tModCnt; k++)
			{
				sprintf(szBuf, ",%d",pep_res.m_AlphaPeptide.m_tModSites[k][0]);
				strTXT += szBuf;			
			}
			

	
			strTXT += "\n";
	
			sprintf(szBuf, "NO%d_Alpha_Modify_Name=%d", rank_i + 1, pep_res.m_AlphaPeptide.m_tModCnt);
			strTXT += szBuf;
			for(size_t k = 0; k < pep_res.m_AlphaPeptide.m_tModCnt; k++)
			{
				sprintf(szBuf, ",%s",m_Condition.m_vSelectedVarMod[pep_res.m_AlphaPeptide.m_tModSites[k][1]].m_strName.c_str());
				strTXT += szBuf;			
			}
	
			strTXT += "\n";
			
			if(pep_res.m_bPair)
			{
				sprintf(szBuf, "NO%d_Beta_SQ=%s\n", rank_i + 1, pep_res.m_BetaPeptide.m_szSequence);
				//sprintf(szBuf, "NO1_SQ=%c.%s.%c\n", pep_res.m_peptide.m_cPrev, pep_res.m_peptide.m_szSequence, pep_res.m_peptide.m_cNext);
				strTXT += szBuf;
				
				sprintf(szBuf, "NO%d_Beta_Proteins=%d", rank_i + 1, pep_res.m_vBetaProteinAC.size());
				strTXT += szBuf;
	
				for(size_t k = 0; k < pep_res.m_vBetaProteinAC.size(); k++)
				{
					sprintf(szBuf, ",%s", pep_res.m_vBetaProteinAC[k].c_str());
					strTXT += szBuf;			
				}
	
				strTXT += "\n";
	
				sprintf(szBuf, "NO%d_Beta_ProteinIDs=%d", rank_i + 1, pep_res.m_vBetaProteinID.size());
				strTXT += szBuf;
				for(size_t k = 0; k < pep_res.m_vBetaProteinID.size(); k++)
				{
					sprintf(szBuf, ",%d",pep_res.m_vBetaProteinID[k]);
					strTXT += szBuf;			
				}	
	
				strTXT += "\n";
	
				// site of the first AA of peptide in protein sequence
				sprintf(szBuf, "NO%d_Beta_ProteinSites=%d", rank_i + 1, pep_res.m_vBetaProteinSite.size());
				strTXT += szBuf;
				for(size_t k = 0; k < pep_res.m_vBetaProteinSite.size(); k++)
				{
					sprintf(szBuf, ",%d",pep_res.m_vBetaProteinSite[k]);
					strTXT += szBuf;			
				}	
				strTXT += "\n";

				sprintf(szBuf, "NO%d_Beta_Modify_Pos=%d", rank_i + 1, pep_res.m_BetaPeptide.m_tModCnt);
				strTXT += szBuf;
				for(size_t k = 0; k < pep_res.m_BetaPeptide.m_tModCnt; k++)
				{
					sprintf(szBuf, ",%d",pep_res.m_BetaPeptide.m_tModSites[k][0]);
					strTXT += szBuf;			
				}
	
				strTXT += "\n";
	
				sprintf(szBuf, "NO%d_Beta_Modify_Name=%d", rank_i + 1, pep_res.m_BetaPeptide.m_tModCnt);
				strTXT += szBuf;
				for(size_t k = 0; k < pep_res.m_BetaPeptide.m_tModCnt; k++)
				{
					sprintf(szBuf, ",%s",m_Condition.m_vSelectedVarMod[pep_res.m_BetaPeptide.m_tModSites[k][1]].m_strName.c_str());
					strTXT += szBuf;			
				}
	
				strTXT += "\n";
				
			}
		}
	}
}

void CXLinkProteinReport::WriteFile(const vector<CXLinkMatchResult> & vResults, const vector<CSpectrum> & vSpectra, string strOutputPath, string strTitle)
{
	string strTXTContent;
	GetLines(vResults,vSpectra,strTXTContent);

	size_t t = m_strpFindFile.find(".pfind");
	if(string::npos == t)
	{
		t = m_strpFindFile.find(".PFIND");
		if(string::npos == t)
		{
			t = m_strpFindFile.find(".");
		}
	}
	int nStart = m_strpFindFile.length();
	while(nStart >= 0 && m_strpFindFile[nStart] != '\\' && m_strpFindFile[nStart] != '/')
	{
		--nStart;
	}
	++nStart;
	if(string::npos != t)
		t -= nStart;
	string strFileHead = m_strpFindFile.substr((size_t)nStart, t);
	if('\"' == strOutputPath[0] && '\"' == strOutputPath[strOutputPath.length() - 1])
	{
		strOutputPath = strOutputPath.substr(1, strOutputPath.length() - 2);
	}
	
	if(SLASH != strOutputPath[strOutputPath.length()-1])
	{
		strOutputPath += SLASH;
	}
	string strFile = strOutputPath + strFileHead + strTitle + "_qry.proteins.txt";
	ofstream of(strFile.c_str());
	of<<strTXTContent;
	of.close();
}




