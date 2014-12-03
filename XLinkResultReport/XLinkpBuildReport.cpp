#include <iostream>
#include <string>
 
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"

#include "XLinkResultReport.h"
#include "XLinkpBuildReport.h"

using namespace std;
using namespace proteomics_sdk;

#define INCLUDE_PROTEIN_SITE
#define INCLUDE_LINKER_ID
#define ONLY_FILTERED_RESULT

CXLinkpBuildReport::CXLinkpBuildReport()
{
	char szBuf[1024] = {0};
	getcwd(szBuf, 1024);
    m_strWorkDir = szBuf;
    if(m_strWorkDir[m_strWorkDir.length() - 1] != SLASH)
    {
    	m_strWorkDir += SLASH;
    }
}

CXLinkpBuildReport::~CXLinkpBuildReport()
{
}

void CXLinkpBuildReport::Init(string strpFindFile,time_t tmStartTime)
{
	m_strpFindFile = strpFindFile;
	m_tmStartTime = tmStartTime;
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

double CXLinkpBuildReport::_GetXLinkerMass(const CXLinkPepResult & pep_res)
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


string CXLinkpBuildReport::_getPepString(const CXLinkPepResult & pep_res)
{
	string strPep;
	char szbuf[10];
	if(pep_res.m_XLink.m_eXLinkType < 0)
		return strPep;
	strPep += pep_res.m_AlphaPeptide.m_szSequence;
	if(pep_res.m_XLink.m_eXLinkType < 1)
		return strPep;
	// fan 2011.12.10 siteno plus 1
//	sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tAlphaSite+1); // fan 2014.1.14 changed again
	sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tAlphaSite);// fan 2013.10.28 changed back. // fan 2014.7.16 changed back again. moved to XLinkPFDFileIO.
	strPep += szbuf;
	if(pep_res.m_XLink.m_eXLinkType < 2)
		return strPep;
	if(pep_res.m_XLink.m_eXLinkType == 2)
	{
		// fan 2011.12.10 siteno plus 1 //changed back
		sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tBetaSite); // fan 2014.1.14 changed again // fan 2014.7.16 changed back again. moved to XLinkPFDFileIO.
		strPep += szbuf;
	}
	else if(pep_res.m_XLink.m_eXLinkType == 3)
	{
		strPep += "-";
		strPep += pep_res.m_BetaPeptide.m_szSequence;
		// fan 2011.12.10 siteno plus 1 //changed back
		sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tBetaSite); // fan 2014.1.14 changed again // fan 2014.7.16 changed back again. moved to XLinkPFDFileIO.
		strPep += szbuf;
		/*
		size_t tSite ;
		if(strcmp(pep_res.m_BetaPeptide.m_szSequence,pep_res.m_AlphaPeptide.m_szSequence)==0 && pep_res.m_XLink.m_tAlphaSite > pep_res.m_XLink.m_tBetaSite)
		{
			strPep = pep_res.m_AlphaPeptide.m_szSequence;
			sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tBetaSite);
			strPep += szbuf;
			strPep += "-";
			strPep += pep_res.m_BetaPeptide.m_szSequence;
			sprintf(szbuf,"(%d)",pep_res.m_XLink.m_tAlphaSite);
			strPep += szbuf;
		}
		*/
	}
	
	return strPep;
}

void CXLinkpBuildReport::GetLines(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra ,string & strTXT)
{
	char szBuf[1024] = {'\0'};
		
	strTXT.clear();

	strTXT += "[Meta]\n";

	strTXT += "Enzyme_List=" + m_Condition.m_strEnzymeListPath + "\n";

	strTXT += "Modify_List=" + m_Condition.m_strModifyPath + "\n";

	strTXT += "AA_List=" + m_Condition.m_strAAListPath + "\n";

	strTXT += "Index_List=" + m_Condition.m_strDBConfPath + "\n";
		
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
	for(size_t k = 0; k < m_Condition.m_vSelectedVarMod.size() - m_Condition.m_vSelectedFixMod.size(); k++)
		strTXT += "," + m_Condition.m_vSelectedVarMod[k].m_strName;
	strTXT += "\n";
		
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

	int nResultCnt = 0;
	for(size_t i = 0 ;i < vResults.size() ; ++i)
	{
#ifdef ONLY_FILTERED_RESULT
		if(vResults[i].m_vPeptideResults.size() > 0)
#endif
		nResultCnt ++ ;
	}
	sprintf(szBuf, "%d", nResultCnt);
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
	
	int nResultId = 0;
	for(size_t i=0; i<vResults.size(); ++i)
	{
#ifdef ONLY_FILTERED_RESULT
		if(vResults[i].m_vPeptideResults.size() <= 0)
			continue;
#endif
		nResultId ++ ;
		sprintf(szBuf, "[Spectrum%d]\n", nResultId);
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

		string strTXT4results = "";
		int nResultCnt = 0;
		
		for(size_t rank_i = 0;rank_i < vResults[i].m_vPeptideResults.size();++rank_i)
		{
			nResultCnt ++;
			const CXLinkPepResult & pep_res = vResults[i].m_vPeptideResults[rank_i];

			sprintf(szBuf, "NO%d_Score=%.5f\n", nResultCnt, pep_res.m_lfScore);
			strTXT4results += szBuf;
	
			sprintf(szBuf, "NO%d_EValue=%E\n", nResultCnt, pep_res.m_lfEvalue);
			strTXT4results += szBuf;
	
			double lfMass1 = 0;
			if(pep_res.m_lfCalc_MH < 0.0001)
			{
				lfMass1 = CXLinkMatchResult::Calc_Theoretical_MH(pep_res, m_Condition.m_bPepMono);
				lfMass1 += _GetXLinkerMass(pep_res);
			}
			else
				lfMass1 = pep_res.m_lfCalc_MH;
			
			sprintf(szBuf, "NO%d_MH=%.5f\n", nResultCnt, lfMass1);
			strTXT4results += szBuf;
			
			string strSQ = "";

			strSQ = _getPepString(pep_res);

#ifdef INCLUDE_LINKER_ID
			sprintf(szBuf,":%d",pep_res.m_XLink.m_nLinkerId);
			strSQ += szBuf;
#endif
			sprintf(szBuf, "NO%d_SQ=%s\n", nResultCnt, strSQ.c_str());
			strTXT4results += szBuf;
	
			// report proteins like this 'pro1-pro2'
			
			vector<string> vProteinsAC;
			vector<string> vProteinId;
			vProteinsAC.clear();
			vProteinId.clear();
			vector<string>().swap(vProteinsAC);
			vector<string>().swap(vProteinId);

			string strTmp;
			size_t tTmp1 = 0,tTmp2 = 0;
			int nSite1,nSite2;
			for(size_t k = 0; k < pep_res.m_vAlphaProteinAC.size(); k++)
			{
				sprintf(szBuf, "%s", pep_res.m_vAlphaProteinAC[k].c_str());
				strTmp = szBuf;
				
				nSite1=nSite2=0;
				if(pep_res.m_vAlphaProteinSite.size()>k)
				{
					nSite1 = pep_res.m_vAlphaProteinSite[k];
					if(pep_res.m_XLink.m_tAlphaSite != -1)
					{
						nSite1 += pep_res.m_XLink.m_tAlphaSite;
					}
					//2012.1.20 fanshengbo 蛋白位点数加一，不从0开始算 //2014.7.19撤销+1，前面已加
//					nSite1++;
				}
				tTmp1 = pep_res.m_vAlphaProteinID[k];
				if(pep_res.m_XLink.m_eXLinkType == 0)
				{
					// here the site is the first ammino acid of the peptide
#ifdef INCLUDE_PROTEIN_SITE
					sprintf(szBuf,"%s(%d)",strTmp.c_str(),nSite1);
#else
					sprintf(szBuf,"%s",strTmp.c_str());
#endif
					vProteinsAC.push_back(szBuf);
					sprintf(szBuf,"%d",tTmp1);
					vProteinId.push_back(szBuf);
				}
				else if(pep_res.m_XLink.m_eXLinkType == 1)
				{
					// here the site is the linking site of the peptide
#ifdef INCLUDE_PROTEIN_SITE
					sprintf(szBuf,"%s(%d)",strTmp.c_str(),nSite1);
#else
					sprintf(szBuf,"%s",strTmp.c_str());
#endif
					vProteinsAC.push_back(szBuf);
					sprintf(szBuf,"%d",tTmp1);
					vProteinId.push_back(szBuf);
				}
				else if(pep_res.m_XLink.m_eXLinkType == 2)
				{
					if(pep_res.m_vAlphaProteinSite.size()>k)
					{
						nSite2 = pep_res.m_vAlphaProteinSite[k];
						if(pep_res.m_XLink.m_tBetaSite != -1)
						{
							nSite2 += pep_res.m_XLink.m_tBetaSite;
						}
						//2012.1.20 fanshengbo 蛋白位点数加一，不从0开始算 //2014.7.19撤销+1，前面已加
//						nSite2++;
					}
#ifdef INCLUDE_PROTEIN_SITE
					sprintf(szBuf,"%s(%d)(%d)",strTmp.c_str(),nSite1,nSite2);
#else
					sprintf(szBuf,"%s",strTmp.c_str());					
#endif
					vProteinsAC.push_back(szBuf);
					sprintf(szBuf,"%d",tTmp1);
					vProteinId.push_back(szBuf);
				}
				else
				{
					for(size_t l = 0; l < pep_res.m_vBetaProteinAC.size(); l++)
					{
						if(pep_res.m_vBetaProteinSite.size()>l)
						{
							nSite2 = pep_res.m_vBetaProteinSite[l];
							if(pep_res.m_XLink.m_tBetaSite != -1)
							{
								nSite2 += pep_res.m_XLink.m_tBetaSite;
							}
							//2012.1.20 fanshengbo 蛋白位点数加一，不从0开始算 //2014.7.19撤销+1，前面已加
//							nSite2++;
						}
#ifdef INCLUDE_PROTEIN_SITE
						sprintf(szBuf,"%s(%d)-%s(%d)",strTmp.c_str(),nSite1,pep_res.m_vBetaProteinAC[l].c_str(),nSite2);
#else
						sprintf(szBuf,"%s-%s",strTmp.c_str(),pep_res.m_vBetaProteinAC[l].c_str());
#endif
						vProteinsAC.push_back(szBuf);
						
						tTmp2 = pep_res.m_vBetaProteinID[l];
						sprintf(szBuf,"%d-%d",tTmp1,tTmp2);				
						vProteinId.push_back(szBuf);
					}
				}
			}
			
			sprintf(szBuf, "NO%d_Proteins=%d", nResultCnt, vProteinsAC.size());
			strTXT4results += szBuf;
			
			for(size_t k = 0;k < vProteinsAC.size() ; ++ k)
			{
				sprintf(szBuf,",%s",vProteinsAC[k].c_str());
				strTXT4results += szBuf;
			}
			
			strTXT4results += "\n";

	
			sprintf(szBuf, "NO%d_ProteinIDs=%d", nResultCnt, vProteinId.size());
			strTXT4results += szBuf;
			
			for(size_t k = 0; k < vProteinId.size(); k++)
			{
				sprintf(szBuf, ",%s",vProteinId[k].c_str());
				strTXT4results += szBuf;			
			}	
			
			strTXT4results += "\n";
						
			int nModCnt = 0;
			nModCnt = pep_res.m_AlphaPeptide.m_tModCnt;
			if(pep_res.m_bPair)
			{
				nModCnt += pep_res.m_BetaPeptide.m_tModCnt;
			}
			
			sprintf(szBuf, "NO%d_Modify_Pos=%d", nResultCnt, nModCnt);
			strTXT4results += szBuf;
			
			for(size_t k = 0; k < pep_res.m_AlphaPeptide.m_tModCnt; k++)
			{
				sprintf(szBuf, ",%d",pep_res.m_AlphaPeptide.m_tModSites[k][0]);
				strTXT4results += szBuf;			
			}
			
			if(pep_res.m_bPair)
			{
				for(size_t k = 0; k < pep_res.m_BetaPeptide.m_tModCnt; k++)
				{
					sprintf(szBuf, ",%d",pep_res.m_AlphaPeptide.m_tLength + pep_res.m_BetaPeptide.m_tModSites[k][0]);
					strTXT4results += szBuf;			
				}
			}
			strTXT4results += "\n";
	
			sprintf(szBuf, "NO%d_Modify_Name=%d", nResultCnt, nModCnt);
			strTXT4results += szBuf;
			for(size_t k = 0; k < pep_res.m_AlphaPeptide.m_tModCnt; k++)
			{
				sprintf(szBuf, ",%s",m_Condition.m_vSelectedVarMod[pep_res.m_AlphaPeptide.m_tModSites[k][1]].m_strName.c_str());
				strTXT4results += szBuf;			
			}
	
			if(pep_res.m_bPair)
			{
				for(size_t k = 0; k < pep_res.m_BetaPeptide.m_tModCnt; k++)
				{
					sprintf(szBuf, ",%s",m_Condition.m_vSelectedVarMod[pep_res.m_BetaPeptide.m_tModSites[k][1]].m_strName.c_str());
					strTXT4results += szBuf;			
				}
			}
			strTXT4results += "\n";
		}

		sprintf(szBuf, "ValidCandidate=%d\n", nResultCnt);
		strTXT += szBuf;
		
		strTXT += strTXT4results;
		
	}
	
}

void CXLinkpBuildReport::_parseline(string strLine,char cSep,vector<string > & vStrs)
{
	vStrs.clear();
	int n = 0;
	char szbuf[1024];
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

bool CXLinkpBuildReport::_getline(FILE * fp, string strTitle,string & strValue)
{
	strValue = "";
	if(fp == NULL)
		return false;
	
	char szbuf[10240];
	if(fgets(szbuf,10240,fp) == NULL)
		return false;
	
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

int CXLinkpBuildReport::_getModifyId(string strModName)
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

void CXLinkpBuildReport::_SetPeptideInfo(CXLinkPepResult & pep_res,string & strPepSeq)
{
	pep_res.m_XLink.m_eXLinkType = XLinkType(0);
	pep_res.m_bPair = false;
	pep_res.m_XLink.m_tAlphaSite = -1;
	pep_res.m_XLink.m_tBetaSite = -1;
	pep_res.m_XLink.m_nLinkerId = 0;
	
	size_t tpos = strPepSeq.find(':');
	size_t tend = strPepSeq.length();
	if(tpos != string::npos)
	{
		tend = tpos;
		string tmpstr = strPepSeq.substr(tpos+1,strPepSeq.length() - tpos - 1);
		pep_res.m_XLink.m_nLinkerId = atoi(tmpstr.c_str());
	}
	
	//NO1_SQ=KWIIPK(0)-KKCK(1):0
	tpos = strPepSeq.find('-');
	string strPepAlpha,strPepBeta;
	if(tpos != string::npos)
	{
		// type : xlink
		pep_res.m_XLink.m_eXLinkType = XLinkType(3);
		pep_res.m_bPair = true;
		strPepAlpha = strPepSeq.substr(0,tpos);
		strPepBeta = strPepSeq.substr(tpos+1,tend-tpos-1);
		size_t t1 = strPepAlpha.find('(');
		size_t t2 = strPepAlpha.find(')');
		string strSeq = strPepAlpha.substr(0,t1);
		pep_res.m_AlphaPeptide.SetPeptideInfor(strSeq.c_str(),strSeq.length(),0,0,0,false);
		string strSiteAlpha = strPepAlpha.substr(t1+1,t2-t1-1);
		pep_res.m_XLink.m_tAlphaSite = atoi(strSiteAlpha.c_str());
		t1 = strPepBeta.find('(');
		t2 = strPepBeta.find(')');
		strSeq = strPepBeta.substr(0,t1);
		pep_res.m_BetaPeptide.SetPeptideInfor(strSeq.c_str(),strSeq.length(),0,0,0,false);
		string strSiteBeta = strPepBeta.substr(t1+1,t2-t1-1);
		pep_res.m_XLink.m_tBetaSite = atoi(strSiteBeta.c_str());
	}
	else
	{
		pep_res.m_bPair = false;
		tpos = strPepSeq.find('(');
		if(tpos == string::npos)
		{
			// type : normal peptide
			//NO1_SQ=KKWIIPK
			pep_res.m_XLink.m_eXLinkType = XLinkType(0);
			string strSeq = strPepSeq.substr(0,tend);
			pep_res.m_AlphaPeptide.SetPeptideInfor(strSeq.c_str(),strSeq.length(),0,0,0,false);
			pep_res.m_XLink.m_tAlphaSite = -1;
			pep_res.m_XLink.m_tBetaSite = -1;
		}
		else
		{
			size_t t1 = tpos;
			size_t t2 = strPepSeq.find('(',t1+1);
			if(t2 == string::npos)
			{
				// type : mono j
				// NO1_SQ=KKWIIPK(0)
				pep_res.m_XLink.m_eXLinkType = XLinkType(1);
				string strSeq = strPepSeq.substr(0,t1);
				pep_res.m_AlphaPeptide.SetPeptideInfor(strSeq.c_str(),strSeq.length(),0,0,0,false);
				t2 = strPepSeq.find(')');
				string strSite = strPepSeq.substr(t1+1,t2-t1-1);
				pep_res.m_XLink.m_tAlphaSite = atoi(strSite.c_str());
				pep_res.m_XLink.m_tBetaSite = -1;
			}
			else
			{
				// NO1_SQ=KKWIIPK(0)(1)
				pep_res.m_XLink.m_eXLinkType = XLinkType(2);
				string strSeq = strPepSeq.substr(0,t1);
				pep_res.m_AlphaPeptide.SetPeptideInfor(strSeq.c_str(),strSeq.length(),0,0,0,false);
				size_t t = strPepSeq.find(')');
				string strSite = strPepSeq.substr(t1 + 1,t-t1-1);
				pep_res.m_XLink.m_tAlphaSite = atoi(strSite.c_str());
				t = strPepSeq.find(')',t2+1);
				strSite = strPepSeq.substr(t2+1,t - t2 -1);
				pep_res.m_XLink.m_tBetaSite = atoi(strSite.c_str());
			}
		}
	}
}

void CXLinkpBuildReport::_SetProteinInfo(CXLinkPepResult & pep_res,string & strProSeq)
{
	size_t tpos = 0;
	size_t tend = strProSeq.length();
	
	tpos = strProSeq.find('-');
	string strProAlpha,strProBeta;
	if(tpos != string::npos)
	{
		// type : xlink
		strProAlpha = strProSeq.substr(0,tpos);
		strProBeta = strProSeq.substr(tpos+1,tend-tpos-1);
		size_t i = 0;
		for(i = 0; i < pep_res.m_vAlphaProteinAC.size(); ++i)
		{
			if(strProAlpha == pep_res.m_vAlphaProteinAC[i])
				break;
		}
		if(i == pep_res.m_vAlphaProteinAC.size())
			pep_res.m_vAlphaProteinAC.push_back(strProAlpha);
		
		i = 0;
		for(i = 0; i < pep_res.m_vBetaProteinAC.size(); ++i)
		{
			if(strProBeta == pep_res.m_vBetaProteinAC[i])
				break;
		}
		if(i == pep_res.m_vBetaProteinAC.size())
			pep_res.m_vBetaProteinAC.push_back(strProBeta);
	}
	else
	{
		pep_res.m_vAlphaProteinAC.push_back(strProSeq);
	}
}

void CXLinkpBuildReport::_SetProteinIDInfo(CXLinkPepResult & pep_res,string & strProSeq)
{
	size_t tpos = 0;
	size_t tend = strProSeq.length();
	
	tpos = strProSeq.find('-');
	string strProAlpha,strProBeta;
	int nPro1 = 0, nPro2 = 0;
	if(tpos != string::npos)
	{
		// type : xlink
		strProAlpha = strProSeq.substr(0,tpos);
		nPro1 = atoi(strProAlpha.c_str());
		strProBeta = strProSeq.substr(tpos+1,tend-tpos-1);
		nPro2 = atoi(strProBeta.c_str());
		size_t i = 0;
		for(i = 0; i < pep_res.m_vAlphaProteinID.size(); ++i)
		{
			if(nPro1 == pep_res.m_vAlphaProteinID[i])
				break;
		}
		if(i == pep_res.m_vAlphaProteinID.size())
			pep_res.m_vAlphaProteinID.push_back(nPro1);
		
		i = 0;
		for(i = 0; i < pep_res.m_vBetaProteinID.size(); ++i)
		{
			if(nPro2 == pep_res.m_vBetaProteinID[i])
				break;
		}
		if(i == pep_res.m_vBetaProteinID.size())
			pep_res.m_vBetaProteinID.push_back(nPro2);
	}
	else
	{
		pep_res.m_vAlphaProteinID.push_back(nPro1);
	}
}

bool CXLinkpBuildReport::LoadFile(string strReportFile, vector<CXLinkMatchResult> & vResults , vector<CSpectrum> & vSpectra )
{
	vResults.clear();
	vSpectra.clear();
	
	char szline[1024];
	char szbuf[1024];
	FILE * fp = fopen(strReportFile.c_str(),"r");
	if(fp == NULL)
	{
		cout << "can't open report file " << strReportFile << endl;
		return false;
	}

	bool bConti = false;
	while((bConti = fgets(szline,1024,fp)))
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
	
	memset(pep_res.m_stMatchInfo.aPepConf, 0, sizeof(int) * 2*MAX_PEPTIDE_LENGTH);
	pep_res.m_stMatchInfo.lfMatchedSpecInt = 0.0;
	pep_res.m_stMatchInfo.lfUnMatchedSpecInt = 0.0;

	for(i = 0 ; i < tSpectraCount ; ++i)
	{
		bool bGood = true;
		string strVal;
		vector<string > vStrs;

		sprintf(szbuf,"[Spectrum%d]",i+1);
		while((bConti = fgets(szline,1024,fp)))
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
			
			sprintf(szbuf,"NO%d_MH=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			pep_res.m_lfCalc_MH = atof(strVal.c_str());

			pep_res.m_XLink.m_eXLinkType = XLinkType(0);
			
			pep_res.m_bPair = false;
			
			pep_res.m_XLink.m_tAlphaSite = -1;
			
			pep_res.m_XLink.m_tBetaSite = -1;
			
			sprintf(szbuf,"NO%d_SQ=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			
			_SetPeptideInfo(pep_res,strVal);
			
			pep_res.m_vAlphaProteinAC.clear();
			pep_res.m_vAlphaProteinID.clear();
			pep_res.m_vBetaProteinAC.clear();
			pep_res.m_vBetaProteinID.clear();
			pep_res.m_vAlphaProteinSite.clear();
			pep_res.m_vBetaProteinSite.clear();
			
			sprintf(szbuf,"NO%d_Proteins=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			_parseline(strVal,',',vStrs);
			for(int k=1;k<vStrs.size();++k)
				_SetProteinInfo(pep_res,vStrs[k]);
			
			sprintf(szbuf,"NO%d_ProteinIDs=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			_parseline(strVal,',',vStrs);
			for(int k=1;k<vStrs.size();++k)
			{
				_SetProteinIDInfo(pep_res,vStrs[k]);
			}
			
			sprintf(szbuf,"NO%d_Modify_Pos=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			_parseline(strVal,',',vStrs);
			pep_res.m_AlphaPeptide.m_tModCnt = atoi(vStrs[0].c_str());
			int nModCnt1 = 0, nModCnt2 = 0;
			vector<bool> vModPep;
			vModPep.clear();
			for(int k=1;k<vStrs.size();++k)
			{
				size_t tSite = atoi(vStrs[k].c_str());
				if(tSite >= pep_res.m_AlphaPeptide.m_tLength)
				{
					nModCnt2 ++;
					tSite -= pep_res.m_AlphaPeptide.m_tLength;
					pep_res.m_BetaPeptide.m_tModSites[nModCnt2-1][0] = tSite;
					// mod on beta peptide
					vModPep.push_back(false);
				}
				else
				{
					nModCnt1 ++;
					pep_res.m_AlphaPeptide.m_tModSites[nModCnt1-1][0] = tSite;
					vModPep.push_back(true);
				}
			}
			pep_res.m_AlphaPeptide.m_tModCnt = nModCnt1;
			pep_res.m_BetaPeptide.m_tModCnt = nModCnt2;
			
			sprintf(szbuf,"NO%d_Modify_Name=",j+1);
			bGood=_getline(fp,szbuf,strVal);
			if(!bGood)
				break;
			_parseline(strVal,',',vStrs);
			nModCnt1 = 0;
			nModCnt2 = 0;
			for(int k=1;k<vStrs.size();++k)
			{
				if(vModPep[k-1])
				{
					nModCnt1 ++;
					pep_res.m_AlphaPeptide.m_tModSites[nModCnt1-1][1] = _getModifyId(vStrs[k]);
				}
				else
				{
					nModCnt2 ++;
					pep_res.m_BetaPeptide.m_tModSites[nModCnt2-1][1] = _getModifyId(vStrs[k]);
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
	
	fclose(fp);
	return (i == tSpectraCount);
	
}
void CXLinkpBuildReport::WriteFile(const vector<CXLinkMatchResult> & vResults , const vector<CSpectrum> & vSpectra , string strOutputPath, string strTitle)
{
	string strTXTContent;
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
	string strFile = strOutputPath + strFileHead + strTitle + ".pbuild.txt";
	
	GetLines(vResults,vSpectra,strTXTContent);
	ofstream of(strFile.c_str());
	of<<strTXTContent;
	of.close();
	
}
void CXLinkpBuildReport::Close(void)
{
	
}

