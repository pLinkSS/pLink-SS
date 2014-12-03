#include <iostream>
#include <cstdio>
#include <string>
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "Trace.h"
#include "Condition.h"
#include "ConditionReader.h"

using namespace proteomics_sdk;
using namespace std;

#ifdef WIN32
const char cSlash = '\\';
#else
const char cSlash = '/';
#endif

namespace proteomics_sdk
{
CConditionReader::CConditionReader(const string & strOption,
		const string & strWorkDir):
	m_strOption(strOption), m_strWorkDir(strWorkDir)
{
}

CConditionReader::~CConditionReader()
{
}

void CConditionReader::Read()
{
   CTrace * pTrace = CTrace::GetInstance();
   pTrace->Info("Loading database settings ...", MODULE_CONDITION);
	_LoadDatabase();

   pTrace->Info("Loading enzyme settings ...", MODULE_CONDITION);
   _LoadEnzyme();

   pTrace->Info("Loading PTM settings ...", MODULE_CONDITION);
   _LoadModify();

   pTrace->Info("Loading XLink settings ...", MODULE_CONDITION);
	_LoadXLinker();
	
   pTrace->Info("Loading ion type settings ...", MODULE_CONDITION);
	_LoadIons();

   pTrace->Info("Loading workflow settings ...", MODULE_CONDITION);
	_LoadFlow();
	
	pTrace->Info("Loading spectra settings ...", MODULE_CONDITION);
	_LoadSpectra();
	
	pTrace->Info("Loading cluster settings ...", MODULE_CONDITION);
	_LoadCluster();

	// _LoadMSearchInfo(); /* abandoned, assured by Jeremy */
	
	pTrace->Info("Loading meta information...", MODULE_CONDITION);
	_LoadMetaDat(); /* add by czhou */

	pTrace->Info("Loading triple peptide parameters...", MODULE_CONDITION);
}

void CConditionReader::Test(const CCondition & c)
{
   CTrace * pTrace = CTrace::GetInstance();
   pTrace->Debug("Selected enzyme name: " + c.m_SelectedEnzyme.m_strName);
}

/* the following members of CCondition has been set
 *     index_content
 *     m_eIndexContent
 *     m_eIndex
 *     m_strIndexPath
 *     m_strDBConfPath
 *     m_strEnzymeListPath
 *     m_vAllEnzyme
 *     m_strModifyPath
 *     m_vAllModification
 *     m_strXLinkerPath
 *     m_vAllXLinker
 *     m_strAAListPath
 *     m_mapAA
 *     m_vSelectedDBName */
void CConditionReader::_LoadDatabase(void)
{
    COptionTool * pOption(NULL);
    pOption = new COptionTool("database", m_strOption.c_str());

    /* index_content */
    string strIdxContent = pOption->GetString("index_content", "PEPTIDE_INDEX");
    IndexContent eIdxCon = PEPTIDE_PAIR_OPEN;
    if ("PEPTIDE_PAIR" == strIdxContent)
    {
    	eIdxCon = PEPTIDE_PAIR;
    }
    else if ("PEPTIDE_PAIR_OPEN" == strIdxContent)
    {
    	eIdxCon = PEPTIDE_PAIR_OPEN;
    	
    	CTrace * pTrace = CTrace::GetInstance();
    	pTrace->Info("Loading simple scoring parameters...\n");
    	pTrace->Info("\t\t  Using open flow for crossLinking identification");
    	_LoadSimpleScore();
    	pTrace->Debug("_LoadSimpleScore() Over");
    }
    m_Condition.m_eIndexContent = eIdxCon;

    /* index_type */
    string strIndexType = pOption->GetString("index_type", "MEMORY_MAPPED");
    IndexType eIndex(MEMORY_MAPPED);
    if("DISK_FILE" == strIndexType)
    {
        eIndex = DISK_FILE;
    }
    else
    {
        eIndex = MEMORY_MAPPED;
    }
    m_Condition.m_eIndex = eIndex;

    /* ----- unused: index_path_flow1 zhangkun czhou */
    string strIndexPath = pOption->GetString("index_path", ".\\task.pindex");
    m_Condition.m_strIndexPath = strIndexPath;
    
    /* path of dbconfig.ini */
    string strTemp = m_strWorkDir + "dbconf.ini";
    string strDBConf = pOption->GetString("index_List", strTemp.c_str());

    /* path of enzyme.ini */
    strTemp = m_strWorkDir + "enzyme.ini";
    string strEnzymeList = pOption->GetString("enzyme_List", strTemp.c_str());

    /* path of modify.ini */
    strTemp = m_strWorkDir + "modify.ini";
    string strModifyList = pOption->GetString("modify_List", strTemp.c_str());
    
    /* path of xlink.ini */
    strTemp = m_strWorkDir + "xlink.ini";
    string strXLinkerList = pOption->GetString("xlink_List", strTemp.c_str());

    /* path of aa.ini */
    strTemp = m_strWorkDir + "aa.ini";
    string strAAList = pOption->GetString("aa_List", strTemp.c_str());

    try
    {
    	m_Condition.m_strDBConfPath = strDBConf;

		CEnzymeConf enzyme_conf(strEnzymeList);
		m_Condition.m_strEnzymeListPath = strEnzymeList;
		m_Condition.m_vAllEnzyme = enzyme_conf.GetEnzymeList();

		CModifyConf modify_conf(strModifyList);
		m_Condition.m_strModifyPath = strModifyList;
		m_Condition.m_vAllModification = modify_conf.GetModifyList();

		if(strXLinkerList != "")
		{
			CXLinkerConf xlinker_conf(strXLinkerList);
			m_Condition.m_strXLinkerPath = strXLinkerList;
			m_Condition.m_vAllXLinker = xlinker_conf.GetXLinkerList();
		}
		else
			m_Condition.m_vAllXLinker.clear();

		CAAConf aa_conf(strAAList);
		m_Condition.m_strAAListPath = strAAList;
		m_Condition.m_mapAA = aa_conf.GetMapAAMass();
    }
    catch(exception & e){
        delete pOption;

        CErrInfo info("CConditionReader", "_LoadDatabase");
        info.Append("ModifyList=" + strModifyList);
        info.Append("EnzymeList=" + strEnzymeList);
        info.Append("AAList=" + strAAList);
        info.Append("DBConfPath=" + strDBConf);
        throw runtime_error(info.Get(e).c_str());
    }
    catch(...)
    {
        delete pOption;

        CErrInfo info("CConditionReader", "_LoadDatabase", "caught an unknown exception.");
        info.Append("ModifyList=" + strModifyList);
        info.Append("EnzymeList=" + strEnzymeList);
        info.Append("AAList=" + strAAList);
        info.Append("DBConfPath=" + strDBConf);
        throw runtime_error(info.Get().c_str());
    }

    /* database name */
    char szBuf[BUFFER_SIZE] = {0};
    m_Condition.m_vSelectedDBName.clear();
    int nTotal = pOption->GetInteger("db_total",1);
    for(int i=0; i<nTotal; ++i)
    {
    	sprintf(szBuf, "db_name%d",i+1);
        string strDB = pOption->GetString(szBuf, "null");
        if("null" == strDB)
        {
            delete pOption;
            CErrInfo info("CConditionReader", "_LoadDatabase",
            		"cannot find the database name.");
            throw runtime_error(info.Get().c_str());
        }
        m_Condition.m_vSelectedDBName.push_back(strDB);
    }
    CTrace * pTrace = CTrace::GetInstance();
    string strLine = "Selected database(s):";
    for(size_t i = 0;i < m_Condition.m_vSelectedDBName.size();++i)
    {
    	strLine += " " + m_Condition.m_vSelectedDBName[i];
    }
    pTrace->Info(strLine);
    delete pOption;
}

void CConditionReader::_LoadEnzyme(void)
{
    COptionTool * pOption(NULL);
    pOption = new COptionTool("enzyme", m_strOption.c_str());

    /* set enzyme */
    string strEnzyme = pOption->GetString("enzyme_name", "null");
    if("null" == strEnzyme)
    {
        delete pOption;
        CErrInfo info("CConditionReader", "_LoadEnzyme",
        		"cannot find the enzyme name.");
        throw runtime_error(info.Get().c_str());
    }

    for(size_t j = 0; j<m_Condition.m_vAllEnzyme.size(); j++)
    {
        if(m_Condition.m_vAllEnzyme[j].m_strName == strEnzyme)
        {
        	m_Condition.m_SelectedEnzyme = m_Condition.m_vAllEnzyme[j];
        	break;
        }
    }

    delete pOption;
}

void CConditionReader::_LoadXLinker(void)
{
    COptionTool * pOption = new COptionTool("xlink", m_strOption.c_str());
    if(pOption == NULL)
    	return ;
    string strKey;

    m_Condition.m_vSelectedXLinker.clear();
    char szBuf[BUFFER_SIZE] = {0};
    int nSelLinkerNum = pOption->GetInteger("linker_total", 1);
    if(nSelLinkerNum > MAX_LINKER_NUM)
    	nSelLinkerNum = MAX_LINKER_NUM;
    for(int i = 0; i < nSelLinkerNum ; ++ i)
    {
    	if(i == 0)
    		sprintf(szBuf, "linker");
    	else
    		sprintf(szBuf, "linker%d",i);
    	
	    string strLinker = pOption->GetString(szBuf, "null");
	    if(strLinker != "null")
	    {
	    	size_t j=0;
	        for(j = 0; j<m_Condition.m_vAllXLinker.size(); j++)
	        {
	            if(m_Condition.m_vAllXLinker[j].m_strName == strLinker)
	            {
	            	m_Condition.m_vSelectedXLinker.push_back(m_Condition.m_vAllXLinker[j]);
	            	break;
	            }
	        }
	        
	        if(j == m_Condition.m_vAllXLinker.size())
	        {
	        	delete pOption;
		        CErrInfo info("CConditionReader", "_LoadXLinker", "cannot find the linker " + strLinker);
		        throw runtime_error(info.Get().c_str());
	        }
	    }
    }
   delete pOption;
}

void CConditionReader::_LoadModify(void)
{
    COptionTool * pOption = new COptionTool("modify", m_strOption.c_str());

    if(pOption == NULL)
    	return ;
    
    vector<CModification> vM = m_Condition.m_vAllModification;

    string strKey;

    //fix modification
    vector<CModification> vM1;
    int nFixTotal = pOption->GetInteger("fix_total",0);
    char szBuf[BUFFER_SIZE] = {0};
    for(int i=0; i<nFixTotal; ++i)
    {
    	sprintf(szBuf, "fix_mod%d",i+1);
        string strFixMod = pOption->GetString(szBuf, "null");
        if(strFixMod != "null")
        {
            for(size_t j = 0; j<vM.size(); j++)
            {
                if(vM[j].m_strName == strFixMod)
                {
                    vM1.push_back(vM[j]);
                }
            }	
        }
        else
        {
            delete pOption;
            CErrInfo info("CConditionReader", "_LoadModify", "cannot find the modification " + string(szBuf));
            throw runtime_error(info.Get().c_str());
        }
    }
    m_Condition.m_vSelectedFixMod = vM1;

    //variable modification
    vector<CModification> vM2;
    int nVarTotal = pOption->GetInteger("var_total",0);
    for(int i=0; i<nVarTotal; ++i)
    {
    	sprintf(szBuf, "var_mod%d",i+1);
    	string strVarMod = pOption->GetString(szBuf, "null");
        if(strVarMod!="null")
        {
            for(size_t j = 0; j<vM.size(); j++)
            {
                if(vM[j].m_strName == strVarMod)
                    vM2.push_back(vM[j]);
            }
        }
        else
         {
             delete pOption;
             CErrInfo info("CConditionReader", "_LoadModify", "cannot find the modification " + string(szBuf));
             throw runtime_error(info.Get().c_str());
         }
    }
    m_Condition.m_vSelectedVarMod = vM2;

    //max_number
    m_Condition.m_nMaxModifyNumber = pOption->GetInteger("max_number", MAX_MODIFY_NUM);
    if(m_Condition.m_nMaxModifyNumber > MAX_MODIFY_NUM)
    	m_Condition.m_nMaxModifyNumber = MAX_MODIFY_NUM;

    delete pOption;
    
}

void CConditionReader::_LoadIons(void)
{
    COptionTool * pOption = new COptionTool("ions", m_strOption.c_str());

    int nMiss= pOption->GetInteger("max_miss", 2);
    m_Condition.m_nMaxMissCleaves = nMiss;
    
    bool bMono = pOption->GetInteger("peptide_mono", 1);
    // modify by emily
    string s;
    
	// multiple mass tolerance windows
    int nTolWndNum = pOption->GetInteger("peptide_tol_total", 1);
    if(nTolWndNum > MAX_PEP_TOL_WND_NUM)
    	nTolWndNum = MAX_PEP_TOL_WND_NUM;
    m_Condition.m_vPepTolWnds.clear();
    for(int i = 0 ;i < nTolWndNum ; ++ i)
    {
    	PEP_TOL_WND stPepTolWnd ;
    	string strTmp ;
    	char szbuf[200] = {0};
    	if(i == 0)
    		sprintf(szbuf,"peptide_tol");
    	else
    		sprintf(szbuf,"peptide_tol%d",i);	
    	strTmp = pOption->GetString(szbuf, "3.0");
    	stPepTolWnd.m_lfPepTol = atof(strTmp.c_str());
    	
    	if(i == 0)
    		sprintf(szbuf,"peptide_tol_type");
    	else
    		sprintf(szbuf,"peptide_tol_type%d",i);
    	strTmp = pOption->GetString(szbuf, "Da");
    	stPepTolWnd.m_strPepTolType = strTmp;
    	
    	if(i == 0)
    		sprintf(szbuf,"peptide_tol_base");
    	else
    		sprintf(szbuf,"peptide_tol_base%d",i);
    	strTmp = pOption->GetString(szbuf, "0.0");
    	stPepTolWnd.m_lfPepTolBase = -atof(strTmp.c_str());
    	
    	if(i == 0)
    		sprintf(szbuf,"peptide_tol_base_type");
    	else
    		sprintf(szbuf,"peptide_tol_base_type%d",i);
    	strTmp = pOption->GetString(szbuf, "Da");
    	stPepTolWnd.m_strPepTolBaseType = strTmp;
    	
    	m_Condition.m_vPepTolWnds.push_back(stPepTolWnd);
    }

    bMono=(bool) pOption->GetInteger("fragment_mono", 1);
    
    s=pOption->GetString("fragment_tol_base", "0.0");
    double lfTolBase = 0;
    sscanf(s.c_str(),"%lf",&lfTolBase);
    s=pOption->GetString("fragment_tol", "1.0");
    double lfTol = 0;
    sscanf(s.c_str(),"%lf",&lfTol);
    string strTolType = pOption->GetString("fragment_tol_type", "Da");
    m_Condition.m_bFragmentMono = bMono;
    m_Condition.m_lfFragmentTolBase = -lfTolBase;
    m_Condition.m_lfFragmentTol = lfTol;
    m_Condition.m_strFragmentTolType = strTolType;
    
    string strTolBaseType = pOption->GetString("fragment_tol_base_type", "Da");
    m_Condition.m_strFragmentTolBaseType = strTolBaseType;

    s=pOption->GetString("min_peptide_tol", "0");
    lfTol = 0;
    sscanf(s.c_str(),"%lf",&lfTol);
    strTolType = pOption->GetString("min_peptide_tol_type", "Da");
    m_Condition.m_lfMinPepTol = lfTol;
    m_Condition.m_strMinPepTolType = strTolType;

	int nIonTypes = pOption->GetInteger("ion_type_total", 0);
	if(nIonTypes > 0)
	{
		char cType = 'b';
		int nCharge = 1;
		int nLossNH3 = 0, nLossH2O = 0;
		double lfOther = 0.0; 
		int nIntraConti = 0;
		int nInterContiWndN = -1;
		int nInterContiWndC = -1;
		int nMaxInterConti = -1;
		int nContainLinker = 0;

		CIonType ion_type;
		m_Condition.m_vIonTypes.clear();
		 
		for(int i = 0;i < nIonTypes;++i)
		{
			char temp[200] = {0};
			sprintf(temp, "ion_type%d", i + 1);
			string strInfor = pOption->GetString(temp, "null");
			if("null" == strInfor)
			{
	            delete pOption;
	            CErrInfo info("CConditionReader", "_LoadIons", "cannot find the iontype " + string(temp));
	            throw runtime_error(info.Get().c_str());
			}
			double lfElement[ELEMENT_NUMBER] = {0};
			if(m_Condition.m_bFragmentMono)
			{
				lfElement[0] = IonMass_Mono_C;
				lfElement[1] = IonMass_Mono_H;
				lfElement[2] = IonMass_Mono_N;
				lfElement[3] = IonMass_Mono_O;
			}
			else
			{
				lfElement[0] = IonMass_Aver_C;
				lfElement[1] = IonMass_Aver_H;
				lfElement[2] = IonMass_Aver_N;
				lfElement[3] = IonMass_Aver_O;
			}
			lfOther = 0;
			istringstream iss(strInfor);
			iss >> cType >> nCharge >> nLossNH3 >> nLossH2O >> lfOther >> nIntraConti >> nInterContiWndN >> nInterContiWndC >> nContainLinker;
			lfOther += m_Condition.m_bFragmentMono ? NH3Mass_Mono * nLossNH3 : NH3Mass_Aver * nLossNH3;
			lfOther += m_Condition.m_bFragmentMono ? H2OMass_Mono * nLossH2O : H2OMass_Aver * nLossH2O;
			
			ion_type.bIntraContinous = (nIntraConti==1);
			ion_type.nInterContinousWnd1 = nInterContiWndN;
			ion_type.nInterContinousWnd2 = nInterContiWndC;
			ion_type.nContainLinker = nContainLinker;
					
			if(nInterContiWndN > nMaxInterConti)
				nMaxInterConti = nInterContiWndN;
			if(nInterContiWndC > nMaxInterConti)
				nMaxInterConti = nInterContiWndC;
			
			ion_type.cSymbol = cType;// added at 2013.11.21
			switch(cType)
			{
			case 'a':
				ion_type.cType = 0;
				ion_type.bNTerm = true;
				lfOther += lfElement[0] + lfElement[3];
				break;
			case 'b':
				ion_type.cType = 0;
				ion_type.bNTerm = true;
				break;
			case 'c':
				ion_type.cType = 0;
				ion_type.bNTerm = true;
				lfOther -= lfElement[2] + 3 * lfElement[1];
				break;
			case 'x':
				ion_type.cType = 1;
				ion_type.bNTerm = false;
				lfOther -= lfElement[0] + lfElement[3] - 2*lfElement[1];
				break;
			case 'y':
				ion_type.cType = 1;
				ion_type.bNTerm = false;
				break;
			case 'z':
				ion_type.cType = 1;
				ion_type.bNTerm = false;
				lfOther += lfElement[2] + 2 * lfElement[1];
				break;
			case 'm':
				ion_type.cType = 2;
				ion_type.bNTerm = false;
				break;
			case 'p':
				ion_type.cType = 3;
				ion_type.bNTerm = false;
				break;
			case 'q':
				ion_type.cType = 4;
				ion_type.bNTerm = false;
				break;
			case 'n':
				// KL: ya-NH3 at the linker site
				ion_type.cType = 5;
				ion_type.bNTerm = false;
				break;
			case 'u'://yP
				ion_type.cType = 6;
				ion_type.bNTerm = false;
				break;
			case 'v'://bP
				ion_type.cType = 7;
				ion_type.bNTerm = true;
				break;
			default:
	            delete pOption;
	            CErrInfo info("CConditionReader", "_LoadIons", string("cannot find the iontype ") + cType);
	            throw runtime_error(info.Get().c_str());
			}
			ion_type.nCharge = nCharge;
			ion_type.nTotalLostVal = (int)(lfOther * MZMULTIPLIER + 0.5);
			m_Condition.m_vIonTypes.push_back(ion_type);
		}
		
		m_Condition.m_tInterContinuousWndNum = nMaxInterConti + 1;
	}
	else
	{
		m_Condition.SetDefaultIons();
	}	
    delete pOption;
    
}

void CConditionReader::_LoadMSearchInfo()
{
	// these variables are for msearch
	COptionTool * pOption = new COptionTool("refined_search", m_strOption.c_str());
	
	// whether save info to path , for the first round of msearch , should set this value to 1
	int nSaveInfo = pOption->GetInteger("saveinfo", 0);
	if(nSaveInfo)
	{
		m_Condition.m_bRSSaveInfo = true;
	}
	else
	{
		m_Condition.m_bRSSaveInfo = false;
	}

	// whether load info from path , for the recursive round of msearch , should set this value to 1
	int nLoadInfo = pOption->GetInteger("loadinfo", 0);
	if(nLoadInfo)
	{
		m_Condition.m_bRSLoadInfo = true;
	}
	else
	{
		m_Condition.m_bRSLoadInfo = false;
	}

	// the path to save the info or load from the info
	char szBuf[256];
	sprintf(szBuf,"temp%c",cSlash);
	string strPath = m_strWorkDir + string(szBuf);
	m_Condition.m_strRSPath = pOption->GetString("path", strPath.c_str());
    if(m_Condition.m_strRSPath[m_Condition.m_strRSPath.length() - 1] != cSlash)
    {
    	m_Condition.m_strRSPath += cSlash;
    }
    // FPR to filter results from the first round of search

    int nFDR =pOption->GetInteger("FDR", 10);
    double lfFDR = (double) nFDR / 100;
    if( 1 < lfFDR)
        lfFDR = 0.05f;
    else if(0.01f > lfFDR)
        lfFDR = 0.01f;
    m_Condition.m_lfRSFDR = lfFDR;
    
    // dbname : protein db index name created in the first round of search
    string strDBName = "temp";
    m_Condition.m_strRSDBName = pOption->GetString("dbname", strDBName.c_str());
    
    // pindex : *.pindex file name used to create protein index in the first round of search
    string strPindexPath = m_strWorkDir + "temp.pindex";
    m_Condition.m_strRSPindexPath = pOption->GetString("pindex_path", strPindexPath.c_str());
    
    // tag : to name all the temp files generated by the first round of search
    m_Condition.m_strRSTag = pOption->GetString("tag", "temp");
    
    
    /*
     * work for refined-search:
     * if loadinfo = 1 then 
     * 	1. load the db created by the first round of search
     *  2. close inner tol and pep gen
     * 
     */
    if(m_Condition.m_bRSLoadInfo)
    {
    	if(m_Condition.m_vSelectedDBName.size())
    	{
    		m_Condition.m_vSelectedDBName.clear();
    	}
    	m_Condition.m_vSelectedDBName.push_back(m_Condition.m_strRSDBName);
    	
    	m_Condition.m_lfMinPepTol = 0;
    	m_Condition.m_tMinScoreNum = 0;
    	
    }

}


// add for xlink open flow
void CConditionReader::_LoadSimpleScore(void)
{
	COptionTool * pOption = new COptionTool("simplescore", m_strOption.c_str());

	CTrace * pTrace = CTrace::GetInstance();
	
	//preproc
    string strPreProc = pOption->GetString("preproc", "PRE_PROC_COASE");
    pTrace->Info("preprocess module for cross-linking open flow : " + strPreProc);
    if("PRE_PROC_XLINK_HCD" == strPreProc)
    {
    	m_Condition.m_eSimplePreProcMethod = PRE_PROC_XLINK_HCD;
    }
    else
    {
    	m_Condition.m_eSimplePreProcMethod = PRE_PROC_XLINK_HCD;
    }
	
    
    //score cutoff
	string strSimpleScoreCutoff;
	strSimpleScoreCutoff = pOption->GetString("simple_score_cutoff", "2.0");
	pTrace->Info("score cutoff for cross-linking open flow : " + strSimpleScoreCutoff);
	m_Condition.m_lfSimpleScoreCutoff = atof(strSimpleScoreCutoff.c_str());
	
	// save pep num
	char szbuf[200];
	m_Condition.m_nSimpleReportPep = pOption->GetInteger("report_peptide_number",500);
	sprintf(szbuf,"save %d candidate peptides for cross-linking open flow \n",m_Condition.m_nSimpleReportPep);
	pTrace->Info(szbuf);
	
	// load ion types
	pTrace->Info("loading ion types for cross-linking open flow");
	int nIonTypes = pOption->GetInteger("ion_type_total", 0);
	if(nIonTypes > 0)
	{
		char cType = 'b';
		int nCharge = 1;
		int nLossNH3 = 0, nLossH2O = 0;
		double lfOther = 0.0; 
		CIonType ion_type;
		m_Condition.m_vSimpleIonTypes.clear();
		for(int i = 0;i < nIonTypes;++i)
		{
			char temp[200] = {0};
			sprintf(temp, "ion_type%d", i + 1);
			string strInfor = pOption->GetString(temp, "null");
			if("null" == strInfor)
			{
	            delete pOption;
	            CErrInfo info("CConditionReader", "_LoadSimpleScore", "cannot find the iontype " + string(temp));
	            throw runtime_error(info.Get().c_str());
			}
			double lfElement[ELEMENT_NUMBER] = {0};
			if(m_Condition.m_bFragmentMono)
			{
				lfElement[0] = IonMass_Mono_C;
				lfElement[1] = IonMass_Mono_H;
				lfElement[2] = IonMass_Mono_N;
				lfElement[3] = IonMass_Mono_O;
			}
			else
			{
				lfElement[0] = IonMass_Aver_C;
				lfElement[1] = IonMass_Aver_H;
				lfElement[2] = IonMass_Aver_N;
				lfElement[3] = IonMass_Aver_O;
			}
			lfOther = 0;
			istringstream iss(strInfor);
			iss >> cType >> nCharge >> nLossNH3 >> nLossH2O >> lfOther;
			lfOther += m_Condition.m_bFragmentMono ? NH3Mass_Mono * nLossNH3 : NH3Mass_Aver * nLossNH3;
			lfOther += m_Condition.m_bFragmentMono ? H2OMass_Mono * nLossH2O : H2OMass_Aver * nLossH2O;
			
			ion_type.bIntraContinous = true;
			ion_type.nInterContinousWnd1 = -1;
			ion_type.nInterContinousWnd2 = -1;
			ion_type.nContainLinker = 0;
			
			switch(cType)
			{
			case 'a':
				ion_type.bNTerm = true;
				ion_type.cType = 0;
				lfOther += lfElement[0] + lfElement[3];
				break;
			case 'b':
				ion_type.bNTerm = true;
				ion_type.cType = 0;
				break;
			case 'c':
				ion_type.bNTerm = true;
				ion_type.cType = 0;
				lfOther -= lfElement[2] + 3 * lfElement[1];
				break;
			case 'x':
				ion_type.bNTerm = false;
				ion_type.cType = 1;
				lfOther -= lfElement[0] + lfElement[3] - 2*lfElement[1];
				break;
			case 'y':
				ion_type.bNTerm = false;
				ion_type.cType = 1;
				break;
			case 'z':
				ion_type.bNTerm = false;
				ion_type.cType = 1;
				lfOther += lfElement[2] + 2 * lfElement[1];
				break;
			default:
	            delete pOption;
	            CErrInfo info("CConditionReader", "_LoadSimpleScore", string("cannot find the iontype ") + cType);
	            throw runtime_error(info.Get().c_str());
			}
			ion_type.nCharge = nCharge;
			ion_type.nTotalLostVal = (int)(lfOther * MZMULTIPLIER + 0.5);
			m_Condition.m_vSimpleIonTypes.push_back(ion_type);
		}
	}
	
	delete pOption;
}

void CConditionReader::_LoadFlow(void)
{
    COptionTool * pOption = new COptionTool("flow", m_strOption.c_str());
    
    //processor_num
    size_t tProcessorNum = pOption->GetInteger("processor_num", 1);
    m_Condition.m_tProcessorNum = tProcessorNum;
    
    // modify by emily
    //min_score_num
    size_t tMinScoreNum = pOption->GetInteger("min_score_num", MIN_SCORE_NUM);
    m_Condition.m_tMinScoreNum = tMinScoreNum;
    
    //max_score_num
    size_t tMaxScoreNum = pOption->GetInteger("max_score_num", MAX_SCORE_NUM);
    m_Condition.m_tMaxScoreNum = tMaxScoreNum;

    //flow
    m_Condition.m_eFlow = FLOW_ALL_IDX;
    string strFlow = pOption->GetString("flow", "FLOW_ALL_IDX");

    if("FLOW_PRO_IDX" == strFlow)
    	m_Condition.m_eFlow = FLOW_PRO_IDX;
    else if("FLOW_XLINK" == strFlow)
    	m_Condition.m_eFlow = FLOW_XLINK;

    //salvo_batch_size
    m_Condition.m_tSalvoBatchSize = pOption->GetInteger("salvo_batch_size", SALVO_BATCH_SIZE);

    if(m_Condition.m_tSalvoBatchSize <= 0)
    	m_Condition.m_tSalvoBatchSize = SALVO_BATCH_SIZE;

    //instrument
    string strInstrument = pOption->GetString("instrument", "UNKNOWN");

    if("ESI-QUAD-TOF" == strInstrument)
    {
        m_Condition.m_eInstrumentType = (INSTRUMENT_ESI_QUAD_TOF);
    }
    else if("MALDI-QUAD-TOF" == strInstrument)
    {
    	m_Condition.m_eInstrumentType = (INSTRUMENT_MALDI_QUAD_TOF);
    }
    else if("ESI-TRAP" == strInstrument)
    {
    	m_Condition.m_eInstrumentType = (INSTRUMENT_ESI_TRAP);
    }
    else if("ESI-QUAD" == strInstrument)
    {
    	m_Condition.m_eInstrumentType = (INSTRUMENT_ESI_QUAD);
    }
    else if("ESI-FTICR" == strInstrument)
    {
    	m_Condition.m_eInstrumentType = (INSTRUMENT_ESI_FTICR);
    }
    else if("ESI-4SECTOR" == strInstrument)
    {
    	m_Condition.m_eInstrumentType = (INSTRUMENT_ESI_4SECTOR);
    }
    else if("MALDI-TOF-PSD" == strInstrument)
    {
    	m_Condition.m_eInstrumentType = (INSTRUMENT_MALDI_TOF_PSD);
    }
    else if("MALDI-TOF-TOF" == strInstrument)
    {
    	m_Condition.m_eInstrumentType = (INSTRUMENT_MALDI_TOF_TOF);
    }
    else if("MALDI-QIT-TOF" == strInstrument)
    {
    	m_Condition.m_eInstrumentType = (INSTRUMENT_MALDI_QIT_TOF);
    }
    else if("FTMS-ECD" == strInstrument)
    {
    	m_Condition.m_eInstrumentType = (INSTRUMENT_FTMS_ECD);
    }
    else
    {
    	m_Condition.m_eInstrumentType = (INSTRUMENT_UNKNOWN);
    }
    
    //preproc
    string strPreProc = pOption->GetString("preproc", "PRE_PROC_DEFAULT");
    if("PRE_PROC_XLINK_HCD" == strPreProc)
    {
    	m_Condition.m_ePreProcMethod = PRE_PROC_XLINK_HCD;
    }
    else
    {
    	m_Condition.m_ePreProcMethod = PRE_PROC_XLINK_HCD;
    }
    //score
    string strScore  = pOption->GetString("score", "SCORE_KSDP");
    ScoreType eScore = SCORE_KSDP;
    if ("SCORE_KSDP" == strScore)
    {
    	eScore = SCORE_KSDP;
    }
    else
    {
    	eScore = SCORE_KSDP;
    }
    m_Condition.m_eScoreMethod = eScore;

    //evaluate
    string strEvaluat = pOption->GetString("evaluat", "EV_DEFAULT");

    EvaluateType eEvaluat = EV_DEFAULT;

    m_Condition.m_eEV = eEvaluat;
    
    
    //gen pep
   string strPepGen = pOption->GetString("pepgen", "PEPGENE_DEFAULT");
   PepGeneType ePepGen = PEPGENE_DEFAULT;

   m_Condition.m_ePepGen = ePepGen;
 
    //add by chihao: different protein inference types
    //modify by czhou
    string strInfer = pOption->GetString("infer_type", "PROTEIN_INFER_DEFAULT");
    
     //max_ev
    string strMaxEV = pOption->GetString("max_ev", "1.0");
    double lfMaxEV = 0;
    if("infinity" == strMaxEV)
        lfMaxEV = 1.7976931348623158e+308;
    else
        sscanf(strMaxEV.c_str(),"%lf",&lfMaxEV);
    if(lfMaxEV<0)
        lfMaxEV=1.0f;

    m_Condition.m_lfMaxEV = lfMaxEV;

    //false_positive_rate
    int nFPR =pOption->GetInteger("false_positive_rate", 5);
    double lfFPR = (double) nFPR / 100;
    if( 1 < lfFPR)
        lfFPR = 0.05f;
    else if(0.01f > lfFPR)
        lfFPR = 0.01f;
    m_Condition.m_lfFPRThreshold = lfFPR;

    //false_positive_sign
    string strFPSign = pOption->GetString("false_positive_sign", "REVERSE_");
    m_Condition.m_strFPSign = strFPSign;

    int nReportPepNum = pOption->GetInteger("report_peptide_number",10);
    m_Condition.m_nReportPep = nReportPepNum;
    //report_protein
    int nReportProNum = pOption->GetInteger("report_protein_number",100);
    m_Condition.m_nReportPro =nReportProNum; 

    string strReportPro  = pOption->GetString("report_protein_type", "PROTEIN_REPORT_HTML");
   if("PROTEIN_INFER_AC" == strReportPro)
    {
    	m_Condition.m_eProReport = PROTEIN_REPORT_TXT;
    }
   else if("PROTEIN_REPORT_CANDIDATE" == strReportPro)
    {
    	m_Condition.m_eProReport = PROTEIN_REPORT_CANDIDATE;
    }
    
   	m_Condition.m_nKSDP_L = pOption->GetInteger("KSDP_L", 5);
   	
   	string strTmp = pOption->GetString("KSDP_G", "0.9");
    sscanf(strTmp.c_str(),"%lf",&m_Condition.m_lfKSDP_G);

    strTmp = pOption->GetString("KSDP_A", "0.5");
    sscanf(strTmp.c_str(),"%lf",&m_Condition.m_lfKSDP_A);
    
    //20100309
    string strLog  = pOption->GetString("log_rank", "LOG_RANK_INFO");
    if("LOG_RANK_DEBUG" == strLog)
    {
    	m_Condition.m_eLogRank = LOG_RANK_DEBUG;
    }
    else
    if("LOG_RANK_ALERT" == strLog)
    {
    	m_Condition.m_eLogRank = LOG_RANK_ALERT;
    }
    else
    {
        if("LOG_RANK_NONE" == strLog)
        {
        	m_Condition.m_eLogRank = LOG_RANK_NONE;
        }
    }
    
    strLog  = pOption->GetString("log_appender", "LOG_APPENDER_CMD");
    if("LOG_RANK_FILE" == strLog)
    {
    	m_Condition.m_eLogAppender = LOG_APPENDER_FILE;
    }
    else
    if("LOG_RANK_CMD_AND_FILE" == strLog)
    {
    	m_Condition.m_eLogAppender = LOG_APPENDER_CMD_AND_FILE;
    }

    strLog  = pOption->GetString("log_module", "MODULE_ALL");    
    const size_t ModuleTypeSize = 17;
    const string s_strModuleType [ModuleTypeSize]= {"MODULE_ALL",
    		"MODULE_SEARCHER",
    		"MODULE_SEARCHENGINE",
    		"MODULE_CONDITION",
    		"MODULE_SPECTRAIO",
    		"MODULE_FLOW",
    		"MODULE_PREPRO",
    		"MODULE_SCORE",
    		"MODULE_MASS2PEP",
    		"MODULE_PROIDX",
    		"MODULE_GENERATOR",
    		"MODULE_EVALUATE",
    		"MODULE_PROINFER",
    		"MODULE_IMPORTER",
    		"MODULE_SPECINDEX",
    		"MODULE_LOADBALANCE",
    		"MODULE_MPISEARCHER"
    };
    strLog  = pOption->GetString("log_module", "MODULE_ALL");
    for(size_t i = 0;i < ModuleTypeSize;++i)
    {
    	if(s_strModuleType[i] == strLog)
    	{
    		m_Condition.m_eModule = (ModuleType)i;
    	}
    }

    //Output path. 20140525
    m_Condition.m_strOutputPath = pOption->GetString("output_path", "NULL");
    if (m_Condition.m_strOutputPath[m_Condition.m_strOutputPath.length() - 1] != cSlash)
    {
    	m_Condition.m_strOutputPath += cSlash;
	}
    CCheck::CheckPath(m_Condition.m_strOutputPath);

    delete pOption;

}

void CConditionReader::_LoadSpectra(void)
{
    COptionTool * pOption = new COptionTool("spectrum", m_strOption.c_str());
//#ifndef WIN32 //20140525 Always read
    m_Condition.m_strSpectraTitle = pOption->GetString("spec_title", "null");
    if("null" == m_Condition.m_strSpectraTitle)
     {
        delete pOption;
        CErrInfo info("CConditionReader", "_LoadSpectra", "m_strSpectraTitle is null");
        throw runtime_error(info.Get().c_str());
     }
//#endif

    //The file path legally has been checked eralier at SearchEngine.cpp, so no longer checked here.
    m_Condition.m_strSpectraPath =  pOption->GetString("spec_path", "null");
    string strType = pOption->GetString("spec_type", "MGF");
    if ("DTA" == strType)
    {
		if (m_Condition.m_strSpectraPath[m_Condition.m_strSpectraPath.length() - 1] != cSlash)
			m_Condition.m_strSpectraPath += cSlash;
	}
    CCheck::CheckPath(m_Condition.m_strSpectraPath);

    delete pOption;
    pOption = NULL;
}

void CConditionReader::_LoadCluster(void)
{
#ifndef WIN32
    COptionTool * pOption = new COptionTool("cluster", m_strOption.c_str());

    m_Condition.m_strSpectraIndexPath = pOption->GetString("block_path", "null"); 
    if("null" == m_Condition.m_strSpectraIndexPath)
     {
        delete pOption;
        CErrInfo info("CConditionReader", "_LoadCluster", "block_path is null");
        throw runtime_error(info.Get().c_str());
     }
    
    //20100309
    CTrace * pTrace = CTrace::GetInstance();
    pTrace->Info("The directory to store spectra index files :" + m_Condition.m_strSpectraIndexPath);

    string strType("SPECTRA_INDEX_BIN_BALANCE"); 
    strType = pOption->GetString("block_format_type", "SPECTRA_INDEX_SIMPLE"); 
    if("SPECTRA_INDEX_LARGE" == strType)
    {
		pTrace->Debug("ConditionReader::SPECTRA_INDEX_LARGE");
      	m_Condition.m_eSpectraIndex = SPECTRA_INDEX_LARGE;
    }    
    else
    {
		pTrace->Debug("ConditionReader::SPECTRA_INDEX_SIMPLE");
     	m_Condition.m_eSpectraIndex = SPECTRA_INDEX_SIMPLE;
    }

    string strLoadBalanceType("LOAD_BALANCE_MASS_DYNAMIC_MUL3"); 
    strLoadBalanceType = pOption->GetString("load_balance_type", "LOAD_BALANCE_MASS_DYNAMIC_MUL3"); 
    if("LOAD_BALANCE_MASS_RANGE" == strLoadBalanceType)
     {

		 pTrace->Debug("LOAD_BALANCE_MASS_RANGE");
        
     	m_Condition.m_eLoadBalance = LOAD_BALANCE_MASS_RANGE;
     }
    else if("LOAD_BALANCE_MASS_SHUFFLE" == strLoadBalanceType)
     {		
		pTrace->Debug("LOAD_BALANCE_MASS_SHUFFLE");

		m_Condition.m_eLoadBalance = LOAD_BALANCE_MASS_SHUFFLE;
     }
    else if("LOAD_BALANCE_MASS_DYNAMIC_MUL2" == strLoadBalanceType)
     {		
		pTrace->Debug("LOAD_BALANCE_MASS_DYNAMIC_MUL2");
		m_Condition.m_eLoadBalance = LOAD_BALANCE_MASS_DYNAMIC_MUL2;
     }
	 else if("LOAD_BALANCE_MASS_DYNAMIC_MUL3" == strLoadBalanceType)
     {	
		pTrace->Debug("LOAD_BALANCE_MASS_DYNAMIC_MUL3");
		m_Condition.m_eLoadBalance = LOAD_BALANCE_MASS_DYNAMIC_MUL3;
     }
	 else if("LOAD_BALANCE_MASS_DYNAMIC_MUL4" == strLoadBalanceType)
     {
		pTrace->Debug("LOAD_BALANCE_MASS_DYNAMIC_MUL4");
		m_Condition.m_eLoadBalance = LOAD_BALANCE_MASS_DYNAMIC_MUL4;
     }
    else if("LOAD_BALANCE_ESTIMATE1" == strLoadBalanceType)
     {		
		pTrace->Debug("LOAD_BALANCE_ESTIMATE1");
		m_Condition.m_eLoadBalance = LOAD_BALANCE_ESTIMATE1;
     }
    else if("LOAD_BALANCE_ESTIMATE2" == strLoadBalanceType)
     {
		pTrace->Debug("LOAD_BALANCE_ESTIMATE2");
		m_Condition.m_eLoadBalance = LOAD_BALANCE_ESTIMATE2;
     }

    m_Condition.m_nMinBlockNum = pOption->GetInteger("processor_num",2);//todo m_nMinBlockNum
    m_Condition.m_nMaxBlockSize = pOption->GetInteger("block_max_num",5000);

    delete pOption;
    pOption = NULL;
    
#endif
}

void CConditionReader::_LoadMetaDat()
{
	try
	{
		CDBConf db(m_Condition.m_strDBConfPath.c_str());
		string strMeta = db.GetMetaName(m_Condition.m_vSelectedDBName[0].c_str(), 
				m_Condition.m_SelectedEnzyme.m_strName.c_str(),
				m_Condition.m_bPepMono);
		string strDBPath = db.GetPath(m_Condition.m_vSelectedDBName[0].c_str(),
				m_Condition.m_SelectedEnzyme.m_strName.c_str());
		
		string strPath = strDBPath + strMeta;
		
		FILE *file = fopen(strPath.c_str(), "rb");
		if (file)
		{
			fseek(file, 0, SEEK_SET);
			size_t tOffset;
			fread((char*)&tOffset, sizeof(size_t), 1, file);
			fread((char*)&m_Condition.m_nCleaveWay, sizeof(size_t), 1, file);
			fseek(file, tOffset - sizeof(bool), SEEK_SET);
			fread((char*)&m_Condition.m_bPep2Pro, sizeof(bool), 1, file);	
			fclose(file);
		}
	}
	catch(exception & e)
	{
        CErrInfo info("CConditionReader", "_LoadMetaDat");
        throw runtime_error(info.Get(e).c_str());
    }
    catch(...)
    {
        CErrInfo info("CConditionReader", "_LoadMetaDat", "caught an unknown exception.");
        throw runtime_error(info.Get().c_str());
    }
}


}

