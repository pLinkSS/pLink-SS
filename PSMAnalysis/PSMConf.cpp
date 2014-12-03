#include <iostream>
//#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "PSMConf.h"

using namespace std;
using namespace proteomics_sdk;

#ifdef WIN32
const char cSlash = '\\';
#else
const char cSlash = '/';
#endif

CPSMConf::CPSMConf()
{
	char szBuf[BUFFER_SIZE] = {0};
	getcwd(szBuf, BUFFER_SIZE);
    m_strWorkDir = szBuf;
    if(m_strWorkDir[m_strWorkDir.length() - 1] != cSlash)
    {
    	m_strWorkDir += cSlash;
    }
    m_strFragmentTolType = "";
}

CPSMConf::~CPSMConf()
{
}

void CPSMConf::Load(string strConfigFile)
{
	m_strOption = strConfigFile;
	
	_LoadPFD();
	_LoadIon();
	_LoadGap();
	_LoadFilter();
}

void CPSMConf::Output()
{
	cout << "pfind : " << m_strpFindFile << endl
	<< "report file : " << m_strReportFile << endl
	<< "outpath : " << m_strOutputPath << endl	
	<< "fragment tol : " << m_lfFragmentTol << endl
	<< "mass scope : " << m_lfMassScope << endl
	<< "spectra list : " << m_strSpectraList << endl
	<< "spectra list tag :" << m_strSpectraListPrefix << endl
	<< "remove noise : " << m_bNoiseRemove << endl;
	
	for(int i=0;i<m_vIonTypes.size();++i)
	{
		cout << "ion type " << i << " n-term : " <<  m_vIonTypes[i].cType 
			<< "	charge : " <<  m_vIonTypes[i].nCharge 
			<< "	lost total : " << m_vIonTypes[i].nTotalLostVal << endl;
	}
	
}


void CPSMConf::_LoadPFD()
{
	COptionTool * pOption(NULL);
    pOption = new COptionTool("pfd", m_strOption.c_str());

    string strIdxContent = pOption->GetString("index_content", "PEPTIDE_INDEX");
    
	m_strpFindFile = pOption->GetString("pFindFile", "test.pfind");
	m_strReportFile = pOption->GetString("ReportFile", "test_qry.proteins.txt");
	m_strInputType = pOption->GetString("InputType", "standard");
	m_strOutputPath = pOption->GetString("OutputPath", ".");
	
	CConditionReader reader(m_strpFindFile, m_strWorkDir);
    try
    {
    	reader.Read();
    }
    catch(exception & e)
    {
    	CErrInfo info("CPSMConf", "_LoadPFD", 
    			"in the function CConditionReader::Read");
    	info.Append("m_strParamFile=" + m_strpFindFile);
    	info.Append("m_strWorkDir=" + m_strWorkDir);
    	throw runtime_error(info.Get(e).c_str());
    }
    catch(...)
    {
    	CErrInfo info("CPSMConf", "_LoadPFD", 
    			"caught an unknown exception in the function CConditionReader::Read");
    	info.Append("m_strParamFile=" + m_strpFindFile);
    	info.Append("m_strWorkDir=" + m_strWorkDir);
    	throw runtime_error(info.Get().c_str());
    }
    
	m_Condition = reader.m_Condition;
    delete pOption;
}

void CPSMConf::_LoadIon()
{
	COptionTool * pOption(NULL);
	pOption = new COptionTool("ion", m_strOption.c_str());

	string strTmp = pOption->GetString("tolerance", "0.5");
	m_lfFragmentTol = atof(strTmp.c_str());
	
	m_strFragmentTolType = pOption->GetString("tolerance_type", "Da");
	
	m_vIonTypes.clear();
	int nIonTypes = pOption->GetInteger("ion_type_total", 0);
	if(nIonTypes > 0)
	{
		char cType = 'b';
		int nCharge = 1;
		int nLossNH3 = 0, nLossH2O = 0;
		double lfOther = 0.0; 
		string strName;
		CIonTypeEx ion_type;
		m_vIonTypes.clear();
		for(int i = 0;i < nIonTypes;++i)
		{
			char temp[200] = {0};
			sprintf(temp, "ion_type%d", i + 1);
			string strInfor = pOption->GetString(temp, "null");
			if("null" == strInfor)
			{
	            delete pOption;
	            CErrInfo info("CPSMConf", "_LoadIon", "cannot find the iontype " + string(temp));
	            throw runtime_error(info.Get().c_str());
			}
			double lfElement[ELEMENT_NUMBER] = {0};
			
			lfElement[0] = IonMass_Mono_C;
			lfElement[1] = IonMass_Mono_H;
			lfElement[2] = IonMass_Mono_N;
			lfElement[3] = IonMass_Mono_O;
			
			lfOther = 0;
			istringstream iss(strInfor);
			iss >> cType >> nCharge >> nLossNH3 >> nLossH2O >> lfOther ;
			ion_type.lfOtherLoss = lfOther;
			lfOther += NH3Mass_Mono * nLossNH3;
			lfOther += H2OMass_Mono * nLossH2O;
			ion_type.nLossH2O = nLossH2O;
			ion_type.nLossNH3 = nLossNH3;
			switch(cType)
			{
			case 'a':
				ion_type.cType = 0;
				lfOther += lfElement[0] + lfElement[3];
				ion_type.lfOtherLoss += lfElement[0] + lfElement[3];
				break;
			case 'b':
				ion_type.cType = 0;
				break;
			case 'c':
				ion_type.cType = 0;
				lfOther -= lfElement[2] + 3 * lfElement[1];
				ion_type.lfOtherLoss -= lfElement[2] + 3 * lfElement[1];
				break;
			case 'x':
				ion_type.cType = 1;
				lfOther -= lfElement[0] + lfElement[3] - 2*lfElement[1];
				ion_type.lfOtherLoss -= lfElement[0] + lfElement[3] - 2*lfElement[1];
				break;
			case 'y':
				ion_type.cType = 1;
				break;
			case 'z':
				ion_type.cType = 1;
				lfOther += lfElement[2] + 2 * lfElement[1];
				ion_type.lfOtherLoss += lfElement[2] + 2 * lfElement[1];
				break;
			case 'm':
				ion_type.cType = 2;
				break;
			case 'p':
				ion_type.cType = 3;
				break;
			case 'q':
				ion_type.cType = 4;
				break;
			case 'u':
				ion_type.cType = 5;
				break;
			case 'v':
				ion_type.cType = 6;
				break;
			default:
	            delete pOption;
	            CErrInfo info("CPSMConf", "_LoadIon", string("cannot find the iontype ") + cType);
	            throw runtime_error(info.Get().c_str());
			}
			ion_type.nCharge = nCharge;
			ion_type.nTotalLostVal = (int)(lfOther * MZMULTIPLIER + 0.5);
			m_vIonTypes.push_back(ion_type);
			
			cout << i+1 << "	" 
				<< int(ion_type.cType) << "	"
				<< ion_type.nCharge << "	"
				<< (ion_type.nTotalLostVal + 0.0)/MZMULTIPLIER << "	" << endl;
		}
	}

	
    delete pOption;
	
}

void CPSMConf::_LoadGap()
{
	COptionTool * pOption(NULL);
	pOption = new COptionTool("gap", m_strOption.c_str());

	string strTmp = pOption->GetString("mass_scope", "0.5");
	m_lfMassScope = atof(strTmp.c_str());
	delete pOption;
}

void CPSMConf::_LoadFilter()
{
	COptionTool * pOption(NULL);
	pOption = new COptionTool("filter", m_strOption.c_str());

	m_strSpectraList = pOption->GetString("spectra_list", "");
	m_strSpectraListPrefix = pOption->GetString("spectra_list_prefix", "");
	
	m_bNoiseRemove = pOption->GetInteger("RemoveNoise", 0);
	
	delete pOption;
}

