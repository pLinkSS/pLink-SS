#include <string>
#include <stdexcept>
#include "../include/sdk.h"
#include "../Mass2PepIndex/Mass2PepIndex.h"
#include "ACEConfigTool.h"
#include "Configure.h"

using namespace std;

CConfigure::CConfigure()
{
}

CConfigure::~CConfigure()
{
}

void CConfigure::Init(string strPIndexName)
{
	this->m_strPIndexName = strPIndexName;
}

void CConfigure::ReadHead()
{
	try
	{
		CACEConfigTool config;
		config.Open(m_strPIndexName);
		
		const string strSection("total");
		m_PIndexHead.nIndexNum = config.GetInteger(strSection, string("db_number"), 1);
		m_PIndexHead.stEnzymeList = config.GetString(strSection, string("enzyme_list"), ".\\enzyme.ini");
		m_PIndexHead.strAAList = config.GetString(strSection, string("aa_list"), ".\\aa.ini");
		m_PIndexHead.strOutpath = config.GetString(strSection, string("out_path"));
		
		config.Close();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("Configure", "ReadHead");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("Configure", "ReadHead","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CConfigure::ReadItems()
{
	try
	{
		CACEConfigTool config;
		config.Open(m_strPIndexName);

		char szTemp[10];
		for(size_t i=1; i<=m_PIndexHead.nIndexNum; ++i)
		{
			PINDEX_ITEM pIndexItem;
			sprintf(szTemp, "db%d", i);
			string strSection(szTemp);
			pIndexItem.bAvrg = config.GetBool(strSection, string("avrg"), true);
			pIndexItem.bMono = config.GetBool(strSection, string("mono"), true);
			pIndexItem.bPep = config.GetBool(strSection, string("pep"), true);
			pIndexItem.bPro = config.GetBool(strSection, string("pro"), true);
			pIndexItem.bMass2Pep = config.GetBool(strSection, string("mass2pep"), false);
			pIndexItem.bPep2Pro = config.GetBool(strSection, string("pep2pro"), false);
			pIndexItem.nEnzymeNum = config.GetInteger(strSection, string("enzyme_number"), 1);
			pIndexItem.nMaxMissSite = config.GetInteger(strSection, string("max_miss_site"), 2);
			pIndexItem.nCleaveWay = config.GetInteger(strSection, string("cleave_way"), 0);
			pIndexItem.nAutoReverse = config.GetInteger(strSection, string("auto_reverse"), 0);
			pIndexItem.strDBName = config.GetString(strSection, string("db_name"), "NULL");
			pIndexItem.strFastaFile = config.GetString(strSection, string("fasta_file"), "NULL");
			pIndexItem.tMassRange = config.GetSizeT(strSection, string("mass_range"), 1000);
			pIndexItem.tMaxMemSize = config.GetSizeT(strSection, string("max_mem_size"), 236870912);
			pIndexItem.tMaxPepLength = config.GetSizeT(strSection, string("max_pep_length"), 100);
			pIndexItem.tMaxPepMass = config.GetSizeT(strSection, string("max_mass"), 6000);
			pIndexItem.tMinPepLength = config.GetSizeT(strSection, string("min_pep_length"), 3);
			pIndexItem.tMinPepMass = config.GetSizeT(strSection, string("min_mass"), 250);
			pIndexItem.tOutTxt = config.GetSizeT(strSection, string("out_txt"), 0);
			pIndexItem.tDatFlow = config.GetSizeT(strSection, string("dat_flow"), 1);
			//cout << "tDatFlow: " << pIndexItem.tDatFlow << endl;
			pIndexItem.tDatsFlow = config.GetSizeT(strSection, string("dats_flow"), 0);
			pIndexItem.bProIndexMap = config.GetSizeT(strSection, string("pro_index_type"), 1);
			pIndexItem.bSaveDebugFile = config.GetBool(strSection, string("save_debug_file"), false);
			pIndexItem.bConsole = config.GetBool(strSection, string("console"), true);
			pIndexItem.tMaxFileSize = config.GetSizeT(strSection, string("pro_file_size"), 500741824);
			//put vector
			char szDigit[10];
			for(int j=1; j<=pIndexItem.nEnzymeNum; ++j)
			{
				sprintf(szDigit, "enzyme%d", j);
				const string strEnzyme(szDigit);
				
				pIndexItem.vtEnzymeNames.push_back(config.GetString(strSection, strEnzyme, "NULL"));
			}
			
			m_PIndexItems.push_back(pIndexItem);
		}
		
		config.Close();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("Configure", "ReadItems");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("Configure", "ReadItems","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}
