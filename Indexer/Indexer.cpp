#include <string>
#include "../include/sdk.h"
#include "../Mass2PepIndex/SPSingleScanIndexCreator.h" 
#include "../Mass2PepIndex/PeptideReaderFactory.h"
#include "../Mass2PepIndex/MetaCreater.h"
#include "../Mass2PepIndex/DiskMass2PepIndexReader.h"

#include "../ProteinIndex/FastaParser.h"
#include "../ProteinIndex/FastaProteinDB.h"

#include "DBConf.h"
#include "Indexer.h"
using namespace ProteinIndex;
using namespace std;
using namespace proteomics_sdk;


CExcuteIndexer::CExcuteIndexer()
{
}
 
CExcuteIndexer::~CExcuteIndexer()
{
}

void _ShowFEndInfo2File(string strFName, double dfTime)
{
//	try
//	{
//		ofstream os("time.txt",ios::app);
//	
//		os << "Time spent on " << strFName << " is:" << (clock() - dfTime)/1000 << " seconds." << endl << endl;
//	}
//	catch(...)
//	{
//		cout << "write time file is error" << endl;
//	}
}

void CExcuteIndexer::Init( ProteinIndex::ProteinWriterType eProteinWriterType, IndexCreatorType ePeptideCreatorType, std::string &strPIndexName)
{	
	if(strPIndexName.empty())
	{
		proteomics_sdk::CErrInfo err_info("CExcuteIndexer", "Init","Can't find config file!");
		throw runtime_error(err_info.Get().c_str());
	}	
	
	CProteinWriterFactory proteinFactory;
	m_proteinIndexerCreator = proteinFactory.GetWriter(eProteinWriterType);
	
	CMass2PepIndexCreatorFactory mass2PepFactory;
	m_peptideIndexerCreator = mass2PepFactory.GetCreator(ePeptideCreatorType);//to make sure which kind of factory


	try
	{
		m_configure.Init(strPIndexName);
		m_configure.ReadHead();
		m_configure.ReadItems();		
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CExcuteIndexer", "Init");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CExcuteIndexer", "Init","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}	
}

void CExcuteIndexer::_WriteConf(PINDEX_HEAD &pIndexHead, PINDEX_ITEM& pIndexItem, size_t bt, bool bMono)
{
	string strMOrA = ".mono";
	if(!bMono)strMOrA = ".avrg";
	string strInfo= "Create PepIndex of " + pIndexItem.strDBName +  "." + pIndexItem.vtEnzymeNames[0] + strMOrA ;

	_ShowFEndInfo2File(strInfo, bt);
	
	CDBConfManage manage;
	manage.Load("dbconf.ini");
	
	DBCONF_ITEM dbConfItem;
	dbConfItem.strDBName = pIndexItem.strDBName.c_str();
	dbConfItem.strEnzymeName = pIndexItem.vtEnzymeNames[0].c_str();	
	dbConfItem.strMetaName = dbConfItem.strDBName + "." + dbConfItem.strEnzymeName + strMOrA + ".meta";
	dbConfItem.strPath = pIndexHead.strOutpath.c_str();
	
	time_t current_time;
	time(&current_time);				
	dbConfItem.strRemark = ctime(&current_time);				
	
	if(manage.CheckItemValidity(dbConfItem)>=0)
	{
		cout << "Warning:\t" << dbConfItem.strMetaName <<"\talready exist, overwrite!"<<endl;
		manage.DelItem(dbConfItem);
		manage.AddItem(dbConfItem);
	}
	else
	{
		manage.AddItem(dbConfItem);
	}
	
	manage.Save();
}

void CExcuteIndexer::_WriteMeta(PINDEX_HEAD &pIndexHead, PINDEX_ITEM& pIndexItem)
{
	CMetaCreater	 m_CMetaCreater;
	m_CMetaCreater.Load(pIndexHead.strOutpath.c_str(), pIndexItem.strDBName.c_str(), pIndexItem.vtEnzymeNames[0].c_str());
	m_CMetaCreater.m_stMetaHead.nCleaveWay = (long)pIndexItem.nCleaveWay;
	m_CMetaCreater.m_stMetaHead.tMultiplier = DEFAULT_MULTPLIER;
	m_CMetaCreater.m_stMetaHead.tPepSQNum = 0;
	m_CMetaCreater.m_stMetaHead.tUniquePepSQNum = 0;
	m_CMetaCreater.m_stMetaHead.tUniqueMassNum = 0;
	m_CMetaCreater.m_stMetaHead.tMinMass = pIndexItem.tMinPepMass * m_CMetaCreater.m_stMetaHead.tMultiplier;
	m_CMetaCreater.m_stMetaHead.tMaxMass = pIndexItem.tMaxPepMass * m_CMetaCreater.m_stMetaHead.tMultiplier;
	m_CMetaCreater.m_stMetaHead.tMinLength = pIndexItem.tMinPepLength;
	m_CMetaCreater.m_stMetaHead.tMaxLength = pIndexItem.tMaxPepLength;
	m_CMetaCreater.m_stMetaHead.tIdxNum = 0;//*************
	m_CMetaCreater.m_stMetaHead.tMissCleaveSites = pIndexItem.nMaxMissSite;
	m_CMetaCreater.m_stMetaHead.bMono = pIndexItem.bMono;
	m_CMetaCreater.m_stMetaHead.bPep2Pro = pIndexItem.bPep2Pro;
	m_CMetaCreater.m_vstMetaItems.clear();
	m_CMetaCreater.WriteMeta();
}

void CExcuteIndexer::Run(void)
{
	double bt = 0;
	
	for(size_t i=0; i < m_configure.m_PIndexItems.size(); ++i)//database 
	{
		PINDEX_HEAD& pIndexHead = m_configure.m_PIndexHead;
		PINDEX_ITEM& pIndexItem = m_configure.m_PIndexItems[i];
		
		bt = clock();

		m_proteinWriter.Init(pIndexHead, pIndexItem);//**********
		m_proteinWriter.Write(m_proteinIndexerCreator);//we should have to restruct the protein Index it is need to construting the decoydata?
				
		string strInfo = "Create Protein Index of " + pIndexItem.strDBName;
		_ShowFEndInfo2File(strInfo, bt);
		
		if(!pIndexItem.bPep)
		{
			_WriteConf(pIndexHead, pIndexItem, bt, true);
			_WriteMeta(pIndexHead, pIndexItem);
			cout << endl << "Do not create peptide index of DB: " << pIndexItem .strDBName << " and enzyme name" << pIndexItem .vtEnzymeNames[0]<< endl;
			continue;
		}
		
		vector<string> vtEnzymeNames = pIndexItem.vtEnzymeNames;
		
		for(long j=0; j< pIndexItem.nEnzymeNum; ++j)//enenzy
		{
			pIndexItem.vtEnzymeNames.clear();
			pIndexItem.vtEnzymeNames.push_back(vtEnzymeNames[j]);
			
			if(pIndexItem.bMono)
			{
				size_t bt = clock();
				
				m_peptideIndexerCreator->Init(pIndexHead, pIndexItem, true);
				m_peptideIndexerCreator->WriteIndex();//write mono
				m_peptideIndexerCreator->Close();
				
				_WriteConf(pIndexHead, pIndexItem, bt, true);
			}
			if(pIndexItem.bAvrg)
			{
				size_t bt = clock();
				
				m_peptideIndexerCreator->Init(pIndexHead, pIndexItem, false);
				m_peptideIndexerCreator->WriteIndex();//write arvg
				m_peptideIndexerCreator->Close();	
				
				_WriteConf(pIndexHead, pIndexItem, bt, false);
			}
			
		}	
	
		pIndexItem.vtEnzymeNames = vtEnzymeNames;
	}
}

void CExcuteIndexer::Tester()
{
	for(size_t i=0; i < m_configure.m_PIndexItems.size(); ++i)//database 
	{				
		vector<string> vtEnzymeNames = m_configure.m_PIndexItems[i].vtEnzymeNames;
		bool	bMono = m_configure.m_PIndexItems[i].bMono;
		bool	bAvrg = m_configure.m_PIndexItems[i].bAvrg;
		
		for(long j=0; j< m_configure.m_PIndexItems[i].nEnzymeNum; ++j)//enenzy
		{
			m_configure.m_PIndexItems[i].vtEnzymeNames.clear();
			m_configure.m_PIndexItems[i].vtEnzymeNames.push_back(vtEnzymeNames[j]);
			
			for(size_t k = 0; k < 2; ++k)//mono or avg
			{
//				size_t bt = clock();
				
				if(!k)//when o == k, we do bMono
				{
					if(!m_configure.m_PIndexItems[i].bMono) continue;
					else m_configure.m_PIndexItems[i].bAvrg = false;// bMono is true, make true bAvrg is false;
				}else
				{
					if(!m_configure.m_PIndexItems[i].bAvrg) continue;
					else m_configure.m_PIndexItems[i].bMono = false;//similar to above
				}			
				
				string strMOrA = ".mono";
				if(!m_configure.m_PIndexItems[i].bMono)strMOrA = ".avrg";
				string strTmp = m_configure.m_PIndexItems[i].strDBName + "." + m_configure.m_PIndexItems[i].vtEnzymeNames[0] + strMOrA;
				string  strMetaFile = strTmp + PEP_META;
				
				ProteinRandomReaderType eType=(ProteinRandomReaderType)(m_configure.m_PIndexItems[i].bProIndexMap);
				CDiskMass2PepIndexReader reader;
				reader.Open(m_configure.m_PIndexHead.strOutpath, strMetaFile, eType);
				
				cout << m_configure.m_PIndexHead.strOutpath << " " << strMetaFile << endl;

				ofstream os("Peptide Number.txt", ios::app);
				os << "all Peptide of " << m_configure.m_PIndexItems[i].strDBName << "is:" << reader.m_pepHandler.m_metaReader.GetMetaHead().tPepSQNum << endl;
				
				size_t tfNum = 0, tsNum = 0;
				for(size_t t = 0; t < reader.m_pepHandler.m_metaReader.GetMetaItems().size(); ++t)
				{
					tfNum += reader.m_pepHandler.m_metaReader.GetMetaItems()[t].tPepSQNum;
					tsNum += reader.m_pepHandler.m_metaReader.GetMetaItems()[t].tUniquePepSQNum;
				}
				os << "after the first time of excluding, the number of peptide of " << m_configure.m_PIndexItems[i].strDBName << "is:" << tfNum << endl;
				os << "after the second time of excluding, the number of peptide of " << m_configure.m_PIndexItems[i].strDBName << "is:" << tsNum << endl;
				
				os.close();

				reader.Close();
				
				if(!k) m_configure.m_PIndexItems[i].bAvrg = bAvrg;//change back of the value
				else m_configure.m_PIndexItems[i].bMono = bMono;
			}
		}
		
		m_configure.m_PIndexItems[i].vtEnzymeNames = vtEnzymeNames;
	}	
}

void CExcuteIndexer::Reader()
{
}

void CExcuteIndexer::Timer()
{
	double bt;

	CFastaProteinDB m_proDB[m_configure.m_PIndexItems.size()];	
	for(size_t i=0; i < m_configure.m_PIndexItems.size(); ++i)//database 
	{
//		PINDEX_HEAD& pIndexHead = m_configure.m_PIndexHead;
//		PINDEX_ITEM& m_configure.m_PIndexItems[i] = m_configure.m_PIndexItems[i];
		
		bt = clock();
		
		
		m_proDB[i].SetDBName(m_configure.m_PIndexItems[i].strDBName);
		m_proDB[i].SetDBPath(m_configure.m_PIndexItems[i].strFastaFile);		
		m_proDB[i].OpenFile();	
		while ( !m_proDB[i].IsEof() ) 
		{
			m_proDB[i].ReadOnePrteinEntry();
		}
		m_proDB[i].CloseFile();
		
		string strInfo= "Read All SQ of " + m_configure.m_PIndexItems[i].strDBName + " in fasta directly";
		_ShowFEndInfo2File(strInfo, bt);
		
		vector<string> vtEnzymeNames = m_configure.m_PIndexItems[i].vtEnzymeNames;
		bool	bMono = m_configure.m_PIndexItems[i].bMono;
		bool	bAvrg = m_configure.m_PIndexItems[i].bAvrg;
		
		for(long j=0; j< m_configure.m_PIndexItems[i].nEnzymeNum; ++j)//enenzy
		{
			m_configure.m_PIndexItems[i].vtEnzymeNames.clear();
			m_configure.m_PIndexItems[i].vtEnzymeNames.push_back(vtEnzymeNames[j]);
			
			for(size_t k = 0; k < 2; ++k)//mono or avg
			{
				if(!k)//when o == k, we do bMono
				{
					if(!m_configure.m_PIndexItems[i].bMono) continue;
					else m_configure.m_PIndexItems[i].bAvrg = false;// bMono is true, make true bAvrg is false;
				}else
				{
					if(!m_configure.m_PIndexItems[i].bAvrg) continue;
					else m_configure.m_PIndexItems[i].bMono = false;//similar to above
				}
	
				bt = clock();
				ProteinIndex::ProteinRandomReaderType eType = (ProteinRandomReaderType)(m_configure.m_PIndexItems[i].bProIndexMap);
				
				ProteinIndex::CProteinRandomReaderFactory ProRandomFactory;		
				CRandomReader *m_proteinExecuter = ProRandomFactory.GetRandomSQReader(eType);
				
				m_proteinExecuter->Open(m_configure.m_PIndexHead.strOutpath,m_configure.m_PIndexItems[i].strDBName);
			
				size_t tProNum = m_proteinExecuter->GetProNum();
				string strProSQ = "";
				
				for(size_t tProID = 0; tProID < tProNum; ++tProID)
				{
					m_proteinExecuter->GetByID( strProSQ, tProID );
				}
				
				string strMOrA = ".mono";
				if(!m_configure.m_PIndexItems[i].bMono)strMOrA = ".avrg";
				strInfo= "Read All SQ of " + m_configure.m_PIndexItems[i].strDBName +  "." + m_configure.m_PIndexItems[i].vtEnzymeNames[0] + strMOrA ;
				if(!m_configure.m_PIndexItems[i].bProIndexMap) strInfo += " By Disk";
				else strInfo += " By MMap";
				_ShowFEndInfo2File(strInfo, bt);			
		
				if(!k) m_configure.m_PIndexItems[i].bAvrg = bAvrg;//change back of the value
				else m_configure.m_PIndexItems[i].bMono = bMono;
			}
		}

	}

}


