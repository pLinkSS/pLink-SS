#include "MetaCreater.h"

using namespace std;
using namespace proteomics_sdk;

namespace Mass2PepIndex
{

CMetaCreater::CMetaCreater()
{
}

CMetaCreater::~CMetaCreater()
{
}

void	CMetaCreater::Load( const std::string strWorkPath, const std::string strDBName, const std::string strEnzymeName)
{
	if(strWorkPath.empty() || strDBName.empty() || strEnzymeName.empty())
	{
		CErrInfo info("CMetaCreater", "Load", "Pointers should not be none!");
		
		info.Append("strWorkPath = " + strWorkPath);
		info.Append("strDBName = " + strDBName);
		info.Append("strEnzymeName = " + (string)(strEnzymeName));
		
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

	try
	{	
		m_strWorkDir = strWorkPath;
		strcpy(m_stMetaHead.cDBName, strDBName.c_str());
		strcpy(m_stMetaHead.cEnzyme, strEnzymeName.c_str());
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaCreater", "Load", "strcpy");

		info.Append("strDBName = " + strDBName);
		info.Append("strEnzymeName = " + (string)(strEnzymeName));
		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaCreater", "Load", "caught an unkown exception");
		throw runtime_error(info.Get().c_str());
	}
}
	
void	CMetaCreater::SetMetaHead(META_HEAD stMetaHead)
{
	m_stMetaHead = stMetaHead;
}

void	CMetaCreater::SetMetaItems(vector<META_ITEM > vstMetaItems)
{
	m_vstMetaItems = vstMetaItems;
}

void	CMetaCreater::WriteMeta()
{
	try
	{
		std::string strMOrA = ".mono";
		
		if(!m_stMetaHead.bMono)
		{
			strMOrA = ".avrg";
		}
		
		m_strMetaName = m_strWorkDir + string(m_stMetaHead.cDBName) + "." + string(m_stMetaHead.cEnzyme)+ strMOrA.c_str() + PEP_META;  
		m_strMetaNameTxt = m_strMetaName + ".txt";
		
		m_fMeta = fopen(m_strMetaName.c_str(), "wb");
		m_fMetaTxt = fopen(m_strMetaNameTxt.c_str(), "w");

		_WriteMetaHead(m_fMeta);
		_WriteMetaItems(m_fMeta);
		
		fclose(m_fMeta);
		fclose(m_fMetaTxt);
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaCreater", "WriteMeta");
		info.Append("m_strMetaName = " + m_strMetaName);		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaCreater", "WriteMeta", "caught an unkown exception");
		throw runtime_error(info.Get().c_str());
	}
}

void	CMetaCreater::_WriteMetaHead(FILE * fMeta)
{
	try
	{
		//bin
		fseek(fMeta, 0, SEEK_SET);
		//fwrite((char *)&m_stMetaHead, _GetMetaHeadSize(), 1, fMeta);
		
		//write one by one
		fwrite((char *)&m_stMetaHead.tOffset, sizeof(size_t), 1, fMeta);
		fwrite((char *)&m_stMetaHead.nCleaveWay, sizeof(size_t), 1, fMeta);
		fwrite((char *)&m_stMetaHead.tMultiplier, sizeof(size_t), 1, fMeta);
		fwrite((char *)&m_stMetaHead.tPepSQNum, sizeof(size_t), 1, fMeta);
		fwrite((char *)&m_stMetaHead.tUniquePepSQNum, sizeof(size_t), 1, fMeta);
		fwrite((char *)&m_stMetaHead.tUniqueMassNum, sizeof(size_t), 1, fMeta);
		fwrite((char *)&m_stMetaHead.tMinMass, sizeof(size_t), 1, fMeta);
		fwrite((char *)&m_stMetaHead.tMaxMass, sizeof(size_t), 1, fMeta);
		fwrite((char *)&m_stMetaHead.tMinLength, sizeof(size_t), 1, fMeta);
		fwrite((char *)&m_stMetaHead.tMaxLength, sizeof(size_t), 1, fMeta);
		fwrite((char *)&m_stMetaHead.cDBName, sizeof(char)*PATH_MAX, 1, fMeta);
		fwrite((char *)&m_stMetaHead.cEnzyme, sizeof(char)*PATH_MAX, 1, fMeta);
		fwrite((char *)&m_stMetaHead.tIdxNum, sizeof(unsigned char), 1, fMeta);
		fwrite((char *)&m_stMetaHead.tMissCleaveSites, sizeof(unsigned char), 1, fMeta);
		fwrite((char *)&m_stMetaHead.bMono, sizeof(bool), 1, fMeta);
		fwrite((char *)&m_stMetaHead.bPep2Pro, sizeof(bool), 1, fMeta);
		
		m_stMetaHead.tOffset = ftell(fMeta);
		
		fseek(fMeta, 0, SEEK_SET);
		fwrite((char *)&m_stMetaHead.tOffset, sizeof(size_t), 1, fMeta);
		
		fseek(fMeta, 0, SEEK_END);
		//txt
		fprintf(m_fMetaTxt, "[Head]\n");
		
		fprintf(m_fMetaTxt, "db_name = %s\n", m_stMetaHead.cDBName);
		fprintf(m_fMetaTxt, "enzyme_name = %s\n", m_stMetaHead.cEnzyme);
		
		fprintf(m_fMetaTxt, "sizeof(head) = %u\n", m_stMetaHead.tOffset);
		fprintf(m_fMetaTxt, "pep_sq_num = %u\n", m_stMetaHead.tPepSQNum);
		fprintf(m_fMetaTxt, "unique_pep_sq_num = %u\n", m_stMetaHead.tUniquePepSQNum);
		
		fprintf(m_fMetaTxt, "min_mass = %u\n", m_stMetaHead.tMinMass/m_stMetaHead.tMultiplier);
		fprintf(m_fMetaTxt, "max_mass = %u\n", m_stMetaHead.tMaxMass/m_stMetaHead.tMultiplier);
		fprintf(m_fMetaTxt, "min_len = %u\n", m_stMetaHead.tMinLength);
		fprintf(m_fMetaTxt, "max_len = %u\n", m_stMetaHead.tMaxLength);
		
		fprintf(m_fMetaTxt, "max_miss_cleave_site = %u\n", m_stMetaHead.tMissCleaveSites);
		fprintf(m_fMetaTxt, "dat_num = %u\n", m_stMetaHead.tIdxNum);
		fprintf(m_fMetaTxt, "multiplier = %u\n", m_stMetaHead.tMultiplier);
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaCreater", "_WriteMetaHead");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaCreater", "_WriteMetaHead", "caught an unkown exception");
		throw runtime_error(info.Get().c_str());
	}
}

void	CMetaCreater::_WriteMetaItems(FILE * fMeta)
{
	try
	{
		//bin
		fseek(fMeta, 0, SEEK_END);
		
		for(size_t t=0; t<m_vstMetaItems.size(); ++t)
		{
			//fwrite((char *)&m_vstMetaItems[t], sizeof(META_ITEM), 1, fMeta);

			fwrite((char *)&m_vstMetaItems[t].tPepSQNum, sizeof(size_t), 1, fMeta);
			fwrite((char *)&m_vstMetaItems[t].tUniquePepSQNum, sizeof(size_t), 1, fMeta);
			fwrite((char *)&m_vstMetaItems[t].tUniqueMassNum, sizeof(size_t), 1, fMeta);
			fwrite((char *)&m_vstMetaItems[t].tMinMass, sizeof(size_t), 1, fMeta);
			fwrite((char *)&m_vstMetaItems[t].tMaxMass, sizeof(size_t), 1, fMeta);
		}
		//txt
		for(size_t t=0; t<m_vstMetaItems.size(); ++t)
		{
			fprintf(m_fMetaTxt, "[Dat%u]\n", t);
			fprintf(m_fMetaTxt, "pep_sq_num = %u\n", m_vstMetaItems[t].tPepSQNum);
			fprintf(m_fMetaTxt, "unique_pep_sq_num = %u\n", m_vstMetaItems[t].tUniquePepSQNum);
			fprintf(m_fMetaTxt, "min_mass = %u\n", m_vstMetaItems[t].tMinMass/m_stMetaHead.tMultiplier);
			fprintf(m_fMetaTxt, "max_mass = %u\n", m_vstMetaItems[t].tMaxMass/m_stMetaHead.tMultiplier);
		}
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaCreater", "_WriteMetaItems");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaCreater", "_WriteMetaItems", "caught an unkown exception");
		throw runtime_error(info.Get().c_str());
	}
}

void	CMetaCreater::WriteMetaItemMass(const int nIndex, const size_t tMassNum)
{
	try
	{
		m_strMetaName = m_strWorkDir + string(m_stMetaHead.cDBName) + "." + string(m_stMetaHead.cEnzyme)+ PEP_META;  
		m_fMeta = fopen(m_strMetaName.c_str(), "r+b");
		
		fseek(m_fMeta, m_stMetaHead.tOffset+nIndex*_GetMetaItemSize()+sizeof(size_t), SEEK_SET);
		printf("index: %u %d", m_stMetaHead.tOffset, nIndex);
		fwrite((char *)&tMassNum, sizeof(size_t), 1, m_fMeta);
		
		fclose(m_fMeta);
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaCreater", "WriteMetaItemMass");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaCreater", "WriteMetaItemMass", "caught an unkown exception");
		throw runtime_error(info.Get().c_str());
	}
}

void    CMetaCreater::WriteTotalMass(const size_t tMassNum)
{
	try
	{
		m_strMetaName = m_strWorkDir + string(m_stMetaHead.cDBName) + "." + string(m_stMetaHead.cEnzyme)+ PEP_META;
		
		m_fMeta = fopen(m_strMetaName.c_str(), "r+b");
		fseek(m_fMeta, 3*sizeof(size_t), SEEK_SET);
		fwrite((char *)&tMassNum, sizeof(size_t), 1, m_fMeta);
		
		fclose(m_fMeta);
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaCreater", "WriteTotalMass");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaCreater", "WriteTotalMass", "caught an unkown exception");
		throw runtime_error(info.Get().c_str());
	}
}
size_t	CMetaCreater::_GetMetaHeadSize(void)const
{
	return sizeof(META_HEAD);
}
size_t	CMetaCreater::_GetMetaItemSize(void)const
{
	return sizeof(META_ITEM);
}

void	CMetaCreater::SetDBName(std::string strDBName)
{
	strcpy(m_stMetaHead.cDBName, strDBName.c_str());
}

void	CMetaCreater::SetEnzymeName(std::string strEnzymeName)
{
	strcpy(m_stMetaHead.cEnzyme, strEnzymeName.c_str());
}
void	CMetaCreater::SetIdxNum(unsigned char tIdxNum)
{
	m_stMetaHead.tIdxNum = tIdxNum;
}

void	CMetaCreater::SetMissCleaves(unsigned char tMissCleaveSites)
{
	m_stMetaHead.tMissCleaveSites = tMissCleaveSites;
}

void	CMetaCreater::SetMultiplier(size_t	tMultiplier)
{
	m_stMetaHead.tMultiplier = tMultiplier;
}

void	CMetaCreater::SetRPepSQNum(size_t	tPepSQNum)
{
	m_stMetaHead.tPepSQNum = tPepSQNum;
}

void	CMetaCreater::SetRPepSQNum(size_t	tPepSQNum, size_t tCurrent)
{
	m_vstMetaItems[tCurrent].tPepSQNum = tPepSQNum;
}
void	CMetaCreater::SetPepSQNum(size_t	tUniquePepSQNum)
{
	m_stMetaHead.tUniquePepSQNum = tUniquePepSQNum;
}

void	CMetaCreater::SetPepSQNum(size_t	tUniquePepSQNum, size_t tCurrent)
{
	m_vstMetaItems[tCurrent].tUniquePepSQNum = tUniquePepSQNum;
}
void	CMetaCreater::SetPepMassNum(size_t	tUniquePepMassNum)
{
	m_stMetaHead.tUniqueMassNum = tUniquePepMassNum;
}
void	CMetaCreater::SetPepMassNum(size_t	tUniquePepMassNum, size_t tCurrent)
{
	m_vstMetaItems[tCurrent].tUniqueMassNum = tUniquePepMassNum;
}
void	CMetaCreater::SetMinPepMass(size_t tMinMass)
{
	m_stMetaHead.tMinMass = tMinMass;
}
void	CMetaCreater::SetMinPepMass(size_t tMinMass, size_t tCurrent)
{
	m_vstMetaItems[tCurrent].tMinMass = tMinMass;
}
void	CMetaCreater::SetMaxPepMass(size_t tMaxMass)
{
	m_stMetaHead.tMaxMass = tMaxMass;
}
void	CMetaCreater::SetMaxPepMass(size_t tMaxMass, size_t tCurrent)
{
	m_vstMetaItems[tCurrent].tMaxMass = tMaxMass;
}

void 	CMetaCreater::SetWorkDir(const char * szDestDir)
{
	try
	{
		string strDestDir( szDestDir );
		string::size_type idx = 0;
		idx = strDestDir.find('\\', idx);
		while (idx != string::npos ) {
			strDestDir[idx] = '/';
			idx = strDestDir.find('\\', idx);
		}
	
		idx = strDestDir.find_last_not_of(" ");
		if (strDestDir[idx] != '/' ) 
		{
			m_strWorkDir = strDestDir.substr(0,idx + 1) + "/";
		}
		else
			m_strWorkDir = strDestDir;
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CMetaCreater", "SetWorkDir");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CMetaCreater", "SetWorkDir", "caught an unkown exception");
		throw runtime_error(info.Get().c_str());
	}
}
}
