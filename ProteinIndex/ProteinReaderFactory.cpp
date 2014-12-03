#include <string>
#include "ProteinReader.h"
#include "ProteinReaderFactory.h"
#include "DiskRandomSQReader.h"
#include "DiskRandomACReader.h"
//#include "MMapRandomSQReader.h"
//#include "MMapRandomACReader.h"
//#include "FiveGramProteinMMapSQReader.h"
#include "FastaSQReader.h"


using namespace std;

namespace ProteinIndex
{

CProteinReaderFactory::CProteinReaderFactory()
{
}

CProteinReaderFactory::~CProteinReaderFactory()
{
}

CProteinReader * CProteinReaderFactory::GetACReader(ProteinDataBaseType eType) const
{
	
	CProteinReader * pReader = 0;
	if(Protein_DB_Index_Disk == eType)
	{
		
	}
	else if(Protein_DB_Index_Memory == eType)
	{
		
		
	}
	else if(Protein_DB_Index_Map == eType)
	{
		
		
	}
	else if(Protein_DB_Fasta == eType)
	{
		
		
	}
	else
	{
//		pReader = new CDiskRandomSQReader();
	}
	return pReader;
}

CProteinReader * CProteinReaderFactory::GetDEReader(ProteinDataBaseType eType) const
{
	
	CProteinReader * pReader = 0;
	if(Protein_DB_Index_Disk == eType)
	{
 
	}
	else if(Protein_DB_Index_Memory == eType)
	{
		
		
	}
	else if(Protein_DB_Index_Map == eType)
	{
		
		
	}
	else if(Protein_DB_Fasta == eType)
	{
		
		
	}
	else
	{
		
		
	}
	return pReader;
}

CProteinReader * CProteinReaderFactory::GetSQReader(ProteinDataBaseType eType) const
{
	CProteinReader * pReader = 0;

	if(Protein_DB_Index_Disk == eType)
	{
//		pReader = new CDiskRandomSQReader();
	}
	else if(Protein_DB_Index_Memory == eType)
	{
		
		
	}
	else if(Protein_DB_Index_Map == eType)
	{
//		pReader = new CMMapRandomSQReader();
		
	}
	else if(Protein_DB_Fasta == eType)
	{
		pReader = new CFastaSQReader();
	}
	else
	{
		
		
	}
	return pReader;
}


//	CRandomReader * CProteinRandomReaderFactory::GetRandomACReader(ProteinLoadType eType) const;
//	CRandomReader * CProteinRandomReaderFactory::GetRandomDEReader(ProteinLoadType eType) const;
CRandomReader * CProteinRandomReaderFactory::GetRandomSQReader(ProteinRandomReaderType eType) const
{

	CRandomReader * pReader = 0;
	return pReader = new CDiskRandomSQReader();
//	if(Protein_Random_FiveGram_MMap == eType)
//	{
////		pReader = new CFiveGramProteinMMapSQReader();
//	}
//	else if(Protein_Random_Index_Disk == eType)
//	{
//		pReader = new CDiskRandomSQReader();
//	}
//	else if(Protein_Random_Index_Memory == eType)
//	{
//		
//		
//	}
//	else if(Protein_Random_Index_Map == eType)
//	{
////		pReader = new CMMapRandomSQReader();
//	}
//	else
//	{
////		pReader = new CMMapRandomSQReader();
//	}
//	
//	return pReader;
}

CRandomReader * CProteinRandomReaderFactory::GetRandomACReader(ProteinRandomReaderType eType) const
{

	CRandomReader * pReader = 0;
	return new CDiskRandomACReader();
//	if(Protein_Random_Index_Disk == eType)
//	{
//		pReader = new CDiskRandomACReader();
//	}
//	else if(Protein_Random_Index_Memory == eType)
//	{
//		
//		
//	}
//	else if(Protein_Random_Index_Map == eType)
//	{
//		pReader = new CMMapRandomACReader();
//	}
//	else
//	{
//		pReader = new CDiskRandomACReader();
//	}
//	
//	return pReader;
}
}
