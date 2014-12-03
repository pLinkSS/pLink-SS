#include "ProteinWriterFactory.h"
#include "DiskProteinWriter.h"
//#include "FiveGramProteinWriter.h" //***************
namespace ProteinIndex
{

CProteinWriterFactory::CProteinWriterFactory()
{
}

CProteinWriterFactory::~CProteinWriterFactory()
{
}
 
CProteinWriter * CProteinWriterFactory::GetWriter(ProteinWriterType eType) const
{
	
	CProteinWriter * pWriter = 0;
	if(Protein_Disk == eType)
	{
		pWriter = new CDiskProteinWriter();
	}
	else if(Protein_DB == eType)
	{
		
		
	}
	else if(Protein_Five_Gram == eType)
	{
//		pWriter = new CFiveGramProteinWriter();//***********
	}
	else
	{
		
	}
	return pWriter;
}

}
