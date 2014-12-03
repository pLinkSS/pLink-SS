#ifndef PROTEINWRITERFACTORY_H_
#define PROTEINWRITERFACTORY_H_
#include "ProteinWriter.h"
namespace ProteinIndex
{

enum ProteinWriterType
{
	Protein_Default = 0,
	Protein_Disk = 1,
	Protein_DB = 2,
	Protein_Five_Gram = 3
};

class CProteinWriterFactory
{
public:
	CProteinWriterFactory();
	CProteinWriter * GetWriter(ProteinWriterType eType) const;
	virtual ~CProteinWriterFactory();
};

}

#endif /*PROTEINWRITERFACTORY_H_*/
