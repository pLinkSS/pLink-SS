#ifndef PEPTIDEREADERFACTORY_H_
#define PEPTIDEREADERFACTORY_H_


#include "../include/sdk.h"

#include "DiskMass2PepIndexReader.h"
//#include "MMapMass2PepIndexReader.h"
#include "Mass2PepIndexReader.h"

namespace Mass2PepIndex
{


enum PeptideFileType
{
	Peptide_FileType_Default = 0,
	Peptide_FileType_Index_Disk = 1,
	Peptide_FileType_Index_Map = 2,
	Peptide_FileType_Txt_Disk = 3,
	Peptide_FileType_Txt_Map = 4
};

class CPeptideReaderFactory
{
public:
	CPeptideReaderFactory(void){}
	virtual~CPeptideReaderFactory(void){}
	
	CMass2PepIndexReader * GetPeptideReader(PeptideFileType eType)const;
};

class CPeptideRandomReaderFactory
{
public:
	CPeptideRandomReaderFactory(void){};
	virtual~CPeptideRandomReaderFactory(void){};
	
	CMass2PepIndexReader * GetPeptideReader(PeptideFileType eType)const;
};
}

#endif /*PEPTIDEREADERFACTORY_H_*/
