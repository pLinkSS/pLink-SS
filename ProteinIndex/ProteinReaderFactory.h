#ifndef PROTEINREADERFACTORY_H_
#define PROTEINREADERFACTORY_H_


namespace ProteinIndex
{

enum ProteinDataBaseType
{
	
	Protein_DB_Index_Disk = 0,
	Protein_DB_Index_Map = 1,
	Protein_DB_Default = 2,
	Protein_DB_Index_Memory = 3,
	Protein_DB_Fasta = 4,
	//Protein_DB_FiveGram_MMap = 5  
//	Protein_XML = 5
};

class CProteinReader;


class CProteinReaderFactory{
public:
	CProteinReaderFactory();
	~CProteinReaderFactory();
	CProteinReader * GetACReader(ProteinDataBaseType eType) const;
	CProteinReader * GetDEReader(ProteinDataBaseType eType) const;
	CProteinReader * GetSQReader(ProteinDataBaseType eType) const;
};

enum ProteinRandomReaderType
{
	Protein_Random_Index_Disk = 0,
	Protein_Random_Index_Map = 1,
	Protein_Random_Default = 2,
	Protein_Random_Index_Memory = 3,
	Protein_Random_FiveGram_MMap = 5
//	Protein_Random_TXT =4,
//	Protein_Random_RDB =5
};

class CRandomReader;

class CProteinRandomReaderFactory
{
public:
	CRandomReader * GetRandomACReader(ProteinRandomReaderType eType) const;
//	CRandomReader * GetRandomDEReader(ProteinLoadType eType) const;
	CRandomReader * GetRandomSQReader(ProteinRandomReaderType eType) const;
};

}

#endif /*PROTEINREADERFACTORY_H_*/
