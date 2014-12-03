
#include "PeptideReaderFactory.h"

using namespace proteomics_sdk;

namespace Mass2PepIndex
{

CMass2PepIndexReader * CPeptideReaderFactory::GetPeptideReader(PeptideFileType eType)const
{
	CMass2PepIndexReader * pReader = NULL;	
	try
	{
		return pReader = new CDiskMass2PepIndexReader();
//
//		
//		if(Peptide_FileType_Index_Disk == eType)
//		{
////			pReader = new CDiskMass2PepIndexReader();
//		}
//		else if(Peptide_FileType_Index_Map == eType)
//		{		
////			pReader = new CMMapPeptideDatReader();
//		}
//		else if(Peptide_FileType_Txt_Disk == eType)
//		{
//			
//			
//		}
//		else if(Peptide_FileType_Txt_Map == eType)
//		{
//			
//			
//		}
//		else
//		{
////			pReader = new CDiskMass2PepIndexReader();
//		}
	}
	catch(runtime_error &e)
	{
		if(pReader)
		{
			delete pReader;
		}
		CErrInfo info("CPeptideReaderFactory", "GetPeptideReader");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		if(pReader)
		{
			delete pReader;
		}
		CErrInfo info("CPeptideReaderFactory", "GetPeptideReader", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
	
	if(!pReader)
	{
		CErrInfo info("CPeptideReaderFactory", "GetPeptideReader", "pReader is null");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
	
	return pReader;
}

CMass2PepIndexReader * CPeptideRandomReaderFactory::GetPeptideReader(PeptideFileType eType)const
{
	CMass2PepIndexReader * pReader = 0;
	try
	{
		return new CDiskMass2PepIndexReader();
//		
//		if(Peptide_FileType_Index_Disk == eType)
//		{
//			pReader = new CDiskMass2PepIndexReader();
//		}
//		else if(Peptide_FileType_Index_Map == eType)
//		{
//			
//			pReader = new CMMapMass2PepIndexReader();
//		}
//		else if(Peptide_FileType_Txt_Disk == eType)
//		{
//			
//			
//		}
//		else if(Peptide_FileType_Txt_Map == eType)
//		{
//			
//			
//		}
//		else
//		{
//			pReader = new CDiskMass2PepIndexReader();
//		}
//		pReader = new CMMapMass2PepIndexReader();
//		pReader = new CDiskMass2PepIndexReader();
	}
	catch(runtime_error &e)
	{
		if(pReader)
		{
			delete pReader;
		}
		CErrInfo info("CPeptideRandomReaderFactory", "GetPeptideReader");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		if(pReader)
		{
			delete pReader;
		}
		CErrInfo info("CPeptideRandomReaderFactory", "GetPeptideReader", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
	if(!pReader)
	{
		CErrInfo info("CPeptideRandomReaderFactory", "GetPeptideReader", "pReader is null");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
	return pReader;
}

}

