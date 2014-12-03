#ifndef MASSRANGEWRITER_H_
#define MASSRANGEWRITER_H_

#include "Mass2PepIndex.h"
//using namespace ProteinIndex;
namespace Mass2PepIndex{
class CMassRangeWriter
{
public:
	CMassRangeWriter();
	virtual ~CMassRangeWriter();
	
	void Init(PINDEX_HEAD stPindexHead, PINDEX_ITEM vstPindexHead);//init member
	void Close();
	
	void WriteIndex(double lfLowMass, double lfHighMass)
	{
		
	}
	void GetMetaITEM()
	{
		
	}
};

}
#endif /*MASSRANGEWRITER_H_*/
