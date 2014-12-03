#ifndef MASSRANGECOUNTER_H_
#define MASSRANGECOUNTER_H_

#include "Mass2PepIndex.h"
//using namespace ProteinIndex;
namespace Mass2PepIndex{
class CMassRangeCounter
{
public:
	CMassRangeCounter();
	virtual ~CMassRangeCounter();
	
	void Init(PINDEX_HEAD stPindexHead, PINDEX_ITEM vstPindexHead);//init member
	void Close();
	
	vector<double> GetMassRange()
	{
		vector<double> vdMassRange;
		return vdMassRange;
	}
};

}
#endif /*MASSRANGECOUNTER_H_*/
