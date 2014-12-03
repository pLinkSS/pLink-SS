#ifndef MPMULITISCANINDEXCREATER_H_
#define MPMULITISCANINDEXCREATER_H_

#include "Mass2PepIndex.h"
#include "Mass2PepIndexCreator.h"

namespace Mass2PepIndex{
class CMPMulitiScanIndexCreator : public CMass2PepIndexCreator
{
public:
	CMPMulitiScanIndexCreator();
	virtual ~CMPMulitiScanIndexCreator();
	
	virtual void Init(const PINDEX_HEAD & pIndexHead, const PINDEX_ITEM & pIndexItem);
	virtual void Close();
	
	virtual void WriteIndex();
	
protected:

	PINDEX_HEAD m_PIndexHead;
	PINDEX_ITEM m_PIndexItem;
	
	vector <double> m_vdfMassRange;
	
	void _ScanMassRange();
	
	void _Excute(double lfLowMass, double lfHighMass);
	void _Dofork();
	
	void _WriteMeta();
	
};

}
#endif /*MPMULITISCANINDEXCREATER_H_*/
