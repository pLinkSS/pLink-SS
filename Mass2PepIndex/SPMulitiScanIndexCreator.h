#ifndef MULITISCANINDEXCREATER_H_
#define MULITISCANINDEXCREATER_H_

#include "Mass2PepIndex.h"
#include "Mass2PepIndexCreator.h"

namespace Mass2PepIndex{
class CSPMulitiScanIndexCreator:public CMass2PepIndexCreator
{
public:
	CSPMulitiScanIndexCreator();
	virtual ~CSPMulitiScanIndexCreator();
	
	virtual void Init(const PINDEX_HEAD &pIndexHead, const PINDEX_ITEM & pIndexItem);
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
#endif /*MULITISCANINDEXCREATER_H_*/
