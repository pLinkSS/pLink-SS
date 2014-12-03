#include "MassRangeCounter.h"
#include "MassRangeWriter.h"

#include "SPMulitiScanIndexCreator.h"

namespace Mass2PepIndex{
CSPMulitiScanIndexCreator::CSPMulitiScanIndexCreator()
{
}

CSPMulitiScanIndexCreator::~CSPMulitiScanIndexCreator()
{
}

void CSPMulitiScanIndexCreator::Close()
{
	
}

void CSPMulitiScanIndexCreator::Init(const PINDEX_HEAD &pIndexHead, const PINDEX_ITEM & pIndexItem) 
{
	m_PIndexHead = pIndexHead;
	m_PIndexItem = pIndexItem;
}

void CSPMulitiScanIndexCreator::WriteIndex()
{

}

void CSPMulitiScanIndexCreator::_ScanMassRange()
{	

}

void CSPMulitiScanIndexCreator::_Excute(double lfLowMass, double lfHighMass)
{

}

void CSPMulitiScanIndexCreator::_Dofork()
{

}

void CSPMulitiScanIndexCreator::_WriteMeta()
{
	
}
}
