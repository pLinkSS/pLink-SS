#ifndef SINGLESCANINDEXCREATER_H_
#define SINGLESCANINDEXCREATER_H_

#include "Mass2PepIndex.h"
#include "SPSPeptideExecuter.h"
#include "Mass2PepIndexCreator.h"

using namespace std;
using namespace ProteinIndex;

namespace Mass2PepIndex{
//typedef long long llong;

class CSPSingleScanIndexCreator : public CMass2PepIndexCreator
{
public:
	CSPSingleScanIndexCreator();
	virtual ~CSPSingleScanIndexCreator();

	virtual void Init(const PINDEX_HEAD &pIndexHead, const PINDEX_ITEM & pIndexItem, bool);
	virtual void Close();
	
	virtual void WriteIndex();


protected:

	PINDEX_HEAD m_PIndexHead;
	PINDEX_ITEM m_PIndexItem;
	
	CSPSPeptideExecuter m_peptideExecuter;

	CRandomReader * m_proteinExecuter;//virrual can not inside of another virtual?

//	vector <PEP_INFO> m_vpMemPep;
};
}
#endif /*SINGLESCANINDEXCREATER_H_*/
