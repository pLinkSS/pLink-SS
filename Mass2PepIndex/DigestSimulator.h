#ifndef DIGESTSIMULATOR_H_
#define DIGESTSIMULATOR_H_

#include<iomanip>
#include "Mass2PepIndex.h"
#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/option.h"
#include "../ProteinIndex/PepCalcFunc.h"
#include "PepFilter.h"

using namespace std;
using namespace proteomics_sdk;



namespace Mass2PepIndex{

class CDigestSimulator
{
public:
	CDigestSimulator(void);
	virtual~CDigestSimulator(void);
	// Clear the contents -- the vector of the cleaved peptides
	void		EmptyContent(void);
	vector<PEP_INFO>&	GetvSubPeptides() ;
	void	DigestPros2(string& strProteinSQ,const proteomics_sdk::CEnzyme& Enzyme, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<PEP_INFO_EX_TAG>	& m_vsDat, size_t tProPos, size_t tProID, size_t tDatNum);
			
	void	SemiDigestPros1(string& strProteinSQ,const proteomics_sdk::CEnzyme& Enzyme, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<PEP_INFO_EX_TAG>	& m_vsDat, size_t tProPos, size_t tProID, size_t tDatNum);

	void	NoneEnzymeDigestPros1(string& strProteinSQ, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<PEP_INFO_EX_TAG>	& m_vsDat, size_t tProPos, size_t tProID, size_t tDatNum);

	void	ComputePepNumAtEachMassBin(string& strProteinSQ, const proteomics_sdk::CEnzyme& Enzyme, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<size_t> & vtPepMassNum, size_t tBitNum);
	void	ComputeSemiPepNumAtEachMassBin(string& strProteinSQ, const proteomics_sdk::CEnzyme& Enzyme, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<size_t> & vtPepMassNum, size_t tBitNum);
	void	ComputeNoneEnzymeAtEachMassBin(string& strProteinSQ, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<size_t> & vtPepMassNum, size_t tBitNum);
protected:
	PEP_INFO			m_stDatBlock;
	
	vector<PEP_INFO>	m_vsDatBlock;
	PEP_INFO_EX			m_vsDatPep2ProBlock;
	PEP_INFO_EX_TAG		m_stDatExTagBlock;
	
	vector<double>		m_vdfPrime;
};

}

#endif /*DIGESTSIMULATOR_H_*/
