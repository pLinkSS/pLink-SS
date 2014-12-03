#ifndef COMMONCLASS_H_
#define COMMONCLASS_H_
namespace proteomics_sdk
{
typedef unsigned char UCHAR;
struct CMzTriple
{
	int nMz;
	int nIonTypeOrder;
	int nPepPosOrder1;
	// for internal ions
	int nPepPosOrder2;
	bool bMatched; // add by fan, 2013.7.5
	// add to count for intra-peptide continuous tag
	bool bNTerm1 ;
	bool bNTerm2 ;
	bool bAddToTag;
	bool bContainLinker;
	int nMassError; //added at 2013.12.9, record the mass difference between theory and experiment Note: has multiple the multiplier
	double dMatchInt; // added at 2013.12.9, record the largest intensity this ion matched.
	int nRKNH; //added at 2014.1.3, record how may R/K/NH3/H in the theroy ion
	double dAAChangeProb;
	double dAAChargeProb;

	static bool Less(const CMzTriple & a, const CMzTriple & b)
	{
		return a.nMz < b.nMz;
	}
	
	void clear()
	{
		nMz = 0;
		nIonTypeOrder = -1;
		nPepPosOrder1 = -1;
		nPepPosOrder2 = -1;
		bNTerm1 = false;
		bNTerm2 = false;
		bAddToTag = false;
		bContainLinker = false;
		bMatched = false;
		nMassError = 100*MZMULTIPLIER; //Notice: for ion mass error large than 100 Da set, will have problem
		dMatchInt = 0.0;
		nRKNH = 0;
		dAAChangeProb = 0;
		dAAChargeProb = 0;
	}

	void MatchInfoClear() // Clear the ion info associated with the match to spectrum
	{
		bMatched = false;
		nMassError = 100*MZMULTIPLIER; //Notice: for ion mass error large than 100 Da set, will have problem
		dMatchInt = 0.0;
		dAAChangeProb = 0;
		dAAChargeProb = 0;
	}
};


const int CO_H = int((IonMass_Mono_C + IonMass_Mono_O - IonMass_Mono_H) * MZMULTIPLIER + 0.5);
const int nhmass_mono_multi = (int)(IonMass_Mono_H*MZMULTIPLIER + 0.5);
const int nomass_mono_multi = (int)(IonMass_Mono_O*MZMULTIPLIER + 0.5);
const int nhmass_mono_multi_3 = (int)(IonMass_Mono_H*MZMULTIPLIER*3 + 0.5);
const int npromass_mono_multi = (int)(IonMass_Proton * MZMULTIPLIER + 0.5);

const int nhmass_avrg_multi = (int)(IonMass_Aver_H*MZMULTIPLIER + 0.5);
const int nomass_avrg_multi = (int)(IonMass_Aver_O*MZMULTIPLIER + 0.5);
const int npmass_multi = (int)(IonMass_Proton*MZMULTIPLIER + 0.5);

}
#endif /*COMMONCLASS_H_*/
