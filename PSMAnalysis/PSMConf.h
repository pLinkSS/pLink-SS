#ifndef PSMCONF_H_
#define PSMCONF_H_

/*
 *	format of *.psm:
 *  [pfd] 
 * 	pFindFile=
	OutputPath=
	TmpFileTotal=
	TmpFileTag=
	SpectraTotal=
 *	[ion]
 * 	tolerance=0.5	// Da ; fragment tolerance
 * 	ion_type_total=4
	ion_type1=b 1 0 0 0
	ion_type2=b 2 0 0 0
	ion_type3=y 1 0 0 0
	ion_type4=y 2 0 0 0
 * 	[gap]
 * 	mass_scope=40 // Da ; mass gap to consider between unmatched and matched peaks
 *	
 */	

//recording the ion type information
struct CIonTypeEx
{
	char cType;
	// 0 : N-term ions
	// 1 : C-term ions
	// 2 : Internal ions
	// 3 : Precursor ions (one peptide) 
	// 4 : Precursor ions (two peptide)
	// 5 : Reporter ions
	int nCharge;
	int nTotalLostVal;
	int nLossNH3;
	int nLossH2O;
	double lfOtherLoss;
	string strName;
};

class CPSMConf
{
public:
	CPSMConf();
	virtual ~CPSMConf();
	void Load(string strConfigFile);
	void Output();
	
	string m_strWorkDir;
	string m_strOption;
	
	string m_strpFindFile;
	string m_strReportFile;
	string m_strInputType;
	string m_strOutputPath;
	
	double m_lfFragmentTol;
	string m_strFragmentTolType;
	vector<CIonTypeEx> m_vIonTypes;
	double m_lfMassScope;
	
	string m_strSpectraList;
	string m_strSpectraListPrefix;
	bool m_bNoiseRemove;
	// condition for pfind
	CCondition m_Condition;
	
protected:
	void _LoadPFD();
	void _LoadIon();
	void _LoadGap();
	void _LoadFilter();
};

#endif /*PSMCONF_H_*/
