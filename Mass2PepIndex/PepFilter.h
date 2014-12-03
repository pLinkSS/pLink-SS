#ifndef CPEPFILTER_H_
#define CPEPFILTER_H_

namespace Mass2PepIndex{
class CPepFilter
{
public:

	CPepFilter(void);
	~CPepFilter(void);

	//Determining whether input Peptide'sequence is satisfied conditions.  true means passed. false means filtered.
	bool IsPassed( size_t tLen = 0, double lfPepMass = 0 ,size_t tMissedSite = 0);

	void SetMaxMissedClvg( size_t tMissedSites );
	void SetMaxPepLength( size_t tMax );
	void SetMinPepLength( size_t tMin );
	void SetMaxPepMass( double lfMaxMass );
	void SetMinPepMass( double lfMinMass );
	void SetMaxPepMassAllowed(double lfMaxMass);
	void SetMinPepMassAllowed(double lfMinMass);

	size_t GetMaxMissedClvg( void );
	size_t GetMaxPepLength( void );
	size_t GetMinPepLength( void );
	double  GetMaxPepMass( void );
	double  GetMinPepMass( void );
	double  GetMaxPepMassAllowed( void );
	double  GetMinPepMassAllowed( void );

protected:

	// Max missed cleaved site, default = 1
	size_t	m_tMaxMissedClvg;

	// The minimal length of peptides to be indexed, default is 5, suggested not blow 4.
	size_t m_tMinPepLength;

	// The maximal length of peptides to be indexed, default is 50, suggested not upper to 100.
	size_t m_tMaxPepLength;

	// The minimal mass value of peptides to be indexed, default is 400, suggested not blow 500.
	double m_lfMinPepMass;

	// The maximal mass value of peptides to be indexed, default is 10000, suggested not upper to 6000.
	double m_lfMaxPepMass;
	
	double m_lfMinPepMassAllowed;
	double m_lfMaxPepMassAllowed;

	bool IsUpperMaxMissClvg(size_t tMissSite);
	bool IsUpperMaxPepLength(size_t tPepLen);
	bool IsUpperMaxPepMass(double lfMass);
	bool IsLowerMinPepLength(size_t tPepLen);
	bool IsLowerMinPepMass(double lfMass);
};
}

#endif
