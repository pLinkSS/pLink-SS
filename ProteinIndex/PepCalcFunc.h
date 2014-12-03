#ifndef CPEPCALCFUNC_H_
#define CPEPCALCFUNC_H_

using namespace std;
using namespace proteomics_sdk;



class CPepCalcFunc
{

public:
	CPepCalcFunc(void);
	CPepCalcFunc(std::string strAAPath);
	~CPepCalcFunc(void);

	void Init(std::string strAAPath);
	// only add all AA residues mass
	void	SumSQResiduesMass(double&	fSum,const std::string& strSQ,bool bMono);
	void	SumSQResiduesMass(size_t&	tSum,const std::string& strSQ,bool bMono);

	// calc the precursor MH+ mass (add 18 +1)
	void	CalcPepMHMass(double&  fMH, const std::string& strSQ,bool bMono);
	void	CalcPepMHMass(size_t& tMH, const std::string& strSQ,bool bMono);

	// The following functions only used for CIonMassIndexProcess for speed.
	void SetPepSQ(std::string& strSQ)
	{
		m_strPepSQ = strSQ;
	}
	
	void SetMapAA();
	void SetAAPath(std::string strAAPath);

	// Get the double-type mass of specified residue
	double GetfAAMass(const char cAA, bool bMono);
	
	// Get the int-type mass of specified residue
	size_t GettAAMass(const char cAA, bool bMono);

protected:
	// map of the  residue's mass
	proteomics_sdk::CMapAAMass m_mapAA;

	std::string m_strPepSQ;
	std::string m_strAAPath;
};

#endif
