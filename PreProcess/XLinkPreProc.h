#ifndef XLINKPREPROC_H_
#define XLINKPREPROC_H_

namespace proteomics_sdk
{

class CSpectrum;
class CCondition;
class CPreProcess;

struct PeakInfo
{
	int nPeakType;
	// 0: common
	// 1: isotop
	// 2: -NH3
	// 3: -H2O
	// 4: MH
	// 5: noise
	int nMainPeakId;
	double lfCalDeltaMz;
	double lfExpDeltaMz;
	int nCharge;
	
	PeakInfo()
	{
		nPeakType = 0;
		nMainPeakId = -1;
		lfCalDeltaMz = 0.0;
		lfExpDeltaMz = 0.0;
		nCharge = 0;
	}
};

class CXLinkPreProc : public CPreProcess
{

public:
	CXLinkPreProc();
	virtual ~CXLinkPreProc();
	virtual void Init(CCondition & condition);
	virtual bool Instrument(InstrumentType eType) const;
	virtual void Run(CSpectrum Input, CSpectrum &Output) ;
	virtual void Close(void);
protected:
	bool _IsMatched(double lfCalMz,double lfExpMz);
	CCondition m_Condition;
};

}

#endif /*XLINKPREPROC_H_*/
