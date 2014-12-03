//============================================================================
// Name        : ProteomicsSDK.cpp
// Author      : Leheng Wang
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <string>
#include <vector>
#include "../include/predefine.h"
#include "ProteomicsSDK.h"

using namespace std;
using namespace proteomics_sdk;

namespace proteomics_sdk {

// add by emily 
// guarantee the correctness of path


void CCheck::CheckPath(string & strPath) {
	if (strPath.find('\"') != string::npos) {
		return;
	}

	bool bMod = false;

	for (size_t i=0; i<strPath.length(); ++i) {
		if ((strPath[i]&0x80) || strPath[i] == ' ' || strPath[i] == '	') {
			bMod = true;
			break;
		}
	}

	if (bMod) {
		strPath = "\"" + strPath + "\"";
	}

}

CErrInfo::CErrInfo(const string &strClass, const string &strMethod,
		const string &strDetail) {
	m_strClass = strClass;
	m_strMethod = strMethod;
	m_strDetail = strDetail;
}

CErrInfo::CErrInfo(const string &strClass, const string &strMethod,
		const string &strDetail, const exception & e) {
	m_strClass = strClass;
	m_strMethod = strMethod;
	m_strDetail = strDetail;
	SetException(e);
}
void CErrInfo::Append(const string &strInfo) {
	if (strInfo.empty())
		return;
	else
		m_strInfo += "\t\t  " + strInfo+"\n";
}

string CErrInfo::Get() const {
	string strError = m_strException;
	strError += "\t  at " + m_strClass + "::" + m_strMethod + "() "
			+ m_strDetail + "\n";
	strError += m_strInfo;
	return strError;
}

string CErrInfo::Get(const exception& e) {
	SetException(e);
	return Get();
}
void CErrInfo::SetException(const exception & e) {
	m_strException = e.what();
}

ofstream& operator<< ( ofstream& os, const CErrInfo& info)
{
	os << endl << "==========================" << endl;
	time_t current_time;
	time(&current_time);
	os << ctime(&current_time) << endl;
	os << info.Get() << endl;
	return os;
}
ostream& operator<< ( ostream& os, const CErrInfo& info)
{
	os << endl << "==========================" << endl;
	time_t current_time;
	time(&current_time);
	os << ctime(&current_time) << endl;
	os << info.Get() << endl;
	return os;
}

CAA::CAA(const CAA & r)
:m_cAA(r.m_cAA),m_strAA(r.m_strAA),m_strAADetail(r.m_strAADetail),
m_lfAvrgMass(r.m_lfAvrgMass),m_lfMonoMass(r.m_lfMonoMass),m_lfMod(r.m_lfMod),
m_tavrgMass(r.m_tavrgMass),m_tmonoMass(r.m_tmonoMass)
{	;}

CAA::CAA()
:m_cAA('A'),m_strAA("Ala"),m_strAADetail("Alanine"),
m_lfAvrgMass(71.0788f),m_lfMonoMass(71.03711f),m_lfMod(0)
{
	m_tavrgMass = size_t(m_lfAvrgMass * nDECIMAL);
	m_tmonoMass = size_t(m_lfMonoMass * nDECIMAL);
}

CAA::CAA( char cAA, string strAA, string strAADe,
		double lfAvg, double lfMono, double lfMod)
:m_cAA(cAA),m_strAA(strAA),m_strAADetail(strAADe),
m_lfAvrgMass(lfAvg),m_lfMonoMass(lfMono), m_lfMod(lfMod)
{
	m_tavrgMass = size_t(m_lfAvrgMass * nDECIMAL);
	m_tmonoMass = size_t(m_lfMonoMass * nDECIMAL);
}

CAA::~CAA(void)
{	m_strAA.clear();m_strAADetail.clear();};

void CAA::InitByString(string& strInit)
{
	istringstream in(strInit);
	in >> m_strAA >> m_strAADetail >> m_lfMonoMass >> m_lfAvrgMass >> m_lfMod;
	m_tavrgMass = size_t(m_lfAvrgMass * nDECIMAL);
	m_tmonoMass = size_t(m_lfMonoMass * nDECIMAL);
}

CMapAAMass::CMapAAMass(const CMapAAMass & m)
{
	m_mapAAMass.clear();
	m_mapAAMass = m.m_mapAAMass;
}

CMapAAMass::CMapAAMass() {
	m_mapAAMass.clear();
	m_mapAAMass['A'] = CAA('A',"Ala", "Alanine", 71.0788f, 71.03711f, 0.0f);
	m_mapAAMass['B'] = CAA('B',"_", "Avg.Asn_and_Asp",114.5962f,114.53494f, 0.0f);
	m_mapAAMass['C'] = CAA('C',"Cys", "Cysteine", 103.1388f, 103.00919f, 0.0f);
	m_mapAAMass['D'] = CAA('D', "Asp", "Aspartic", 115.0886f, 115.02694f, 0.0f);
	m_mapAAMass['E'] = CAA('E', "Glu", "Glutamic", 129.1155f, 129.04259f, 0.0f);
	m_mapAAMass['F'] = CAA('F', "Phe", "Phenylalanine", 147.1766f, 147.06841f, 0.0f);
	m_mapAAMass['G'] = CAA('G', "Gly", "Glycine", 57.0519f, 57.02146f, 0.0f);
	m_mapAAMass['H'] = CAA('H', "His", "Histidine", 137.1411f, 137.05891f, 0.0f);
	m_mapAAMass['I'] = CAA('I', "Ile", "Isoleucine", 113.1594f, 113.08406f, 0.0f);
	m_mapAAMass['J'] = CAA('J', "_", "_", 0.0f, 0.0f, 0.0f);
	m_mapAAMass['K'] = CAA('K', "Lys", "Lysine", 128.1741f, 128.09496f, 0.0f);
	m_mapAAMass['L'] = CAA('L', "Leu", "Leucine", 113.1594f, 113.08406f, 0.0f);
	m_mapAAMass['M'] = CAA('M', "Met", "Methionine", 131.1926f, 131.04049f, 0.0f);
	m_mapAAMass['N'] = CAA('N', "Asn" , "Asparagine", 114.1038f, 114.04293f, 0.0f);
	m_mapAAMass['O'] = CAA('O', "_", "Ornithine", 114.1472f, 114.07931f, 0.0f);
	m_mapAAMass['P'] = CAA('P', "Pro", "Proline", 97.1167f, 97.05276f, 0.0f);
	m_mapAAMass['Q'] = CAA('Q', "Gln", "Glutamine", 128.1307f, 128.05858f, 0.0f);
	m_mapAAMass['R'] = CAA('R', "Arg", "Arginine", 156.1875f, 156.10111f, 0.0f);
	m_mapAAMass['S'] = CAA('S', "Ser", "Serine", 87.0782f, 87.03203f, 0.0f);
	m_mapAAMass['T'] = CAA('T', "Thr", "Threonine", 101.1051f, 101.04768f, 0.0f);
	m_mapAAMass['U'] = CAA('U', "_", "_", 0.0f, 0.0f, 0.0f);
	m_mapAAMass['V'] = CAA('V', "Val", "Valine", 99.1326f, 99.06841f, 0.0f);
	m_mapAAMass['W'] = CAA('W', "Trp", "Tryptophan", 186.2132f, 186.07931f, 0.0f);
	m_mapAAMass['X'] = CAA('X', "_", "Leu_or_Ile", 113.1594f, 113.08406f, 0.0f);
	m_mapAAMass['Y'] = CAA('Y', "Tyr", "Tyrosine", 163.1760f, 163.06333f, 0.0f);
	m_mapAAMass['Z'] = CAA('Z', "_", "Avg.Gln_and_Glu", 128.6231f, 128.55059f, 0.0f);
}
void CMapAAMass::SetAAMap(char cAA, CAA& aa) {m_mapAAMass[cAA] = aa;};
CAA& CMapAAMass::GetAA(char cAA) {return m_mapAAMass[cAA];};

CXLinker::CXLinker():m_bHomo(true),m_lfMonoMass_dif(0.0f),m_lfAvrgMass_dif(0.0f),m_lfMLMonoMass_dif(0.0f),m_lfMLAvrgMass_dif(0.0f)
{
	m_strAlphaAA.reserve(10);
	m_strBetaAA.reserve(10);
}

CXLinker::~CXLinker() {}

void CXLinker::InitByString(string & strInit)
{
	istringstream in(strInit);
	in >> m_strAlphaAA >> m_strBetaAA >> m_lfMonoMass_dif >> m_lfAvrgMass_dif >> m_lfMLMonoMass_dif >> m_lfMLAvrgMass_dif;

	// todo : define homo 
	if(m_strAlphaAA == m_strBetaAA)
	{
		m_bHomo = true;
	}
	else
	{
		m_bHomo = false;
	}

	// N-terminal of protein is '['
	// C-terminal of protein is ']'
	// N-terminal of peptide is '('
	// C-terminal of peptide is ')'
	// all ammino is '*' 

	if(m_strAlphaAA == "*")
	{
		m_strAlphaAA = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	}
	if(m_strBetaAA == "*")
	{
		m_strBetaAA = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	}
}

CModification::CModification() :m_lfMonoMass_dif(0.0f),m_lfAvrgMass_dif(0.0f), m_tNLSize(0)
{
	m_strAA.reserve(10);
}

CModification::~CModification() {};

void CModification::InitByString(string& strInit)
{
	istringstream in(strInit);
	string strModType;
	in >> m_strAA >> strModType >> m_lfMonoMass_dif >> m_lfAvrgMass_dif >> m_tNLSize;
	if("PEP_N" == strModType)
	m_eModType = MT_PEP_NTERM;
	else if("PEP_C" == strModType)
	m_eModType = MT_PEP_CTERM;
	else if("PRO_N" == strModType)
	m_eModType = MT_PRO_NTERM;
	else if("PRO_C" == strModType)
	m_eModType = MT_PRO_CTERM;
	else m_eModType = MT_NORMAL;
	if(m_strAA == "[ALL]" || m_strAA == "[all]")
	{
		m_strAA = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	}
	double lfMonoMass, lfAvrgMass;
	m_vlfMonoNeutralLoss_dif.clear();
	m_vlfAvrgNeutralLoss_dif.clear();
	for(size_t i = 0;i < m_tNLSize;++i)
	{
		in >> lfMonoMass >> lfAvrgMass;
		m_vlfMonoNeutralLoss_dif.push_back(lfMonoMass);
		m_vlfAvrgNeutralLoss_dif.push_back(lfAvrgMass);
	}
}

CPeak::CPeak(double mz, double intensity, int nmz) {lfMz=mz;lfIntensity=intensity;nMz = nmz;nCharge=0;};

CPeak::CPeak() {lfMz=0;lfIntensity=0;nMz = 0;nCharge=0;};

CPeak& CPeak::operator=(const CPeak& rhs)
{
	if(&rhs != this)
	{
		lfMz = rhs.lfMz;
		lfIntensity = rhs.lfIntensity;
		nMz = rhs.nMz;
		nCharge = rhs.nCharge;
	}
	return *this;
}

CEnzyme::CEnzyme(): m_bNTerm(false),m_strName("Trypsin"),m_strCleave("KR"),m_strNotCleave("P") {};
CEnzyme::~CEnzyme() {};

CEnzyme& CEnzyme::operator=(const CEnzyme& rhs)
{
	if( &rhs != this )
	{
		m_bNTerm = rhs.m_bNTerm; // .GetIsNTerm();
		m_strName = rhs.m_strName; //.GetName();
		m_strCleave = rhs.m_strCleave; //.GetCleaveString();
		m_strNotCleave= rhs.m_strNotCleave; //.GetNotCleav();
	}
	return *this;
}

bool CEnzyme::GetIsNTerm() const {return m_bNTerm;};
string CEnzyme::GetCleaveString()const {return m_strCleave;};
string CEnzyme::GetName()const {return m_strName;};
string CEnzyme::GetNotCleave()const {return m_strNotCleave;};

void CEnzyme::SetIsNTerm(bool bNTerm) {m_bNTerm = bNTerm;};
void CEnzyme::SetCleaveString(const char * szCleave) {m_strCleave = szCleave;};
void CEnzyme::SetName(const char * szName) {m_strName = szName;};
void CEnzyme::SetNotCleave(const char * szNoClv) {m_strNotCleave = szNoClv;};

bool CEnzyme::empty(void) const {return (m_strName.empty() || m_strCleave.empty());};
void CEnzyme::InitByString(string& strInit)
{
	istringstream isIn(strInit);
	char cTemp;
	isIn >> m_strCleave >> m_strNotCleave >> cTemp;
	m_bNTerm = (cTemp == 'N');
}

CSimplePeptide::CSimplePeptide():m_tLength(0), m_lfMass(0.0), m_tModCnt(0),
		m_tFixedModCnt(0), m_cMissCleave(0), m_cPrev('*'), m_cNext('*'),m_cEnd(0)
{
	memset(m_szSequence, 0, sizeof(m_szSequence));
};
CSimplePeptide::CSimplePeptide(const CSimplePeptide & peptide)
{
	*this = peptide;
}

double CSimplePeptide::Calc_Theoretical_MH(const CSimplePeptide & pep, bool bPepMono)
{
	if(bPepMono)
		return pep.m_lfMass + IonMass_Mono_O + IonMass_Mono_H * 2.0 + IonMass_Proton;
	else
		return pep.m_lfMass + IonMass_Aver_O + IonMass_Aver_H * 3.0;
}

double CSimplePeptide::Calc_Theoretical_M(const CSimplePeptide & pep)
{
	return pep.m_lfMass + IonMass_Mono_O + IonMass_Mono_H * 2.0;
}

CSimplePeptide &CSimplePeptide::operator=(const CSimplePeptide & peptide)
{
	strcpy(m_szSequence, peptide.m_szSequence);
	m_tLength = peptide.m_tLength;
	m_lfMass = peptide.m_lfMass;
	m_tModCnt = peptide.m_tModCnt;
	m_cMissCleave = peptide.m_cMissCleave;
	m_cPrev = peptide.m_cPrev;
	m_cNext = peptide.m_cNext;
	m_cEnd = peptide.m_cEnd;
	memcpy(&m_tModSites[0][0], &peptide.m_tModSites[0][0], sizeof(m_tModSites));
	m_tFixedModCnt = peptide.m_tFixedModCnt;
	return *this;
}

void CSimplePeptide::SetPeptideInfor(const char * pszSequence, size_t tLength, double lfMass, unsigned char cMissCleave, unsigned char cEnd, bool bCopyNCTerm)
{
	if (tLength > MAX_PEPTIDE_LENGTH) // added at 2013.12.8 by fan, for long sequence
	{
		tLength = MAX_PEPTIDE_LENGTH; // only copy the MAX_PEPTIDE_LENGTH char.
	}
	if(bCopyNCTerm)
	{
		tLength -= 4;
		memcpy(m_szSequence, pszSequence + 2, tLength);
		m_szSequence[tLength] = 0;
		m_tLength = tLength;
		m_lfMass = lfMass;
		m_cMissCleave = cMissCleave;
		m_cPrev = pszSequence[0];
		m_cNext = pszSequence[tLength + 3];
		m_cEnd = cEnd;
		m_tModCnt = 0;
		m_tFixedModCnt = 0;
	}
	else
	{
		memcpy(m_szSequence, pszSequence, tLength); //changed from strcpy to memcpy at 2013.12.8
		m_szSequence[tLength] = 0;
		m_tLength = tLength;
		m_lfMass = lfMass;
		m_cMissCleave = cMissCleave;
		m_cEnd = cEnd;
		m_cPrev = '*';//added at 2013.11.19
		m_cNext = '*';
		m_tModCnt = 0;//added at 2013.11.21
		m_tFixedModCnt = 0;
	}
}

bool CSimplePeptide::EqualPeptide(const CSimplePeptide & peptide)
{
	if(strcmp(peptide.m_szSequence, m_szSequence) != 0
			|| peptide.m_tModCnt != m_tModCnt)
	return false;
	for(size_t i = 0;i < m_tModCnt;++i)
	{
		if(peptide.m_tModSites[i][0] != m_tModSites[i][0]
				|| peptide.m_tModSites[i][1] != m_tModSites[i][1])
		{
			return false;
		}
	}
	return true;
}

bool CSimplePeptide::IsProNTerm()
{
	//return '-' == m_cPrev;
	return (m_cEnd == 1)||(m_cEnd == 3);
}

bool CSimplePeptide::IsProNTerm() const
{
	return (m_cEnd == 1)||(m_cEnd == 3);
}
bool CSimplePeptide::IsProCTerm()
{
	//return '-' == m_cNext;
	return (m_cEnd == 2)||(m_cEnd == 3);
}
bool CSimplePeptide::IsProCTerm() const
{
	return (m_cEnd == 2)||(m_cEnd == 3);
}

CXLinkItem & CXLinkItem::operator=(const CXLinkItem & xlinkitem)
{
	setXLinkItemInfo(xlinkitem.m_eXLinkType, xlinkitem.m_tAlphaSite, xlinkitem.m_tBetaSite, xlinkitem.m_nLinkerId);
	return *this;
}

void CXLinkItem::setXLinkItemInfo(const XLinkType &eXLinkType, const int nAlphaSite,
		const int nBetaSite, const int nLinkerId)
{
	m_eXLinkType = eXLinkType;
	m_tAlphaSite = nAlphaSite;
	m_tBetaSite = nBetaSite;
	m_nLinkerId = nLinkerId;
}

CXLinkPepResult::CXLinkPepResult()
{
	m_bPair = false;
	m_lfScore = 0.0;
	m_lfCalc_MH = 0.0;
	m_lfEvalue = 1.0;
	m_bEV = false;
	m_lfFPR = 0;

	m_lfAlphaEvalue = 1.0;
	m_lfBetaEvalue = 1.0;
}

bool CXLinkPepResult::Score_Greater(const CXLinkPepResult & a, const CXLinkPepResult & b)
{
	return a.m_lfScore> b.m_lfScore;
}

void CXLinkPepResult::SwapAlphaBeta()
{
	CSimplePeptide tempPep(m_AlphaPeptide);
	m_AlphaPeptide = m_BetaPeptide;
	m_BetaPeptide = tempPep;
	int nSite = m_XLink.m_tAlphaSite;
	m_XLink.m_tAlphaSite = m_XLink.m_tBetaSite;
	m_XLink.m_tBetaSite = nSite;
	m_vAlphaProteinAC.swap(m_vBetaProteinAC);
	m_vAlphaProteinID.swap(m_vBetaProteinID);
	m_vAlphaProteinSite.swap(m_vBetaProteinSite);
}

CXLinkPepResult & CXLinkPepResult::operator=(const CXLinkPepResult & pepResult)
{
	m_AlphaPeptide = pepResult.m_AlphaPeptide;
	m_BetaPeptide = pepResult.m_BetaPeptide;
	m_XLink = pepResult.m_XLink;

	m_bEV = pepResult.m_bEV;
	m_bPair = pepResult.m_bPair;
	m_lfCalc_MH = pepResult.m_lfCalc_MH;
	m_lfEvalue = pepResult.m_lfEvalue;
	m_lfFPR = pepResult.m_lfFPR;
	m_lfScore = pepResult.m_lfScore;
	m_stMatchInfo = pepResult.m_stMatchInfo;
	m_vAlphaProteinAC = pepResult.m_vAlphaProteinAC;
	m_vAlphaProteinID = pepResult.m_vAlphaProteinID;
	m_vAlphaProteinSite = pepResult.m_vAlphaProteinSite;
	m_vBetaProteinAC = pepResult.m_vBetaProteinAC;
	m_vBetaProteinID = pepResult.m_vBetaProteinID;
	m_vBetaProteinSite = pepResult.m_vBetaProteinSite;

	m_lfAlphaEvalue = pepResult.m_lfAlphaEvalue;
	m_lfBetaEvalue = pepResult.m_lfBetaEvalue;

	return *this;
}

CXLinkOpenPepResult::CXLinkOpenPepResult()
{
	m_lfScore = 0.0;
	m_lfOpenMass = 0.0;
	m_nLinkSite = -1;
	// todo:
	m_nLinkerId = 0;
}
void CXLinkOpenPepResult::resetPeptideInfo(string strSq, double lfMass,
			unsigned char cMissCleave, unsigned char cEnd, bool bCopyNCTerm)
{
	m_peptide.SetPeptideInfor(strSq.c_str(), strSq.length(), lfMass,
			cMissCleave, cEnd, bCopyNCTerm);
	m_peptide.m_tModCnt = 0;
	m_lfScore = 0.0;
	m_lfOpenMass = 0.0;
	m_nLinkSite = -1;
	m_nLinkerId = -1;
}

CXLinkOpenMiddlePepResult::CXLinkOpenMiddlePepResult()
{
	m_lfScore = 0.0;
	m_lfOpenMass = 0.0;
	m_nLinkSite = -1;
	// todo:
	m_nLinkerId = 0;

	m_lfSidePepMass = 0.0;
	m_nSidePepSite = -1;
}
void CXLinkOpenMiddlePepResult::resetPeptideInfo(string strSq, double lfMass,
			unsigned char cMissCleave, unsigned char cEnd, bool bCopyNCTerm)
{
	m_peptide.SetPeptideInfor(strSq.c_str(), strSq.length(), lfMass,
			cMissCleave, cEnd, bCopyNCTerm);
	m_peptide.m_tModCnt = 0;
	m_lfScore = 0.0;
	m_lfOpenMass = 0.0;
	m_nLinkSite = -1;
	m_nLinkerId = -1;

	m_lfSidePepMass = 0.0;
	m_nSidePepSite = -1;
}

void CXLinkOpenMiddlePepResult::SetPep(CXLinkOpenPepResult pep)
{
	this->m_peptide = pep.m_peptide;
	m_lfScore = pep.m_lfScore;
	m_nLinkSite = pep.m_nLinkSite;
	m_nLinkerId = pep.m_nLinkerId;
	m_lfOtherMinMass = pep.m_lfOtherMinMass;
	m_lfOtherMaxMass = pep.m_lfOtherMaxMass;

	m_lfOpenMass = pep.m_lfOpenMass;
	m_lfSidePepMass = 0.0;
	m_nSidePepSite = -1;
}

CPeptideResult::CPeptideResult()
{
	m_lfEvalue = m_lfScore = 0;
}
bool CPeptideResult::Score_Greater(const CPeptideResult & a, const CPeptideResult & b)
{
	return a.m_lfScore> b.m_lfScore;
}

/*
 * refined-search needed
 * store evalue coefficients 
 */

EVCoefficient::EVCoefficient()
{
	lfCoef0 = 0;
	lfCoef1 = 0;
}

int EVCoefficient::GetSize()
{
	return 2*sizeof(double);
}

bool EVCoefficient::WriteToFile(FILE * fp)
{
	if(fp == NULL)
	return false;

	if(1 != fwrite(&lfCoef0,sizeof(double),1,fp))
	return false;

	if(1 != fwrite(&lfCoef1,sizeof(double),1,fp))
	return false;

	return true;
}

bool EVCoefficient::ReadFromFile(FILE * fp)
{
	if(fp == NULL)
	return false;

	if(1 != fread(&lfCoef0,sizeof(double),1,fp))
	return false;

	if(1 != fread(&lfCoef1,sizeof(double),1,fp))
	return false;

	return true;
}

CSpectrum::CSpectrum():m_nCharge(0), m_lfIntensity(0.0f), m_lfMH(0.0f), m_lfMZ(0.0f), m_strFilePath(""), m_pPeaks(NULL), m_tPeaksNum(0),m_tScanNo(0), m_lfRetentionTime(0.0f), m_strActivationType(""), m_tPrecursorScanNo(0), m_strInstrumentType("") {};
CSpectrum::~CSpectrum() {clear();};

void CSpectrum::CopyBasicItems(const CSpectrum & spec)
{
	m_nCharge = spec.m_nCharge;
	m_lfIntensity = spec.m_lfIntensity;
	m_lfSqtMaxInten = spec.m_lfSqtMaxInten;
	m_lfMH = spec.m_lfMH;
	m_lfMZ = spec.m_lfMZ;
	m_strFilePath = spec.m_strFilePath;
	m_vHash = spec.m_vHash;
	m_tScanNo = spec.m_tScanNo;
	m_lfRetentionTime = spec.m_lfRetentionTime;
	m_strActivationType = spec.m_strActivationType;
	m_tPrecursorScanNo = spec.m_tPrecursorScanNo;
	m_strInstrumentType = spec.m_strInstrumentType;
}

CSpectrum & CSpectrum::operator = (const CSpectrum & spec)
{
	CopyBasicItems(spec);
	clear();
	m_tPeaksNum = spec.m_tPeaksNum;
	if(0 != m_tPeaksNum)
	{
		m_pPeaks = new CPeak[m_tPeaksNum];
		memcpy(m_pPeaks, spec.m_pPeaks, m_tPeaksNum * sizeof(CPeak));
	}
	return *this;
}

CSpectrum::CSpectrum(const CSpectrum & spec)
{
	CopyBasicItems(spec);
	m_tPeaksNum = spec.m_tPeaksNum;
	if(0 != m_tPeaksNum)
	{
		m_pPeaks = new CPeak[m_tPeaksNum];
		memcpy(m_pPeaks, spec.m_pPeaks, m_tPeaksNum * sizeof(CPeak));
	}
}

bool CSpectrum::intesity_greater( const CPeak & elem1, const CPeak & elem2 ) 
{
    return elem1.lfIntensity> elem2.lfIntensity;
}

bool CSpectrum::intesity_lesser( const CPeak & elem1, const CPeak & elem2 ) 
{
    return elem1.lfIntensity < elem2.lfIntensity;
}

bool CSpectrum::mz_lesser(const CPeak & elem1, const CPeak & elem2 ) 
{
    return elem1.lfMz < elem2.lfMz;
}

void CSpectrum::clear(void)
{
	if(m_pPeaks != NULL)
	{
		if(0 != m_tPeaksNum)
		delete [] m_pPeaks;

		m_pPeaks = NULL;
		m_tPeaksNum = 0;
	}
}

void CSpectrum::ReAlloc(size_t tNum)
{
	if(m_tPeaksNum> 0)
	{
		CPeak * p = new CPeak[m_tPeaksNum];
		memcpy(p, m_pPeaks, sizeof(CPeak) * m_tPeaksNum);
		size_t t = m_tPeaksNum;
		clear();
		m_pPeaks = new CPeak[tNum];
		m_tPeaksNum = t;
		memcpy(m_pPeaks, p, sizeof(CPeak) * m_tPeaksNum);
		delete [] p;
	}
	else
	{
		if(NULL != m_pPeaks)
		{
			delete [] m_pPeaks;
		}

		if(tNum>0)
		m_pPeaks = new CPeak[tNum];
	}
}

void CSpectrum::Create_Hash_Table()
{
	m_vHash.clear();
	if(0 == m_tPeaksNum)
	return;
	m_vHash.reserve(10000);
	size_t tTemp1 = 0, tTemp2 = m_pPeaks[0].nMz / MZMULTIPLIER;

	for(size_t i = 0;i <= tTemp2;++i)
	{
		m_vHash.push_back(0);
	}

	for(size_t i = 1;i < m_tPeaksNum;++i)
	{
		if((tTemp1 = size_t(1.0 * m_pPeaks[i].nMz / MZMULTIPLIER))> tTemp2)
		{
			for(size_t j = tTemp2 + 1;j < tTemp1;++j)
			{
				m_vHash.push_back(i - 1);
			}
			m_vHash.push_back(i);
			tTemp2 = tTemp1;
		}
	}
}

CXLinkMatchResult::CXLinkMatchResult()
{
	m_tScore = 0;
	m_tRealCandidate = 0;
}

void CXLinkMatchResult::remove_invalid(void)
{
	vector<CXLinkPepResult> vTemp;
	size_t i = 0;
	for(i = 0;i < m_vPeptideResults.size();++i)
	{
		if(m_vPeptideResults[i].m_bEV)
		{
			vTemp.push_back(m_vPeptideResults[i]);
		}
	}
	vTemp.swap(m_vPeptideResults);
}

double CXLinkMatchResult::Calc_Theoretical_MH(const CXLinkPepResult & pep_res, bool bPepMono)
{
	if(pep_res.m_XLink.m_eXLinkType == X_LINK)
	{
		if(bPepMono)
		return pep_res.m_AlphaPeptide.m_lfMass + pep_res.m_BetaPeptide.m_lfMass + 2*(IonMass_Mono_O + IonMass_Mono_H * 2.0) + IonMass_Proton;
		else
		return pep_res.m_AlphaPeptide.m_lfMass + pep_res.m_BetaPeptide.m_lfMass + 2*(IonMass_Aver_O + IonMass_Aver_H * 2.0) + IonMass_Aver_H;
	}
	else
	{
		if(bPepMono)
		return pep_res.m_AlphaPeptide.m_lfMass + IonMass_Mono_O + IonMass_Mono_H * 2.0 + IonMass_Proton;
		else
		return pep_res.m_AlphaPeptide.m_lfMass + IonMass_Aver_O + IonMass_Aver_H * 3.0;
	}

}

CXLinkMatchResult & CXLinkMatchResult::operator=(const CXLinkMatchResult & matchResult)
{
	m_tRealCandidate = matchResult.m_tRealCandidate;
	m_tScore = matchResult.m_tScore;
	m_vPeptideResults = matchResult.m_vPeptideResults;
	m_vlfScores = matchResult.m_vlfScores;
	return *this;
}

CProtein::CProtein() {}
CProtein::~CProtein() {}

CProtein::CProtein(string strAC, string strDE, string strSQ)
:m_strAC(strAC),m_strDE(strDE), m_strSQ(strSQ) {}

CProtein::CProtein(const CProtein & pro)
:m_nID(pro.m_nID), m_strAC(pro.m_strAC),m_strDE(pro.m_strDE), m_strSQ(pro.m_strSQ) {}

CProtein CProtein::Reverse(void)
{
	CProtein RevPro;
	RevPro.m_strAC = "REVERSE_" + m_strAC;
	RevPro.m_strDE = m_strDE;
	RevPro.m_strSQ.assign(m_strSQ.rbegin(), m_strSQ.rend());
	return RevPro;
}

void CCleaver::Cleave_specific(const string & strSQ)
{
	m_viCleaves.clear();
	m_vlfMasses.clear();

	int len = (int)strSQ.length();
	int i = 0;

	m_viCleaves.push_back(-1);

	double lfMass = 0.0;

	for(i = 0;i < len;++i)
	{
		lfMass += m_lfMass[strSQ[i] - 'A'];
		if(m_enzyme.m_strCleave.find(strSQ[i]) != string::npos)
		{
			m_viCleaves.push_back(i);
			m_vlfMasses.push_back(lfMass);
			lfMass = 0.0;
		}

	}
	if(len - 1 - m_viCleaves.back()> 0)
	{
		m_viCleaves.push_back(len - 1);
		m_vlfMasses.push_back(lfMass);
	}

}

void CCleaver::SetEnzyme(const CEnzyme & enzyme)
{
	m_enzyme = enzyme;
}

void CCleaver::InitValue(double pMass[], const CEnzyme & enzyme)
{
	memcpy(m_lfMass, pMass, sizeof(double) * 26);
	SetEnzyme(enzyme);
}

void CAssignedProtein::EmptyContents() {m_vpContainPep.clear();vector< pair<const CPeptideResult * ,size_t> >().swap(m_vpContainPep);};

double CAssignedProtein::GetTotalScore(void) const
{
	double total = 0;
	for(vector< pair<const CPeptideResult * ,size_t> >::const_iterator it = m_vpContainPep.begin(); it != m_vpContainPep.end(); ++it)
	total += it->first->m_lfScore;
	return total;
}

bool CAssignedProtein::Total_Score_Greater(const CAssignedProtein& elem1,CAssignedProtein& elem2)
{
	return ( elem1.GetTotalScore()> elem2.GetTotalScore());
}

double CAssignedProtein::GetTotalEValue(void) const
{
	double total = 1;
	for(vector< pair<const CPeptideResult * ,size_t> >::const_iterator it = m_vpContainPep.begin(); it != m_vpContainPep.end(); ++it)
	total *= it->first->m_lfEvalue;
	return total;
}

bool CAssignedProtein::Total_EValue_Lesser(const CAssignedProtein& elem1,const CAssignedProtein& elem2)
{
	return ( elem1.GetTotalEValue() < elem2.GetTotalEValue());
}

double CAssignedProtein::GetTotalFPR(void) const
{
	double total = 1;
	for(vector< pair<const CPeptideResult * ,size_t> >::const_iterator it = m_vpContainPep.begin(); it != m_vpContainPep.end(); ++it)
	total *= (*it).first->m_lfFPR;
	return total;
}

double CAssignedProtein::GetPepAverRank(void) const
{
	double total = 0;
	for(vector< pair<const CPeptideResult * ,size_t> >::const_iterator it = m_vpContainPep.begin(); it != m_vpContainPep.end(); ++it)
	total += it->second + 1;
	total /= m_vpContainPep.size();
	return total;
}
#ifndef IN_VC_PROGRAM_

bool CFileFind::Open(const string & strFilePath)
{
	m_pDir = opendir(strFilePath.c_str());
	if (0 == m_pDir)
	{
		return false;
	}
	return true;
}
bool CFileFind::GetNextFile()
{
	return (0 != (m_pent=readdir(m_pDir)));

}
void CFileFind::Close()
{
	closedir(m_pDir);
	m_pDir = 0;
	m_pent = 0;

}
#endif

CSearchState::CSearchState()
{
	m_lfStartTime = clock();
	time(&m_tStartTime);
	m_nTotalSpectra = 0;
	m_nCurrentSpectra = 0;
	m_bSet = false;
}

//add for xlink output;
void CSearchState::XlinkConsole()
{
#ifdef WIN32
	system("cls");
#else
	system("clear");
#endif
	const int WIDTH = 60;
	for(int i = 0;i < WIDTH;++i)
	{
		printf("%c", '*');
	}
	printf("\n");
	//output the condition
	//printf("%s\n", m_strCondition.c_str());
	//output the time
	printf("Start Time: %s\n", ctime(&m_tStartTime));
	printf("Time cost: %.3lf\n", (clock() - m_lfStartTime) / CLOCKS_PER_SEC);
	//output the spectra
	printf("Total Spectra: %d, Completed: %d\n", m_nTotalSpectra,
			m_nCurrentSpectra);
	for(int i = 0;i < WIDTH;++i)
	{
		printf("%c", '-');
	}
	printf("\n");
	//output the processor
	for(size_t i = 0;i < m_vThreadState.size();++i)
	{
		printf("Thread%d: current DB ratio(%f) total DB ratio(%f)\n",
				(int)m_vThreadState[i].m_tID, m_vThreadState[i].m_lfCurrentDbRatio,
				m_vThreadState[i].m_lfTotalDbRatio);
	}
	for(int i = 0;i < WIDTH;++i)
	{
		printf("%c", '*');
	}
	printf("\n");
}

void CSearchState::GetProgress(string & strInfo)
{
	strInfo = "";
	char szbuf[1024];

	sprintf(szbuf,"Searching...\n\tStart Time: %s", ctime(&m_tStartTime));
	strInfo += szbuf;

	sprintf(szbuf,"\tTime cost: %.3lf\n", (clock() - m_lfStartTime) / CLOCKS_PER_SEC);
	strInfo += szbuf;

	//output the spectra
	int nCompleted = m_nCurrentSpectra;

	//output the processor
	for(size_t i = 0;i < m_vThreadState.size();++i)
	{
		sprintf(szbuf,"\tThread%d : %d / %d\n",
				(int)m_vThreadState[i].m_tID, m_vThreadState[i].m_nCurrentSpectra,
				m_vThreadState[i].m_nTotalSpectra);
		nCompleted += m_vThreadState[i].m_nCurrentSpectra;
		strInfo += szbuf;
	}

	sprintf(szbuf,"\tTotal : %d / %d\n", nCompleted, m_nTotalSpectra);
	strInfo += szbuf;
}

void CSearchState::SetTotalSpectra(int nTotal)
{
	while(m_bSet) {}
	m_bSet = true;
	m_nTotalSpectra = nTotal;
	m_bSet = false;
}
void CSearchState::SetCurrentSpectra(int nCurrent)
{
	while(m_bSet) {}
	m_bSet = true;
	m_nCurrentSpectra = nCurrent;
	m_bSet = false;
}
void CSearchState::SetThreadTotal(size_t tID, int nTotal)
{
	while(m_bSet) {}
	m_bSet = true;
	m_vThreadState[tID].m_nTotalSpectra = nTotal;
	m_bSet = false;
}
void CSearchState::SetThreadCurrent(size_t tID, int nCurrent)
{
	while(m_bSet) {}
	m_bSet = true;
	m_vThreadState[tID].m_nCurrentSpectra = nCurrent;
	m_bSet = false;
}
// set for xlink output
void CSearchState::XlinkSetThreadCurrent(size_t tID, float lfCur, float lfTotal)
{
	while(m_bSet) {}
	m_bSet = true;
	m_vThreadState[tID].m_lfCurrentDbRatio = lfCur;
	m_vThreadState[tID].m_lfTotalDbRatio = lfTotal;
	m_bSet = false;
}

} // end of namespace proteomics_sdk
