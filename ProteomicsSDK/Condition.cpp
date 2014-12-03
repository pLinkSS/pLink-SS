#include <string>
#include <vector>
#include "../include/predefine.h"
#include "ProteomicsSDK.h"
#include "Condition.h"
using namespace std;
using namespace proteomics_sdk;

namespace proteomics_sdk {

// modify by emily
CCondition::CCondition(void) :
		m_bPepMono(true), m_bFragmentMono(true),
		m_lfFragmentTolBase(0), m_lfFragmentTol(0),
		m_eScoreMethod(SCORE_KSDP),
		m_eInstrumentType(INSTRUMENT_UNKNOWN),
		m_ePreProcMethod(PRE_PROC_XLINK_HCD),
		m_eSimplePreProcMethod(PRE_PROC_XLINK_HCD),
		m_nReportPep(10), m_nSimpleReportPep(500), m_nReportPro(2000),
		m_nMaxMissCleaves(2),
		m_eFlow(FLOW_ALL_IDX), m_eEV(EV_DEFAULT), m_ePepGen(PEPGENE_DEFAULT),
		m_lfFPRThreshold(0.01), m_strFPSign("REVERSE_"),
		m_lfMinPepTol(0), m_strMinPepTolType("Da"),
		m_nKSDP_L(0), m_lfKSDP_G(0.0f), m_lfKSDP_A(0.0f),
		m_nMaxModifyNumber(MAX_MODIFY_NUM),
		m_lfSimpleScoreCutoff(2.0),
		m_eIndex(MEMORY_MAPPED),
		m_lfMaxEV(1.0f),
		m_tSalvoBatchSize(SALVO_BATCH_SIZE),
		m_tProcessorNum(1),
		m_tMinScoreNum(MIN_SCORE_NUM), m_tMaxScoreNum(MAX_SCORE_NUM),
		m_eInferType(PROTEIN_INFER_DEFAULT),
		m_nMinBlockNum(2), m_nMaxBlockSize(5000),
		m_eSpectraIndex(SPECTRA_INDEX_SIMPLE),
		m_eLoadBalance(LOAD_BALANCE_MASS_DYNAMIC_MUL3),
		m_tInterContinuousWndNum(0),
		m_bRSSaveInfo(false), m_bRSLoadInfo(false),
		m_strRSPath("temp\\"), m_lfRSFDR(0.1),
		m_strRSPindexPath("temp.pindex"), m_strRSDBName("temp"),
		m_strRSTag("temp"),
		m_eModule(MODULE_ALL),
		m_eLogRank(LOG_RANK_INFO), m_eLogAppender(LOG_APPENDER_CMD),
		m_strSumoSeq("MLGG"), m_strSumoSite('G'), m_strFixSumoSite('K')
{
}

bool CCondition::ValidateElement(void) const {
	return (!m_vIonTypes.empty() || !m_vAllDataBaseName.empty()
			|| !m_vAllEnzyme.empty() || !m_vAllModification.empty()
			|| m_mapAA.m_mapAAMass.size() != 27);
}

bool CCondition::ValidateAll(void) const {
	return (!ValidateElement() || !m_vSelectedDBName.empty()
			|| !m_SelectedEnzyme.empty() || !m_strFragmentTolType.empty());
}

void CCondition::clear(void) {
	m_vIonTypes.clear();
	m_vSimpleIonTypes.clear();
	m_vAllDataBaseName.clear();
	m_vAllEnzyme.clear();
	m_vAllModification.clear();
	m_vSelectedDBName.clear();
	m_vSelectedFixMod.clear();
	m_vSelectedVarMod.clear();
	m_mapAA.m_mapAAMass.clear();
	// pfind-xlink
	m_vAllXLinker.clear();
	m_vPepTolWnds.clear();
	m_vSelectedXLinker.clear();
}

void CCondition::SetDefaultIons(void) {
	m_vIonTypes.clear();
	m_vSimpleIonTypes.clear();
	CIonType iontype;
	for (int i = 0; i < ELEMENT_NUMBER; ++i)
		iontype.nLost[i] = 0;
	iontype.bNTerm = true;
	iontype.nCharge = 1;
	iontype.nTotalLostVal = 0;
	m_vIonTypes.push_back(iontype);
	m_vSimpleIonTypes.push_back(iontype);
	iontype.nCharge = 2;
	m_vIonTypes.push_back(iontype);
	m_vSimpleIonTypes.push_back(iontype);
	iontype.bNTerm = false;
	m_vIonTypes.push_back(iontype);
	m_vSimpleIonTypes.push_back(iontype);
	iontype.nCharge = 1;
	m_vIonTypes.push_back(iontype);
	m_vSimpleIonTypes.push_back(iontype);
}

double CCondition::GetExpMass(double lfMH, int nCharge) const {
	double lfExpMass = 0;
	if (m_bPepMono)
		lfExpMass = lfMH + IonMass_Proton * (nCharge - 1)
				- IonMass_Proton * nCharge;
	//todo : lfExpMass = lfMH + IonMass_Mono_H * (nCharge-1) - IonMass_Proton* nCharge;

	else
		lfExpMass = lfMH - IonMass_Aver_H;

	return lfExpMass;
}

double CCondition::GetMZ(double lfMH, int nCharge) const {
	double lfMZ = 0;
	if (m_bPepMono)
		lfMZ = (lfMH + (nCharge - 1) * IonMass_Proton) / nCharge;
	else
		lfMZ = (lfMH + (nCharge - 1) * IonMass_Aver_H) / nCharge;
	return lfMZ;
}
;
double CCondition::GetNeutraMass(double lfExpMass) const {
	double lfPepMass = 0;

	if (m_bPepMono)
		lfPepMass = lfExpMass - IonMass_Mono_O - IonMass_Mono_H * 2.0;
	else
		lfPepMass = lfExpMass - IonMass_Aver_O - IonMass_Aver_H * 2.0;

	return lfPepMass;
}

// modify by emily
double CCondition::GetPepTolBase(double lfMZ, int nCharge) const {
	return GetPepTolBase(0, lfMZ, nCharge);
}

double CCondition::GetPepTolBase(int nWndId, double lfMZ, int nCharge) const {
	double lfTolBase = 0;

	if (nWndId >= 0 && nWndId < m_vPepTolWnds.size()) {
		if (0 == (m_vPepTolWnds[nWndId].m_strPepTolBaseType.compare("%")))
			lfTolBase = lfMZ * m_vPepTolWnds[nWndId].m_lfPepTolBase * 0.01
					* nCharge;
		else if (0
				== (m_vPepTolWnds[nWndId].m_strPepTolBaseType.compare("mmu")))
			lfTolBase = m_vPepTolWnds[nWndId].m_lfPepTolBase * 0.001;
		else if (0
				== (m_vPepTolWnds[nWndId].m_strPepTolBaseType.compare("ppm")))
			lfTolBase = lfMZ * m_vPepTolWnds[nWndId].m_lfPepTolBase * 0.000001
					* nCharge;
		else if (0 == (m_vPepTolWnds[nWndId].m_strPepTolBaseType.compare("Da")))
			lfTolBase = m_vPepTolWnds[nWndId].m_lfPepTolBase;
	}
	return lfTolBase;
}

double CCondition::GetPepTol(double lfMZ, int nCharge) const {
	return GetPepTol(0, lfMZ, nCharge);
}

double CCondition::GetPepTol(int nWndId, double lfMZ, int nCharge) const {
	double lfTol = m_lfMinPepTol;
	if (nWndId >= 0 && nWndId < m_vPepTolWnds.size()) {
		if (0 == (m_vPepTolWnds[nWndId].m_strPepTolType.compare("%")))
			lfTol = lfMZ * m_vPepTolWnds[nWndId].m_lfPepTol * 0.01 * nCharge;
		else if (0 == (m_vPepTolWnds[nWndId].m_strPepTolType.compare("mmu")))
			lfTol = m_vPepTolWnds[nWndId].m_lfPepTol * 0.001;
		else if (0 == (m_vPepTolWnds[nWndId].m_strPepTolType.compare("ppm")))
			lfTol = lfMZ * m_vPepTolWnds[nWndId].m_lfPepTol * 0.000001
					* nCharge;
		else if (0 == (m_vPepTolWnds[nWndId].m_strPepTolType.compare("Da")))
			lfTol = m_vPepTolWnds[nWndId].m_lfPepTol;
	}
	return lfTol;
}

double CCondition::GetMinPepTol(double lfMZ, int nCharge) const {
	double lfTol = m_lfMinPepTol;

	if (0 < lfTol) {
		if (0 == (m_strMinPepTolType.compare("%")))
			lfTol = lfMZ * m_lfMinPepTol * 0.01 * nCharge;
		else if (0 == (m_strMinPepTolType.compare("mmu")))
			lfTol = m_lfMinPepTol * 0.001;
		else if (0 == (m_strMinPepTolType.compare("ppm")))
			lfTol = lfMZ * m_lfMinPepTol * 0.000001 * nCharge;
		else if (0 == (m_strMinPepTolType.compare("Da")))
			lfTol = m_lfMinPepTol;
	} else {
		lfTol = MIN_PEPTIDE_TOL;
	}

	return lfTol;
}

void CCondition::GetMHMassBorder(double lfMH, int nCharge, double &lfLower,
		double &lfUpper) const {
	//double lfExpMass = GetExpMass(lfMH, nCharge);

	double lfMZ = GetMZ(lfMH, nCharge);

	double lfTol = GetPepTol(lfMZ, nCharge);

	double lfTolBase = GetPepTolBase(lfMZ, nCharge);

	double lfMinTol = GetMinPepTol(lfMZ, nCharge);

	if (lfTol < lfMinTol) {
		lfLower = lfMH - lfMinTol + lfTolBase;
		lfUpper = lfMH + lfMinTol + lfTolBase;
	} else {
		lfLower = lfMH - lfTol + lfTolBase;
		lfUpper = lfMH + lfTol + lfTolBase;
	}
}

void CCondition::GetMHMassBorder(int nWnd, double lfMH, int nCharge,
		double &lfLower, double &lfUpper) const {
	// double lfExpMass = GetExpMass(lfMH, nCharge);

	double lfMZ = GetMZ(lfMH, nCharge);

	double lfTol = GetPepTol(nWnd, lfMZ, nCharge);

	double lfTolBase = GetPepTolBase(nWnd, lfMZ, nCharge);

	double lfMinTol = GetMinPepTol(lfMZ, nCharge);

	if (lfTol < lfMinTol) {
		lfLower = lfMH - lfMinTol + lfTolBase;
		lfUpper = lfMH + lfMinTol + lfTolBase;
	} else {
		lfLower = lfMH - lfTol + lfTolBase;
		lfUpper = lfMH + lfTol + lfTolBase;
	}
}

void CCondition::GetMassBorder(double lfMH, int nCharge, double &lfLower,
		double &lfUpper) const {
	// modify by emily
	double lfExpMass = GetExpMass(lfMH, nCharge);

	double lfPepMass = GetNeutraMass(lfExpMass);

	double lfMZ = GetMZ(lfMH, nCharge);

	double lfTol = GetPepTol(lfMZ, nCharge);

	double lfTolBase = GetPepTolBase(lfMZ, nCharge);

	double lfMinTol = GetMinPepTol(lfMZ, nCharge);

	if (lfTol < lfMinTol) {
		lfLower = lfPepMass - lfMinTol + lfTolBase;
		lfUpper = lfPepMass + lfMinTol + lfTolBase;
	} else {
		lfLower = lfPepMass - lfTol + lfTolBase;
		lfUpper = lfPepMass + lfTol + lfTolBase;
	}
}

void CCondition::GetMassBorder(int nWnd, double lfMH, int nCharge,
		double &lfLower, double &lfUpper) const {
	// modify by emily
	double lfExpMass = GetExpMass(lfMH, nCharge);

	double lfPepMass = GetNeutraMass(lfExpMass);

	double lfMZ = GetMZ(lfMH, nCharge);

	double lfTol = GetPepTol(nWnd, lfMZ, nCharge);

	double lfTolBase = GetPepTolBase(nWnd, lfMZ, nCharge);

	double lfMinTol = GetMinPepTol(lfMZ, nCharge);

	if (lfTol < lfMinTol) {
		lfLower = lfPepMass - lfMinTol + lfTolBase;
		lfUpper = lfPepMass + lfMinTol + lfTolBase;
	} else {
		lfLower = lfPepMass - lfTol + lfTolBase;
		lfUpper = lfPepMass + lfTol + lfTolBase;
	}
}

double CCondition::GetLeastNegativeMod() const {
	double lfNegativeMod = 0;

	for (size_t i = 0; i < m_vSelectedFixMod.size(); ++i) {
		if (m_bFragmentMono) {
			if (0 > m_vSelectedFixMod[i].m_lfMonoMass_dif)
				if (lfNegativeMod > m_vSelectedFixMod[i].m_lfMonoMass_dif)
					lfNegativeMod = m_vSelectedFixMod[i].m_lfMonoMass_dif;
		} else {
			if (0 > m_vSelectedFixMod[i].m_lfAvrgMass_dif)
				if (lfNegativeMod > m_vSelectedFixMod[i].m_lfAvrgMass_dif)
					lfNegativeMod = m_vSelectedFixMod[i].m_lfAvrgMass_dif;
		}
	}

	for (size_t i = 0; i < m_vSelectedVarMod.size(); ++i) {
		if (m_bFragmentMono) {
			if (0 > m_vSelectedVarMod[i].m_lfMonoMass_dif)
				if (lfNegativeMod > m_vSelectedVarMod[i].m_lfMonoMass_dif)
					lfNegativeMod = m_vSelectedVarMod[i].m_lfMonoMass_dif;
		} else {
			if (0 > m_vSelectedVarMod[i].m_lfAvrgMass_dif)
				if (lfNegativeMod > m_vSelectedVarMod[i].m_lfAvrgMass_dif)
					lfNegativeMod = m_vSelectedVarMod[i].m_lfAvrgMass_dif;
		}
	}

	return lfNegativeMod;
}

double CCondition::GetMaxMassDif() const {
	double lfMaxModMass = 0;

	for (size_t i = 0; i < m_vSelectedFixMod.size(); ++i) {
		if (m_bFragmentMono) {
			if (lfMaxModMass < m_vSelectedFixMod[i].m_lfMonoMass_dif)
				lfMaxModMass = m_vSelectedFixMod[i].m_lfMonoMass_dif;
		} else {
			if (lfMaxModMass < m_vSelectedFixMod[i].m_lfAvrgMass_dif)
				lfMaxModMass = m_vSelectedFixMod[i].m_lfAvrgMass_dif;
		}
	}

	for (size_t i = 0; i < m_vSelectedVarMod.size(); ++i) {
		if (m_bFragmentMono) {
			if (lfMaxModMass < m_vSelectedVarMod[i].m_lfMonoMass_dif)
				lfMaxModMass = m_vSelectedVarMod[i].m_lfMonoMass_dif;
		} else {
			if (lfMaxModMass < m_vSelectedVarMod[i].m_lfAvrgMass_dif)
				lfMaxModMass = m_vSelectedVarMod[i].m_lfAvrgMass_dif;
		}
	}

	return lfMaxModMass;
}

string CCondition::GetInstrumentType(void) {
	if (INSTRUMENT_ESI_QUAD_TOF == m_eInstrumentType)
		return "INSTRUMENT_ESI_QUAD_TOF";
	else if (INSTRUMENT_ESI_TRAP == m_eInstrumentType)
		return "INSTRUMENT_ESI_TRAP";
	else if (INSTRUMENT_ESI_QUAD == m_eInstrumentType)
		return "INSTRUMENT_ESI_QUAD";
	else if (INSTRUMENT_ESI_FTICR == m_eInstrumentType)
		return "INSTRUMENT_ESI_FTICR";
	else if (INSTRUMENT_ESI_4SECTOR == m_eInstrumentType)
		return "INSTRUMENT_ESI_4SECTOR";
	else if (INSTRUMENT_MALDI_QUAD_TOF == m_eInstrumentType)
		return "INSTRUMENT_MALDI_QUAD_TOF";
	else if (INSTRUMENT_MALDI_TOF_PSD == m_eInstrumentType)
		return "INSTRUMENT_MALDI_TOF_PSD";
	else if (INSTRUMENT_MALDI_TOF_TOF == m_eInstrumentType)
		return "INSTRUMENT_MALDI_TOF_TOF";
	else if (INSTRUMENT_MALDI_QIT_TOF == m_eInstrumentType)
		return "INSTRUMENT_MALDI_QIT_TOF";
	else if (INSTRUMENT_FTMS_ECD == m_eInstrumentType)
		return "INSTRUMENT_FTMS_ECD";
	else
		return "INSTRUMENT_UNKNOWN";
}

}
