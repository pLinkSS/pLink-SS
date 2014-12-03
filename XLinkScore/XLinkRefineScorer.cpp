#include <iostream>
#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "XLinkRefineScorer.h"
using namespace std;


#define MAX_INNER_ION_LENGTH 4

#define PEAKMAX 4000





namespace proteomics_sdk
{

void CXLinkRefineScorer::Initialize(const CCondition & condition)
{

	_SetCondition(condition);

	m_pbMatched = new char[m_nMaxPeaksNum];

	size_t tNLSize = 0;
	for (size_t i = 0; i < m_Condition.m_vSelectedVarMod.size(); ++i)
	{
		tNLSize +=
				m_Condition.m_vSelectedVarMod[i].m_vlfAvrgNeutralLoss_dif.size();
	}
	m_ionmz.clear();
	m_ionmz.resize(
			10 * MAX_IONTYPE_NUM * 2 * MAX_PEPTIDE_LENGTH * (tNLSize + 1));

	memset(m_FullInterContiWnd, 0, sizeof(int) * MAX_IONTYPE_NUM);
	for (size_t i = 0; i < m_Condition.m_vIonTypes.size(); ++i)
	{
		const CIonType & IonType = m_Condition.m_vIonTypes[i];
		if (IonType.nInterContinousWnd1 >= 0)
		{
			if (IonType.cType == 2)
				m_FullInterContiWnd[IonType.nInterContinousWnd1] += l_win;
			else
				m_FullInterContiWnd[IonType.nInterContinousWnd1]++;
		}
		if (IonType.nInterContinousWnd2 >= 0)
		{
			if (IonType.cType == 2)
				m_FullInterContiWnd[IonType.nInterContinousWnd2] += l_win;
			else
				m_FullInterContiWnd[IonType.nInterContinousWnd2]++;
		}
	}


}
void CXLinkRefineScorer::Close(void)
{
	if (m_pbMatched)
		delete[] m_pbMatched;
	m_pbMatched = NULL;
	/*
	 if(m_pPeakWeight)
	 delete [] m_pPeakWeight;
	 m_pPeakWeight = NULL;
	 */

}
void CXLinkRefineScorer::_SetCondition(const CCondition & cond)
{
	m_Condition = cond;
	if (0 == m_Condition.m_nKSDP_L)
		m_Condition.m_nKSDP_L = 5;
	if (0 == m_Condition.m_lfKSDP_G)
		m_Condition.m_lfKSDP_G = 0.9f;
	if (0 == m_Condition.m_lfKSDP_A)
		m_Condition.m_lfKSDP_A = 0.5f;
	/*Ô¤¼ÆËã³£Êý*/
	l_win = m_Condition.m_nKSDP_L;
	alfa = m_Condition.m_lfKSDP_A;
	gamma = m_Condition.m_lfKSDP_G;

	for (int i = 0; i < (MAX_PEPTIDE_LENGTH << 2); ++i)
	{
		exp_table[i] = exp((-1) * gamma * i);
	}

	l1 = (int) floor((l_win - 1) / 2.0);
	l2 = (int) ceil((l_win - 1) / 2.0);




	if (0 == (m_Condition.m_strFragmentTolType.compare("%")))
	{
		m_cTolType = 0;
		m_lfTolMultiplier = m_Condition.m_lfFragmentTol * 0.01;
	}
	else if (0 == (m_Condition.m_strFragmentTolType.compare("mmu")))
	{
		m_cTolType = 1;
		m_lfTolMultiplier = m_Condition.m_lfFragmentTol * MMU_MULTIPLIER;
	}
	else if (0 == (m_Condition.m_strFragmentTolType.compare("ppm")))
	{
		m_cTolType = 2;
		m_lfTolMultiplier = m_Condition.m_lfFragmentTol * 0.000001;
	}
	else if (0 == (m_Condition.m_strFragmentTolType.compare("Da")))
	{
		m_cTolType = 3;
		m_lfTolMultiplier = m_Condition.m_lfFragmentTol * MZMULTIPLIER;
	}
	else
	{
		m_cTolType = 3;
		m_lfTolMultiplier = m_Condition.m_lfFragmentTol * MZMULTIPLIER;
	}


	if (0 == (m_Condition.m_strFragmentTolBaseType.compare("%")))
	{
		m_cTolBaseType = 0;
		m_lfTolBaseMultiplier = m_Condition.m_lfFragmentTolBase * 0.01;
	}
	else if (0 == (m_Condition.m_strFragmentTolBaseType.compare("mmu")))
	{
		m_cTolBaseType = 1;
		m_lfTolBaseMultiplier = m_Condition.m_lfFragmentTolBase
				* MMU_MULTIPLIER;
	}
	else if (0 == (m_Condition.m_strFragmentTolBaseType.compare("ppm")))
	{
		m_cTolBaseType = 2;
		m_lfTolBaseMultiplier = m_Condition.m_lfFragmentTolBase * 0.000001;
	}
	else if (0 == (m_Condition.m_strFragmentTolBaseType.compare("Da")))
	{
		m_cTolBaseType = 3;
		m_lfTolBaseMultiplier = m_Condition.m_lfFragmentTolBase * MZMULTIPLIER;
	}
	else
	{
		m_cTolBaseType = 3;
		m_lfTolBaseMultiplier = m_Condition.m_lfFragmentTolBase * MZMULTIPLIER;
	}


	for (size_t i = 0; i < m_Condition.m_vSelectedFixMod.size(); ++i)
	{
		m_Condition.m_vSelectedVarMod.push_back(
				m_Condition.m_vSelectedFixMod[i]);
	}


	CAAConf aa(m_Condition.m_strAAListPath);
	CMapAAMass mapAAMass = aa.GetMapAAMass();
	for (int i = 0; i < 26; ++i)
	{
		if (m_Condition.m_bFragmentMono)
			m_nAAMass[i] = (int) (mapAAMass.m_mapAAMass[i + 'A'].m_lfMonoMass
					* MZMULTIPLIER);
		else
			m_nAAMass[i] = (int) (mapAAMass.m_mapAAMass[i + 'A'].m_lfAvrgMass
					* MZMULTIPLIER);
	}
}
void CXLinkRefineScorer::SetSpectrum(CSpectrum & spec)
{

	m_pSpectrum = &spec;
	/*
	 if(m_pPeakWeight && (m_nMaxPeaksNum < m_pSpectrum->m_tPeaksNum))
	 {
	 delete [] m_pPeakWeight;
	 m_pPeakWeight = NULL;
	 m_pPeakWeight = new double[m_pSpectrum->m_tPeaksNum];
	 }
	 memset(m_pPeakWeight, 1, sizeof(double) * m_pSpectrum->m_tPeaksNum);
	 double lfBasePeak = 0.0;
	 for(size_t i = 0;i < m_pSpectrum->m_tPeaksNum ; ++ i)
	 {
	 if(m_pSpectrum->m_pPeaks[i].lfIntensity > lfBasePeak)
	 lfBasePeak = m_pSpectrum->m_pPeaks[i].lfIntensity;
	 }
	 for(size_t i = 0;i < m_pSpectrum->m_tPeaksNum ; ++ i)
	 {
	 m_pPeakWeight[i] = m_pSpectrum->m_pPeaks[i].lfIntensity/lfBasePeak;
	 }
	 */
}
void CXLinkRefineScorer::SetPeptide(CXLinkPepResult & pep)
{
	m_pPeptide = &pep;
	m_bComputeMz = false;
}

void CXLinkRefineScorer::_PrintMZ()
{
	if (strcmp(m_pPeptide->m_AlphaPeptide.m_szSequence, "SYIGCFAKQSYSD"))
		return;





	cout << "TypePos1\tPos2\tCharge\tMz" << endl;
	for (int i = 0; i < m_tionmzSize; i++)
	{
		CIonType IonType = m_Condition.m_vIonTypes[m_ionmz[i].nIonTypeOrder];
		cout << IonType.cSymbol << m_ionmz[i].nPepPosOrder1 << "\t"
				<< m_ionmz[i].nPepPosOrder2 << "\t" << IonType.nCharge << "\t"
				<< m_ionmz[i].nMz << endl;
	}
}

double CXLinkRefineScorer::Score()
{

	if (0 == m_pSpectrum->m_tPeaksNum)
	{
		return 0.0;
	}
	if (!m_bComputeMz)
	{
		_ComputeMZ();
		m_bComputeMz = true;
	}



	_Match();
	double Fen = 0.0;
	int Isn = 0;
	double s = 0;
	int win = 0;
	int k = 0;
	double Kp1 = 0, Kp2 = 0;
	double Ki = 0;
	double lfwin = 0.0;
	size_t Pl1 = m_pPeptide->m_AlphaPeptide.m_tLength;

	/*
	 for (size_t i=0;i<m_Condition.m_vIonTypes.size();++i)
	 {
	 if((m_pSpectrum->m_nCharge > 1 && m_Condition.m_vIonTypes[i].nCharge >= m_pSpectrum->m_nCharge)
	 || (m_pSpectrum->m_nCharge == 1 && m_Condition.m_vIonTypes[i].nCharge > 1))
	 continue;

	 if(!m_Condition.m_vIonTypes[i].bIntraContinous)
	 continue;

	 ++ionstypenum;

	 win=0;
	 int nUB = 0,nLB = 0;
	 for (int j=0;(j<=l2)&&(j< (int)m_pPeptide->m_AlphaPeptide.m_tLength -1);++j)
	 {
	 nUB ++ ;
	 if (m_iontag_pep[i][j]>=0)
	 ++win;
	 }
	 kernel_pep1+=exp_table[nUB - win];

	 for (size_t j=1;j< m_pPeptide->m_AlphaPeptide.m_tLength -1; ++j)
	 {
	 k=j+l2;

	 if (k< (int)m_pPeptide->m_AlphaPeptide.m_tLength -1)
	 {
	 nUB = k;
	 if (m_iontag_pep[i][k]>=0)
	 ++win;
	 }
	 else
	 {
	 nUB = (int)m_pPeptide->m_AlphaPeptide.m_tLength -2;
	 }

	 k=j-l1-1;
	 if (k>=0)
	 {
	 nLB = k;
	 if (m_iontag_pep[i][k]>=0)
	 --win;
	 }
	 else
	 {
	 nLB = -1;
	 }
	 kernel_pep1+=exp_table[nUB - nLB - win];
	 }

	 if(m_pPeptide->m_bPair)
	 {
	 win=0;
	 nLB = 0;
	 nUB = 0;
	 for (int j=0;(j<=l2)&&(j< (int)m_pPeptide->m_BetaPeptide.m_tLength -1);++j)
	 {
	 nUB ++ ;
	 if (m_iontag_pep[i][tPepLen1 + j]>=0)
	 ++win;
	 }
	 kernel_pep2+=exp_table[nUB - win];

	 for (size_t j=1;j< m_pPeptide->m_BetaPeptide.m_tLength -1; ++j)
	 {
	 k=j+l2;
	 if (k< (int)m_pPeptide->m_BetaPeptide.m_tLength -1)
	 {
	 nUB = k;
	 if (m_iontag_pep[i][tPepLen1 + k]>=0)
	 ++win;
	 }
	 else
	 {
	 nUB = (int)m_pPeptide->m_BetaPeptide.m_tLength -2;
	 }
	 k=j-l1-1;
	 if (k>=0)
	 {
	 nLB = k;
	 if (m_iontag_pep[i][tPepLen1 + k]>=0)
	 --win;
	 }
	 else
	 {
	 nLB = -1;
	 }
	 kernel_pep2+=exp_table[nUB - nLB - win];
	 }
	 }
	 }
	 */

	int Iic = 0;
	for (size_t i = 0; i < m_Condition.m_vIonTypes.size(); ++i)
	{
		if ((m_pSpectrum->m_nCharge > 1
				&& m_Condition.m_vIonTypes[i].nCharge >= m_pSpectrum->m_nCharge)
				|| (m_pSpectrum->m_nCharge == 1
						&& m_Condition.m_vIonTypes[i].nCharge > 1))
			continue;

		if (!m_Condition.m_vIonTypes[i].bIntraContinous)
			continue;

		++Isn;

		const CIonType & It = m_Condition.m_vIonTypes[i];

		int nBegin = 0, nEnd = 0;
		if (It.nContainLinker == 0 || m_pPeptide->m_XLink.m_eXLinkType < 3)
		{
			nBegin = 0;
			nEnd = int(m_pPeptide->m_AlphaPeptide.m_tLength - 1);
		}
		else
		{
			if (It.cType == 0)
			{
				if (It.nContainLinker == 1)
				{

					nBegin = 0;
					nEnd = int(m_pPeptide->m_XLink.m_tAlphaSite);
				}
				else
				{

					nBegin = int(m_pPeptide->m_XLink.m_tAlphaSite);
					nEnd = int(m_pPeptide->m_AlphaPeptide.m_tLength - 1);
				}
			}
			else
			{
				if (It.nContainLinker == 1)
				{

					nBegin = 0;
					nEnd = int(
							m_pPeptide->m_AlphaPeptide.m_tLength - 1
									- m_pPeptide->m_XLink.m_tAlphaSite);
				}
				else
				{

					nBegin = int(
							m_pPeptide->m_AlphaPeptide.m_tLength - 1
									- m_pPeptide->m_XLink.m_tAlphaSite);
					nEnd = int(m_pPeptide->m_AlphaPeptide.m_tLength - 1);
				}
			}
		}
		nBegin = 0;
		nEnd = int(m_pPeptide->m_AlphaPeptide.m_tLength - 1);

		win = 0;
		Iic += nEnd - nBegin;

		for (int j = nBegin; (j <= nBegin + l2) && (j < nEnd); ++j)
			if (m_iontag_pep[i][j] < 0)
				++win;
		Kp1 += exp_table[win];

		for (int j = nBegin + 1; j < nEnd; ++j)
		{
			k = j + l2;
			if (k < nEnd)
				if (m_iontag_pep[i][k] < 0)
					++win;
			k = j - l1 - 1;
			if (k >= nBegin)
				if (m_iontag_pep[i][k] < 0)
					--win;
			Kp1 += exp_table[win];
		}

		if (m_pPeptide->m_bPair)
		{
			nBegin = 0;
			nEnd = 0;
			if (It.nContainLinker == 0)
			{
				nBegin = 0;
				nEnd = int(m_pPeptide->m_BetaPeptide.m_tLength - 1);
			}
			else
			{
				if (It.cType == 0)
				{
					if (It.nContainLinker == 1)
					{

						nBegin = 0;
						nEnd = int(m_pPeptide->m_XLink.m_tBetaSite);
					}
					else
					{

						nBegin = int(m_pPeptide->m_XLink.m_tBetaSite);
						nEnd = int(m_pPeptide->m_BetaPeptide.m_tLength - 1);
					}
				}
				else
				{
					if (It.nContainLinker == 1)
					{

						nBegin = 0;
						nEnd = int(
								m_pPeptide->m_BetaPeptide.m_tLength - 1
										- m_pPeptide->m_XLink.m_tBetaSite);
					}
					else
					{

						nBegin = int(
								m_pPeptide->m_BetaPeptide.m_tLength - 1
										- m_pPeptide->m_XLink.m_tBetaSite);
						nEnd = int(m_pPeptide->m_BetaPeptide.m_tLength - 1);
					}
				}
			}
			nBegin = 0;
			nEnd = int(m_pPeptide->m_BetaPeptide.m_tLength - 1);

			Iic += (nEnd - nBegin);
			win = 0;
			for (int j = nBegin; (j <= nBegin + l2) && (j < nEnd); ++j)
				if (m_iontag_pep[i][Pl1 + j] < 0)
					++win;
			Kp2 += exp_table[win];

			for (int j = nBegin + 1; j < nEnd; ++j)
			{
				k = j + l2;
				if (k < nEnd)
					if (m_iontag_pep[i][Pl1 + k] < 0)
						++win;
				k = j - l1 - 1;
				if (k >= nBegin)
					if (m_iontag_pep[i][Pl1 + k] < 0)
						--win;
				Kp2 += exp_table[win];
			}
		}
	}


	int Imc = 0;
	for (size_t i = 0; i < m_Condition.m_tInterContinuousWndNum; ++i)
	{
		for (size_t j = 0; j < m_pPeptide->m_AlphaPeptide.m_tLength; ++j)
		{

			Imc = 10 - m_InterContiWnd[i][j];
			if (Imc < 0)
				Imc = 0;

			Ki += exp_table[Imc];


		}

		if (m_pPeptide->m_bPair)
		{
			for (size_t j = 0; j < m_pPeptide->m_BetaPeptide.m_tLength; ++j)
			{

				Imc = 10 - m_InterContiWnd[i][Pl1 + j];
				if (Imc < 0)
					Imc = 0;

				Ki += exp_table[Imc];

			}
		}

	}

	s = 0;


	for (size_t i = 0; i < m_pSpectrum->m_tPeaksNum; ++i)
	{
		if (m_pbMatched[i])
		{

			s += m_pSpectrum->m_pPeaks[i].lfIntensity;
		}
		/*
		 else if(m_pbMatched[i] == 2)
		 {



		 double weight = 1.5;

		 s+= weight*m_pSpectrum->m_pPeaks[i].lfIntensity;
		 }
		 */
	}

	s = s / m_pSpectrum->m_lfSqtMaxInten * 100;

	int Pl = 0;
	Pl += m_pPeptide->m_AlphaPeptide.m_tLength - 1;
	if (m_pPeptide->m_bPair)
	{
		Pl += m_pPeptide->m_BetaPeptide.m_tLength - 1;
	}

	double Mk = 0.0;

	Mk = 0;

#ifndef USEAVERAGE
	if (Kp1 != 0 && Kp2 != 0)
	{
		Mk = 4 / (1 / Kp1 + 1 / Kp2);
	}
	else
	{
		Mk = Kp1 + Kp2;
	}
#else
	Mk = Kp1 + Kp2;
#endif


	Fen = s * pow(Mk + Ki, (double) alfa)
			/ (Iic/*ionstypenum*(nPepLen)*/);

	/*
	 fprintf(m_fp,"%s(%d)-%s(%d)	%f\n",
	 m_pPeptide->m_AlphaPeptide.m_szSequence,
	 m_pPeptide->m_XLink.m_tAlphaSite,
	 m_pPeptide->m_BetaPeptide.m_szSequence,
	 m_pPeptide->m_XLink.m_tBetaSite,
	 lfScore);
	 */

	/*
	 cout << "kernel_pep1 = " << kernel_pep1 << endl
	 <<  "kernel_pep2 = " << kernel_pep2 << endl
	 << "mean_kernel = " << lfMeanKernel << endl
	 << "kernel_inter = " << kernel_inter << endl
	 << "s = " << s << endl;
	 */

	_CalMatchInfo();

#ifdef _DEBUG2
	char szBuf[1024];

	cout << "m_pSpectrum->m_lfSqtMaxInten=" << m_pSpectrum->m_lfSqtMaxInten << endl;

	m_strDebug = "";
	sprintf(szBuf,"%s %s",m_pSpectrum->m_strFilePath.c_str(),m_pPeptide->m_AlphaPeptide.m_szSequence);
	m_strDebug += szBuf;

	if(m_pPeptide->m_bPair)
	{
		sprintf(szBuf,"-%s",m_pPeptide->m_BetaPeptide.m_szSequence);
		m_strDebug += szBuf;
	}

	sprintf(szBuf,"\nspectra matched : \n");
	m_strDebug += szBuf;

	for(int i=0;i<m_pSpectrum->m_tPeaksNum;++i)
	{
		if ( m_pbMatched[i] )
		{
			sprintf(szBuf,"%.5f	%.5f\n",m_pSpectrum->m_pPeaks[i].lfMz,m_pSpectrum->m_pPeaks[i].lfIntensity);
			m_strDebug += szBuf;
		}
	}

	sprintf(szBuf,"pep matched : \n");
	m_strDebug += szBuf;

	sprintf(szBuf,"=================================\n");
	m_strDebug += szBuf;

	for(int i=0;i<m_Condition.m_vIonTypes.size();++i)
	{
		for(int j=0;j<m_pPeptide->m_AlphaPeptide.m_tLength - 1;++j)
		{
			sprintf(szBuf,"%d	",m_iontag_pep[i][j]);

			m_strDebug += szBuf;
		}
		sprintf(szBuf,"\n");
		m_strDebug += szBuf;
	}
	sprintf(szBuf,"\n=================================\n");
	m_strDebug += szBuf;

	if(m_pPeptide->m_bPair)
	{
		sprintf(szBuf,"=================================\n");
		m_strDebug += szBuf;

		for(int i=0;i<m_Condition.m_vIonTypes.size();++i)
		{
			for(int j=0;j<m_pPeptide->m_BetaPeptide.m_tLength - 1;++j)
			{
				sprintf(szBuf,"%d	",m_iontag_pep[i][m_pPeptide->m_AlphaPeptide.m_tLength+j]);

				m_strDebug += szBuf;
			}
			sprintf(szBuf,"\n");
			m_strDebug += szBuf;
		}
		sprintf(szBuf,"\n=================================\n");
		m_strDebug += szBuf;
	}

	sprintf(szBuf,"kernel_pep1 = %.5f	kernel_pep2 = %.5f	s = %.5f	lfScore = %.5f\n",Kp1,Kp2,s,Fen);
	m_strDebug += szBuf;

	_CalMatchInfo();

	m_strDebug += string("Pep1 confidence : \n");
	for(int i=0;i<m_pPeptide->m_AlphaPeptide.m_tLength;++i)
	{
		sprintf(szBuf,"%d	",m_pPeptide->m_stMatchInfo.aPepConf[i]);
		m_strDebug += szBuf;
	}
	m_strDebug += "\n";

	if(m_pPeptide->m_bPair)
	{
		m_strDebug += "Pep2 confidence : \n";
		for(int i=0;i<m_pPeptide->m_BetaPeptide.m_tLength;++i)
		{
			sprintf(szBuf,"%d	",m_pPeptide->m_stMatchInfo.aPepConf[m_pPeptide->m_AlphaPeptide.m_tLength + i]);
			m_strDebug += szBuf;
		}
		m_strDebug += "\n";
	}

	sprintf(szBuf,"matched - unmatched : %.5f - %.5f\n",m_pPeptide->m_stMatchInfo.lfMatchedSpecInt,m_pPeptide->m_stMatchInfo.lfUnMatchedSpecInt);
	m_strDebug += szBuf;
	cout << m_strDebug << endl;
	char tmpch;
	cin >> tmpch;
#endif

	return Fen;

}

double CXLinkRefineScorer::i2Score()
{
	if (0 == m_pSpectrum->m_tPeaksNum)
	{
		return 0.0;
	}
	if (!m_bComputeMz)
	{
		_ComputeMZ();
		m_bComputeMz = true;
	}
	if (!m_tionmzSize)
		return 0;

	_Match();

	double s = 0;
	for (size_t i = 0; i < m_pSpectrum->m_tPeaksNum; ++i)
	{
		if (m_pbMatched[i])
			s += m_pSpectrum->m_pPeaks[i].lfIntensity;
	}

	double nTotal = 0;
	for (size_t i = 0; i < m_pSpectrum->m_tPeaksNum; ++i)
	{
		nTotal += m_pSpectrum->m_pPeaks[i].lfIntensity;
	}
	if (nTotal <= 0)
	{
		cout << "Total spectrum intensity not legal.";
		return 0;
	}

	double peakratio = s / nTotal;

	double nTotalIonMatched = 0;

	for (int i = 0; i < m_tionmzSize; i++)
		if (m_ionmz[i].bMatched)
			nTotalIonMatched++;

	double ionratio = nTotalIonMatched / m_tionmzSize;

	return peakratio * ionratio;
}

double CXLinkRefineScorer::i2HyperGeoMetric()
{
	if (0 == m_pSpectrum->m_tPeaksNum)
	{
		return 0.0;
	}
	if (!m_bComputeMz)
	{
		_ComputeMZ();
		m_bComputeMz = true;
	}
	if (!m_tionmzSize)
		return 0;

	_Match();

	double s = 0;
	for (size_t i = 0; i < m_pSpectrum->m_tPeaksNum; ++i)
	{
		if (m_pbMatched[i])
			s += m_pSpectrum->m_pPeaks[i].lfIntensity;
	}

	double nTotal = 0;
	for (size_t i = 0; i < m_pSpectrum->m_tPeaksNum; ++i)
	{
		nTotal += m_pSpectrum->m_pPeaks[i].lfIntensity;
	}
	if (nTotal <= 0)
	{
		cout << "Total spectrum intensity not legal.";
		return 0;
	}

	double peakratio = s / nTotal;

	double nTotalIonMatched = 0;
	for (int i = 0; i < m_tionmzSize; i++)
		if (m_ionmz[i].bMatched)
			nTotalIonMatched++;

	double ionratio = CalcuCombine(m_pSpectrum->m_tPeaksNum, nTotalIonMatched)
			* CalcuCombine(PEAKMAX - m_pSpectrum->m_tPeaksNum,
					m_tionmzSize - nTotalIonMatched)
			/ (double) CalcuCombine(PEAKMAX, m_tionmzSize);

	return peakratio * ionratio;
}

double CXLinkRefineScorer::CalcuCombine(int a, int b)
{
	if (a < b || a <= 0 || b <= 0)
		return 0;

	double result = 1;

	for (int i = 0; i < b; i++)
	{
		result *= (a - i);
		result /= (i + 1);
	}

	return result;
}

void CXLinkRefineScorer::_GetMassBorder(int nMz, int nChg, int &nMin, int &nMax)
{
	double lfTol = m_lfTolMultiplier;
	double lfTolBase = m_lfTolBaseMultiplier;
	if ((m_cTolType & 0x01) == 0)
	{
		lfTol *= nMz;
	}
	if ((m_cTolBaseType & 0x01) == 0)
	{
		lfTolBase *= nMz;
	}
	nMin = int(nMz - lfTol + lfTolBase);
	nMax = int(nMz + lfTol + lfTolBase);
}

void CXLinkRefineScorer::_Match()
{
	const vector<CIonType> & vIonTypes = m_Condition.m_vIonTypes;
	memset(&m_iontag_pep[0], -1, sizeof(m_iontag_pep[0]) * vIonTypes.size());

	memset(&m_InterContiWnd[0], 0,
			sizeof(m_InterContiWnd[0]) * vIonTypes.size());

	if (m_pbMatched && (m_nMaxPeaksNum < m_pSpectrum->m_tPeaksNum))
	{
		delete[] m_pbMatched;
		m_pbMatched = NULL;
		m_nMaxPeaksNum = m_pSpectrum->m_tPeaksNum;
		m_pbMatched = new char[m_pSpectrum->m_tPeaksNum];
	}

	memset(m_pbMatched, 0, sizeof(char) * m_pSpectrum->m_tPeaksNum);

	int nPepLen1 = m_pPeptide->m_AlphaPeptide.m_tLength;

	for (size_t pos_pep = 0; pos_pep < m_tionmzSize; ++pos_pep)
	{
		if ((m_pSpectrum->m_nCharge == 1
				&& vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge
						> m_pSpectrum->m_nCharge)
				|| (m_pSpectrum->m_nCharge > 1
						&& vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge
								>= m_pSpectrum->m_nCharge))
		{
			continue;
		}

		int nMin, nMax;
		_GetMassBorder(m_ionmz[pos_pep].nMz,
				vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge, nMin, nMax);
		int nTemp = nMin / MZMULTIPLIER;

		if (nTemp >= (int) m_pSpectrum->m_vHash.size())
			continue;

		if (nTemp < 0)
			continue;

		size_t pos_spec = m_pSpectrum->m_vHash[nTemp];

		while (pos_spec < m_pSpectrum->m_tPeaksNum
				&& m_pSpectrum->m_pPeaks[pos_spec].nMz < nMin)
		{
			++pos_spec;

		}
		while (pos_spec < m_pSpectrum->m_tPeaksNum
				&& m_pSpectrum->m_pPeaks[pos_spec].nMz <= nMax)
		{
			if (m_pSpectrum->m_pPeaks[pos_spec].nCharge == 0
					|| m_pSpectrum->m_pPeaks[pos_spec].nCharge
							== vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge)
			{
				/*
				 if(m_ionmz[pos_pep].nIonTypeOrder == 0)
				 {
				 double u = -4.1249e-005;
				 double v = -0.0362;
				 double a =  -0.2363;
				 double b = 0.2880;

				 double lfExpMz = m_pSpectrum->m_pPeaks[pos_spec].lfMz;
				 double lfCalMz = (m_ionmz[pos_pep].nMz+0.0)/MZMULTIPLIER;
				 double lfTol = lfExpMz - lfCalMz;
				 double lfIntensity = m_pSpectrum->m_pPeaks[pos_spec].lfIntensity;
				 double zvalue = (lfTol - (u*lfExpMz+v))/(b*pow(lfIntensity,a));
				 if(zvalue > 3 || zvalue < -3)
				 {
				 ++pos_spec;
				 continue;
				 }

				 }
				 */

				if (m_ionmz[pos_pep].bContainLinker)
					m_pbMatched[pos_spec] = 2;
				else
					m_pbMatched[pos_spec] = 1;

				m_ionmz[pos_pep].bMatched = 1;

				if (m_ionmz[pos_pep].bAddToTag)
				{
					m_iontag_pep[m_ionmz[pos_pep].nIonTypeOrder][m_ionmz[pos_pep].nPepPosOrder1] =
							pos_spec;

				}

				if (vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nInterContinousWnd1
						>= 0
						|| vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nInterContinousWnd2
								>= 0)
				{
					int nWndN =
							vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nInterContinousWnd1;
					int nWndC =
							vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nInterContinousWnd2;
					int nPos = -1;
					if (m_ionmz[pos_pep].nPepPosOrder1 >= 0)
					{
						nPos = m_ionmz[pos_pep].nPepPosOrder1;
						if (nWndN >= 0 && m_ionmz[pos_pep].bNTerm1)
						{
							m_InterContiWnd[nWndN][nPos]++;
						}
						else if (nWndC >= 0 && !m_ionmz[pos_pep].bNTerm1)
						{
							m_InterContiWnd[nWndC][nPos]++;
						}
					}
					if (m_ionmz[pos_pep].nPepPosOrder2 >= 0)
					{
						nPos = m_ionmz[pos_pep].nPepPosOrder2;
						if (nWndN >= 0 && m_ionmz[pos_pep].bNTerm2)
						{
							m_InterContiWnd[nWndN][nPos]++;
						}
						else if (nWndC >= 0 && !m_ionmz[pos_pep].bNTerm2)
						{
							m_InterContiWnd[nWndC][nPos]++;
						}
					}
				}

			}
			++pos_spec;
		}
	}

}

void CXLinkRefineScorer::_CalMatchInfo()
{

	m_pPeptide->m_stMatchInfo.lfMatchedSpecInt = 0.0;
	m_pPeptide->m_stMatchInfo.lfUnMatchedSpecInt = 0.0;
	memset(m_pPeptide->m_stMatchInfo.aPepConf, 0,
			sizeof(int) * 2 * MAX_PEPTIDE_LENGTH);


	int aCleavConf[MAX_PEPTIDE_LENGTH + 1];
	memset(aCleavConf, 0, sizeof(int) * (1 + MAX_PEPTIDE_LENGTH));
	aCleavConf[0] = 1;
	aCleavConf[m_pPeptide->m_AlphaPeptide.m_tLength] = 1;

	for (int i = 0; i < (int) m_Condition.m_vIonTypes.size(); ++i)
	{
		const CIonType & IonType = m_Condition.m_vIonTypes[i];

		if (!IonType.bIntraContinous)
			continue;

		for (size_t j = 0; j < m_pPeptide->m_AlphaPeptide.m_tLength - 1; ++j)
		{
			if (m_iontag_pep[i][j] < 0)
			{
				continue;
			}
			if (IonType.cType == 0)
			{
				aCleavConf[j + 1]++;
			}
			else if (IonType.cType == 1)
			{
				aCleavConf[m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j]++;
			}
		}
	}

	for (size_t i = 0; i < m_pPeptide->m_AlphaPeptide.m_tLength; ++i)
	{
		m_pPeptide->m_stMatchInfo.aPepConf[i] =
				(aCleavConf[i] <= aCleavConf[i + 1]) ?
						aCleavConf[i] : aCleavConf[i + 1];
	}

	if (m_pPeptide->m_bPair)
	{
		memset(aCleavConf, 0, sizeof(int) * (MAX_PEPTIDE_LENGTH + 1));
		aCleavConf[0] = 1;
		aCleavConf[m_pPeptide->m_BetaPeptide.m_tLength] = 1;

		int nPepLen1 = m_pPeptide->m_AlphaPeptide.m_tLength;
		for (int i = 0; i < (int) m_Condition.m_vIonTypes.size(); ++i)
		{
			const CIonType & IonType = m_Condition.m_vIonTypes[i];

			if (!IonType.bIntraContinous)
				continue;

			for (size_t j = 0; j < m_pPeptide->m_BetaPeptide.m_tLength - 1; ++j)
			{
				if (m_iontag_pep[i][nPepLen1 + j] < 0)
				{
					continue;
				}
				if (IonType.cType == 0)
				{
					aCleavConf[j + 1]++;
				}
				else if (IonType.cType == 1)
				{
					aCleavConf[m_pPeptide->m_BetaPeptide.m_tLength - 1 - j]++;
				}
			}
		}

		for (size_t i = 0; i < m_pPeptide->m_BetaPeptide.m_tLength; ++i)
		{
			m_pPeptide->m_stMatchInfo.aPepConf[m_pPeptide->m_AlphaPeptide.m_tLength
					+ i] =
					(aCleavConf[i] <= aCleavConf[i + 1]) ?
							aCleavConf[i] : aCleavConf[i + 1];
		}
	}


	for (size_t i = 0; i < m_pSpectrum->m_tPeaksNum; ++i)
	{
		if (m_pbMatched[i])
		{
			m_pPeptide->m_stMatchInfo.lfMatchedSpecInt +=
					m_pSpectrum->m_pPeaks[i].lfIntensity;
		}
		else
		{
			m_pPeptide->m_stMatchInfo.lfUnMatchedSpecInt +=
					m_pSpectrum->m_pPeaks[i].lfIntensity;
		}
	}

}

void CXLinkRefineScorer::_ComputeMZ()
{



	int nPepLen1 = m_pPeptide->m_AlphaPeptide.m_tLength;
	int nPepLen2 = 0;

	int nPepAAMass[2 * MAX_PEPTIDE_LENGTH] =
	{ 0 };

	if (m_pPeptide->m_bPair)
	{
		nPepLen2 = m_pPeptide->m_BetaPeptide.m_tLength;
	}
	else
	{
		nPepLen2 = 0;
	}

	int i;
	for (int i = 0; i < nPepLen1; ++i)
	{
		nPepAAMass[i] = m_nAAMass[m_pPeptide->m_AlphaPeptide.m_szSequence[i]
				- 65];
	}

	for (i = 0; i < nPepLen2; ++i)
	{
		nPepAAMass[nPepLen1 + i] =
				m_nAAMass[m_pPeptide->m_BetaPeptide.m_szSequence[i] - 65];
	}


	if (m_Condition.m_bFragmentMono)
	{
		for (size_t i = 0; i < m_pPeptide->m_AlphaPeptide.m_tModCnt; ++i)
		{
			nPepAAMass[m_pPeptide->m_AlphaPeptide.m_tModSites[i][0]] +=
					(int) (m_Condition.m_vSelectedVarMod[m_pPeptide->m_AlphaPeptide.m_tModSites[i][1]].m_lfMonoMass_dif
							* MZMULTIPLIER);
		}

		if (m_pPeptide->m_bPair)
		{
			for (size_t i = 0; i < m_pPeptide->m_BetaPeptide.m_tModCnt; ++i)
			{
				nPepAAMass[nPepLen1
						+ m_pPeptide->m_BetaPeptide.m_tModSites[i][0]] +=
						(int) (m_Condition.m_vSelectedVarMod[m_pPeptide->m_BetaPeptide.m_tModSites[i][1]].m_lfMonoMass_dif
								* MZMULTIPLIER);
			}
		}
	}
	else
	{
		for (size_t i = 0; i < m_pPeptide->m_AlphaPeptide.m_tModCnt; ++i)
		{
			nPepAAMass[m_pPeptide->m_AlphaPeptide.m_tModSites[i][0]] +=
					(int) (m_Condition.m_vSelectedVarMod[m_pPeptide->m_AlphaPeptide.m_tModSites[i][1]].m_lfAvrgMass_dif
							* MZMULTIPLIER);
		}

		if (m_pPeptide->m_bPair)
		{
			for (size_t i = 0; i < m_pPeptide->m_BetaPeptide.m_tModCnt; ++i)
			{
				nPepAAMass[nPepLen1
						+ m_pPeptide->m_BetaPeptide.m_tModSites[i][0]] +=
						(int) (m_Condition.m_vSelectedVarMod[m_pPeptide->m_BetaPeptide.m_tModSites[i][1]].m_lfAvrgMass_dif
								* MZMULTIPLIER);
			}
		}

	}


	int bPepAAMass[2 * MAX_PEPTIDE_LENGTH] =
	{ 0 };
	int yPepAAMass[2 * MAX_PEPTIDE_LENGTH] =
	{ 0 };

	int nLinkerMass = 0;

	if (m_pPeptide->m_XLink.m_eXLinkType == 1)
	{
		if (m_Condition.m_bFragmentMono)
		{
			nLinkerMass =
					(int) (m_Condition.m_vSelectedXLinker[m_pPeptide->m_XLink.m_nLinkerId].m_lfMLMonoMass_dif
							* MZMULTIPLIER);
			nPepAAMass[m_pPeptide->m_XLink.m_tAlphaSite] += nLinkerMass;
		}
		else
		{
			nLinkerMass =
					(int) (m_Condition.m_vSelectedXLinker[m_pPeptide->m_XLink.m_nLinkerId].m_lfMLAvrgMass_dif
							* MZMULTIPLIER);
			nPepAAMass[m_pPeptide->m_XLink.m_tAlphaSite] += nLinkerMass;
		}
	}
	else if (m_pPeptide->m_XLink.m_eXLinkType > 1)
	{
		if (m_Condition.m_bFragmentMono)
			nLinkerMass =
					(int) (m_Condition.m_vSelectedXLinker[m_pPeptide->m_XLink.m_nLinkerId].m_lfMonoMass_dif
							* MZMULTIPLIER);
		else
			nLinkerMass =
					(int) (m_Condition.m_vSelectedXLinker[m_pPeptide->m_XLink.m_nLinkerId].m_lfAvrgMass_dif
							* MZMULTIPLIER);
	}

	int nTmpMass = 0;
	int nPepMass1 = 0, nPepMass2 = 0;

	nPepMass1 = 0;
	nPepMass2 = 0;

	for (i = 0; i < nPepLen1; ++i)
	{
		nPepMass1 += nPepAAMass[i];
	}

	if (m_Condition.m_bFragmentMono)
		nPepMass1 += 2 * nhmass_mono_multi + nomass_mono_multi;
	else
		nPepMass1 += 2 * nhmass_avrg_multi + nomass_avrg_multi;

	for (i = 0; i < nPepLen2; ++i)
	{
		nPepMass2 += nPepAAMass[nPepLen1 + i];
	}

	if (m_Condition.m_bFragmentMono)
		nPepMass2 += 2 * nhmass_mono_multi + nomass_mono_multi;
	else
		nPepMass2 += 2 * nhmass_avrg_multi + nomass_avrg_multi;

	int nNTermDeltaMass = npmass_multi;
	int nCTermDeltaMass;
	int nPMass = npmass_multi;
	if (m_Condition.m_bFragmentMono)
	{
		nCTermDeltaMass = 2 * nhmass_mono_multi + nomass_mono_multi
				+ npmass_multi;

	}
	else
	{
		nCTermDeltaMass = 2 * nhmass_avrg_multi + nomass_avrg_multi
				+ npmass_multi;

	}

	if (m_pPeptide->m_XLink.m_eXLinkType < 2)
	{

		bPepAAMass[0] = (int) (nPepAAMass[0] + nNTermDeltaMass);
		yPepAAMass[0] = (int) (nPepAAMass[nPepLen1 - 1] + nCTermDeltaMass);

		for (int i = 1; i < nPepLen1 - 1; ++i)
		{
			bPepAAMass[i] = bPepAAMass[i - 1] + nPepAAMass[i];
		}

		for (int i = 1, j = nPepLen1 - 2; j > 0; ++i, --j)
		{
			yPepAAMass[i] = yPepAAMass[i - 1] + nPepAAMass[j];
		}

	}
	else if (m_pPeptide->m_XLink.m_eXLinkType == 2)
	{



		nTmpMass = int(nNTermDeltaMass);
		for (i = 0; i < m_pPeptide->m_XLink.m_tAlphaSite; ++i)
		{
			nTmpMass += nPepAAMass[i];
			bPepAAMass[i] = nTmpMass;
		}

		for (; i < m_pPeptide->m_XLink.m_tBetaSite; ++i)
		{
			nTmpMass += nPepAAMass[i];
			bPepAAMass[i] = nPepMass1;
		}
		nTmpMass += nLinkerMass;
		for (; i < nPepLen1; ++i)
		{
			nTmpMass += nPepAAMass[i];
			bPepAAMass[i] = nTmpMass;
		}


		nTmpMass = int(nCTermDeltaMass);
		for (i = nPepLen1 - 1; i > m_pPeptide->m_XLink.m_tBetaSite; --i)
		{
			nTmpMass += nPepAAMass[i];
			yPepAAMass[nPepLen1 - i - 1] = nTmpMass;
		}

		for (; i > m_pPeptide->m_XLink.m_tAlphaSite; --i)
		{
			nTmpMass += nPepAAMass[i];
			yPepAAMass[nPepLen1 - i - 1] = nPepMass1;
		}
		nTmpMass += nLinkerMass;
		for (; i > 0; --i)
		{
			nTmpMass += nPepAAMass[i];
			yPepAAMass[nPepLen1 - i - 1] = nTmpMass;
		}
	}
	else
	{



		nTmpMass = int(nNTermDeltaMass);
		for (i = 0; i < nPepLen1; ++i)
		{
			nTmpMass += nPepAAMass[i];
			if (i == m_pPeptide->m_XLink.m_tAlphaSite)
			{
				nTmpMass += nPepMass2 + nLinkerMass;
			}
			bPepAAMass[i] = nTmpMass;
		}

		nTmpMass = int(nCTermDeltaMass);
		for (i = nPepLen1 - 1; i > 0; --i)
		{
			nTmpMass += nPepAAMass[i];
			if (i == m_pPeptide->m_XLink.m_tAlphaSite)
			{
				nTmpMass += nPepMass2 + nLinkerMass;
			}

			yPepAAMass[nPepLen1 - 1 - i] = nTmpMass;
		}


		nTmpMass = int(nNTermDeltaMass);
		for (i = 0; i < nPepLen2; ++i)
		{
			nTmpMass += nPepAAMass[nPepLen1 + i];
			if (i == m_pPeptide->m_XLink.m_tBetaSite)
			{
				nTmpMass += nLinkerMass + nPepMass1;
			}
			bPepAAMass[nPepLen1 + i] = nTmpMass;
		}

		nTmpMass = int(nCTermDeltaMass);
		for (i = nPepLen2 - 1; i > 0; --i)
		{
			nTmpMass += nPepAAMass[nPepLen1 + i];
			if (i == m_pPeptide->m_XLink.m_tBetaSite)
			{
				nTmpMass += nLinkerMass + nPepMass1;
			}
			yPepAAMass[nPepLen1 + nPepLen2 - 1 - i] = nTmpMass;
		}
	}

	CMzTriple triple;
	triple.clear();
	m_tionmzSize = 0;
	m_ionmz.clear();
	bool bContainLinker = false;

	int j;
	for (i = 0; i < (int) m_Condition.m_vIonTypes.size(); ++i)
	{
		const CIonType & IonType = m_Condition.m_vIonTypes[i];
		if (IonType.cType >= 2)
			continue;
		bContainLinker = false;

		for (j = 0; j < (int) m_pPeptide->m_AlphaPeptide.m_tLength - 1; ++j)
		{
			bContainLinker = false;
			if (m_pPeptide->m_XLink.m_eXLinkType == 2)
			{
				if (IonType.cType == 0
						&& (j >= m_pPeptide->m_XLink.m_tAlphaSite
								&& j < m_pPeptide->m_XLink.m_tBetaSite))
				{
					continue;
				}
				else if (IonType.cType == 1
						&& (m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j
								> m_pPeptide->m_XLink.m_tAlphaSite
								&& m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j
										<= m_pPeptide->m_XLink.m_tBetaSite))
				{
					continue;
				}
			}
#ifdef ONLY_CONSIDER_FRAGMENTS_WITH_LINKER

			if(m_pPeptide->m_XLink.m_eXLinkType == 3)
			{
				if(IonType.cType==0 && j < m_pPeptide->m_XLink.m_tAlphaSite)
				{
					continue;
				}
				else if(IonType.cType==1 && m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j > m_pPeptide->m_XLink.m_tAlphaSite)
				{
					continue;
				}
			}
#endif


			if (m_pPeptide->m_XLink.m_eXLinkType == 3)
			{
				if (IonType.cType == 0)
				{
					if (j >= m_pPeptide->m_XLink.m_tAlphaSite)
					{

						bContainLinker = true;
						if (IonType.nContainLinker == 1)

							continue;
					}
					else
					{

						bContainLinker = false;
						if (IonType.nContainLinker == 2)

							continue;
					}
				}
				else if (IonType.cType == 1)
				{
					if (m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j
							<= m_pPeptide->m_XLink.m_tAlphaSite)
					{

						bContainLinker = true;
						if (IonType.nContainLinker == 1)

							continue;
					}
					else
					{

						bContainLinker = false;
						if (IonType.nContainLinker == 2)

							continue;
					}
				}
			}

			triple.clear();
			triple.bContainLinker = bContainLinker;
			triple.nIonTypeOrder = i;
			triple.bAddToTag = IonType.bIntraContinous;
			triple.nPepPosOrder1 = j;
			triple.nPepPosOrder2 = -1;

			if (IonType.cType == 0)
				triple.bNTerm1 = true;
			else
				triple.bNTerm1 = false;

			triple.bNTerm2 = false;
			triple.nMz = IonType.cType == 0 ? bPepAAMass[j] : yPepAAMass[j];
			triple.nMz += (IonType.nCharge - 1) * npmass_multi
					- IonType.nTotalLostVal;
			triple.nMz /= IonType.nCharge;

			m_ionmz.push_back(triple);

			++m_tionmzSize;

		}
		for (size_t k = 0; k < m_pPeptide->m_AlphaPeptide.m_tModCnt; ++k)
		{
			for (size_t ww = 0;
					ww
							< m_Condition.m_vSelectedVarMod[m_pPeptide->m_AlphaPeptide.m_tModSites[k][1]].m_tNLSize;
					++ww)
			{
				double lfLoss =
						m_Condition.m_vSelectedVarMod[m_pPeptide->m_AlphaPeptide.m_tModSites[k][1]].m_vlfMonoNeutralLoss_dif[ww];
				lfLoss *= MZMULTIPLIER;
				triple.clear();

				for (int j = m_pPeptide->m_AlphaPeptide.m_tModSites[k][0];
						j < (int) m_pPeptide->m_AlphaPeptide.m_tLength - 1; ++j)
				{
					bContainLinker = false;
					if (m_pPeptide->m_XLink.m_eXLinkType == 2)
					{
						if (IonType.cType == 0
								&& (j >= m_pPeptide->m_XLink.m_tAlphaSite
										&& j < m_pPeptide->m_XLink.m_tBetaSite))
						{
							continue;
						}
						else if (IonType.cType == 1
								&& (m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j
										> m_pPeptide->m_XLink.m_tAlphaSite
										&& m_pPeptide->m_AlphaPeptide.m_tLength
												- 1 - j
												<= m_pPeptide->m_XLink.m_tBetaSite))
						{
							continue;
						}
					}
#ifdef ONLY_CONSIDER_FRAGMENTS_WITH_LINKER

					if(m_pPeptide->m_XLink.m_eXLinkType == 3)
					{
						if(IonType.cType==0 && j < m_pPeptide->m_XLink.m_tAlphaSite)
						{
							continue;
						}
						else if(IonType.cType==1 && m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j > m_pPeptide->m_XLink.m_tAlphaSite)
						{
							continue;
						}
					}
#endif

					if (m_pPeptide->m_XLink.m_eXLinkType == 3)
					{
						if (IonType.cType == 0)
						{
							if (j >= m_pPeptide->m_XLink.m_tAlphaSite)
							{

								bContainLinker = true;
								if (IonType.nContainLinker == 1)

									continue;
							}
							else
							{

								bContainLinker = false;
								if (IonType.nContainLinker == 2)

									continue;
							}
						}
						else if (IonType.cType == 1)
						{
							if (m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j
									<= m_pPeptide->m_XLink.m_tAlphaSite)
							{

								bContainLinker = true;
								if (IonType.nContainLinker == 1)

									continue;
							}
							else
							{

								bContainLinker = false;
								if (IonType.nContainLinker == 2)

									continue;
							}
						}
					}
					triple.bContainLinker = bContainLinker;
					triple.nIonTypeOrder = i;
					triple.bAddToTag = false;
					triple.nPepPosOrder1 = j;
					triple.nPepPosOrder2 = -1;

					if (IonType.cType == 0)
						triple.bNTerm1 = true;
					else
						triple.bNTerm1 = false;

					triple.bNTerm2 = false;
					triple.nMz =
							IonType.cType == 0 ? bPepAAMass[j] : yPepAAMass[j];
					triple.nMz = int(triple.nMz - lfLoss);
					triple.nMz += (IonType.nCharge - 1) * npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					if (triple.nMz < 0)
						continue;
					triple.nMz /= IonType.nCharge;
					m_ionmz.push_back(triple);
					m_tionmzSize++;
				}
			}
		}

		if (m_pPeptide->m_bPair)
		{
			for (j = 0; j < (int) m_pPeptide->m_BetaPeptide.m_tLength - 1; ++j)
			{
#ifdef ONLY_CONSIDER_FRAGMENTS_WITH_LINKER

				if(m_pPeptide->m_XLink.m_eXLinkType == 3)
				{
					if(IonType.cType==0 && j < m_pPeptide->m_XLink.m_tBetaSite)
					{
						continue;
					}
					else if(IonType.cType==1 && m_pPeptide->m_BetaPeptide.m_tLength - 1 - j > m_pPeptide->m_XLink.m_tBetaSite)
					{
						continue;
					}
				}
#endif


				bContainLinker = false;
				if (m_pPeptide->m_XLink.m_eXLinkType == 3)
				{
					if (IonType.cType == 0)
					{
						if (j >= m_pPeptide->m_XLink.m_tBetaSite)
						{

							bContainLinker = true;
							if (IonType.nContainLinker == 1)

								continue;
						}
						else
						{

							bContainLinker = false;
							if (IonType.nContainLinker == 2)

								continue;
						}
					}
					else if (IonType.cType == 1)
					{
						if (m_pPeptide->m_BetaPeptide.m_tLength - 1 - j
								<= m_pPeptide->m_XLink.m_tBetaSite)
						{

							bContainLinker = true;
							if (IonType.nContainLinker == 1)

								continue;
						}
						else
						{

							bContainLinker = false;
							if (IonType.nContainLinker == 2)

								continue;
						}
					}
				}

				triple.clear();
				triple.bContainLinker = bContainLinker;
				triple.nIonTypeOrder = i;
				triple.bAddToTag = IonType.bIntraContinous;
				triple.nPepPosOrder1 = nPepLen1 + j;
				triple.nPepPosOrder2 = -1;

				if (IonType.cType == 0)
				{
					triple.bNTerm1 = true;
				}
				else if (IonType.cType == 1)
				{
					triple.bNTerm1 = false;
				}

				triple.bNTerm2 = false;
				triple.nMz =
						IonType.cType == 0 ?
								bPepAAMass[nPepLen1 + j] :
								yPepAAMass[nPepLen1 + j];
				triple.nMz += (IonType.nCharge - 1) * npmass_multi
						- IonType.nTotalLostVal;
				triple.nMz /= IonType.nCharge;
				++m_tionmzSize;

				m_ionmz.push_back(triple);

			}
			for (size_t k = 0; k < m_pPeptide->m_BetaPeptide.m_tModCnt; ++k)
			{
				for (size_t ww = 0;
						ww
								< m_Condition.m_vSelectedVarMod[m_pPeptide->m_BetaPeptide.m_tModSites[k][1]].m_tNLSize;
						++ww)
				{
					double lfLoss =
							m_Condition.m_vSelectedVarMod[m_pPeptide->m_BetaPeptide.m_tModSites[k][1]].m_vlfMonoNeutralLoss_dif[ww];
					lfLoss *= MZMULTIPLIER;

					triple.clear();
					for (int j = m_pPeptide->m_BetaPeptide.m_tModSites[k][0];
							j < (int) m_pPeptide->m_BetaPeptide.m_tLength - 1;
							++j)
					{
#ifdef ONLY_CONSIDER_FRAGMENTS_WITH_LINKER

						if(m_pPeptide->m_XLink.m_eXLinkType == 3)
						{
							if(IonType.cType==0 && j < m_pPeptide->m_XLink.m_tBetaSite)
							{
								continue;
							}
							else if(IonType.cType==1 && m_pPeptide->m_BetaPeptide.m_tLength - 1 - j > m_pPeptide->m_XLink.m_tBetaSite)
							{
								continue;
							}
						}
#endif

						bContainLinker = false;
						if (m_pPeptide->m_XLink.m_eXLinkType == 3)
						{
							if (IonType.cType == 0)
							{
								if (j >= m_pPeptide->m_XLink.m_tBetaSite)
								{

									bContainLinker = true;
									if (IonType.nContainLinker == 1)

										continue;
								}
								else
								{

									bContainLinker = false;
									if (IonType.nContainLinker == 2)

										continue;
								}
							}
							else if (IonType.cType == 1)
							{
								if (m_pPeptide->m_BetaPeptide.m_tLength - 1 - j
										<= m_pPeptide->m_XLink.m_tBetaSite)
								{

									bContainLinker = true;
									if (IonType.nContainLinker == 1)

										continue;
								}
								else
								{

									bContainLinker = false;
									if (IonType.nContainLinker == 2)

										continue;
								}
							}
						}
						triple.bContainLinker = bContainLinker;
						triple.nIonTypeOrder = i;
						triple.bAddToTag = false;
						triple.nPepPosOrder1 = nPepLen1 + j;
						triple.nPepPosOrder2 = -1;

						if (IonType.cType == 0)
						{
							triple.bNTerm1 = true;
						}
						else if (IonType.cType == 1)
						{
							triple.bNTerm1 = false;
						}

						triple.bNTerm2 = false;
						triple.nMz =
								IonType.cType == 0 ?
										bPepAAMass[nPepLen1 + j] :
										yPepAAMass[nPepLen1 + j];
						triple.nMz = int(triple.nMz - lfLoss);
						triple.nMz += (IonType.nCharge - 1) * npmass_multi;
						triple.nMz -= IonType.nTotalLostVal;
						if (triple.nMz < 0)
							continue;
						triple.nMz /= IonType.nCharge;
						m_ionmz.push_back(triple);
						m_tionmzSize++;
					}
				}
			}

		}
	}


	triple.clear();

	triple.bContainLinker = false;
	triple.bNTerm1 = true;
	triple.bNTerm2 = false;
	triple.bAddToTag = false;
	for (size_t i = 0; i < m_Condition.m_vIonTypes.size(); ++i)
	{
		const CIonType & IonType = m_Condition.m_vIonTypes[i];
		if (IonType.cType == 3)
		{
			triple.nIonTypeOrder = i;
			triple.nMz = nPepMass1 + (IonType.nCharge) * npmass_multi;
			triple.nMz -= IonType.nTotalLostVal;
			triple.nMz /= IonType.nCharge;
			triple.nPepPosOrder1 = nPepLen1 - 1;
			triple.nPepPosOrder2 = -1;
			m_tionmzSize++;
			m_ionmz.push_back(triple);
			if (nPepLen2 > 0)
			{
				triple.nIonTypeOrder = i;
				triple.nMz = nPepMass2 + (IonType.nCharge) * npmass_multi;
				triple.nMz -= IonType.nTotalLostVal;
				triple.nMz /= IonType.nCharge;
				triple.nPepPosOrder1 = nPepLen1 + nPepLen2 - 1;
				triple.nPepPosOrder2 = -1;
				m_tionmzSize++;
				m_ionmz.push_back(triple);
			}
		}
	}


	vector<int> vKLIonOrder;
	vKLIonOrder.clear();
	for (int i = 0; i < (int) m_Condition.m_vIonTypes.size(); ++i)
	{
		const CIonType & IonType = m_Condition.m_vIonTypes[i];
		if (IonType.cType == 5)
		{
			vKLIonOrder.push_back(i);
		}
	}
	if (vKLIonOrder.size() > 0)
	{
		if (m_pPeptide->m_XLink.m_eXLinkType == 3)
		{

			int nSum = nPepMass2 + nPepAAMass[m_pPeptide->m_XLink.m_tAlphaSite]
					+ nLinkerMass + npmass_multi;
			triple.clear();
			triple.nPepPosOrder1 = nPepLen1 - 1
					- m_pPeptide->m_XLink.m_tAlphaSite;
			triple.nPepPosOrder2 = m_pPeptide->m_XLink.m_tAlphaSite;

			triple.bNTerm1 = false;
			triple.bNTerm2 = true;
			triple.bAddToTag = false;
			for (size_t w = 0; w < vKLIonOrder.size(); ++w)
			{
				const CIonType & IonType =
						m_Condition.m_vIonTypes[vKLIonOrder[w]];
				triple.nIonTypeOrder = vKLIonOrder[w];
				triple.nMz = nSum + (IonType.nCharge - 1) * npmass_multi;
				triple.nMz -= IonType.nTotalLostVal;
				triple.nMz /= IonType.nCharge;
				m_tionmzSize++;
				m_ionmz.push_back(triple);
			}


			nSum = nPepMass1
					+ nPepAAMass[nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite]
					+ nLinkerMass + npmass_multi;
			triple.clear();
			triple.nPepPosOrder1 = nPepLen1 + nPepLen2 - 1
					- m_pPeptide->m_XLink.m_tBetaSite;
			triple.nPepPosOrder2 = nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite;

			triple.bNTerm1 = false;
			triple.bNTerm2 = true;
			triple.bAddToTag = false;
			for (size_t w = 0; w < vKLIonOrder.size(); ++w)
			{
				const CIonType & IonType =
						m_Condition.m_vIonTypes[vKLIonOrder[w]];
				triple.nIonTypeOrder = vKLIonOrder[w];
				triple.nMz = nSum + (IonType.nCharge - 1) * npmass_multi;
				triple.nMz -= IonType.nTotalLostVal;
				triple.nMz /= IonType.nCharge;
				m_tionmzSize++;
				m_ionmz.push_back(triple);
			}
		}
	}




	triple.clear();
	triple.nPepPosOrder1 = -1;
	triple.nPepPosOrder2 = -1;
	triple.bNTerm1 = false;
	triple.bNTerm2 = false;
	triple.bAddToTag = false;
	triple.bContainLinker = false;

	for (int m = 0; m < m_Condition.m_vIonTypes.size(); ++m)
	{
		const CIonType & IonType = m_Condition.m_vIonTypes[m];
		if (IonType.cType == 6)
		{
			triple.nIonTypeOrder = m;

			if (m_pPeptide->m_XLink.m_eXLinkType == 3)
			{
				for (int i = 1; i < m_pPeptide->m_XLink.m_tAlphaSite; i++)
				{
					triple.nMz = 0;


					for (int j = i; j < nPepLen1; j++)
					{
						triple.nMz += nPepAAMass[j];

					}
					triple.nPepPosOrder1 = nPepLen1 - i - 1;
					triple.nPepPosOrder2 = nPepLen1;
					triple.nMz += 2 * nhmass_mono_multi + nomass_mono_multi;

					triple.nMz += npmass_multi;
					triple.bContainLinker = true;

					triple.nMz += (IonType.nCharge - 1) * npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					triple.nMz /= IonType.nCharge;
					m_tionmzSize++;
					m_ionmz.push_back(triple);
				}
				for (int i = 1; i < m_pPeptide->m_XLink.m_tBetaSite; i++)
				{
					triple.nMz = 0;


					for (int j = i; j < nPepLen2; j++)
					{
						triple.nMz += nPepAAMass[nPepLen1 + j];

					}
					triple.nPepPosOrder1 = nPepLen1 + nPepLen2 - i - 1;
					triple.nPepPosOrder2 = nPepLen1 + nPepLen2;
					triple.nMz += 2 * nhmass_mono_multi + nomass_mono_multi;

					triple.nMz += npmass_multi;
					triple.bContainLinker = true;

					triple.nMz += (IonType.nCharge - 1) * npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					triple.nMz /= IonType.nCharge;
					m_tionmzSize++;
					m_ionmz.push_back(triple);
				}
			}
		}
	}


	triple.clear();
	triple.nPepPosOrder1 = -1;
	triple.nPepPosOrder2 = -1;
	triple.bNTerm1 = false;
	triple.bNTerm2 = false;
	triple.bAddToTag = false;
	triple.bContainLinker = false;

	for (int m = 0; m < m_Condition.m_vIonTypes.size(); ++m)
	{
		const CIonType & IonType = m_Condition.m_vIonTypes[m];
		if (IonType.cType == 7)
		{
			triple.nIonTypeOrder = m;

			if (m_pPeptide->m_XLink.m_eXLinkType == 3)
			{
				for (int i = m_pPeptide->m_XLink.m_tAlphaSite - 1;
						i < nPepLen1 - 1; i++)
				{
					triple.nMz = 0;


					for (int j = 0; j <= i; j++)
					{
						triple.nMz += nPepAAMass[j];

					}
					triple.nPepPosOrder1 = -1;
					triple.nPepPosOrder2 = i;

					triple.nMz += npmass_multi;
					triple.bContainLinker = true;

					triple.nMz += (IonType.nCharge - 1) * npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					triple.nMz /= IonType.nCharge;
					m_tionmzSize++;
					m_ionmz.push_back(triple);






				}
				for (int i = m_pPeptide->m_XLink.m_tBetaSite - 1;
						i < nPepLen2 - 1; i++)
				{
					triple.nMz = 0;


					for (int j = 0; j <= i; j++)
					{
						triple.nMz += nPepAAMass[nPepLen1 + j];

					}
					triple.nPepPosOrder1 = nPepLen1;
					triple.nPepPosOrder2 = i + nPepLen1;

					triple.nMz += npmass_multi;
					triple.bContainLinker = true;

					triple.nMz += (IonType.nCharge - 1) * npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					triple.nMz /= IonType.nCharge;
					m_tionmzSize++;
					m_ionmz.push_back(triple);
				}
			}
		}
	}


	vector<int> vInternalIonOrder;
	vInternalIonOrder.clear();
	for (int i = 0; i < (int) m_Condition.m_vIonTypes.size(); ++i)
	{
		const CIonType & IonType = m_Condition.m_vIonTypes[i];
		if (IonType.cType == 2)
		{
			vInternalIonOrder.push_back(i);
		}
	}
	if (vInternalIonOrder.size() <= 0)
		return;

	int nSum = 0;
	int nPepLen = nPepLen1 + nPepLen2;
	int nAAnum = 0;
	bContainLinker = false;
	for (int i = 0; i < nPepLen; ++i)
	{
		for (int j = i; j < nPepLen; ++j)
		{



			nSum = 0;
			bool bNTerm1 = false;
			bool bNTerm2 = false;
			bContainLinker = false;
			nAAnum = 0;
			triple.clear();
			if (i < nPepLen1)
			{
				if (j < nPepLen1)
				{
					bNTerm1 = false;
					bNTerm2 = true;
					nAAnum = 0;
					bContainLinker = false;
					if (i == 0 || j == nPepLen1 - 1)
						continue;


					nSum = 0;
					for (int k = i; k <= j; ++k)
					{
						nSum += nPepAAMass[k];
						nAAnum++;
					}

					nSum += npmass_multi;

					if (m_pPeptide->m_XLink.m_eXLinkType == 2)
					{
						if (i <= m_pPeptide->m_XLink.m_tAlphaSite
								&& j >= m_pPeptide->m_XLink.m_tBetaSite)
						{
							nSum += nLinkerMass;
						}
						else if (i <= m_pPeptide->m_XLink.m_tAlphaSite
								&& j < m_pPeptide->m_XLink.m_tBetaSite)
						{
							continue;
						}
						else if (i > m_pPeptide->m_XLink.m_tAlphaSite
								&& j >= m_pPeptide->m_XLink.m_tBetaSite)
						{
							continue;
						}

					}
					else if (m_pPeptide->m_XLink.m_eXLinkType == 3)
					{
						if (i <= m_pPeptide->m_XLink.m_tAlphaSite
								&& j >= m_pPeptide->m_XLink.m_tAlphaSite)
						{
							nSum += (nLinkerMass + nPepMass2);
							nAAnum += nPepLen2;
							bContainLinker = true;
						}
					}

					if (nAAnum <= MAX_INNER_ION_LENGTH)
					{
						triple.nPepPosOrder1 = nPepLen1 - 1 - i;
						triple.nPepPosOrder2 = j;

						triple.bNTerm1 = bNTerm1;
						triple.bNTerm2 = bNTerm2;
						triple.bAddToTag = false;
						for (size_t w = 0; w < vInternalIonOrder.size(); ++w)
						{
							const CIonType & IonType =
									m_Condition.m_vIonTypes[vInternalIonOrder[w]];

							if (IonType.nContainLinker == 0
									|| (IonType.nContainLinker == 1
											&& bContainLinker == false)
									|| (IonType.nContainLinker == 2
											&& bContainLinker == true))
							{
								triple.nIonTypeOrder = vInternalIonOrder[w];
								triple.bContainLinker = bContainLinker;
								triple.nMz = nSum
										+ (IonType.nCharge - 1) * npmass_multi;
								triple.nMz -= IonType.nTotalLostVal;
								triple.nMz /= IonType.nCharge;
								m_tionmzSize++;
								m_ionmz.push_back(triple);
							}
						}
					}
				}
				else
				{
					if (i <= m_pPeptide->m_XLink.m_tAlphaSite)
					{
						if (j <= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite)
						{


							bNTerm1 = false;
							bNTerm2 = false;
							nAAnum = 0;
							bContainLinker = false;

							if (i != 0 && j != nPepLen1)
							{

								nSum = 0;
								for (int k = i; k < nPepLen1; ++k)
								{
									nSum += nPepAAMass[k];
									nAAnum++;
								}
								for (int k = j; k < nPepLen; ++k)
								{
									nSum += nPepAAMass[k];
									nAAnum++;
								}

								nSum += 2
										* (2 * nhmass_mono_multi
												+ nomass_mono_multi)
										+ npmass_multi;
								nSum += nLinkerMass;
								bContainLinker = true;

								triple.nPepPosOrder1 = nPepLen1 - 1 - i;
								triple.nPepPosOrder2 = nPepLen2 - 1
										- (j - nPepLen1) + nPepLen1;

								triple.bNTerm1 = bNTerm1;
								triple.bNTerm2 = bNTerm2;
								triple.bAddToTag = false;

								/* todo :
								 * don't know how to calculate the ion mass
								 *
								 *
								 if(nAAnum <= MAX_INNER_ION_LENGTH)
								 {
								 for(int w = 0 ; w < vInternalIonOrder.size(); ++ w)
								 {
								 const CIonType & IonType = m_psmconf.m_vIonTypes[vInternalIonOrder[w]];

								 if(IonType.nContainLinker == 0
								 || (IonType.nContainLinker == 1 && bContainLinker == false)
								 || (IonType.nContainLinker == 2 && bContainLinker == true ))
								 {
								 triple.nIonTypeOrder = vInternalIonOrder[w];
								 triple.bContainLinker = bContainLinker;
								 triple.nMz = nSum + (IonType.nCharge - 1)*npmass_multi;
								 triple.nMz -= IonType.nTotalLostVal;
								 triple.nMz /= IonType.nCharge;
								 m_tionmzSize++;
								 m_ionmz.push_back(triple);
								 }
								 }
								 }
								 */
							}
						}
					}
					if (j >= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite)
					{
						bNTerm1 = false;
						bNTerm2 = true;
						nAAnum = 0;
						bContainLinker = false;

						if (i != 0 && j != nPepLen - 1)
						{

							nSum = 0;
							for (int k = i; k <= j; ++k)
							{
								nSum += nPepAAMass[k];
								nAAnum++;
							}

							nSum += 2 * nhmass_mono_multi + nomass_mono_multi
									+ npmass_multi;

							nSum += nLinkerMass;
							bContainLinker = true;
							triple.nPepPosOrder1 = nPepLen1 - 1 - i;
							triple.nPepPosOrder2 = (j - nPepLen1) + nPepLen1;

							triple.bNTerm1 = bNTerm1;
							triple.bNTerm2 = bNTerm2;
							triple.bAddToTag = false;

							if (nAAnum <= MAX_INNER_ION_LENGTH)
							{
								for (size_t w = 0; w < vInternalIonOrder.size();
										++w)
								{
									const CIonType & IonType =
											m_Condition.m_vIonTypes[vInternalIonOrder[w]];

									if (IonType.nContainLinker == 0
											|| (IonType.nContainLinker == 1
													&& bContainLinker == false)
											|| (IonType.nContainLinker == 2
													&& bContainLinker == true))
									{
										triple.nIonTypeOrder =
												vInternalIonOrder[w];
										triple.bContainLinker = bContainLinker;
										triple.nMz = nSum
												+ (IonType.nCharge - 1)
														* npmass_multi;
										triple.nMz -= IonType.nTotalLostVal;
										triple.nMz /= IonType.nCharge;
										m_tionmzSize++;
										m_ionmz.push_back(triple);
									}
								}
							}

						}

					}
					if (i >= m_pPeptide->m_XLink.m_tAlphaSite)
					{
						if (j <= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite)
						{
							bNTerm1 = true;
							bNTerm2 = false;
							nAAnum = 0;
							bContainLinker = false;
							if (i != nPepLen1 - 1 && j != nPepLen1)
							{

								nSum = 0;
								for (int k = 0; k <= i; ++k)
								{
									nSum += nPepAAMass[k];
									nAAnum++;
								}
								for (int k = j; k < nPepLen; ++k)
								{
									nSum += nPepAAMass[k];
									nAAnum++;
								}

								nSum += (2 * nhmass_mono_multi
										+ nomass_mono_multi) + npmass_multi;

								nSum += nLinkerMass;
								bContainLinker = true;
								triple.nPepPosOrder1 = i;
								triple.nPepPosOrder2 = nPepLen2 - 1
										- (j - nPepLen1) + nPepLen1;
								triple.bNTerm1 = bNTerm1;
								triple.bNTerm2 = bNTerm2;
								triple.bAddToTag = false;

								if (nAAnum <= MAX_INNER_ION_LENGTH)
								{
									for (size_t w = 0;
											w < vInternalIonOrder.size(); ++w)
									{
										const CIonType & IonType =
												m_Condition.m_vIonTypes[vInternalIonOrder[w]];
										if (IonType.nContainLinker == 0
												|| (IonType.nContainLinker == 1
														&& bContainLinker
																== false)
												|| (IonType.nContainLinker == 2
														&& bContainLinker
																== true))
										{
											triple.nIonTypeOrder =
													vInternalIonOrder[w];
											triple.bContainLinker =
													bContainLinker;
											triple.nMz = nSum
													+ (IonType.nCharge - 1)
															* npmass_multi;
											triple.nMz -= IonType.nTotalLostVal;
											triple.nMz /= IonType.nCharge;
											m_tionmzSize++;
											m_ionmz.push_back(triple);
										}
									}
								}

							}

						}
						if (j >= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite)
						{

							bNTerm1 = true;
							bNTerm2 = true;
							nAAnum = 0;
							bContainLinker = false;
							if (i != nPepLen1 - 1 && j != nPepLen - 1)
							{

								nSum = 0;
								for (int k = 0; k <= i; ++k)
								{
									nSum += nPepAAMass[k];
									nAAnum++;
								}
								for (int k = nPepLen1; k <= j; ++k)
								{
									nSum += nPepAAMass[k];
									nAAnum++;
								}

								nSum += npmass_multi;

								if (i == nPepLen1 - 1)
								{
									nSum += 2 * nhmass_mono_multi
											+ nomass_mono_multi;
								}
								if (j == nPepLen - 1)
								{
									nSum += 2 * nhmass_mono_multi
											+ nomass_mono_multi;
								}

								nSum += nLinkerMass;
								bContainLinker = true;
								/* todo :
								 * don't know how to calculate the ion mass
								 *
								 triple.nPepPosOrder1 = i;
								 triple.nPepPosOrder2 = (j - nPepLen1) + nPepLen1;
								 triple.bNTerm1 = bNTerm1;
								 triple.bNTerm2 = bNTerm2;
								 triple.bAddToTag = false;

								 if(nAAnum <= 4)
								 {
								 for(int w = 0 ; w < vInternalIonOrder.size(); ++ w)
								 {
								 const CIonType & IonType = m_psmconf.m_vIonTypes[vInternalIonOrder[w]];
								 if(IonType.nContainLinker == 0
								 || (IonType.nContainLinker == 1 && bContainLinker == false)
								 || (IonType.nContainLinker == 2 && bContainLinker == true ))
								 {
								 triple.nIonTypeOrder = vInternalIonOrder[w];
								 triple.bContainLinker = bContainLinker;
								 triple.nMz = nSum + (IonType.nCharge - 1)*npmass_multi;
								 triple.nMz -= IonType.nTotalLostVal;
								 triple.nMz /= IonType.nCharge;
								 m_tionmzSize++;
								 m_ionmz.push_back(triple);
								 }
								 }
								 }
								 */
							}

						}
					}

				}
			}
			else
			{
				bNTerm1 = false;
				bNTerm2 = true;
				nAAnum = 0;
				bContainLinker = false;
				if (i == nPepLen1 || j == nPepLen - 1)
					continue;

				nSum = 0;
				for (int k = i; k <= j; ++k)
				{
					nSum += nPepAAMass[k];
					nAAnum++;
				}

				nSum += npmass_multi;

				if ((i <= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite)
						&& (j >= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite))
				{
					nSum += (nLinkerMass + nPepMass1);
					nAAnum += nPepLen1;
					bContainLinker = true;
				}

				triple.nPepPosOrder1 = nPepLen2 - 1 - (i - nPepLen1) + nPepLen1;
				triple.nPepPosOrder2 = (j - nPepLen1) + nPepLen1;

				triple.bNTerm1 = bNTerm1;
				triple.bNTerm2 = bNTerm2;
				triple.bAddToTag = false;

				if (nAAnum <= MAX_INNER_ION_LENGTH)
				{
					for (size_t w = 0; w < vInternalIonOrder.size(); ++w)
					{
						const CIonType & IonType =
								m_Condition.m_vIonTypes[vInternalIonOrder[w]];

						if (IonType.nContainLinker == 0
								|| (IonType.nContainLinker == 1
										&& bContainLinker == false)
								|| (IonType.nContainLinker == 2
										&& bContainLinker == true))
						{
							triple.nIonTypeOrder = vInternalIonOrder[w];
							triple.bContainLinker = bContainLinker;
							triple.nMz = nSum
									+ (IonType.nCharge - 1) * npmass_multi;
							triple.nMz -= IonType.nTotalLostVal;
							triple.nMz /= IonType.nCharge;
							m_tionmzSize++;
							m_ionmz.push_back(triple);
						}
					}
				}
			}
		}
	}

}

}
