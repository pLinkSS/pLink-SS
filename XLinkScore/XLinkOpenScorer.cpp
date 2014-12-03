#include <iostream>
#include <set>
#include <vector>
#include <map>
#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/interface.h"
#include "../include/option.h"
#include <iomanip>
#include "XLinkOpenScorer.h"
using namespace std;

#define MINUS_INFINITE -10000

#define debug
#ifdef debug
#define PR(x) cout<<setprecision(6)<<"\t"<<#x<<": "<<x;
#define PRL(x) cout<<setprecision(6)<<"\t"<<#x<<": "<<x<<endl;
#else
#define PR(x)
#define PRL(x)
#endif

namespace proteomics_sdk
{
double CXLinkOpenScorer::m_lfLogNPerm[100000] = {0.0};
double CXLinkOpenScorer::m_lfLogTagLenCDF[120][21]= {0.0};
double CXLinkOpenScorer::m_lfLogStdNormCDF[101] = {
		 -1.5064998e+001,
		 -1.4551183e+001,
		 -1.4047029e+001,
		 -1.3552525e+001,
		 -1.3067660e+001,
		 -1.2592420e+001,
		 -1.2126791e+001,
		 -1.1670761e+001,
		 -1.1224313e+001,
		 -1.0787432e+001,
		 -1.0360101e+001,
		 -9.9423044e+000,
		 -9.5340221e+000,
		 -9.1352354e+000,
		 -8.7459236e+000,
		 -8.3660653e+000,
		 -7.9956375e+000,
		 -7.6346161e+000,
		 -7.2829755e+000,
		 -6.9406885e+000,
		 -6.6077262e+000,
		 -6.2840582e+000,
		 -5.9696520e+000,
		 -5.6644732e+000,
		 -5.3684849e+000,
		 -5.0816483e+000,
		 -4.8039217e+000,
		 -4.5352608e+000,
		 -4.2756184e+000,
		 -4.0249442e+000,
		 -3.7831843e+000,
		 -3.5502813e+000,
		 -3.3261738e+000,
		 -3.1107961e+000,
		 -2.9040780e+000,
		 -2.7059444e+000,
		 -2.5163149e+000,
		 -2.3351033e+000,
		 -2.1622175e+000,
		 -1.9975588e+000,
		 -1.8410216e+000,
		 -1.6924928e+000,
		 -1.5518513e+000,
		 -1.4189678e+000,
		 -1.2937038e+000,
		 -1.1759118e+000,
		 -1.0654340e+000,
		 -9.6210282e-001,
		 -8.6573952e-001,
		 -7.7615459e-001,
		 -6.9314718e-001,
		 -6.1650501e-001,
		 -5.4600435e-001,
		 -4.8141016e-001,
		 -4.2247637e-001,
		 -3.6894642e-001,
		 -3.2055397e-001,
		 -2.7702394e-001,
		 -2.3807370e-001,
		 -2.0341461e-001,
		 -1.7275378e-001,
		 -1.4579608e-001,
		 -1.2224636e-001,
		 -1.0181180e-001,
		 -8.4204403e-002,
		 -6.9143456e-002,
		 -5.6357984e-002,
		 -4.5589029e-002,
		 -3.6591704e-002,
		 -2.9136948e-002,
		 -2.3012909e-002,
		 -1.8025916e-002,
		 -1.4001006e-002,
		 -1.0782028e-002,
		 -8.2313205e-003,
		 -6.2290255e-003,
		 -4.6720852e-003,
		 -3.4729977e-003,
		 -2.5584002e-003,
		 -1.8675561e-003,
		 -1.3508100e-003,
		 -9.6807164e-004,
		 -6.8737413e-004,
		 -4.8354103e-004,
		 -3.3698604e-004,
		 -2.3265614e-004,
		 -1.5912125e-004,
		 -1.0780554e-004,
		 -7.2350661e-005,
		 -4.8097501e-005,
		 -3.1671743e-005,
		 -2.0657720e-005,
		 -1.3345838e-005,
		 -8.5399419e-006,
		 -5.4125586e-006,
		 -3.3976789e-006,
		 -2.1124569e-006,
		 -1.3008083e-006,
		 -7.9332847e-007,
		 -4.7918339e-007,
		 -2.8665161e-007
};

double CXLinkOpenScorer::m_lfLogGamaCDF[100] =
{
		/*
		 -1.1618134e+001,
		 -9.1202390e+000,
		 -7.6891650e+000,
		 -6.6947100e+000,
		 -5.9393725e+000,
		 -5.3351777e+000,
		 -4.8351952e+000,
		 -4.4114100e+000,
		 -4.0457503e+000,
		 -3.7258755e+000,
		 -3.4429822e+000,
		 -3.1905693e+000,
		 -2.9636985e+000,
		 -2.7585306e+000,
		 -2.5720209e+000,
		 -2.4017143e+000,
		 -2.2456017e+000,
		 -2.1020184e+000,
		 -1.9695689e+000,
		 -1.8470724e+000,
		 -1.7335201e+000,
		 -1.6280438e+000,
		 -1.5298903e+000,
		 -1.4384024e+000,
		 -1.3530029e+000,
		 -1.2731821e+000,
		 -1.1984877e+000,
		 -1.1285163e+000,
		 -1.0629068e+000,
		 -1.0013343e+000,
		 -9.4350593e-001,
		 -8.8915612e-001,
		 -8.3804378e-001,
		 -7.8994916e-001,
		 -7.4467143e-001,
		 -7.0202660e-001,
		 -6.6184563e-001,
		 -6.2397290e-001,
		 -5.8826482e-001,
		 -5.5458859e-001,
		 -5.2282123e-001,
		 -4.9284858e-001,
		 -4.6456450e-001,
		 -4.3787019e-001,
		 -4.1267349e-001,
		 -3.8888834e-001,
		 -3.6643429e-001,
		 -3.4523599e-001,
		 -3.2522283e-001,
		 -3.0632858e-001,
		 -2.8849099e-001,
		 -2.7165158e-001,
		 -2.5575530e-001,
		 -2.4075032e-001,
		 -2.2658780e-001,
		 -2.1322169e-001,
		 -2.0060853e-001,
		 -1.8870728e-001,
		 -1.7747921e-001,
		 -1.6688768e-001,
		 -1.5689806e-001,
		 -1.4747762e-001,
		 -1.3859536e-001,
		 -1.3022197e-001,
		 -1.2232969e-001,
		 -1.1489224e-001,
		 -1.0788474e-001,
		 -1.0128362e-001,
		 -9.5066552e-002,
		 -8.9212386e-002,
		 -8.3701087e-002,
		 -7.8513674e-002,
		 -7.3632162e-002,
		 -6.9039515e-002,
		 -6.4719593e-002,
		 -6.0657108e-002,
		 -5.6837578e-002,
		 -5.3247289e-002,
		 -4.9873254e-002,
		 -4.6703173e-002,
		 -4.3725406e-002,
		 -4.0928934e-002,
		 -3.8303326e-002,
		 -3.5838715e-002,
		 -3.3525768e-002,
		 -3.1355655e-002,
		 -2.9320028e-002,
		 -2.7410995e-002,
		 -2.5621097e-002,
		 -2.3943285e-002,
		 -2.2370902e-002,
		 -2.0897659e-002,
		 -1.9517620e-002,
		 -1.8225180e-002,
		 -1.7015052e-002,
		 -1.5882246e-002,
		 -1.4822058e-002,
		 -1.3830051e-002,
		 -1.2902044e-002,
		 -1.2034097e-002,
		 */


		 -4.4219216e+001,
		 -3.2149175e+001,
		 -2.5371026e+001,
		 -2.0758557e+001,
		 -1.7331973e+001,
		 -1.4654785e+001,
		 -1.2494092e+001,
		 -1.0710814e+001,
		 -9.2151832e+000,
		 -7.9458485e+000,
		 -6.8589907e+000,
		 -5.9221954e+000,
		 -5.1107846e+000,
		 -4.4055127e+000,
		 -3.7910571e+000,
		 -3.2549983e+000,
		 -2.7871084e+000,
		 -2.3788434e+000,
		 -2.0229733e+000,
		 -1.7133074e+000,
		 -1.4444870e+000,
		 -1.2118275e+000,
		 -1.0111966e+000,
		 -8.3891962e-001,
		 -6.9170693e-001,
		 -5.6659636e-001,
		 -4.6090962e-001,
		 -3.7221857e-001,
		 -2.9831982e-001,
		 -2.3721563e-001,
		 -1.8709952e-001,
		 -1.4634536e-001,
		 -1.1349844e-001,
		 -8.7267675e-002,
		 -6.6517785e-002,
		 -5.0261063e-002,
		 -3.7648115e-002,
		 -2.7957534e-002,
		 -2.0584551e-002,
		 -1.5028930e-002,
		 -1.0882489e-002,
		 -7.8166558e-003,
		 -5.5704853e-003,
		 -3.9394577e-003,
		 -2.7653288e-003,
		 -1.9271765e-003,
		 -1.3337011e-003,
		 -9.1675874e-004,
		 -6.2604699e-004,
		 -4.2482269e-004,
		 -2.8651563e-004,
		 -1.9209531e-004,
		 -1.2805540e-004,
		 -8.4893273e-005,
		 -5.5978502e-005,
		 -3.6721209e-005,
		 -2.3968029e-005,
		 -1.5568120e-005,
		 -1.0064553e-005,
		 -6.4769464e-006,
		 -4.1497701e-006,
		 -2.6473629e-006,
		 -1.6818760e-006,
		 -1.0641885e-006,
		 -6.7071505e-007,
		 -4.2111540e-007,
		 -2.6342349e-007,
		 -1.6418860e-007,
		 -1.0197875e-007,
		 -6.3124185e-008,
		 -3.8944046e-008,
		 -2.3948795e-008,
		 -1.4681156e-008,
		 -8.9723492e-009,
		 -5.4670791e-009,
		 -3.3215480e-009,
		 -2.0123023e-009,
		 -1.2157486e-009,
		 -7.3252338e-010,
		 -4.4020465e-010,
		 -2.6385805e-010,
		 -1.5775936e-010,
		 -9.4092290e-011,
		 -5.5984994e-011,
		 -3.3233194e-011,
		 -1.9682367e-011,
		 -1.1630807e-011,
		 -6.8579586e-012,
		 -4.0349946e-012,
		 -2.3691049e-012,
		 -1.3881118e-012,
		 -8.1179508e-013,
		 -4.7373216e-013,
		 -2.7600144e-013,
		 -1.6042723e-013,
		 -9.3147712e-014,
		 -5.3956839e-014,
		 -3.1197267e-014,
		 -1.7985613e-014,
		 -1.0325074e-014,


		 /*
		 -4.7649148e+001,
		 -3.4951637e+001,
		 -2.7799305e+001,
		 -2.2916487e+001,
		 -1.9276569e+001,
		 -1.6422217e+001,
		 -1.4109428e+001,
		 -1.2192505e+001,
		 -1.0577426e+001,
		 -9.1999473e+000,
		 -8.0142017e+000,
		 -6.9862753e+000,
		 -6.0903661e+000,
		 -5.3063694e+000,
		 -4.6182973e+000,
		 -4.0132096e+000,
		 -3.4804685e+000,
		 -3.0112071e+000,
		 -2.5979411e+000,
		 -2.2342805e+000,
		 -1.9147118e+000,
		 -1.6344311e+000,
		 -1.3892148e+000,
		 -1.1753187e+000,
		 -9.8939805e-001,
		 -8.2844529e-001,
		 -6.8974002e-001,
		 -5.7081010e-001,
		 -4.6940106e-001,
		 -3.8345249e-001,
		 -3.1107994e-001,
		 -2.5056131e-001,
		 -2.0032660e-001,
		 -1.5895012e-001,
		 -1.2514427e-001,
		 -9.7754208e-002,
		 -7.5752617e-002,
		 -5.8234143e-002,
		 -4.4409140e-002,
		 -3.3596513e-002,
		 -2.5215660e-002,
		 -1.8777643e-002,
		 -1.3875805e-002,
		 -1.0176140e-002,
		 -7.4077281e-003,
		 -5.3535194e-003,
		 -3.8417213e-003,
		 -2.7379602e-003,
		 -1.9383289e-003,
		 -1.3633636e-003,
		 -9.5293698e-004,
		 -6.6201571e-004,
		 -4.5719945e-004,
		 -3.1394678e-004,
		 -2.1438607e-004,
		 -1.4561366e-004,
		 -9.8388613e-005,
		 -6.6144575e-005,
		 -4.4250395e-005,
		 -2.9463058e-005,
		 -1.9527108e-005,
		 -1.2884176e-005,
		 -8.4642956e-006,
		 -5.5372603e-006,
		 -3.6076245e-006,
		 -2.3411012e-006,
		 -1.5133514e-006,
		 -9.7459988e-007,
		 -6.2535185e-007,
		 -3.9983129e-007,
		 -2.5475603e-007,
		 -1.6177339e-007,
		 -1.0239120e-007,
		 -6.4599288e-008,
		 -4.0629143e-008,
		 -2.5475714e-008,
		 -1.5926730e-008,
		 -9.9281867e-009,
		 -6.1714437e-009,
		 -3.8256611e-009,
		 -2.3651432e-009,
		 -1.4583675e-009,
		 -8.9693619e-010,
		 -5.5025839e-010,
		 -3.3674963e-010,
		 -2.0559221e-010,
		 -1.2522416e-010,
		 -7.6098239e-011,
		 -4.6140980e-011,
		 -2.7915559e-011,
		 -1.6852852e-011,
		 -1.0152879e-011,
		 -6.1038952e-012,
		 -3.6622927e-012,
		 -2.1930235e-012,
		 -1.3106183e-012,
		 -7.8181905e-013,
		 -4.6551651e-013,
		 -2.7666758e-013,
		 -1.6409096e-013,
		*/

};

CXLinkOpenScorer::CXLinkOpenScorer():m_pPeptide(NULL),m_pSpectrum(NULL),m_pnMatched(NULL),m_pPeakRank(NULL),m_nMaxPeaksNum(5000), m_ionmz(NULL)
{
	m_tIonNum = 0;
	_InitLogNPerm();
	_InitLogTagLenCDF();

	m_bComputeMz = false;
	m_bMiddlePepUsed = false;
	
	m_tionmzSize = 0;
	m_nPepLen = 0;
	m_nMatchedNo = 0;
}

CXLinkOpenScorer::~CXLinkOpenScorer()
{
	Close();
	/*
	if(m_fp)
		fclose(m_fp);
	*/
}

bool SPEC_RANK_LESS(const pair<double,size_t> & rank1, const pair<double,size_t> & rank2)
{
	return rank1.first > rank2.first; 
}


void CXLinkOpenScorer::_SetPeakRank()
{
	if(m_pPeakRank)
		delete []m_pPeakRank;
		
	m_pPeakRank = new size_t[m_pSpectrum->m_tPeaksNum];
	vector<pair<double,size_t> > vRanks;
	
	for(size_t i = 0; i < m_pSpectrum->m_tPeaksNum ; ++i)
	{
		vRanks.push_back(pair<double,size_t>(m_pSpectrum->m_pPeaks[i].lfIntensity,i));
	}
	
	sort(vRanks.begin(),vRanks.end(),SPEC_RANK_LESS);
	for(size_t i = 0; i < m_pSpectrum->m_tPeaksNum ; ++i)
	{
		m_pPeakRank[vRanks[i].second] = i;
	}
}

void CXLinkOpenScorer::_InitLogNPerm()
{
	for(int i = 0 ;i < 100000 ; ++ i )
		CXLinkOpenScorer::m_lfLogNPerm[i] = 0.0;
	
	double lfSum = 0.0;
	for(int i = 1 ;i < 100000 ; ++ i)
	{
		lfSum += log(i);
		CXLinkOpenScorer::m_lfLogNPerm[i] = lfSum;
	}
}

void CXLinkOpenScorer::_InitLogTagLenCDF()
{
	for(int i = 0 ;i < 120 ; ++ i)
	{
		for(int j = 0 ;j < 21; ++ j)
			m_lfLogTagLenCDF[i][j] =  0;
	}
	
	FILE * fp = fopen("logtaglen.txt","r");
	
	if(fp == NULL)
		return ;
	
	char szbuf[1024];
	int nLine = 0;
	int nNumCnt = 0;
	while(fgets(szbuf,1024,fp))
	{
		if(szbuf[strlen(szbuf)-1] == 0x0a)
			szbuf[strlen(szbuf)-1] = 0 ;
		if(szbuf[strlen(szbuf)-1] == 0x0d)
			szbuf[strlen(szbuf)-1] = 0 ;
		
		int nSet = -1;
		bool bNum = false;
		nNumCnt = 0;
		
		for(int i = 0 ;szbuf[i] && i<1024; ++ i)
		{
			if(szbuf[i] == ' ' || szbuf[i] == '	')
			{
				if(bNum == false)
				{
					nSet = i;
					continue;
				}
				szbuf[i] = 0;
				
				double lftmp = atof(&szbuf[nSet+1]);
				
				if(lftmp > 0)
					lftmp = 0;

				if(nLine < 120 && nNumCnt < 21)
				{
					m_lfLogTagLenCDF[nLine][nNumCnt] = lftmp;
					nNumCnt ++ ;
				}
				
				if(nNumCnt >= 21)
					break;

				nSet = i;
				bNum = false;
			}
			else
				bNum = true;
		}
		if(bNum)
		{
			double lftmp = atof(&szbuf[nSet+1]);
			if(lftmp > 0)
				lftmp = 0;

			if(nLine < 120 && nNumCnt < 21)
			{
				m_lfLogTagLenCDF[nLine][nNumCnt] = lftmp;
				nNumCnt ++ ;
			}
		}
		nLine++;
		if(nLine >= 120)
			break;
	}
	fclose(fp);
	

	
}
void CXLinkOpenScorer::Initialize(const CCondition & condition)
{
	_SetCondition(condition);
	m_pnMatched = new int[m_nMaxPeaksNum];
	size_t tNLSize = 0;
	for(size_t i = 0;i < m_Condition.m_vSelectedVarMod.size();++i)
	{
		tNLSize += m_Condition.m_vSelectedVarMod[i].m_vlfAvrgNeutralLoss_dif.size();
	}
	m_ionmz = new CMzTriple[MAX_IONTYPE_NUM * MAX_PEPTIDE_LENGTH * (tNLSize + 1)];
	
	m_bMiddlePepUsed = false;
}

void CXLinkOpenScorer::Close(void)
{
	if(m_pnMatched)
		delete [] m_pnMatched;
	m_pnMatched = NULL;
	if(m_ionmz)
		delete [] m_ionmz;
	m_ionmz = NULL;
	if(m_pPeakRank)
		delete [] m_pPeakRank;
	m_pPeakRank = NULL;
}

void CXLinkOpenScorer::SetSpectrum(CSpectrum & spec)
{
	m_pSpectrum = &spec;
	_SetPeakRank();
}

void CXLinkOpenScorer::SetPeptide(CXLinkOpenPepResult & pep)
{
	m_pPeptide = &pep;
	m_bComputeMz = false;
	m_bMiddlePepUsed = false;
	m_nPepLen = int(m_pPeptide->m_peptide.m_tLength);
}

void CXLinkOpenScorer::SetPeptide(CXLinkOpenMiddlePepResult & pep)
{
	m_pMidPeptide = pep;
	m_bComputeMz = false;
	m_bMiddlePepUsed = true;
	m_nPepLen = int(m_pMidPeptide.m_peptide.m_tLength);
}

double CXLinkOpenScorer::Score()
{
	if(0 == m_pSpectrum->m_tPeaksNum)
	{
		return 0.0;
	}

	if(!m_bComputeMz)
	{
		_ComputeMZ();
		m_bComputeMz = true;
	}
	_Match();
	
	double F = 0.0;
	
	double K = 0.5; 

	int N,n,l,x,n0,x1;
	N = n = l = x = x1 = n0 = 0;

	N = int((m_pSpectrum->m_lfMH)/K);
	
	n = m_tIonNum;
	
	n0 = 10;
	if(n0 > n)
		n0 = n;
	
	int Pp = 0;
	int Bi = -1;
	double Pw = 0.0;
	vector<double> Puw ;
	Puw.clear();
	
	for(size_t i = 0; i < m_pSpectrum->m_tPeaksNum ; ++i )
	{
		if(int(m_pSpectrum->m_pPeaks[i].lfMz/K)>Bi)
		{
			l++;
			if(m_pnMatched[i]>=0)
			{
				x++;
				Pp += m_pPeakRank[i];
				
				double lfTmpTol ;
				lfTmpTol = ((_ShiftMass(m_pnMatched[i])+0.0)/MZMULTIPLIER - m_pSpectrum->m_pPeaks[i].lfMz);
				Puw.push_back(lfTmpTol);
				Pw += lfTmpTol;
				if(m_pPeakRank[i] < n0)
				{
					x1 ++;
				}
			}
			Bi = int(m_pSpectrum->m_pPeaks[i].lfMz/K);
		}
	}

	double Ji = 0.0;

	double Lpc = 0.0;
	double Lp = 0.0; 
	Lp = _CalProb(N, n, l, x);
	
	if(Lp == 0)
		return 0;
	
	Lpc = Lp;

	for(int i=x+1;i<=l;++i)
	{
		Lp = _CalProb(N,n,l,i);
		
		if(Lp != MINUS_INFINITE)
		{
			if(Lpc == MINUS_INFINITE)
				Lpc = Lp;
			else
				Lpc = log(exp(Lpc) + exp(Lp));
		}
	}

	double Lpr = 0.0;
	double Apr = (x == 0 || l < 10 ? 0.5 : (Pp/x)/(l+0.0));
	double Sig = ( x == 0 ? 0.1238 : 0.2887/sqrt(x) );
	
	int Lsnc = int(((Apr - 0.5)/Sig - (-5))/(0.1));
	if(Lsnc < 0)
		Lsnc = 0;
	else if(Lsnc > 100)
		Lsnc = 100;
	Lpr = m_lfLogStdNormCDF[Lsnc]; 

	double Lptl = 0.0;
	double R = GetMaxTagLength();

	int Pl = m_nPepLen - 1;

	if(Pl < 1)
		Pl = 1;
	if(Pl > 120)
		Pl = 120;
	Lptl = m_lfLogTagLenCDF[Pl-1][int(R/0.05)];

	double Lpl = 0.0;
	double Apt = (x == 0 ? 0 : (Pw/x));
	double Pts = 0.0;
	
	for(size_t i = 0;i < Puw.size() ; ++i)
	{
		Pts += (Puw[i] - Apt)*(Puw[i] - Apt);
	}
	if(x > 5)
	{
		Pts /= (x-1);
		Pts = sqrt(Pts);
		 
		int Lgc = int((Pts - 0.01)/0.01);
		if(Lgc < 0)
			Lgc = 0;
		else if(Lgc > 99)
			Lgc = 99;
		
		Lpl = m_lfLogGamaCDF[Lgc];
	}

	Ji = Lpc + Lpr + Lptl;
	
	F = - Ji;
	return F ; 
}


double CXLinkOpenScorer::_CalProb(int N, int n, int n0, int l, int x, int x1)
{
	if( l >= N )

		return 0;
	if(n >= N)

		return 0;
	double lfTmpLogProb = 0 , lfLogProb = 0;
	lfTmpLogProb = _C(n0,x1);
	lfLogProb += lfTmpLogProb;
	
	lfTmpLogProb = _C(n - n0 ,x - x1);
	if(lfTmpLogProb < 0)
		return MINUS_INFINITE;

	lfLogProb += lfTmpLogProb;
	lfTmpLogProb = _C(N-n,l-x);
	
	if(lfTmpLogProb < 0)
		return MINUS_INFINITE;
	lfLogProb += lfTmpLogProb;
	lfTmpLogProb = _C(N,l);
	lfLogProb -= lfTmpLogProb;
	return lfLogProb; 
}

double CXLinkOpenScorer::_CalProb(int N, int n, int l, int x)
{
	if( l >= N )
		return 0;
	if(n >= N)
		return 0;
	double lfTmpLogProb = 0 , lfLogProb = 0;
	lfTmpLogProb = _C(n,x);
	
	if(lfTmpLogProb < 0)
		return MINUS_INFINITE;
	lfLogProb += lfTmpLogProb;
	lfTmpLogProb = _C(N-n,l-x);
	
	if(lfTmpLogProb < 0)
		return MINUS_INFINITE;
	lfLogProb += lfTmpLogProb;
	lfTmpLogProb = _C(N,l);
	lfLogProb -= lfTmpLogProb;
	return lfLogProb; 
}

double CXLinkOpenScorer::_C(int D,int U)
{
	if(D < U)
		return -1;
	else
	{
		double lftmp;
		lftmp = CXLinkOpenScorer::m_lfLogNPerm[D] - CXLinkOpenScorer::m_lfLogNPerm[U] - CXLinkOpenScorer::m_lfLogNPerm[D-U];
		return lftmp;
	}
		
}

double CXLinkOpenScorer::GetMaxTagLength()
{
	int nMaxContiLen = 0;
	int nContiLen = 0;
	for(int i = 0;i < (int)m_Condition.m_vSimpleIonTypes.size();++i)
	{
		const CIonType & IonType = m_Condition.m_vSimpleIonTypes[i];
		if(IonType.cType >= 2)
			continue;
		
		nContiLen = 0;
		for(size_t j = 0;j < m_nPepLen - 1 ; ++j)
		{
			if(m_iontag_pep[i][j]==false)
			{
				if(nContiLen > nMaxContiLen)
					nMaxContiLen = nContiLen;
				nContiLen = 0;
			}
			else
			{
				nContiLen ++ ;
			}
		}
		if(nContiLen > nMaxContiLen)
			nMaxContiLen = nContiLen;
	}
	return (nMaxContiLen+0.0)/(m_nPepLen - 1.0);
}

void CXLinkOpenScorer::_Match()
{
	const vector<CIonType> & It = m_Condition.m_vSimpleIonTypes;
	memset(&m_iontag_pep[0], false, sizeof(m_iontag_pep[0]) * It.size());
	
	if(m_pnMatched && (m_nMaxPeaksNum < m_pSpectrum->m_tPeaksNum) )
	{
		delete [] m_pnMatched;
		m_pnMatched = NULL;
		m_nMaxPeaksNum = m_pSpectrum->m_tPeaksNum;
		m_pnMatched = new int[m_pSpectrum->m_tPeaksNum];
	}
	memset(m_pnMatched, -1, sizeof(int) * m_pSpectrum->m_tPeaksNum);
	
	int Pl = m_nPepLen;
	
	m_tIonNum = 0;
	m_nMatchedNo = 0;
	for(size_t Ppe = 0;
		Ppe < m_tionmzSize;++Ppe)
	{
		if((m_pSpectrum->m_nCharge == 1 && It[m_ionmz[Ppe].nIonTypeOrder].nCharge > m_pSpectrum->m_nCharge)
				 || (m_pSpectrum->m_nCharge > 1 && It[m_ionmz[Ppe].nIonTypeOrder].nCharge >= m_pSpectrum->m_nCharge))
		{
			continue;
		}
		m_tIonNum++;
		int nMz = _ShiftMass(Ppe);
		

		int nMin, nMax;
		_GetMassBorder(nMz, It[m_ionmz[Ppe].nIonTypeOrder].nCharge, nMin, nMax);
		int nTemp = nMin / MZMULTIPLIER;
		
		if(nTemp >= (int)m_pSpectrum->m_vHash.size() || nTemp < 0)
			continue;
		size_t pos_spec = m_pSpectrum->m_vHash[nTemp];

		while(pos_spec < m_pSpectrum->m_tPeaksNum && 
				m_pSpectrum->m_pPeaks[pos_spec].nMz < nMin)
		{
			++pos_spec;

		}

		bool bF = true;
		while(pos_spec < m_pSpectrum->m_tPeaksNum
				&& m_pSpectrum->m_pPeaks[pos_spec].nMz <= nMax)
		{
			if(m_pSpectrum->m_pPeaks[pos_spec].nCharge == 0 || m_pSpectrum->m_pPeaks[pos_spec].nCharge == It[m_ionmz[Ppe].nIonTypeOrder].nCharge)
			{
				m_pnMatched[pos_spec] = int(Ppe);
				if(m_ionmz[Ppe].bAddToTag)
				{
					m_iontag_pep[m_ionmz[Ppe].nIonTypeOrder][m_ionmz[Ppe].nPepPosOrder1] = true;	
				}
				if (bF)
				{
					++m_nMatchedNo;
					bF = false;
				}
			}
			++pos_spec;
		}
	}
}

int CXLinkOpenScorer::_ShiftMass(size_t tIonId)
{
	if(tIonId >= m_tionmzSize || tIonId < 0)
		return 0;
	
	if (!m_bMiddlePepUsed)
		if(!m_pPeptide)
			return 0;
	
	int nPepLen = m_nPepLen;
	int nLinkSite = -1;
	double lfOpenMass = 0.0;
	if (!m_bMiddlePepUsed)
	{
		nLinkSite = m_pPeptide->m_nLinkSite;
		lfOpenMass = m_pPeptide->m_lfOpenMass;
	}
	else
	{
		nLinkSite = m_pMidPeptide.m_nLinkSite;
		lfOpenMass = m_pMidPeptide.m_lfOpenMass;
	}

	const CMzTriple & triple = m_ionmz[tIonId];
	const CIonType & IonType = m_Condition.m_vSimpleIonTypes[triple.nIonTypeOrder];
	
	if(nLinkSite == -1 || lfOpenMass == 0.0)
		return triple.nMz;
	
	int nMz;
	nMz = triple.nMz;
	
	if(IonType.bNTerm)
	{
		if(triple.nPepPosOrder1 >= nLinkSite)
		{
			nMz *= IonType.nCharge;
			nMz += int(lfOpenMass*MZMULTIPLIER);
			nMz /= IonType.nCharge;
		}
	}
	else
	{
		if(nPepLen - 1 - triple.nPepPosOrder1 <= nLinkSite)
		{
			nMz *= IonType.nCharge;
			nMz += lfOpenMass*MZMULTIPLIER;
			nMz /= IonType.nCharge;
		}
	}
	return nMz;
}

void CXLinkOpenScorer::_ComputeMZ()
{
	size_t Pl = 0;
	string S = "";
	size_t Mc = 0;
	vector < vector <int> > Ms;

	Pl = m_nPepLen;
	if (!m_bMiddlePepUsed)
	{
		for (size_t i = 0; i < Pl; ++i)
			S.push_back(m_pPeptide->m_peptide.m_szSequence[i]);
		Mc = m_pPeptide->m_peptide.m_tModCnt;
		for (size_t i = 0; i < Mc; ++i)
		{
			vector <int> vModTemp;
			vModTemp.push_back(m_pPeptide->m_peptide.m_tModSites[i][0]);
			vModTemp.push_back(m_pPeptide->m_peptide.m_tModSites[i][1]);
			Ms.push_back(vModTemp);
		}
	}
	else
	{
		for (size_t i = 0; i < Pl; ++i)
			S.push_back(m_pMidPeptide.m_peptide.m_szSequence[i]);
		Mc = m_pMidPeptide.m_peptide.m_tModCnt;
		for (size_t i = 0; i < Mc; ++i)
		{
			vector <int> vModTemp;
			vModTemp.push_back(m_pMidPeptide.m_peptide.m_tModSites[i][0]);
			vModTemp.push_back(m_pMidPeptide.m_peptide.m_tModSites[i][1]);
			Ms.push_back(vModTemp);
		}
	}

	int Pam[MAX_PEPTIDE_LENGTH] = {0};
	
	int i;
	for(int i=0;i<Pl;++i)
	{
		Pam[i]=m_nAAMass[S[i]-65];
	}
	
	for (size_t i=0;i<Mc;++i)
	{
		Pam[Ms[i][0]] += (int)(m_Condition.m_vSelectedVarMod[Ms[i][1]].m_lfMonoMass_dif*MZMULTIPLIER);
	}
	
	if (m_bMiddlePepUsed)
		Pam[m_pMidPeptide.m_nSidePepSite] += (int)(m_pMidPeptide.m_lfSidePepMass * MZMULTIPLIER);

	int Pam1[MAX_PEPTIDE_LENGTH] = {0};
	int Pam2[MAX_PEPTIDE_LENGTH] = {0};
	
	int Pm = 0;
	for(i = 0;i<Pl;++i)
	{
		Pm += Pam[i];
	}
	
	int Tdm1 = npmass_multi;
	int Tdm2 ;

	Tdm2 = 2*nhmass_mono_multi + nomass_mono_multi + npmass_multi;

	Pam1[0] = (int)(Pam[0] + Tdm1);
	Pam2[0] = (int)(Pam[Pl - 1] + Tdm2);

	for(int i = 1;i < Pl;++i)
	{
		Pam1[i] = Pam1[i - 1] + Pam[i];
	}
	
	for(int i = 1,j = Pl - 2;j > 0;++i,--j)
	{
		Pam2[i] = Pam2[i - 1] + Pam[j];
	}
	
	CMzTriple Tri;
	m_tionmzSize = 0;
	
	int j;
	for(i = 0;i < (int)m_Condition.m_vSimpleIonTypes.size();++i)
	{
		const CIonType & It = m_Condition.m_vSimpleIonTypes[i];
		for(j = 0;j < (int)Pl - 1;++j)
		{
			m_ionmz[m_tionmzSize].nIonTypeOrder = i;
			m_ionmz[m_tionmzSize].bAddToTag = true;
			m_ionmz[m_tionmzSize].nPepPosOrder1 = j;
			m_ionmz[m_tionmzSize].nMz = It.bNTerm?Pam1[j]:Pam2[j];
			m_ionmz[m_tionmzSize].nMz += (It.nCharge - 1) * npmass_multi - It.nTotalLostVal;
			m_ionmz[m_tionmzSize].nMz /= It.nCharge;
			++m_tionmzSize;
		}
		
		for(size_t k = 0;k < Mc;++k)
		{
			for(size_t ww = 0;ww < m_Condition.m_vSelectedVarMod[Ms[k][1]].m_tNLSize;++ww)
			{
				double lfLoss = m_Condition.m_vSelectedVarMod[Ms[k][1]].m_vlfMonoNeutralLoss_dif[ww];
				lfLoss *= MZMULTIPLIER;
				
				for(int j = Ms[k][0];j < (int)Pl - 1;++j)
				{
					Tri.nIonTypeOrder = i;
					Tri.bAddToTag = false;
					Tri.nPepPosOrder1 = j;
					Tri.nMz = It.bNTerm?Pam1[j]:Pam2[j];
					Tri.nMz = int(Tri.nMz-lfLoss);
					Tri.nMz += (It.nCharge - 1) * npmass_multi;
					Tri.nMz -= It.nTotalLostVal;
					if(Tri.nMz < 0)
						continue;
					Tri.nMz /= It.nCharge;
					m_ionmz[m_tionmzSize++] = Tri;
				}
			}
		}
	}
}

void CXLinkOpenScorer::_GetMassBorder(int nMz, int nChg, int &nMin, int &nMax)
{
	double lfTol = m_lfTolMultiplier;
	double lfTolBase = m_lfTolBaseMultiplier;
	if((m_cTolType & 0x01)== 0)
	{
		lfTol *= nMz;
	}
	if((m_cTolBaseType & 0x01)== 0)
	{
		lfTolBase *= nMz;
	}
	nMin = int(nMz - lfTol + lfTolBase);
	nMax = int(nMz + lfTol + lfTolBase);
}

void CXLinkOpenScorer::_SetCondition(const CCondition & cond)
{
	m_Condition = cond;
	
	if ( 0 == (m_Condition.m_strFragmentTolType.compare("%"))) 
	{
		m_cTolType = 0;
		m_lfTolMultiplier = m_Condition.m_lfFragmentTol * 0.01;
	}
	else if ( 0 == (m_Condition.m_strFragmentTolType.compare("mmu")))
	{
		m_cTolType = 1;
		m_lfTolMultiplier = m_Condition.m_lfFragmentTol * MMU_MULTIPLIER;
	}
	else if ( 0 == (m_Condition.m_strFragmentTolType.compare("ppm"))) 
	{
		m_cTolType = 2;
		m_lfTolMultiplier = m_Condition.m_lfFragmentTol * 0.000001;
	}
	else if ( 0 == (m_Condition.m_strFragmentTolType.compare("Da"))) 
	{
		m_cTolType = 3;
		m_lfTolMultiplier = m_Condition.m_lfFragmentTol * MZMULTIPLIER;
	}
	else
	{
		m_cTolType = 3;
		m_lfTolMultiplier = m_Condition.m_lfFragmentTol * MZMULTIPLIER;
	}
	
	if ( 0 == (m_Condition.m_strFragmentTolBaseType.compare("%"))) 
	{
		m_cTolBaseType = 0;
		m_lfTolBaseMultiplier = m_Condition.m_lfFragmentTolBase * 0.01;
	}
	else if ( 0 == (m_Condition.m_strFragmentTolBaseType.compare("mmu")))
	{
		m_cTolBaseType = 1;
		m_lfTolBaseMultiplier = m_Condition.m_lfFragmentTolBase * MMU_MULTIPLIER;
	}
	else if ( 0 == (m_Condition.m_strFragmentTolBaseType.compare("ppm"))) 
	{
		m_cTolBaseType = 2;
		m_lfTolBaseMultiplier = m_Condition.m_lfFragmentTolBase * 0.000001;
	}
	else if ( 0 == (m_Condition.m_strFragmentTolBaseType.compare("Da"))) 
	{
		m_cTolBaseType = 3;
		m_lfTolBaseMultiplier = m_Condition.m_lfFragmentTolBase * MZMULTIPLIER;
	}
	else
	{
		m_cTolBaseType = 3;
		m_lfTolBaseMultiplier = m_Condition.m_lfFragmentTolBase * MZMULTIPLIER;
	}
	
	for(size_t i = 0;i < m_Condition.m_vSelectedFixMod.size();++i)
	{
		m_Condition.m_vSelectedVarMod.push_back(m_Condition.m_vSelectedFixMod[i]);
	}
	CAAConf aa(m_Condition.m_strAAListPath);
	CMapAAMass mapAAMass = aa.GetMapAAMass();
	for(int i = 0;i < 26;++i)
	{
		if(m_Condition.m_bFragmentMono)
			m_nAAMass[i] = (int)(mapAAMass.m_mapAAMass[i + 'A'].m_lfMonoMass * MZMULTIPLIER);
		else
			m_nAAMass[i] = (int)(mapAAMass.m_mapAAMass[i + 'A'].m_lfAvrgMass * MZMULTIPLIER);
	}
}

void CXLinkOpenScorer::_CalMatchInfo()
{
	
}

size_t CXLinkOpenScorer::GetMatchNo()
{
	return m_nMatchedNo;
}

}
