#include <iostream>
//#include "../include/sdk.h"
#include "../include/predefine.h"
#include "../include/interface.h"
#include "../include/option.h"
#include "PSMConf.h"
#include "PSMatch.h"
using namespace std;

//#define _DEBUG2

#define MINUS_INFINITE -10000

double CPSMatch::m_lfLogNPerm[100000] = {0.0};
double CPSMatch::m_lfLogTagLenCDF[120][21]= {0.0};
double CPSMatch::m_lfLogStdNormCDF[101] = {
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

double CPSMatch::m_lfLogGamaCDF[100] =
{
		
		 // for etd 
		 // gama parameter : 3.7093 0.1071
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
		
		// for etd (correct) 
		// gama parameter : 18.3999    0.0138
		/*
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
		 */
		
		 // for cid : 
		 // gama parameter : 19.2796    0.0142
		
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
		
};

                    
CPSMatch::CPSMatch()
{
	m_pbMatched = NULL;
	
	m_pPeptide = NULL;
	m_pSpectrum = NULL;
	m_debugFp = NULL;
	
	m_SingleIonmz = NULL;
	m_tSingleIonMzSize = 0;
	m_pSingleMatched = NULL;
	m_lfSingleMatched = NULL;
	m_pPeakRank = NULL;
	_InitLogNPerm();
	_InitLogTagLenCDF();
}

CPSMatch::~CPSMatch()
{
	Close();
}

bool MzAscend(const pair<double,int> & elem1, const pair<double,int> & elem2)
{
	return elem1.first < elem2.first;
}

void CPSMatch::Preprocess()
{
	vector<pair<double,int> > vPeakMz;
	vPeakMz.clear();
	size_t tPeakNum = m_pSpectrum->m_tPeaksNum;
	int nCharge = m_pSpectrum->m_nCharge;
	for(size_t i = 0 ;i < tPeakNum ; ++i)
	{
		vPeakMz.push_back(pair<double,int>(m_pSpectrum->m_pPeaks[i].lfMz,i));
	}
	
	// ion sorting
	sort(vPeakMz.begin(),vPeakMz.end(),MzAscend);

	// noise peak removal
	vector<double> vNoiseMass;
	vNoiseMass.push_back(108);
	vNoiseMass.push_back(153);
	vNoiseMass.push_back(200);
	int nCurNoiseId = 0;
	for(size_t i = 0 ;i < tPeakNum && nCurNoiseId < vNoiseMass.size(); ++ i)
	{
		if(vPeakMz[i].first < vNoiseMass[nCurNoiseId])
			continue;
		else if(vPeakMz[i].first > vNoiseMass[nCurNoiseId]+1)
		{
			while(vPeakMz[i].first > vNoiseMass[nCurNoiseId]+1)
			{
				nCurNoiseId++;
				if(nCurNoiseId >= vNoiseMass.size())
					break;
			}
		}
		else
		{
			if(int(vPeakMz[i].first) == vNoiseMass[nCurNoiseId])
			{
				m_vPeakInfo[vPeakMz[i].second].nPeakType = 5;
				m_vPeakInfo[vPeakMz[i].second].nCharge = 0;
			}
		}
	}
	
	// MH removal
	vector<double > vMHMass;
	for(int i = nCharge ;i >= 1 ; -- i)
	{
		vMHMass.push_back((m_pSpectrum->m_lfMH + (i-1)*IonMass_Proton)/i);
	}
	
	int nCurChargeId = 0;
	for(size_t i = 0 ;i <tPeakNum ; ++ i)
	{
		if(vPeakMz[i].first < vMHMass[nCurChargeId]-1)
			continue;
		else if(vPeakMz[i].first > vMHMass[nCurChargeId]+1)
		{
			while(vPeakMz[i].first > vMHMass[nCurChargeId]+1 && nCurChargeId < nCharge-1)
			{
				nCurChargeId++;
			}
		}
		else
		{
			if(_IsMatched(vMHMass[nCurChargeId],vPeakMz[i].first))
			{
				m_vPeakInfo[vPeakMz[i].second].nPeakType = 4;
				m_vPeakInfo[vPeakMz[i].second].nCharge = nCharge - nCurChargeId;
			}
		}
	}
	
	
	// gap detection
	// step1 : for isotop ions
	for(size_t i = 0 ;i < tPeakNum ; ++ i)
	{
		for(size_t j = i + 1 ; j < tPeakNum ; ++ j)
		{
			double lfDeltaMz = vPeakMz[j].first - vPeakMz[i].first;
			if(lfDeltaMz > 1.5)
				break;
			for(int k =nCharge ; k >= 1 ; --k)
			{
				if(_IsMatched(vPeakMz[i].first+IonMass_Mono_H/k,vPeakMz[j].first))
				{
					double lfRatio = m_pSpectrum->m_pPeaks[vPeakMz[j].second].lfIntensity / m_pSpectrum->m_pPeaks[vPeakMz[i].second].lfIntensity;
					if(lfRatio <= 2)
					{
						if(m_vPeakInfo[vPeakMz[j].second].nPeakType != 0)
							// peak j has been identified
							continue;
						
						if(m_vPeakInfo[vPeakMz[i].second].nPeakType == 2 || m_vPeakInfo[vPeakMz[i].second].nPeakType == 3 || m_vPeakInfo[vPeakMz[i].second].nPeakType == 5)
							// peak i has been identified with -NH3 or -H2O or noise 
							continue;
						
						if(m_vPeakInfo[vPeakMz[j].second].nCharge != 0 && m_vPeakInfo[vPeakMz[j].second].nCharge != k)
							// peak j has been identified with different charge state
							continue;
						
						if(m_vPeakInfo[vPeakMz[i].second].nCharge != 0 && m_vPeakInfo[vPeakMz[i].second].nCharge != k)
							// peak i has been identified with different charge state
							continue;
						
						// peak j is a isotope of peak i
						m_vPeakInfo[vPeakMz[j].second].nPeakType = 1;
						m_vPeakInfo[vPeakMz[j].second].nMainPeakId = i;
						m_vPeakInfo[vPeakMz[j].second].nCharge = k;
						m_vPeakInfo[vPeakMz[i].second].nCharge = k;
						break;
					}
				}
			}
		}
	}
	// step2 : for -NH3 or -H2O ions
	for(int i = tPeakNum - 1 ;i >= 0 ; --i)
	{
		for(int j = i - 1 ; j >= 0 ; -- j)
		{
			double lfDeltaMz = vPeakMz[i].first - vPeakMz[j].first;
			
			if(lfDeltaMz > H2OMass_Aver + 1.5)
				break;
			
			for(int k = nCharge ; k >= 1 ; --k)
			{
				if(_IsMatched(vPeakMz[j].first+NH3Mass_Mono/k,vPeakMz[i].first))
				{
					double lfRatio = m_pSpectrum->m_pPeaks[vPeakMz[j].second].lfIntensity / m_pSpectrum->m_pPeaks[vPeakMz[i].second].lfIntensity;
					if(lfRatio < 1.5)
					{
						if(m_vPeakInfo[vPeakMz[j].second].nPeakType != 0)
							// peak j has been identified
							continue;
						
						if(m_vPeakInfo[vPeakMz[i].second].nPeakType != 0 && m_vPeakInfo[vPeakMz[i].second].nPeakType != 4)
							// peak i has been identified 
							continue;
						
						if(m_vPeakInfo[vPeakMz[j].second].nCharge != 0 && m_vPeakInfo[vPeakMz[j].second].nCharge != k)
							// peak j has been identified with different charge state
							continue;
						
						if(m_vPeakInfo[vPeakMz[i].second].nCharge != 0 && m_vPeakInfo[vPeakMz[i].second].nCharge != k)
							// peak i has been identified with different charge state
							continue;
						
						// peak j is a isotope of peak i
						m_vPeakInfo[vPeakMz[j].second].nPeakType = 2;
						m_vPeakInfo[vPeakMz[j].second].nMainPeakId = i;
						m_vPeakInfo[vPeakMz[j].second].nCharge = k;
						m_vPeakInfo[vPeakMz[i].second].nCharge = k;
						break;
					}
				}
				else if(_IsMatched(vPeakMz[j].first+H2OMass_Mono/k,vPeakMz[i].first))
				{
					double lfRatio = m_pSpectrum->m_pPeaks[vPeakMz[j].second].lfIntensity / m_pSpectrum->m_pPeaks[vPeakMz[i].second].lfIntensity;
					if(lfRatio < 1.5)
					{
						if(m_vPeakInfo[vPeakMz[j].second].nPeakType != 0)
							// peak j has been identified
							continue;
						
						if(m_vPeakInfo[vPeakMz[i].second].nPeakType != 0 && m_vPeakInfo[vPeakMz[i].second].nPeakType != 4)
							// peak i has been identified 
							continue;
						
						if(m_vPeakInfo[vPeakMz[j].second].nCharge != 0 && m_vPeakInfo[vPeakMz[j].second].nCharge != k)
							// peak j has been identified with different charge state
							continue;
						
						if(m_vPeakInfo[vPeakMz[i].second].nCharge != 0 && m_vPeakInfo[vPeakMz[i].second].nCharge != k)
							// peak i has been identified with different charge state
							continue;
						
						// peak j is a isotope of peak i
						m_vPeakInfo[vPeakMz[j].second].nPeakType = 3;
						m_vPeakInfo[vPeakMz[j].second].nMainPeakId = i;
						m_vPeakInfo[vPeakMz[j].second].nCharge = k;
						m_vPeakInfo[vPeakMz[i].second].nCharge = k;
						break;
					}
				}
			}
		}
	}
	
}

bool CPSMatch::_IsMatched(double lfCalMz,double lfExpMz)
{
	double lfTol = m_psmconf.m_lfFragmentTol;
	if(m_psmconf.m_strFragmentTolType == "ppm")
	{
		lfTol *= 0.000001*lfExpMz;
	}
	else if(m_psmconf.m_strFragmentTolType == "%")
	{
		lfTol *= 0.01*lfExpMz;
	}
	else if(m_psmconf.m_strFragmentTolType == "mmu")
	{
		lfTol *= 0.001*lfExpMz;
	}
	if(lfCalMz <= lfExpMz + lfTol && lfCalMz >= lfExpMz - lfTol)
		return true;
	else
		return false;
}

string CPSMatch::GetIonDescription(int IonId)
{
	if(IonId < 0 || IonId >= int(m_tionmzSize))
		return "";
	
	string strDes = "";
	char szbuf[256];
	
	bool bPepA = true;
	if(m_ionmz[IonId].nPepPosOrder1<m_pPeptide->m_AlphaPeptide.m_tLength)
		bPepA = true;
	else
		bPepA = false;
	
	double lfLostVal = (m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].nTotalLostVal + 0.0)/MZMULTIPLIER;
	int nPepPos;
	
	if(m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].cType == 0)
	{
		if(bPepA)
			strDes += "N-";
		else
			strDes += "n-";
		
		if(lfLostVal == -17.0264 || lfLostVal == 0.984)
		{
		
			if(bPepA)
				strDes += "C";
			else
				strDes += "c";
		}
		else if(lfLostVal == -16.0186 || lfLostVal == 1.0078 || lfLostVal == 1.9918 || lfLostVal == 19.0184)
		{
			
			if(bPepA)
				strDes += "C*";
			else
				strDes += "c*";
		}
		else if(lfLostVal == 0 || lfLostVal == 17.0265 || lfLostVal == 18.0106 || lfLostVal == 35.0371)
		{
		
			if(bPepA)
				strDes += "B";
			else
				strDes += "b";
		}
		else if(lfLostVal == 27.9949 || lfLostVal == 45.0215 || lfLostVal == 46.0055 || lfLostVal == 63.032)
		{
		
			if(bPepA)
				strDes += "A";
			else
				strDes += "a";
		}
		else
		{
			sprintf(szbuf,"[%f]",lfLostVal);
			strDes += szbuf;
		}
	}
	else if(m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].cType == 1)
	{
		if(bPepA)
			strDes += "C-";
		else
			strDes += "c-";
								

		if(lfLostVal == -25.9792 || lfLostVal == -8.9526 || lfLostVal == -7.9686|| lfLostVal == 9.0578)
		{
			if(bPepA)
				strDes += "X";
			else
				strDes += "x";
		}
		else if(lfLostVal == 0 || lfLostVal == 17.0265 || lfLostVal == 18.0106|| lfLostVal == 35.0371)
		{
			if(bPepA)
				strDes += "Y";
			else
				strDes += "y";
		}
		else if(lfLostVal == 15.0109 || lfLostVal == 32.0374 || lfLostVal == 33.0215|| lfLostVal == 50.048)
		{
			// z+H
			if(bPepA)
				strDes += "Z*";
			else
				strDes += "z*";	
		}
		else if(lfLostVal == 16.0187|| lfLostVal == 33.0453|| lfLostVal == 34.0293|| lfLostVal == 51.0558)
		{
			if(bPepA)
				strDes += "Z";
			else
				strDes += "z";
		}
		else
		{
			sprintf(szbuf,"[%f]",lfLostVal);
			strDes += szbuf;
		}
	}
	else if(m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].cType == 2)
	{
		
		strDes += "M-";
		
		bool bPepA1,bPepA2;
		bPepA1 = false;
		bPepA2 = false;
		if(m_ionmz[IonId].nPepPosOrder1 < m_pPeptide->m_AlphaPeptide.m_tLength)
			bPepA1 = true;
		else
			bPepA1 = false;
		
		if(m_ionmz[IonId].nPepPosOrder2 < m_pPeptide->m_AlphaPeptide.m_tLength)
			bPepA2 = true;
		else
			bPepA2 =false;

		if(lfLostVal == -43.0057 || lfLostVal ==-25.9792 || lfLostVal == -24.9951|| lfLostVal == -7.9686)
		{
			// xc
			if(bPepA1)
			{
				if(m_ionmz[IonId].bNTerm1)
					strDes += "C";
				else
					strDes += "X";
				
			}
			else
			{
				if(m_ionmz[IonId].bNTerm1)
					strDes += "c";
				else
					strDes += "x";
			}
		
			
			if(m_ionmz[IonId].nPepPosOrder1>=m_pPeptide->m_AlphaPeptide.m_tLength)
				nPepPos = m_ionmz[IonId].nPepPosOrder1 - m_pPeptide->m_AlphaPeptide.m_tLength;
			else
				nPepPos = m_ionmz[IonId].nPepPosOrder1;
			sprintf(szbuf,"%d",nPepPos);
			strDes += szbuf;
			
			if(bPepA2)
			{
				if(m_ionmz[IonId].bNTerm2)
					strDes += "C";
				else
					strDes += "X";
			}
			else
			{
				if(m_ionmz[IonId].bNTerm2)
					strDes += "c";
				else
					strDes += "x";
			}
	
			if(m_ionmz[IonId].nPepPosOrder2>=m_pPeptide->m_AlphaPeptide.m_tLength)
				nPepPos = m_ionmz[IonId].nPepPosOrder2 - m_pPeptide->m_AlphaPeptide.m_tLength;
			else
				nPepPos = m_ionmz[IonId].nPepPosOrder2;
			sprintf(szbuf,"%d",nPepPos);
			strDes += szbuf;

		}
		else if(lfLostVal == -17.0264 || lfLostVal == 0.9841  )
		{

			// yc
			if(bPepA1)
			{
				if(m_ionmz[IonId].bNTerm1)
					strDes += "C";
				else
					strDes += "Y";
				
			}
			else
			{
				if(m_ionmz[IonId].bNTerm1)
					strDes += "c";
				else
					strDes += "y";
			}
		
			
			if(m_ionmz[IonId].nPepPosOrder1>=m_pPeptide->m_AlphaPeptide.m_tLength)
				nPepPos = m_ionmz[IonId].nPepPosOrder1 - m_pPeptide->m_AlphaPeptide.m_tLength;
			else
				nPepPos = m_ionmz[IonId].nPepPosOrder1;
			sprintf(szbuf,"%d",nPepPos);
			strDes += szbuf;
			
			if(bPepA2)
			{
				if(m_ionmz[IonId].bNTerm2)
					strDes += "C";
				else
					strDes += "Y";
			}
			else
			{
				if(m_ionmz[IonId].bNTerm2)
					strDes += "c";
				else
					strDes += "y";
			}
	
			if(m_ionmz[IonId].nPepPosOrder2>=m_pPeptide->m_AlphaPeptide.m_tLength)
				nPepPos = m_ionmz[IonId].nPepPosOrder2 - m_pPeptide->m_AlphaPeptide.m_tLength;
			else
				nPepPos = m_ionmz[IonId].nPepPosOrder2;
			sprintf(szbuf,"%d",nPepPos);
			strDes += szbuf;
			
		}
		else if(lfLostVal == -1.0077 || lfLostVal == 16.0187 || lfLostVal == 17.0027|| lfLostVal == 34.0293)
		{
			// zc
			if(bPepA1)
			{
				if(m_ionmz[IonId].bNTerm1)
					strDes += "C";
				else
					strDes += "Z";
				
			}
			else
			{
				if(m_ionmz[IonId].bNTerm1)
					strDes += "c";
				else
					strDes += "z";
			}
		
			
			if(m_ionmz[IonId].nPepPosOrder1>=m_pPeptide->m_AlphaPeptide.m_tLength)
				nPepPos = m_ionmz[IonId].nPepPosOrder1 - m_pPeptide->m_AlphaPeptide.m_tLength;
			else
				nPepPos = m_ionmz[IonId].nPepPosOrder1;
			sprintf(szbuf,"%d",nPepPos);
			strDes += szbuf;
			
			if(bPepA2)
			{
				if(m_ionmz[IonId].bNTerm2)
					strDes += "C";
				else
					strDes += "Z";
			}
			else
			{
				if(m_ionmz[IonId].bNTerm2)
					strDes += "c";
				else
					strDes += "z";
			}
	
			if(m_ionmz[IonId].nPepPosOrder2>=m_pPeptide->m_AlphaPeptide.m_tLength)
				nPepPos = m_ionmz[IonId].nPepPosOrder2 - m_pPeptide->m_AlphaPeptide.m_tLength;
			else
				nPepPos = m_ionmz[IonId].nPepPosOrder2;
			sprintf(szbuf,"%d",nPepPos);
			strDes += szbuf;
		
		}
		else if(lfLostVal == 0|| lfLostVal == 17.0265|| lfLostVal ==18.0106 || lfLostVal == 35.0371)
		{
			// yb
		
			if(bPepA1)
			{
				if(m_ionmz[IonId].bNTerm1)
					strDes += "B";
				else
					strDes += "Y";
				
			}
			else
			{
				if(m_ionmz[IonId].bNTerm1)
					strDes += "b";
				else
					strDes += "y";
			}
		
			
			if(m_ionmz[IonId].nPepPosOrder1>=m_pPeptide->m_AlphaPeptide.m_tLength)
				nPepPos = m_ionmz[IonId].nPepPosOrder1 - m_pPeptide->m_AlphaPeptide.m_tLength;
			else
				nPepPos = m_ionmz[IonId].nPepPosOrder1;
			sprintf(szbuf,"%d",nPepPos);
			strDes += szbuf;
			
			if(bPepA2)
			{
				if(m_ionmz[IonId].bNTerm2)
					strDes += "B";
				else
					strDes += "Y";
			}
			else
			{
				if(m_ionmz[IonId].bNTerm2)
					strDes += "b";
				else
					strDes += "y";
			}
	
			if(m_ionmz[IonId].nPepPosOrder2>=m_pPeptide->m_AlphaPeptide.m_tLength)
				nPepPos = m_ionmz[IonId].nPepPosOrder2 - m_pPeptide->m_AlphaPeptide.m_tLength;
			else
				nPepPos = m_ionmz[IonId].nPepPosOrder2;
			sprintf(szbuf,"%d",nPepPos);
			strDes += szbuf;
			
		}
		else if(lfLostVal ==27.9949 || lfLostVal == 45.0214|| lfLostVal == 46.0055 || lfLostVal ==63.032 )
		{
			// ya
		
			if(bPepA1)
			{
				if(m_ionmz[IonId].bNTerm1)
					strDes += "A";
				else
					strDes += "Y";
				
			}
			else
			{
				if(m_ionmz[IonId].bNTerm1)
					strDes += "a";
				else
					strDes += "y";
			}
		
			
			if(m_ionmz[IonId].nPepPosOrder1>=m_pPeptide->m_AlphaPeptide.m_tLength)
				nPepPos = m_ionmz[IonId].nPepPosOrder1 - m_pPeptide->m_AlphaPeptide.m_tLength;
			else
				nPepPos = m_ionmz[IonId].nPepPosOrder1;
			sprintf(szbuf,"%d",nPepPos);
			strDes += szbuf;
			
			if(bPepA2)
			{
				if(m_ionmz[IonId].bNTerm2)
					strDes += "A";
				else
					strDes += "Y";
			}
			else
			{
				if(m_ionmz[IonId].bNTerm2)
					strDes += "a";
				else
					strDes += "y";
			}
	
			if(m_ionmz[IonId].nPepPosOrder2>=m_pPeptide->m_AlphaPeptide.m_tLength)
				nPepPos = m_ionmz[IonId].nPepPosOrder2 - m_pPeptide->m_AlphaPeptide.m_tLength;
			else
				nPepPos = m_ionmz[IonId].nPepPosOrder2;
			sprintf(szbuf,"%d",nPepPos);
			strDes += szbuf;
		}
		else
		{
			sprintf(szbuf,"[%f]",lfLostVal);
			strDes += szbuf;
		}
	}
	else if(m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].cType == 3)
	{
		if(m_ionmz[IonId].nAAnum == m_pPeptide->m_AlphaPeptide.m_tLength)
		{
			sprintf(szbuf,"P-LA-[%f]",lfLostVal);
		}
		else
		{
			if(m_ionmz[IonId].nAAnum == m_pPeptide->m_BetaPeptide.m_tLength)
			{
				sprintf(szbuf,"P-LB-[%f]",lfLostVal);
			}
		}
		strDes += szbuf;
	}
	else if(m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].cType == 4)
	{
		sprintf(szbuf,"Q-MH-");
		strDes += szbuf;
	}
	else if(m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].cType == 5)
	{
		sprintf(szbuf,"P-");
		strDes += szbuf;
		if(m_ionmz[IonId].nPepPosOrder1 < m_pPeptide->m_AlphaPeptide.m_tLength)
		{
			strDes += "Y";
		}
		else
		{
			strDes += "y";
		}
	
		if(m_ionmz[IonId].nPepPosOrder1>=m_pPeptide->m_AlphaPeptide.m_tLength)
			nPepPos = m_ionmz[IonId].nPepPosOrder1 - m_pPeptide->m_AlphaPeptide.m_tLength;
		else
			nPepPos = m_ionmz[IonId].nPepPosOrder1;
		sprintf(szbuf,"%d",nPepPos);
		strDes += szbuf;
		if(lfLostVal)
		{
			sprintf(szbuf,"-[%f]",lfLostVal);
			strDes += szbuf;
		}
	}
	else if(m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].cType == 6)
	{
		sprintf(szbuf,"P-");
		strDes += szbuf;

		if(m_ionmz[IonId].nPepPosOrder2 < m_pPeptide->m_AlphaPeptide.m_tLength)
		{
			strDes += "B";
		}
		else
		{
			strDes += "b";
		}
		if(m_ionmz[IonId].nPepPosOrder2>=m_pPeptide->m_AlphaPeptide.m_tLength)
			nPepPos = m_ionmz[IonId].nPepPosOrder2 - m_pPeptide->m_AlphaPeptide.m_tLength;
		else
			nPepPos = m_ionmz[IonId].nPepPosOrder2;
		sprintf(szbuf,"%d",nPepPos);
		strDes += szbuf;
		if(lfLostVal)
		{
			sprintf(szbuf,"-[%f]",lfLostVal);
			strDes += szbuf;
		}
	}
	if(m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].cType < 2)
	{
		if(m_ionmz[IonId].nPepPosOrder1>=m_pPeptide->m_AlphaPeptide.m_tLength)
			nPepPos = m_ionmz[IonId].nPepPosOrder1 - m_pPeptide->m_AlphaPeptide.m_tLength;
		else
			nPepPos = m_ionmz[IonId].nPepPosOrder1;
		sprintf(szbuf,"%d",nPepPos);
		strDes += szbuf;
	}
	
	for(int i=0;i<m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].nCharge;++i)
	{
		strDes += '+';
	}
	
	if(m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].nLossH2O > 0)
	{
		sprintf(szbuf,"-Loss-%d-H2O",m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].nLossH2O);
		strDes += szbuf;
	}
	if(m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].nLossNH3 > 0)
	{
		sprintf(szbuf,"-Loss-%d-NH3",m_psmconf.m_vIonTypes[m_ionmz[IonId].nIonTypeOrder].nLossNH3);
		strDes += szbuf;
	}
	
	return strDes;
}

void CPSMatch::Initialize(const CPSMConf & psmconf, FILE * fp)
{
	SetCondition(psmconf);
	m_fp = fp; 
	
	m_debugFp = fopen("test.psm.debug.output","w");
	
	if(m_pbMatched)
		delete [] m_pbMatched;
	m_pbMatched = NULL;
	if(m_pSingleMatched)
		delete [] m_pSingleMatched;
	
	if(m_lfSingleMatched)
			delete [] m_lfSingleMatched;
		
	if(m_pPeakRank)
		delete [] m_pPeakRank;
	
	this->m_vMatchedIons.clear();
	
	size_t tNLSize = 0;
	for(size_t i = 0;i < m_psmconf.m_Condition.m_vSelectedVarMod.size();++i)
	{
		tNLSize += m_psmconf.m_Condition.m_vSelectedVarMod[i].m_vlfAvrgNeutralLoss_dif.size();
	}
	
	m_ionmz.clear();
	m_ionmz.resize(10*MAX_IONTYPE_NUM * 2 * MAX_PEPTIDE_LENGTH * (tNLSize + 1));
	//m_ionmz = new CMzTriple[MAX_IONTYPE_NUM * 2 * MAX_PEPTIDE_LENGTH * (tNLSize + 1)];
	if(m_SingleIonmz)
		delete [] m_SingleIonmz;
	m_SingleIonmz = new CMzTriple[10*MAX_IONTYPE_NUM * 2 * MAX_PEPTIDE_LENGTH * (tNLSize + 1)];
	
	m_vIonTypeMatchInfo.clear();
	m_vIonTypeMatchInfo.resize(10*MAX_IONTYPE_NUM);
	
	m_vXlinkIonTypeMatchInfo.clear();
	m_vXlinkIonTypeMatchInfo.resize(10*MAX_IONTYPE_NUM);
	
	m_vCommonIonTypeMatchInfo.clear();
	m_vCommonIonTypeMatchInfo.resize(10*MAX_IONTYPE_NUM);
	
	m_vPeakInfo.clear();

}

void CPSMatch::_InitLogNPerm()
{
	for(int i = 0 ;i < 100000 ; ++ i )
		CPSMatch::m_lfLogNPerm[i] = 0.0;
	
	double lfSum = 0.0;
	for(int i = 1 ;i < 100000 ; ++ i)
	{
		lfSum += log(i);
		CPSMatch::m_lfLogNPerm[i] = lfSum;
	}
}

void CPSMatch::_InitLogTagLenCDF()
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
	/*
	for(int i = 0 ;i < 120 ; ++ i)
	{
		for(int j = 0 ;j < 21; ++ j)
			cout << m_lfLogTagLenCDF[i][j] << "	";
		cout << endl;
		char tmpch;
		cin >> tmpch;
	}
	*/	
	
}

void CPSMatch::Close(void)
{
	if(m_pbMatched)
		delete [] m_pbMatched;
	m_pbMatched = NULL;
	
	fclose(m_debugFp);
}

void CPSMatch::SetCondition(const CPSMConf & psmconf)
{
	m_psmconf = psmconf;
	//预计算误差扩大倍数
	// modify by emily
	// for fragment tolerance
	if ( 0 == (m_psmconf.m_strFragmentTolType.compare("%"))) 
	{
		m_cTolType = 0;
		m_lfTolMultiplier = m_psmconf.m_lfFragmentTol * 0.01;
	}
	else if ( 0 == (m_psmconf.m_strFragmentTolType.compare("mmu")))
	{
		m_cTolType = 1;
		m_lfTolMultiplier = m_psmconf.m_lfFragmentTol * MMU_MULTIPLIER;
	}
	else if ( 0 == (m_psmconf.m_strFragmentTolType.compare("ppm"))) 
	{
		m_cTolType = 2;
		m_lfTolMultiplier = m_psmconf.m_lfFragmentTol * 0.000001;
	}
	else if ( 0 == (m_psmconf.m_strFragmentTolType.compare("Da"))) 
	{
		m_cTolType = 3;
		m_lfTolMultiplier = m_psmconf.m_lfFragmentTol * MZMULTIPLIER;
	}
	else
	{
		m_cTolType = 3;
		m_lfTolMultiplier = m_psmconf.m_lfFragmentTol * MZMULTIPLIER;
	}
		
	
	//append the modification information
	for(size_t i = 0;i < m_psmconf.m_Condition.m_vSelectedFixMod.size();++i)
	{
		m_psmconf.m_Condition.m_vSelectedVarMod.push_back(m_psmconf.m_Condition.m_vSelectedFixMod[i]);
	}

	//init the aa mass;
	CAAConf aa(m_psmconf.m_Condition.m_strAAListPath);
	CMapAAMass mapAAMass = aa.GetMapAAMass();
	for(int i = 0;i < 26;++i)
	{
		m_nAAMass[i] = (int)(mapAAMass.m_mapAAMass[i + 'A'].m_lfMonoMass * MZMULTIPLIER);
	}
}

void CPSMatch::SetSpectrum(CSpectrum & spec)
{
	m_pSpectrum = &spec;
	_SetPeakRank();
	m_vPeakInfo.clear();
	for(size_t i = 0 ;i < m_pSpectrum->m_tPeaksNum ; ++i)
	{
		m_vPeakInfo.push_back(PeakInfo());
	}
}

void CPSMatch::SetPeptide(CXLinkPepResult & pep)
{
	m_pPeptide = &pep;
}

void CPSMatch::_GetMassBorder(int nMz, int nChg, int &nMin, int &nMax)
{
	double lfTol = m_lfTolMultiplier;
	
	if((m_cTolType & 0x01)== 0)
	{
		lfTol *= nMz ;
	}
	
	nMin = int(nMz - lfTol);
	nMax = int(nMz + lfTol);
}

double CPSMatch::_GetIonTol(size_t tSpecId,size_t tIonId)
{
	double lfExpMz = m_pSpectrum->m_pPeaks[tSpecId].lfMz;
	double lfCalMz = (m_ionmz[m_vMatchedIons[tSpecId][tIonId]].nMz + 0.0)/MZMULTIPLIER;
	
	double lfTol;
	if ( 0 == m_cTolType) 
	{
		lfTol =  (lfCalMz - lfExpMz)/lfCalMz * 100;
	}
	else if ( 1 == m_cTolType)
	{
		lfTol =  (lfCalMz - lfExpMz)/lfCalMz * 1000;
	}
	else if ( 2 == m_cTolType) 
	{
		lfTol = (lfCalMz - lfExpMz)/lfCalMz * 1000000;
	}
	else
	{
		lfTol = (lfCalMz - lfExpMz);
	}
	if(lfTol < 0)
		return -lfTol;
	else
		return lfTol;
			
}

bool CPSMatch::_HasIsotope(size_t tSpecId, size_t tIonId)
{
	if(tSpecId >= m_pSpectrum->m_tPeaksNum)
		return false;
	double lfCurMz = m_pSpectrum->m_pPeaks[tSpecId].lfMz;
	double lfTol = m_psmconf.m_lfFragmentTol;
	int nCharge = m_psmconf.m_vIonTypes[m_ionmz[m_vMatchedIons[tSpecId][tIonId]].nIonTypeOrder].nCharge;
	double lfCalMz = lfCurMz + (IonMass_Mono_H)/(nCharge+0.0);

	if ( 0 == m_cTolType) 
	{
		lfTol /= 100;
		lfTol *= lfCalMz;
	}
	else if (1 == m_cTolType)
	{
		lfTol /= 1000;
		lfTol *= lfCalMz;
	}
	else if (2 == m_cTolType) 
	{
		lfTol /= 1000000;
		lfTol *= lfCalMz;
	}
	
	double lfMzLB = lfCalMz - lfTol;
	double lfMzUB = lfCalMz + lfTol;

	
	size_t i;
	for(i = tSpecId + 1; i < m_pSpectrum->m_tPeaksNum && m_pSpectrum->m_pPeaks[i].lfMz >= lfMzLB; ++i);
	if(i >= m_pSpectrum->m_tPeaksNum)
		return false;
	if(m_pSpectrum->m_pPeaks[i].lfMz <= lfMzUB)
	{
		return true;
	}
	else
		return false;
}


double CPSMatch::GetMaxTagLength(int nPep)
{
	// 对普通离子计算谱峰连续性
	/*
	int nMaxContiLen = 0;
	int nContiLen = 0;
	for(int i = 0;i < (int)m_psmconf.m_vIonTypes.size();++i)
	{
		const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[i];
		if(IonType.cType >= 2)
			continue;
		
		nContiLen = 0;
		if(nPep == 0)
		{
			// for alpha peptide
			for(size_t j = 0;j < m_pPeptide->m_AlphaPeptide.m_tLength - 1 ; ++j)
			{
				if(m_single_iontag_pep1[i][j]==false)
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
		else
		{
			// for beta peptide
			nContiLen = 0;
			if(m_pPeptide->m_bPair)
			{
				for(size_t j = 0;j < m_pPeptide->m_BetaPeptide.m_tLength - 1 ; ++j)
				{
					if(m_single_iontag_pep2[i][j]==false)
					{
						if(nContiLen > nMaxContiLen)
							nMaxContiLen = nContiLen;
						nContiLen = 0;
					}
					else
					{
						nContiLen ++;
					}
				}
				if(nContiLen > nMaxContiLen)
					nMaxContiLen = nContiLen;
			}
		}
	}
	
	if(nPep == 0)
	{
		return (nMaxContiLen+0.0)/(m_pPeptide->m_AlphaPeptide.m_tLength - 1.0);
	}
	else
	{
		if(m_pPeptide->m_bPair)
			return (nMaxContiLen+0.0)/(m_pPeptide->m_BetaPeptide.m_tLength - 1.0);
	}
	*/
}

void CPSMatch::GetPSMMatchInfo()
{
	// 统计每种类型的离子的重要性
	// 初始化
	// 倒数第二种： KL型：y_a_-NH3
	// 特征： 
	//1. iontypeorder<=215 && iontypeorder>=210，bContainLinker=true, bSameSite=true
	//2. iontypeorder<=59 && iontypeorder>=54,bContainLinker=true, pepOrder-nterm = 0
	// 最后一种：La/Lb
	// 特征：
	//1. iontypeorder<=299 && iontypeorder >= 288 
	//2. iontypeorder<=365 && iontypeorder >= 360
	
	int nLastIonTypeOrder1 = m_psmconf.m_vIonTypes.size();
	int nLastIonTypeOrder2 = nLastIonTypeOrder1 + 1;
	
	// 遍历所有匹配的离子
	bool bExist = false;
	for(size_t i = 0;i < m_vIonMatchInfo.size(); ++ i)
	{
		struct IonMatchInfoEx & matchIon = m_vIonMatchInfo[i];
		
		if(matchIon.bContainLinker)
		{
			m_vXlinkIonTypeMatchInfo[matchIon.nIonTypeOrder].nMatchIonCount++;
			m_vXlinkIonTypeMatchInfo[matchIon.nIonTypeOrder].lfMatchAvrgItensity += matchIon.lfIntensity;
			m_vXlinkIonTypeMatchInfo[matchIon.nIonTypeOrder].lfMatchAvrgTolerance += matchIon.lfTol;
			m_vXlinkIonTypeMatchInfo[matchIon.nIonTypeOrder].lfAvrgAAnum += matchIon.nAAnum;
		}
		else
		{
			m_vCommonIonTypeMatchInfo[matchIon.nIonTypeOrder].nMatchIonCount++;
			m_vCommonIonTypeMatchInfo[matchIon.nIonTypeOrder].lfMatchAvrgItensity += matchIon.lfIntensity;
			m_vCommonIonTypeMatchInfo[matchIon.nIonTypeOrder].lfMatchAvrgTolerance += matchIon.lfTol;
			m_vCommonIonTypeMatchInfo[matchIon.nIonTypeOrder].lfAvrgAAnum += matchIon.nAAnum;
		}
		
		m_vIonTypeMatchInfo[matchIon.nIonTypeOrder].nMatchIonCount++;
		m_vIonTypeMatchInfo[matchIon.nIonTypeOrder].lfMatchAvrgItensity += matchIon.lfIntensity;
		m_vIonTypeMatchInfo[matchIon.nIonTypeOrder].lfMatchAvrgTolerance += matchIon.lfTol;
		m_vIonTypeMatchInfo[matchIon.nIonTypeOrder].lfAvrgAAnum += matchIon.nAAnum;

		if( (matchIon.nIonTypeOrder >= 210 && matchIon.nIonTypeOrder <= 215 && matchIon.bContainLinker &&  matchIon.bSameSite) 
				|| (matchIon.nIonTypeOrder >= 54 && matchIon.nIonTypeOrder <= 59 && matchIon.bContainLinker && matchIon.nPepPosOrder1 == 0) )
		{
			m_vIonTypeMatchInfo[nLastIonTypeOrder1].nMatchIonCount++;
			m_vIonTypeMatchInfo[nLastIonTypeOrder1].lfMatchAvrgItensity += matchIon.lfIntensity;
			m_vIonTypeMatchInfo[nLastIonTypeOrder1].lfMatchAvrgTolerance += matchIon.lfTol;
			m_vIonTypeMatchInfo[nLastIonTypeOrder1].lfAvrgAAnum += matchIon.nAAnum;
		}

		if( (matchIon.nIonTypeOrder >= 288 && matchIon.nIonTypeOrder <= 299) 
				|| (matchIon.nIonTypeOrder >= 360 && matchIon.nIonTypeOrder <= 365))
		{
			m_vIonTypeMatchInfo[nLastIonTypeOrder2].nMatchIonCount++;
			m_vIonTypeMatchInfo[nLastIonTypeOrder2].lfMatchAvrgItensity += matchIon.lfIntensity;
			m_vIonTypeMatchInfo[nLastIonTypeOrder2].lfMatchAvrgTolerance += matchIon.lfTol;
			m_vIonTypeMatchInfo[nLastIonTypeOrder2].lfAvrgAAnum += matchIon.nAAnum;
		}
		
		if(matchIon.nCharge == m_pSpectrum->m_nCharge)
		{
			if(matchIon.nIonType == 0 || matchIon.nIonType == 1)
			{
				fprintf(m_debugFp,"%s\t",m_pSpectrum->m_strFilePath.c_str());
				fprintf(m_debugFp,"%d\t%d\t%d\t%f\t%d\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t\n",
						matchIon.nIonTypeOrder + 1,
						matchIon.nIonType,
						matchIon.bContainLinker,
						matchIon.lfIntensity,
						matchIon.nRank,
						matchIon.lfCalMz,
						matchIon.lfExpMz,
						matchIon.lfTol,
						matchIon.nCharge,
						matchIon.nLossNH3,
						matchIon.nLossH2O,
						matchIon.nPepPosOrder1,
						matchIon.nPepPosOrder2,
						matchIon.bSameSite,
						matchIon.lfOtherLoss,
						matchIon.lfExpTotalLoss);
			}
		}
	}
	
	// 计算所有信息
	int nTotalMatchIonCount = m_vIonMatchInfo.size();
	for(size_t i = 0;i <= nLastIonTypeOrder2; ++i)
	{
		m_vIonTypeMatchInfo[i].nIonTypeOrder = i;
		if(m_vIonTypeMatchInfo[i].nMatchIonCount > 0)
		{
			m_vIonTypeMatchInfo[i].lfMatchCountRatio = (m_vIonTypeMatchInfo[i].nMatchIonCount+0.0)/nTotalMatchIonCount;
			m_vIonTypeMatchInfo[i].lfMatchGainRatio = (m_vIonTypeMatchInfo[i].nMatchIonCount+0.0)/m_vIonTypeMatchInfo[i].nTheoIonCount;
			m_vIonTypeMatchInfo[i].lfMatchAvrgItensity /= m_vIonTypeMatchInfo[i].nMatchIonCount ;
			m_vIonTypeMatchInfo[i].lfMatchAvrgTolerance /= m_vIonTypeMatchInfo[i].nMatchIonCount ;
			m_vIonTypeMatchInfo[i].lfAvrgAAnum /= m_vIonTypeMatchInfo[i].nMatchIonCount;
		}
		
		m_vIonTypeMatchInfo[i].lfSignificance = m_vIonTypeMatchInfo[i].lfMatchCountRatio*m_vIonTypeMatchInfo[i].lfMatchGainRatio*m_vIonTypeMatchInfo[i].lfMatchAvrgItensity;
	}
	
	// for xlink ion type
	for(size_t i = 0;i < m_psmconf.m_vIonTypes.size(); ++i)
	{
		m_vXlinkIonTypeMatchInfo[i].nIonTypeOrder = i;
		if(m_vXlinkIonTypeMatchInfo[i].nMatchIonCount > 0)
		{
			m_vXlinkIonTypeMatchInfo[i].lfMatchCountRatio = (m_vXlinkIonTypeMatchInfo[i].nMatchIonCount+0.0)/nTotalMatchIonCount;
			m_vXlinkIonTypeMatchInfo[i].lfMatchGainRatio = (m_vXlinkIonTypeMatchInfo[i].nMatchIonCount+0.0)/m_vXlinkIonTypeMatchInfo[i].nTheoIonCount;
			m_vXlinkIonTypeMatchInfo[i].lfMatchAvrgItensity /= m_vXlinkIonTypeMatchInfo[i].nMatchIonCount ;
			m_vXlinkIonTypeMatchInfo[i].lfMatchAvrgTolerance /= m_vXlinkIonTypeMatchInfo[i].nMatchIonCount ;
			m_vXlinkIonTypeMatchInfo[i].lfAvrgAAnum /= m_vXlinkIonTypeMatchInfo[i].nMatchIonCount;
		}
		
		m_vXlinkIonTypeMatchInfo[i].lfSignificance = m_vXlinkIonTypeMatchInfo[i].lfMatchCountRatio*m_vXlinkIonTypeMatchInfo[i].lfMatchGainRatio*m_vXlinkIonTypeMatchInfo[i].lfMatchAvrgItensity;
	}
	
	// for common ion type
	for(size_t i = 0;i < m_psmconf.m_vIonTypes.size(); ++i)
	{
		m_vCommonIonTypeMatchInfo[i].nIonTypeOrder = i;
		if(m_vCommonIonTypeMatchInfo[i].nMatchIonCount > 0)
		{
			m_vCommonIonTypeMatchInfo[i].lfMatchCountRatio = (m_vCommonIonTypeMatchInfo[i].nMatchIonCount+0.0)/nTotalMatchIonCount;
			m_vCommonIonTypeMatchInfo[i].lfMatchGainRatio = (m_vCommonIonTypeMatchInfo[i].nMatchIonCount+0.0)/m_vCommonIonTypeMatchInfo[i].nTheoIonCount;
			m_vCommonIonTypeMatchInfo[i].lfMatchAvrgItensity /= m_vCommonIonTypeMatchInfo[i].nMatchIonCount ;
			m_vCommonIonTypeMatchInfo[i].lfMatchAvrgTolerance /= m_vCommonIonTypeMatchInfo[i].nMatchIonCount ;
			m_vCommonIonTypeMatchInfo[i].lfAvrgAAnum /= m_vCommonIonTypeMatchInfo[i].nMatchIonCount;
		}
		
		m_vCommonIonTypeMatchInfo[i].lfSignificance = m_vCommonIonTypeMatchInfo[i].lfMatchCountRatio*m_vCommonIonTypeMatchInfo[i].lfMatchGainRatio*m_vCommonIonTypeMatchInfo[i].lfMatchAvrgItensity;
	}

	// 对b/y离子计算谱峰连续性
	for(int i = 0;i < (int)m_psmconf.m_vIonTypes.size();++i)
	{
		const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[i];
		if(IonType.cType >= 2)
			continue;
		
		int nContiLen = 0;
		int nMaxContiLen = 0;
		// for alpha peptide
		for(size_t j = 0;j < m_pPeptide->m_AlphaPeptide.m_tLength - 1 ; ++j)
		{
			if(m_iontag_pep1[i][j]<0)
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
		
		m_vIonTypeMatchInfo[i].lfMatchContiRatio = (nMaxContiLen+0.0) / (m_pPeptide->m_AlphaPeptide.m_tLength - 1);
		
		// for beta peptide
		nContiLen = 0;
		nMaxContiLen = 0;
		if(m_pPeptide->m_bPair)
		{
			for(size_t j = 0;j < m_pPeptide->m_BetaPeptide.m_tLength - 1 ; ++j)
			{
				if(m_iontag_pep2[i][j]<0)
				{
					if(nContiLen > nMaxContiLen)
						nMaxContiLen = nContiLen;
					nContiLen = 0;
				}
				else
				{
					nContiLen ++;
				}
			}
			if(nContiLen > nMaxContiLen)
				nMaxContiLen = nContiLen;
			
			m_vIonTypeMatchInfo[i].lfMatchContiRatio += (nMaxContiLen+0.0) / (m_pPeptide->m_BetaPeptide.m_tLength - 1);
			m_vIonTypeMatchInfo[i].lfMatchContiRatio /= 2;
		}
	}
	
	// for xlink ion and common ion
	for(int i = 0;i < (int)m_psmconf.m_vIonTypes.size();++i)
	{
		const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[i];
		if(IonType.cType >= 2)
			continue;
		
		int nContiLen1 = 0;
		int nMaxContiLen1 = 0;
		int nContiLen2 = 0;
		int nMaxContiLen2 = 0;
		int nLen1 = 0;
		int nLen2 = 0;
		size_t tSite = 0;
		if(IonType.cType == 0)
		{
			// N-term
			tSite = m_pPeptide->m_XLink.m_tAlphaSite;
		}
		else
		{
			// C-term
			tSite = m_pPeptide->m_AlphaPeptide.m_tLength - 1 - m_pPeptide->m_XLink.m_tAlphaSite;
		}
		// for alpha peptide

		for(size_t j = 0;j < m_pPeptide->m_AlphaPeptide.m_tLength - 1 ; ++j)
		{
			if(j < tSite)
			{
				// common ion
				nLen1 ++;
				if(m_iontag_pep1[i][j]<0)
				{
					if(nContiLen1 > nMaxContiLen1)
						nMaxContiLen1 = nContiLen1;
					nContiLen1 = 0;
				}
				else
				{
					nContiLen1 ++ ;
				}
			}
			else
			{
				// xlink ion
				nLen2 ++ ;
				if(m_iontag_pep1[i][j]<0)
				{
					if(nContiLen2 > nMaxContiLen2)
						nMaxContiLen2 = nContiLen2;
					nContiLen2 = 0;
				}
				else
				{
					nContiLen2 ++ ;
				}
			}
		}
		if(nContiLen1 > nMaxContiLen1)
			nMaxContiLen1 = nContiLen1;
		if(nContiLen2 > nMaxContiLen2)
			nMaxContiLen2 = nContiLen2;
		
		if(nLen1 > 0)
			m_vCommonIonTypeMatchInfo[i].lfMatchContiRatio = (nMaxContiLen1+0.0) / (nLen1);
		if(nLen2 > 0)
			m_vXlinkIonTypeMatchInfo[i].lfMatchContiRatio = (nMaxContiLen2+0.0) / (nLen2);
		
		// for beta peptide
		
		nContiLen1 = 0;
		nMaxContiLen1 = 0;
		nContiLen2 = 0;
		nMaxContiLen2 = 0;
		nLen1 = 0;
		nLen2 = 0;
		if(IonType.cType == 0)
		{
			// N-term
			tSite = m_pPeptide->m_XLink.m_tBetaSite;
		}
		else
		{
			// C-term
			tSite = m_pPeptide->m_BetaPeptide.m_tLength - 1 - m_pPeptide->m_XLink.m_tBetaSite;
		}
		if(m_pPeptide->m_bPair)
		{
			for(size_t j = 0;j < m_pPeptide->m_BetaPeptide.m_tLength - 1 ; ++j)
			{
				if(j < tSite)
				{
					// common
					nLen1 ++;
					if(m_iontag_pep2[i][j]<0)
					{
						if(nContiLen1 > nMaxContiLen1)
							nMaxContiLen1 = nContiLen1;
						nContiLen1 = 0;
					}
					else
					{
						nContiLen1 ++;
					}
				}
				else
				{
					// xlink
					nLen2 ++;
					if(m_iontag_pep2[i][j]<0)
					{
						if(nContiLen2 > nMaxContiLen2)
							nMaxContiLen2 = nContiLen2;
						nContiLen2 = 0;
					}
					else
					{
						nContiLen2 ++;
					}
				}
			}
			if(nContiLen1 > nMaxContiLen1)
				nMaxContiLen1 = nContiLen1;
			if(nContiLen2 > nMaxContiLen2)
				nMaxContiLen2 = nContiLen2;

			if(nLen1 > 0)
				m_vCommonIonTypeMatchInfo[i].lfMatchContiRatio += (nMaxContiLen1+0.0) / (nLen1);
			if(nLen2 > 0)
				m_vXlinkIonTypeMatchInfo[i].lfMatchContiRatio += (nMaxContiLen2+0.0) / (nLen2);
			m_vCommonIonTypeMatchInfo[i].lfMatchContiRatio /= 2;
			m_vXlinkIonTypeMatchInfo[i].lfMatchContiRatio /= 2;
		}
	}
	
}

void CPSMatch::_CalMatchInfo()
{
	// initialize
	m_pPeptide->m_stMatchInfo.lfMatchedSpecInt = 0.0;
	m_pPeptide->m_stMatchInfo.lfUnMatchedSpecInt = 0.0;
	memset(m_pPeptide->m_stMatchInfo.aPepConf,0,sizeof(int)*2*MAX_PEPTIDE_LENGTH);
	
	// calculate the confidence for each amino of the peptide
	int aCleavConf[MAX_PEPTIDE_LENGTH+1];
	memset(aCleavConf,0,sizeof(int)*(1+MAX_PEPTIDE_LENGTH));
	aCleavConf[0] = 1;
	aCleavConf[m_pPeptide->m_AlphaPeptide.m_tLength] = 1;
	
	for(int i = 0;i < (int)m_psmconf.m_vIonTypes.size();++i)
	{
		const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[i];
		for(int j = 0;j < m_pPeptide->m_AlphaPeptide.m_tLength - 1 ; ++j)
		{
			if(m_iontag_pep1[i][j]<0)
			{
				continue;
			}
			if(IonType.cType == 0)
			{
				aCleavConf[j+1] ++ ;
			}
			else if(IonType.cType == 1)
			{
				aCleavConf[m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j] ++;
			}
		}
	}
	
	for(int i=0;i < m_pPeptide->m_AlphaPeptide.m_tLength ; ++i)
	{
		m_pPeptide->m_stMatchInfo.aPepConf[i] = (aCleavConf[i] <= aCleavConf[i+1])?aCleavConf[i]:aCleavConf[i+1];
	}
	
	if(m_pPeptide->m_bPair)
	{
		memset(aCleavConf,0,sizeof(int)*(MAX_PEPTIDE_LENGTH+1));
		aCleavConf[0] = 1;
		aCleavConf[m_pPeptide->m_BetaPeptide.m_tLength] = 1;
		
		for(int i = 0;i < (int)m_psmconf.m_vIonTypes.size();++i)
		{
			const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[i];
			for(int j = 0;j < m_pPeptide->m_BetaPeptide.m_tLength - 1 ; ++j)
			{
				if(m_iontag_pep2[i][j]<0)
				{
					continue;
				}
				if(IonType.cType == 0)
				{
					aCleavConf[j+1] ++ ;
				}
				else if(IonType.cType == 1)
				{
					aCleavConf[m_pPeptide->m_BetaPeptide.m_tLength - 1 - j] ++;
				}
			}
		}
		
		for(int i=0;i < m_pPeptide->m_BetaPeptide.m_tLength ; ++i)
		{
			m_pPeptide->m_stMatchInfo.aPepConf[m_pPeptide->m_AlphaPeptide.m_tLength + i] = (aCleavConf[i] <= aCleavConf[i+1])?aCleavConf[i]:aCleavConf[i+1];
		}
	}

	// calculate the intensity sum for matched and unmatched ions 
	for(int i=0;i<m_pSpectrum->m_tPeaksNum;++i)
	{
		if ( m_pbMatched[i] )
		{
			m_pPeptide->m_stMatchInfo.lfMatchedSpecInt += m_pSpectrum->m_pPeaks[i].lfIntensity;
		}
		else
		{
			m_pPeptide->m_stMatchInfo.lfUnMatchedSpecInt += m_pSpectrum->m_pPeaks[i].lfIntensity;
		}
	}
	
}
double CPSMatch::GetMatchOdd(int nPepId)
{
	if(nPepId < 0 || nPepId > 1)
		return 0.0;

	double lfBinWidth = 0.5; // bin width is 0.5 Da

	int N,n,l,x,n0,x1;
	N = n = l = x = x1 = n0 = 0;

	N = int((m_pSpectrum->m_lfMH)/lfBinWidth);
	
	if(nPepId == 0)
	{
		n = m_tSingleIonNum1;
	}
	else
	{
		n = m_tSingleIonNum2;
	}
	
	int nLookfor = 0;
	if(nPepId == 0)
		nLookfor = 1;
	else
		nLookfor = -1;
	
	int nCurBinId = -1;
	
	// main peak is the top ten peaks
	n0 = 10;
	if(n0 > n)
		n0 = n;
	
	int nPeakRank = 0;
	double lfPeakTol = 0.0;
	vector<double> vPeakTol;
	vPeakTol.clear();
	
	for(size_t i = 0; i < m_pSpectrum->m_tPeaksNum ; ++i )
	{
		if(int(m_pSpectrum->m_pPeaks[i].lfMz/lfBinWidth)>nCurBinId)
		{
			l++;
			if(m_pSingleMatched[i] == nLookfor || m_pSingleMatched[i] == 2)
			{
				x++;
				nPeakRank += m_pPeakRank[i];
				lfPeakTol += m_lfSingleMatched[i];
				vPeakTol.push_back(m_lfSingleMatched[i]);
				if(m_pPeakRank[i] < n0)
				{
					x1 ++;
				}
			}
			nCurBinId = int(m_pSpectrum->m_pPeaks[i].lfMz/lfBinWidth);
		}
	}
	
	double lfOdd = 0.0;
	double lfLogProb = 0.0; 
	
	//lfLogProb = _CalProb(N, n, n0, l, x,x1);
	lfLogProb = _CalProb(N, n, l, x);
	
	if(lfLogProb == 0)
		return 0;
	
	
	// mu = 0.5578 sigma = 0.1238
	double lfAvgPeakRank = (x == 0 || l < 10 ? 0.5 : (nPeakRank/x)/(l+0.0));
	double lfsigma = ( x == 0 ? 0.1238 : 0.2887/sqrt(x) );

	int nLogStdNormCdfId = int(((lfAvgPeakRank - 0.5)/lfsigma - (-5))/(0.1));
	if(nLogStdNormCdfId < 0)
		nLogStdNormCdfId = 0;
	else if(nLogStdNormCdfId > 100)
		nLogStdNormCdfId = 100;
	
	double lfAvgPeakTol = (x == 0 ? 0 : (lfPeakTol/x));
	double lfPeakTolStd = 0.0;
	double lfLogProbTol = 0.0;
	
	for(size_t i = 0;i < vPeakTol.size() ; ++i)
	{
		lfPeakTolStd += (vPeakTol[i] - lfAvgPeakTol)*(vPeakTol[i] - lfAvgPeakTol);
	}
	if(x > 5)
	{
		lfPeakTolStd /= (x-1);
		lfPeakTolStd = sqrt(lfPeakTolStd);
		 
		int nLogGamaCdfId = int((lfPeakTolStd - 0.01)/0.01);
		if(nLogGamaCdfId < 0)
			nLogGamaCdfId = 0;
		else if(nLogGamaCdfId > 99)
			nLogGamaCdfId = 99;
		
		lfLogProbTol = m_lfLogGamaCDF[nLogGamaCdfId];
	}

	
	lfOdd = lfLogProb ;

	/*
	cout << "N = " << N << endl
	<< "n = " << n << endl
	<< "x = " << x << endl
	<< "l = " << l << endl
	<< "log prob = " << lfLogProb << endl;
	char tmpch;
	cin >> tmpch;
	*/
	
	for(int i=x+1;i<=l;++i)
	{
		lfLogProb = _CalProb(N,n,l,i);
		//cout << "lfOdd = " << lfOdd << endl
		//<< i << " : " << lfLogProb << endl;
		
		if(lfLogProb != MINUS_INFINITE)
		{
			if(lfOdd == MINUS_INFINITE)
				lfOdd = lfLogProb;
			else
				lfOdd = log(exp(lfOdd) + exp(lfLogProb));
		}
	}
	
	// for log tag length
	
	double lfRatio = GetMaxTagLength(nPepId);
	int nPepLen = 0;
	if(nPepId == 0)
	{
		nPepLen = int(m_pPeptide->m_AlphaPeptide.m_tLength) - 1;
	}
	else
	{
		nPepLen = int(m_pPeptide->m_BetaPeptide.m_tLength) - 1;
	}
	if(nPepLen < 1)
		nPepLen = 1;
	if(nPepLen > 120)
		nPepLen = 120;
	double lfLogProbTagLen = m_lfLogTagLenCDF[nPepLen-1][int(lfRatio/0.05)];
	
	/*
	cout << m_pSpectrum->m_strFilePath << endl;
	cout << "taglen = " << lfRatio * nPepLen << endl
		<< "peplen = " << nPepLen<< endl
		<< "ratio = " << lfRatio << endl
		<< "int(lfRatio/0.05)=" << int(lfRatio/0.05) << endl
		<< "prob taglen =" << m_lfLogTagLenCDF[nPepLen-1][int(lfRatio/0.05)] << endl
		<< "prob(taglen) = " << lfLogProbTagLen << endl;
	char tmpch;
	cin >> tmpch;
	*/
	/*	
		cout << "N = " << N << endl
		<< "n = " << n << endl
		<< "x = " << x << endl
		<< "l = " << l << endl
		<< "avrg peak rank =" << lfAvgPeakRank << endl
		<< "taglen = " << lfRatio*nPepLen << endl
		<< "nPepLen = " << nPepLen << endl
 		<< "log prob of tag length = " << lfLogProbTagLen << endl
		<< "log prob final = " <<  lfOdd << endl;
	*/
	
	return lfOdd + m_lfLogStdNormCDF[nLogStdNormCdfId] + lfLogProbTagLen;// + lfLogProbTol;
}

/*
double CPSMatch::GetMatchOdd(int nPepId)
{
	// return matching probability of one peptide under the condition of random match
	if(nPepId < 0 || nPepId > 1)
		return 0.0;
	
	int nLookfor = 0;
	if(nPepId == 0)
		nLookfor = 1;
	else
		nLookfor = -1;
	
	double lfBinWidth = 0.5; // bin width is 0.5 Da

	int N,n,l,x;
	N = n = l = x = 0;
	double lfMinMz, lfMaxMz;
	lfMinMz = m_pSpectrum->m_pPeaks[0].lfMz;
	lfMaxMz = m_pSpectrum->m_pPeaks[m_pSpectrum->m_tPeaksNum-1].lfMz;
	N = int(lfMaxMz/lfBinWidth) - int(lfMinMz/lfBinWidth) + 1;
	
	if(nPepId == 0)
	{
		l = m_tSingleIonNum1;
	}
	else
	{
		l = m_tSingleIonNum2;
	}

	int nCurBinId = -1;
	for(size_t i = 0; i < m_pSpectrum->m_tPeaksNum ; ++i )
	{
		if(int(m_pSpectrum->m_pPeaks[i].lfMz/lfBinWidth)>nCurBinId)
		{
			n++;
			if(m_pSingleMatched[i] == nLookfor || m_pSingleMatched[i] == 2)
				x++;
			nCurBinId = int(m_pSpectrum->m_pPeaks[i].lfMz/lfBinWidth);
		}
	}
	
	//return (x+0.0)/(l+0.0);
	
	// p = C(n,x)*C(N-n,l-x)/C(N,l) 
	// odd = log(p)
	
	//N = 5; n = 4; x = 0 ; l = 4;
	double lfOdd = 0.0;
	double lfLogProb = 0.0; 
	lfLogProb = _CalProb(N, n, l, x);
	// LogProb > 0 (Prob is minus infinity)means impossible incidents
	// LogProb = 0 (Prob is 1)	means inevitable incidents
	
	
	cout << "N = " << N << endl
		<< "n = " << n << endl
		<< "x = " << x << endl
		<< "l = " << l << endl
		<< "log prob = " << lfLogProb << endl;
	char tmpch;
	cin >> tmpch;
	
	
	if(lfLogProb == 0)
		return 0;
	lfOdd = lfLogProb;
	
	for(int i=x+1;i<=l;++i)
	{
		lfLogProb = _CalProb(N,n,l,i);
		//cout << "lfOdd = " << lfOdd << endl
		//<< i << " : " << lfLogProb << endl;
		
		if(lfLogProb != MINUS_INFINITE)
		{
			if(lfOdd == MINUS_INFINITE)
				lfOdd = lfLogProb;
			else
				lfOdd = log(exp(lfOdd) + exp(lfLogProb));
		}
	}
	return lfOdd;
}
*/

/*
double CPSMatch::_CalProb(int N, int n, int l, int x)
{
	//伯松分布
	if( l >= N )
		return 0;
	if(n >= N)
		return 0;
	double mu = (n+0.0)/(N+0.0)*l;
	return (x*log(mu)-mu-m_lfLogNPerm[x]);
}
*/

double CPSMatch::_CalProb(int N, int n, int n0, int l, int x, int x1)
{
	// 超几何分布
	// return log(Prob)
	if( l >= N )
		// inevitable
		return 0;
	if(n >= N)
		// inevitable
		return 0;
	double lfTmpLogProb = 0 , lfLogProb = 0;
	lfTmpLogProb = _C(n0,x1);
	lfLogProb += lfTmpLogProb;
	
	lfTmpLogProb = _C(n - n0 ,x - x1);
	if(lfTmpLogProb < 0)
		// impossible
		return MINUS_INFINITE;

	lfLogProb += lfTmpLogProb;
	lfTmpLogProb = _C(N-n,l-x);
	
	if(lfTmpLogProb < 0)
		// impossible
		return MINUS_INFINITE;
	lfLogProb += lfTmpLogProb;
	lfTmpLogProb = _C(N,l);
	lfLogProb -= lfTmpLogProb;
	return lfLogProb; 
}


double CPSMatch::_CalProb(int N, int n, int l, int x)
{
	// 超几何分布
	// return log(Prob)
	if( l >= N )
		// inevitable
		return 0;
	if(n >= N)
		// inevitable
		return 0;
	double lfTmpLogProb = 0 , lfLogProb = 0;
	lfTmpLogProb = _C(n,x);
	
	if(lfTmpLogProb < 0)
		// impossible
		return MINUS_INFINITE;
	lfLogProb += lfTmpLogProb;
	lfTmpLogProb = _C(N-n,l-x);
	
	if(lfTmpLogProb < 0)
		// impossible
		return MINUS_INFINITE;
	lfLogProb += lfTmpLogProb;
	lfTmpLogProb = _C(N,l);
	lfLogProb -= lfTmpLogProb;
	return lfLogProb; 
}


double CPSMatch::_C(int D,int U)
{
	if(D < U)
		// impossible
		return -1;
	else
	{
		double lftmp;
		lftmp = CPSMatch::m_lfLogNPerm[D] - CPSMatch::m_lfLogNPerm[U] - CPSMatch::m_lfLogNPerm[D-U];
		//cout << "C(" << D << "," << U << ")=" << lftmp << endl;
		return lftmp;
	}
		
}

void CPSMatch::Output()
{
	if(!m_fp)
		return;
	
	_CalMatchInfo();
	
	string strTXT;
	char szbuf[1024];
	
	sprintf(szbuf,"spectrum=%s\n",m_pSpectrum->m_strFilePath.c_str());
	strTXT += szbuf;
	sprintf(szbuf,"MH=%f\n",m_pSpectrum->m_lfMH);
	strTXT += szbuf;
	
	//sprintf(szbuf,"Max Intensity=%f\n",m_pSpectrum->m_lfSqtMaxInten);
	//strTXT += szbuf;
	
	sprintf(szbuf,"sequence=%s",m_pPeptide->m_AlphaPeptide.m_szSequence);
	strTXT += szbuf;
	
	if(m_pPeptide->m_bPair)
	{
		sprintf(szbuf,"-%s",m_pPeptide->m_BetaPeptide.m_szSequence);
		strTXT += szbuf;
	}
	strTXT += "\n";
	
	
	sprintf(szbuf,"site=%d-%d\n",m_pPeptide->m_XLink.m_tAlphaSite,m_pPeptide->m_XLink.m_tBetaSite);
	strTXT += szbuf;
	
	/*
	for(int j=0;j<m_pPeptide->m_AlphaPeptide.m_tModCnt;++j)
	{
		sprintf(szbuf,"%d,%d ",m_pPeptide->m_AlphaPeptide.m_tModSites[j][0],m_pPeptide->m_AlphaPeptide.m_tModSites[j][1]+1);
		strTXT += szbuf;
	}
	
	if(m_pPeptide->m_bPair)
	{
		sprintf(szbuf,"|	");
		strTXT += szbuf;
		for(int j=0;j<m_pPeptide->m_BetaPeptide.m_tModCnt;++j)
		{
			sprintf(szbuf,"%d,%d ",m_pPeptide->m_BetaPeptide.m_tModSites[j][0],m_pPeptide->m_BetaPeptide.m_tModSites[j][1]+1);
			strTXT += szbuf;
		}
	}
	strTXT += "\n";
	*/
	
	// add match info
	strTXT += "matchinfo=";
	for(size_t k=0;k<m_pPeptide->m_AlphaPeptide.m_tLength;++k)
	{
		sprintf(szbuf,"%d",m_pPeptide->m_stMatchInfo.aPepConf[k]);
		strTXT += szbuf;
	}

	if(m_pPeptide->m_bPair)
	{
		strTXT += "-";
		for(size_t k=0;k<m_pPeptide->m_BetaPeptide.m_tLength;++k)
		{
			sprintf(szbuf,"%d",m_pPeptide->m_stMatchInfo.aPepConf[m_pPeptide->m_AlphaPeptide.m_tLength + k]);
			strTXT += szbuf;
		}
	}
	strTXT += "\n";
	
	for(int i=0;i<m_pSpectrum->m_tPeaksNum;++i)
	{
		
		sprintf(szbuf,"%.5f	%.5f	",m_pSpectrum->m_pPeaks[i].lfIntensity,m_pSpectrum->m_pPeaks[i].lfMz);
		strTXT += szbuf;
		
		if(m_pbMatched[i])
		{
			for(int j=0;j<this->m_vMatchedIons[i].size();++j)
			{
				strTXT += GetIonDescription(m_vMatchedIons[i][j]);
				if(m_ionmz[m_vMatchedIons[i][j]].bContainLinker)
				{
					strTXT += "[x]";
				}
				else
				{
					strTXT += "[s]";
				}
				strTXT +="	";
			}
		}
		else
		{
			double lfGap;
			
			/*
			if(m_pSpectrum->m_pPeaks[i].lfIntensity >= 0.1)
			{
				double lfMz = (m_lfPepMass1 + 2*IonMass_Proton)/2;
				fprintf(m_debugFp,"%f\n",m_pSpectrum->m_pPeaks[i].lfMz - lfMz);
				lfMz = (m_lfPepMass2 + 2*IonMass_Proton)/2;
				fprintf(m_debugFp,"%f\n",m_pSpectrum->m_pPeaks[i].lfMz - lfMz);
			}
			*/
			
			for(int j=i-1;j>=0;--j)
			{
				lfGap = m_pSpectrum->m_pPeaks[i].lfMz - m_pSpectrum->m_pPeaks[j].lfMz ; 
				
				if(lfGap > m_psmconf.m_lfMassScope)
					break;
				
				
				if(m_pbMatched[j])
				{
					sprintf(szbuf,"[%s]+[%f]	",GetIonDescription(m_vMatchedIons[j][0]).c_str(),lfGap);
					strTXT += szbuf;
					
					//fprintf(m_debugFp,"%f\n",lfGap);
					break;
				}
			}
			for(int j=i+1;j<m_pSpectrum->m_tPeaksNum;++j)
			{
				lfGap = m_pSpectrum->m_pPeaks[j].lfMz - m_pSpectrum->m_pPeaks[i].lfMz ; 
				
				if(lfGap > m_psmconf.m_lfMassScope)
					break;
				
				if(m_pbMatched[j])
				{
					sprintf(szbuf,"[%s]-[%f]	",GetIonDescription(m_vMatchedIons[j][0]).c_str(),lfGap);
					strTXT += szbuf;
					
					//fprintf(m_debugFp,"%f\n",lfGap);
					break;
				}
			}
		}
	
		strTXT += "\n";
	}
	
	fprintf(m_fp,"%s",strTXT.c_str());
	

	/*
	fprintf(m_debugFp,"spectrum=%s\n",m_pSpectrum->m_strFilePath.c_str());	
	fprintf(m_debugFp,"sequence=%s(%d)-%s(%d)\n"
			,m_pPeptide->m_AlphaPeptide.m_szSequence
			,m_pPeptide->m_XLink.m_tAlphaSite
			,m_pPeptide->m_BetaPeptide.m_szSequence
			,m_pPeptide->m_XLink.m_tBetaSite);
	*/
	
	//fprintf(m_debugFp,"%s\t%f\n",m_pSpectrum->m_strFilePath.c_str(),m_pSpectrum->m_lfSqtMaxInten);
	
	struct IonMatchInfoEx stIonMatchInfo;
	
	
	int nTotal = 0;
	int nMatch = 0;
	int nMatchTop10 = 0;
	int nMatch10Percent = 0;
	int n10PercentTotal = 0;
	
	//fprintf(m_debugFp,"UnMatchIons:\n");
	
	for(size_t i = 0; i < m_pSpectrum->m_tPeaksNum ; ++i )
	{
		nTotal ++ ;
		if(m_pSpectrum->m_pPeaks[i].lfIntensity >= 0.1)
			n10PercentTotal++;
		if(m_pbMatched[i])
		{
			nMatch++;
			if(m_pPeakRank[i]<10)
				nMatchTop10++;
			if(m_pSpectrum->m_pPeaks[i].lfIntensity >= 0.1)
				nMatch10Percent++;
		}
		else
		{
			//fprintf(m_debugFp,"%d	%f	%f\n",m_pPeakRank[i],m_pSpectrum->m_pPeaks[i].lfMz,m_pSpectrum->m_pPeaks[i].lfIntensity);
			if(m_pPeakRank[i]<10)
			{
				//fprintf(m_debugFp,"%d	%f	%f\n",m_pPeakRank[i],m_pSpectrum->m_pPeaks[i].lfMz,m_pSpectrum->m_pPeaks[i].lfIntensity);
			}
		}
		
	}
	
	int nRemovedPeak = 0;
	int nErrorRemovedPeak = 0;
	m_lfNRErrorRatio = 0;
	fprintf(m_debugFp,"%s\n",m_pSpectrum->m_strFilePath.c_str());
	for(size_t i = 0 ;i < m_pSpectrum->m_tPeaksNum ; ++i)
	{
		if(m_vPeakInfo[i].nPeakType == 0)
		{
			// selected peak
			
		}
		else
		{
			// removed peak
			nRemovedPeak ++ ;
			if(m_pbMatched[i])
			{
				nErrorRemovedPeak ++ ;
				fprintf(m_debugFp,"%f	%f	%d	",m_pSpectrum->m_pPeaks[i].lfMz,m_pSpectrum->m_pPeaks[i].lfIntensity,m_vPeakInfo[i].nPeakType);
				for(int j=0;j<this->m_vMatchedIons[i].size();++j)
				{
					string strIon = GetIonDescription(m_vMatchedIons[i][j]);
					fprintf(m_debugFp,"%s	",strIon.c_str());
				}
				fprintf(m_debugFp,"\n");
			}
		}
	}
	
	fprintf(m_debugFp,"%d/%d\n",nErrorRemovedPeak,nRemovedPeak);
	m_lfNRErrorRatio = (nErrorRemovedPeak + 0.0)/nRemovedPeak ; 
	
	/*
	fprintf(m_debugFp,"MatchIons :\n");
	fprintf(m_debugFp,"id\tTypeId\tClassId\tContainLinker\tIntensity\tRank\tCalculatedMz\tExperimentMz\tTolerence\tChargeState\tLossNH3\tLossH2O\tPepPos1\tPepPos2\tSameSite\tOtherLoss\tExpLoss\t\n");
	*/
	
	for(size_t i  = 0;i < m_vIonMatchInfo.size() ; ++ i)
	{
		stIonMatchInfo = m_vIonMatchInfo[i];
		
		/*
		fprintf(m_debugFp,"%d\t%d\t%d\t%d\t%f\t%d\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t\n",
				i,
				stIonMatchInfo.nIonTypeOrder,
				stIonMatchInfo.nIonType,
				stIonMatchInfo.bContainLinker,
				stIonMatchInfo.lfIntensity,
				stIonMatchInfo.nRank,
				stIonMatchInfo.lfCalMz,
				stIonMatchInfo.lfExpMz,
				stIonMatchInfo.lfTol,
				stIonMatchInfo.nCharge,
				stIonMatchInfo.nLossNH3,
				stIonMatchInfo.nLossH2O,
				stIonMatchInfo.nPepPosOrder1,
				stIonMatchInfo.nPepPosOrder2,
				stIonMatchInfo.bSameSite,
				stIonMatchInfo.lfOtherLoss,
				stIonMatchInfo.lfExpTotalLoss);
		*/
		//fprintf(m_debugFp,"%f\t%f\t%f\n",stIonMatchInfo.lfTol,stIonMatchInfo.lfIntensity,stIonMatchInfo.lfExpMz);
		
	}
	//fprintf(m_debugFp,"Top10/Match/Total = %d/%d/%d\n",nMatchTop10,nMatch,nTotal);
	
	m_lfTop10Ratio = (nMatchTop10+0.0)/10;
	if(n10PercentTotal > 0)
		m_lf10PercentRatio = (nMatch10Percent+0.0)/n10PercentTotal;
	else
		m_lf10PercentRatio = 0.0;
	if(nTotal > 0)
		m_lfTotalRatio = (nMatch+0.0)/nTotal;
	else
		m_lfTotalRatio = 0.0;
	//fprintf(m_debugFp,"%f	%f	%f\n\n",m_lfTop10Ratio,m_lf10PercentRatio,m_lfTotalRatio);
	

}

bool SPEC_RANK_LESS(const pair<double,size_t> & rank1, const pair<double,size_t> & rank2)
{
	return rank1.first > rank2.first; 
}


void CPSMatch::_SetPeakRank()
{
	if(m_pPeakRank)
		delete m_pPeakRank;
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

void CPSMatch::Match()
{
	const vector<CIonTypeEx> & vIonTypes = m_psmconf.m_vIonTypes;
	memset(&m_iontag_pep1[0], -1, sizeof(m_iontag_pep1[0]) * vIonTypes.size());
	memset(&m_iontag_pep2[0], -1, sizeof(m_iontag_pep2[0]) * vIonTypes.size());
	/*
	memset(&m_single_iontag_pep1[0], false, sizeof(m_single_iontag_pep1[0]) * vIonTypes.size());
	memset(&m_single_iontag_pep2[0], false, sizeof(m_single_iontag_pep2[0]) * vIonTypes.size());
	*/
	m_vIonMatchInfo.clear();
	if(m_pbMatched)
	{
		delete [] m_pbMatched;
		m_pbMatched = NULL;
	}
	m_pbMatched = new bool[m_pSpectrum->m_tPeaksNum];
	
	if(m_pSingleMatched)
	{
		delete [] m_pSingleMatched;
		m_pSingleMatched = NULL;
	}
	
	if(m_lfSingleMatched)
	{
		delete [] m_lfSingleMatched;
		m_lfSingleMatched = NULL;
	}
	
	m_pSingleMatched = new int[m_pSpectrum->m_tPeaksNum];
	m_lfSingleMatched = new double[m_pSpectrum->m_tPeaksNum];
	
	m_vMatchedIons.clear();
	vector<vector<int> >().swap(m_vMatchedIons);
	
	m_vMatchedIons.resize(m_pSpectrum->m_tPeaksNum);
	for(int i=0;i<m_pSpectrum->m_tPeaksNum;++i)
		m_vMatchedIons.push_back(vector<int>());
	
	memset(m_pbMatched, 0, sizeof(bool) * m_pSpectrum->m_tPeaksNum);
	memset(m_pSingleMatched, 0 , sizeof(int) * m_pSpectrum->m_tPeaksNum);
	memset(m_lfSingleMatched, -1 , sizeof(double) * m_pSpectrum->m_tPeaksNum);
	
	int nPepLen1 = m_pPeptide->m_AlphaPeptide.m_tLength;
	m_vIonMatchInfo.clear();
	
	for(size_t i = 0;i <= m_psmconf.m_vIonTypes.size()+1;++i)
	{
		m_vIonTypeMatchInfo[i].clear();
		m_vCommonIonTypeMatchInfo[i].clear();
		m_vXlinkIonTypeMatchInfo[i].clear();
	}

	for(size_t pos_pep = 0;
		pos_pep < m_tionmzSize;++pos_pep)
	{
		// todo : modify
		if((m_pSpectrum->m_nCharge == 1 && vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge > m_pSpectrum->m_nCharge)
				 || (m_pSpectrum->m_nCharge > 1 && vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].cType != 4 && vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge >= m_pSpectrum->m_nCharge)
				 || (m_pSpectrum->m_nCharge > 1 && vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].cType == 4 && vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge > m_pSpectrum->m_nCharge))
		{

			continue;
		}
		
		int nLastIonTypeOrder1 = m_psmconf.m_vIonTypes.size();
		int nLastIonTypeOrder2 = nLastIonTypeOrder1 + 1;
		
		m_vIonTypeMatchInfo[m_ionmz[pos_pep].nIonTypeOrder].nTheoIonCount++;
		if(m_ionmz[pos_pep].bContainLinker)
			m_vXlinkIonTypeMatchInfo[m_ionmz[pos_pep].nIonTypeOrder].nTheoIonCount++;
		else
			m_vCommonIonTypeMatchInfo[m_ionmz[pos_pep].nIonTypeOrder].nTheoIonCount++;
	
		if((m_ionmz[pos_pep].nIonTypeOrder >= 210 && m_ionmz[pos_pep].nIonTypeOrder <= 215 && m_ionmz[pos_pep].bContainLinker &&  m_ionmz[pos_pep].bSameSite)
				|| (m_ionmz[pos_pep].nIonTypeOrder >= 54 && m_ionmz[pos_pep].nIonTypeOrder <= 59 && m_ionmz[pos_pep].bContainLinker && m_ionmz[pos_pep].nPepPosOrder1 == 0) )
		{
			m_vIonTypeMatchInfo[nLastIonTypeOrder1].nTheoIonCount++;
		}
		
		if( (m_ionmz[pos_pep].nIonTypeOrder >= 288 && m_ionmz[pos_pep].nIonTypeOrder <= 299) 
						|| (m_ionmz[pos_pep].nIonTypeOrder >= 360 && m_ionmz[pos_pep].nIonTypeOrder <= 365))
		{
			m_vIonTypeMatchInfo[nLastIonTypeOrder2].nTheoIonCount++;
		}


		//pos_pep指向一个合适电荷的离子峰
		int nMin, nMax;
		_GetMassBorder(m_ionmz[pos_pep].nMz, vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge, nMin, nMax);
		int nTemp = nMin / MZMULTIPLIER;
		
		if(nTemp >= (int)m_pSpectrum->m_vHash.size() || nTemp < 0)
			continue;
		size_t pos_spec = m_pSpectrum->m_vHash[nTemp];

		while(pos_spec < m_pSpectrum->m_tPeaksNum && 
				m_pSpectrum->m_pPeaks[pos_spec].nMz < nMin)
		{
			++pos_spec;

		}
		
	
		double lfMaxIntensity = 0;
		int nSpecId = -1;
		while(pos_spec < m_pSpectrum->m_tPeaksNum
				&& m_pSpectrum->m_pPeaks[pos_spec].nMz <= nMax)
		{
			if(m_pSpectrum->m_pPeaks[pos_spec].nCharge == 0 || m_pSpectrum->m_pPeaks[pos_spec].nCharge == vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge)
			{
				if(m_pSpectrum->m_pPeaks[pos_spec].lfIntensity > lfMaxIntensity)
				{
					lfMaxIntensity = m_pSpectrum->m_pPeaks[pos_spec].lfIntensity;
					nSpecId = pos_spec;
				}
			}
			++pos_spec;
		}
		
		if(nSpecId >= 0)
		{
			m_pbMatched[nSpecId] = true;
			this->m_vMatchedIons[nSpecId].push_back(pos_pep);
		
			
			// one match
			/*
			struct IonMatchInfoEx stIonMatchInfo;
			stIonMatchInfo.clear();
	
			stIonMatchInfo.nIonTypeOrder = m_ionmz[pos_pep].nIonTypeOrder;
			stIonMatchInfo.nIonType = vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].cType;
			stIonMatchInfo.nPepPosOrder1 = m_ionmz[pos_pep].nPepPosOrder1;
			stIonMatchInfo.nPepPosOrder2 = m_ionmz[pos_pep].nPepPosOrder2;
			stIonMatchInfo.nCharge = vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge;
			stIonMatchInfo.lfExpMz = m_pSpectrum->m_pPeaks[nSpecId].lfMz;
			stIonMatchInfo.lfCalMz = (m_ionmz[pos_pep].nMz+0.0)/MZMULTIPLIER;
			stIonMatchInfo.lfTol =stIonMatchInfo.lfExpMz - stIonMatchInfo.lfCalMz;
			stIonMatchInfo.lfIntensity = m_pSpectrum->m_pPeaks[nSpecId].lfIntensity;
			stIonMatchInfo.nRank = m_pPeakRank[nSpecId];
			stIonMatchInfo.nLossH2O = vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nLossH2O;
			stIonMatchInfo.nLossNH3 = vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nLossNH3;
			stIonMatchInfo.bContainLinker = m_ionmz[pos_pep].bContainLinker;
			stIonMatchInfo.bSameSite = m_ionmz[pos_pep].bSameSite;
			stIonMatchInfo.lfOtherLoss = vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].lfOtherLoss;
			stIonMatchInfo.lfExpTotalLoss = (vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nTotalLostVal+0.0)/MZMULTIPLIER - stIonMatchInfo.lfTol;
			stIonMatchInfo.nAAnum = m_ionmz[pos_pep].nAAnum;
			
			m_vIonMatchInfo.push_back(stIonMatchInfo);
			
			if(m_ionmz[pos_pep].bAddToTag)
			{
				if( m_ionmz[pos_pep].nPepPosOrder1 >= nPepLen1 )
				{
					m_iontag_pep2[m_ionmz[pos_pep].nIonTypeOrder][m_ionmz[pos_pep].nPepPosOrder1 - nPepLen1] = nSpecId;
				}
				else
				{
					m_iontag_pep1[m_ionmz[pos_pep].nIonTypeOrder][m_ionmz[pos_pep].nPepPosOrder1] = nSpecId;	
				}
			}
			*/
			
		}
	}
	
	//char tmpch;
	//cin >> tmpch;
	// remove redundant matched ions
	
	vector<size_t> vMatchedIons;
	for(size_t i=0;i<m_pSpectrum->m_tPeaksNum;++i)
	{
		bool bSimpleType = false;
		bool bOtherType = false;
		vMatchedIons.clear();
		for(size_t j = 0 ;j <m_vMatchedIons[i].size() ; ++ j )
		{
			size_t pos_pep = m_vMatchedIons[i][j];
			
			if(vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].cType < 2)
			{
				bSimpleType = true;
				vMatchedIons.push_back(pos_pep);
			}
			else
			{
				bOtherType = true;
			}
		}
		
		if(bSimpleType && bOtherType)
		{
			m_vMatchedIons[i].clear();
			for(size_t j = 0 ;j <vMatchedIons.size() ; ++ j )
			{
				m_vMatchedIons[i].push_back(vMatchedIons[j]);
			}
		}
		
		for(size_t j = 0 ;j <m_vMatchedIons[i].size() ; ++ j )
		{
			size_t pos_pep = m_vMatchedIons[i][j];
			struct IonMatchInfoEx stIonMatchInfo;
			stIonMatchInfo.clear();
	
			stIonMatchInfo.nIonTypeOrder = m_ionmz[pos_pep].nIonTypeOrder;
			stIonMatchInfo.nIonType = vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].cType;
			stIonMatchInfo.nPepPosOrder1 = m_ionmz[pos_pep].nPepPosOrder1;
			stIonMatchInfo.nPepPosOrder2 = m_ionmz[pos_pep].nPepPosOrder2;
			stIonMatchInfo.nCharge = vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nCharge;
			stIonMatchInfo.lfExpMz = m_pSpectrum->m_pPeaks[i].lfMz;
			stIonMatchInfo.lfCalMz = (m_ionmz[pos_pep].nMz+0.0)/MZMULTIPLIER;
			stIonMatchInfo.lfTol =stIonMatchInfo.lfExpMz - stIonMatchInfo.lfCalMz;
			stIonMatchInfo.lfIntensity = m_pSpectrum->m_pPeaks[i].lfIntensity;
			stIonMatchInfo.nRank = m_pPeakRank[i];
			stIonMatchInfo.nLossH2O = vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nLossH2O;
			stIonMatchInfo.nLossNH3 = vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nLossNH3;
			stIonMatchInfo.bContainLinker = m_ionmz[pos_pep].bContainLinker;
			stIonMatchInfo.bSameSite = m_ionmz[pos_pep].bSameSite;
			stIonMatchInfo.lfOtherLoss = vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].lfOtherLoss;
			stIonMatchInfo.lfExpTotalLoss = (vIonTypes[m_ionmz[pos_pep].nIonTypeOrder].nTotalLostVal+0.0)/MZMULTIPLIER - stIonMatchInfo.lfTol;
			stIonMatchInfo.nAAnum = m_ionmz[pos_pep].nAAnum;
			
			m_vIonMatchInfo.push_back(stIonMatchInfo);
			
			if(m_ionmz[pos_pep].bAddToTag)
			{
				if( m_ionmz[pos_pep].nPepPosOrder1 >= nPepLen1 )
				{
					m_iontag_pep2[m_ionmz[pos_pep].nIonTypeOrder][m_ionmz[pos_pep].nPepPosOrder1 - nPepLen1] = i;
				}
				else
				{
					m_iontag_pep1[m_ionmz[pos_pep].nIonTypeOrder][m_ionmz[pos_pep].nPepPosOrder1] = i;	
				}
			}
		}
	}
	

	m_tSingleIonNum1 = 0;
	m_tSingleIonNum2 = 0;
	
	for(size_t pos_pep = 0;
		pos_pep < m_tSingleIonMzSize;++pos_pep)
	{
		if((m_pSpectrum->m_nCharge == 1 && vIonTypes[m_SingleIonmz[pos_pep].nIonTypeOrder].nCharge > m_pSpectrum->m_nCharge)
				 || (m_pSpectrum->m_nCharge > 1 && vIonTypes[m_SingleIonmz[pos_pep].nIonTypeOrder].nCharge >= m_pSpectrum->m_nCharge))
		{
			continue;
		}
		//pos_pep指向一个合适电荷的离子峰
		
		if(m_SingleIonmz[pos_pep].nPepPosOrder1 >= nPepLen1)
			m_tSingleIonNum2 ++;
		else
			m_tSingleIonNum1 ++;
		
		int nMin, nMax;
		_GetMassBorder(m_SingleIonmz[pos_pep].nMz, vIonTypes[m_SingleIonmz[pos_pep].nIonTypeOrder].nCharge, nMin, nMax);
		
		/*
		if(m_SingleIonmz[pos_pep].nPepPosOrder < nPepLen1)
		{
			cout <<  "iontype :" << int(m_SingleIonmz[pos_pep].nIonTypeOrder) << endl
				<< "peppos :" << int(m_SingleIonmz[pos_pep].nPepPosOrder) << endl
				<< "mz :" << m_SingleIonmz[pos_pep].nMz << endl;
		}
		*/
		
		int nTemp = nMin / MZMULTIPLIER;
		
		if(nTemp >= (int)m_pSpectrum->m_vHash.size() || nTemp < 0)
			continue;
		size_t pos_spec = m_pSpectrum->m_vHash[nTemp];

		while(pos_spec < m_pSpectrum->m_tPeaksNum && 
				m_pSpectrum->m_pPeaks[pos_spec].nMz < nMin)
		{
			++pos_spec;

		}
		
		while(pos_spec < m_pSpectrum->m_tPeaksNum
				&& m_pSpectrum->m_pPeaks[pos_spec].nMz <= nMax)
		{
			
			m_lfSingleMatched[pos_spec] = ( m_SingleIonmz[pos_pep].nMz - m_pSpectrum->m_pPeaks[pos_spec].nMz + 0.0)/MZMULTIPLIER;
				
			if(m_SingleIonmz[pos_pep].nPepPosOrder1 >= nPepLen1)
			{
				if(m_pSingleMatched[pos_spec] == 0)
					m_pSingleMatched[pos_spec] = -1;
				else if(m_pSingleMatched[pos_spec] == 1)
					m_pSingleMatched[pos_spec] = 2;
					
			}
			else
			{
				if(m_pSingleMatched[pos_spec] == 0)
					m_pSingleMatched[pos_spec] = 1;
				else if(m_pSingleMatched[pos_spec] == -1)
					m_pSingleMatched[pos_spec] = 2;
			}
			
			/*
			if(m_SingleIonmz[pos_pep].bAddToTag)
			{
				if( m_SingleIonmz[pos_pep].nPepPosOrder1 >= nPepLen1 )
				{
					m_single_iontag_pep2[m_SingleIonmz[pos_pep].nIonTypeOrder][m_SingleIonmz[pos_pep].nPepPosOrder1 - nPepLen1] = true;
				}
				else
				{
					m_single_iontag_pep1[m_SingleIonmz[pos_pep].nIonTypeOrder][m_SingleIonmz[pos_pep].nPepPosOrder1] = true;	
				}
			}
			*/
			++pos_spec;
		}
		
	}
	
}

void CPSMatch::ComputeMZ()
{
	int nPepLen1 = m_pPeptide->m_AlphaPeptide.m_tLength;
	int nPepLen2 = 0;
	
	int nPepAAMass[2*MAX_PEPTIDE_LENGTH] = {0};
	
	if(m_pPeptide->m_bPair)
	{
		nPepLen2 = m_pPeptide->m_BetaPeptide.m_tLength;
	}
	else
	{
		nPepLen2 = 0;
	}
	
	int i;
	for(int i=0;i<nPepLen1;++i)
	{
		nPepAAMass[i]=m_nAAMass[m_pPeptide->m_AlphaPeptide.m_szSequence[i]-65];
	}
	
	for(i=0 ;i<nPepLen2; ++i)
	{
		nPepAAMass[nPepLen1 + i]=m_nAAMass[m_pPeptide->m_BetaPeptide.m_szSequence[i]-65];		
	}
	
	// modification
	for (size_t i=0;i<m_pPeptide->m_AlphaPeptide.m_tModCnt;++i)
	{
		nPepAAMass[m_pPeptide->m_AlphaPeptide.m_tModSites[i][0]] += (int)(m_psmconf.m_Condition.m_vSelectedVarMod[m_pPeptide->m_AlphaPeptide.m_tModSites[i][1]].m_lfMonoMass_dif*MZMULTIPLIER);
	}
	
	if(m_pPeptide->m_bPair)
	{
		for (size_t i=0;i<m_pPeptide->m_BetaPeptide.m_tModCnt;++i)
		{
			nPepAAMass[nPepLen1 + m_pPeptide->m_BetaPeptide.m_tModSites[i][0]] += (int)(m_psmconf.m_Condition.m_vSelectedVarMod[m_pPeptide->m_BetaPeptide.m_tModSites[i][1]].m_lfMonoMass_dif*MZMULTIPLIER);
		}
	}
	

	//compute b and y ions masses
	int bPepAAMass[2 * MAX_PEPTIDE_LENGTH] = {0};
	int yPepAAMass[2 * MAX_PEPTIDE_LENGTH] = {0};
	
	int mPepAAMass[2 * MAX_PEPTIDE_LENGTH][2 * MAX_PEPTIDE_LENGTH] = {0};
	
	// for single peptide
	int bSinglePepAAMass[2 * MAX_PEPTIDE_LENGTH] = {0};
	int ySinglePepAAMass[2 * MAX_PEPTIDE_LENGTH] = {0};
		
	if(m_pPeptide->m_XLink.m_eXLinkType == 1)
	{
		nPepAAMass[m_pPeptide->m_XLink.m_tAlphaSite] += (int)(m_psmconf.m_Condition.m_vSelectedXLinker[m_pPeptide->m_XLink.m_nLinkerId].m_lfMLMonoMass_dif*MZMULTIPLIER);
	}
	
	int nLinkerMass = 0;
	nLinkerMass = (int)(m_psmconf.m_Condition.m_vSelectedXLinker[m_pPeptide->m_XLink.m_nLinkerId].m_lfMonoMass_dif*MZMULTIPLIER);
	
	int nTmpMass = 0;
	int nPepMass1 = 0 , nPepMass2 = 0;
	
	nPepMass1 = 0;
	nPepMass2 = 0;
	
	for(i = 0;i<nPepLen1;++i)
	{
		nPepMass1 += nPepAAMass[i];
	}
	
	nPepMass1 += 2*nhmass_mono_multi + nomass_mono_multi ;
	
	for(i=0;i<nPepLen2;++i)
	{
		nPepMass2 += nPepAAMass[nPepLen1 + i];
	}
	
	nPepMass2 += 2*nhmass_mono_multi + nomass_mono_multi ;
	
	m_lfPepMass1 = (nPepMass1+0.0)/MZMULTIPLIER;
	m_lfPepMass2 = (nPepMass2+0.0)/MZMULTIPLIER;
	
	int nNTermDeltaMass = npmass_multi;
	int nCTermDeltaMass ;

	nCTermDeltaMass = 2*nhmass_mono_multi + nomass_mono_multi + npmass_multi;
	
	// for single peptide
	int nDeltaMass1 = 0;
	int nDeltaMass2 = 0;
	// delta mass 1 including mass of peptide2 plus water plus linker mass  
	nDeltaMass1 = int(m_pSpectrum->m_lfMH*MZMULTIPLIER) - nPepMass1 - npmass_multi;
	nDeltaMass2 = int(m_pSpectrum->m_lfMH*MZMULTIPLIER) - nPepMass2 - npmass_multi;

	// todo : avoid the influence of precursor ion tolerance
	
	if(m_pPeptide->m_XLink.m_eXLinkType < 2)
	{
		// for none link or mono link	
		bPepAAMass[0] = (int)(nPepAAMass[0] + nNTermDeltaMass);
		yPepAAMass[0] = (int)(nPepAAMass[nPepLen1 - 1] + nCTermDeltaMass);

		for(int i = 1;i < nPepLen1;++i)
		{
			bPepAAMass[i] = bPepAAMass[i - 1] + nPepAAMass[i];
		}
		
		for(int i = 1,j = nPepLen1 - 2;j > 0;++i,--j)
		{
			yPepAAMass[i] = yPepAAMass[i - 1] + nPepAAMass[j];
		}
		
		for(i=0;i<nPepLen1;++i)
		{
			bSinglePepAAMass[i] = bPepAAMass[i];
			ySinglePepAAMass[i] = yPepAAMass[i];
		}

	}
	else if(m_pPeptide->m_XLink.m_eXLinkType == 2)
	{
		// for loop link

		// b ions
		nTmpMass = int(nNTermDeltaMass);
		for(i = 0;i < m_pPeptide->m_XLink.m_tAlphaSite ;++i)
		{
			nTmpMass += nPepAAMass[i];
			bPepAAMass[i] = nTmpMass;
		}

		for(;i < m_pPeptide->m_XLink.m_tBetaSite ;++i)
		{
			nTmpMass += nPepAAMass[i];
			bPepAAMass[i] = nPepMass1;
		}
		nTmpMass += nLinkerMass;
		for(;i<nPepLen1;++i)
		{
			nTmpMass += nPepAAMass[i];
			bPepAAMass[i] = nTmpMass;
		}
		
		// y ions
		nTmpMass = int(nCTermDeltaMass);
		for(i = nPepLen1 - 1;i > m_pPeptide->m_XLink.m_tBetaSite ;--i)
		{
			nTmpMass += nPepAAMass[i];
			yPepAAMass[nPepLen1 - i - 1] = nTmpMass;
		}
		
		for(;i > m_pPeptide->m_XLink.m_tAlphaSite ;--i)
		{
			nTmpMass += nPepAAMass[i];
			yPepAAMass[nPepLen1 - i - 1] = nPepMass1;
		}
		nTmpMass += nLinkerMass;
		for(;i > 0 ;--i)
		{
			nTmpMass += nPepAAMass[i];
			yPepAAMass[nPepLen1 - i - 1] = nTmpMass;
		}
		
		
		for(i=0;i<nPepLen1;++i)
		{
			bSinglePepAAMass[i] = bPepAAMass[i];
			ySinglePepAAMass[i] = yPepAAMass[i];
		}
	}
	else
	{
		// for x-link
		
		// b ions of alpha pep
		nTmpMass = int(nNTermDeltaMass);
		
		int nTmpMass1 = int(nNTermDeltaMass);
		for(i = 0;i < nPepLen1 ;++i)
		{
			nTmpMass += nPepAAMass[i];
			nTmpMass1 += nPepAAMass[i];
			if(i == m_pPeptide->m_XLink.m_tAlphaSite)
			{
				nTmpMass += nPepMass2 + nLinkerMass;
				nTmpMass1 += nDeltaMass1;
			}
			
			bPepAAMass[i] = nTmpMass;
			bSinglePepAAMass[i] = nTmpMass1; 
			
		}
		// y ions of alpha pep
		nTmpMass = int(nCTermDeltaMass);
		nTmpMass1 = int(nCTermDeltaMass);
		for(i = nPepLen1 - 1 ; i > 0 ; --i)
		{
			nTmpMass += nPepAAMass[i];
			nTmpMass1 += nPepAAMass[i];
			if(i == m_pPeptide->m_XLink.m_tAlphaSite)
			{
				nTmpMass += nPepMass2 + nLinkerMass;
				nTmpMass1 += nDeltaMass1;
			}
			
			yPepAAMass[nPepLen1 - 1 - i] = nTmpMass;
			ySinglePepAAMass[nPepLen1 - 1 - i] = nTmpMass1;
		}
		
		// b ions of beta pep
		nTmpMass = int(nNTermDeltaMass);
		nTmpMass1 = int(nNTermDeltaMass);
		
		for(i = 0; i< nPepLen2 ; ++i)
		{
			nTmpMass += nPepAAMass[nPepLen1 + i];
			nTmpMass1 += nPepAAMass[nPepLen1 + i];
			
			if(i == m_pPeptide->m_XLink.m_tBetaSite)
			{
				nTmpMass += nLinkerMass + nPepMass1;
				nTmpMass1 += nDeltaMass2;
			}
			bPepAAMass[nPepLen1 + i] = nTmpMass;
			bSinglePepAAMass[nPepLen1 + i] = nTmpMass1;
			
		}
		// y ions of beta pep
		nTmpMass = int(nCTermDeltaMass);
		nTmpMass1 = int(nCTermDeltaMass);
		
		for(i = nPepLen2 - 1 ;i > 0 ; -- i)
		{
			nTmpMass += nPepAAMass[nPepLen1 + i];
			nTmpMass1 += nPepAAMass[nPepLen1 + i];
			
			if(i == m_pPeptide->m_XLink.m_tBetaSite)
			{
				nTmpMass += nLinkerMass + nPepMass1;
				nTmpMass1 += nDeltaMass2;
			}
			yPepAAMass[nPepLen1 + nPepLen2 - 1 - i] = nTmpMass;
			// bug for 2011.3.6
			ySinglePepAAMass[nPepLen1 + nPepLen2 - 1 - i] = nTmpMass1;
		}
	}

	CMzTriple triple;
	triple.clear();
	m_tionmzSize = 0;
	m_tSingleIonMzSize = 0;
	m_ionmz.clear();

	int j;
	for(i = 0;i < (int)m_psmconf.m_vIonTypes.size();++i)
	{
		const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[i];
		if(IonType.cType >= 2)
			continue;
		
		// only for n-term or c-term ions
		for(j = 0;j < (int)m_pPeptide->m_AlphaPeptide.m_tLength - 1;++j)
		{
			if(m_pPeptide->m_XLink.m_eXLinkType == 2)
			{	
				if(IonType.cType==0 && (j >= m_pPeptide->m_XLink.m_tAlphaSite && j < m_pPeptide->m_XLink.m_tBetaSite))
				{
					continue;
				}
				else if(IonType.cType==1 && (m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j > m_pPeptide->m_XLink.m_tAlphaSite && m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j <= m_pPeptide->m_XLink.m_tBetaSite))
				{
					continue;
				}
			}
			
			triple.clear();
			triple.nIonTypeOrder = i;
			triple.bAddToTag = true;
			triple.nPepPosOrder1 = j;
			triple.nPepPosOrder2 = -1;
			
			triple.bContainLinker = false;
			triple.bSameSite = false;
			
			if(IonType.cType == 0)
			{
				if(m_pPeptide->m_XLink.m_eXLinkType == 2 && j >= m_pPeptide->m_XLink.m_tBetaSite)
				{
					triple.bContainLinker = true;
					triple.nAAnum = j + 1;
				}
				else if(m_pPeptide->m_XLink.m_eXLinkType == 1 && j >= m_pPeptide->m_XLink.m_tAlphaSite)
				{
					triple.bContainLinker = true;
					triple.nAAnum = j + 1;
				}
				else if(m_pPeptide->m_XLink.m_eXLinkType == 3 && j >= m_pPeptide->m_XLink.m_tAlphaSite)
				{
					triple.bContainLinker = true;
					triple.nAAnum = j + 1 + nPepLen2;
				}
				else
				{
					triple.bContainLinker = false;
					triple.nAAnum = j + 1;
				}
			}
			else if(IonType.cType == 1)
			{
				if(m_pPeptide->m_XLink.m_eXLinkType == 2 && m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j <= m_pPeptide->m_XLink.m_tAlphaSite)
				{
					triple.bContainLinker = true;
					triple.nAAnum = j + 1;
				}
				else if(m_pPeptide->m_XLink.m_eXLinkType == 1 && m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j <= m_pPeptide->m_XLink.m_tAlphaSite)
				{
					triple.bContainLinker = true;
					triple.nAAnum = j + 1;
				}
				else if(m_pPeptide->m_XLink.m_eXLinkType == 3 && m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j <= m_pPeptide->m_XLink.m_tAlphaSite)
				{
					triple.bContainLinker = true;
					triple.nAAnum = j + 1 + nPepLen2 ;
				}
				else
				{
					triple.bContainLinker = false;
					triple.nAAnum = j + 1;
				}
			}
			
			triple.bNTerm1 = false;
			triple.bNTerm2 = false;
			triple.nMz = IonType.cType == 0?bPepAAMass[j]:yPepAAMass[j];
			triple.nMz += (IonType.nCharge - 1) * npmass_multi - IonType.nTotalLostVal;
			triple.nMz /= IonType.nCharge;
			
			m_ionmz.push_back(triple);
			
			++m_tionmzSize;

			m_SingleIonmz[m_tSingleIonMzSize].nIonTypeOrder = i;
			m_SingleIonmz[m_tSingleIonMzSize].bAddToTag = true;
			m_SingleIonmz[m_tSingleIonMzSize].nPepPosOrder1 = j;
			m_SingleIonmz[m_tSingleIonMzSize].nPepPosOrder2 = -1;
			m_SingleIonmz[m_tSingleIonMzSize].bNTerm1 = false;
			m_SingleIonmz[m_tSingleIonMzSize].bNTerm2 = false;
			m_SingleIonmz[m_tSingleIonMzSize].nMz = IonType.cType == 0 ?bSinglePepAAMass[j]:ySinglePepAAMass[j];
			m_SingleIonmz[m_tSingleIonMzSize].nMz += (IonType.nCharge - 1) * npmass_multi - IonType.nTotalLostVal;
			m_SingleIonmz[m_tSingleIonMzSize].nMz /= IonType.nCharge;
			++m_tSingleIonMzSize;
			
			
		}
		for(size_t k = 0;k < m_pPeptide->m_AlphaPeptide.m_tModCnt;++k)
		{
			for(size_t ww = 0;ww < m_psmconf.m_Condition.m_vSelectedVarMod[m_pPeptide->m_AlphaPeptide.m_tModSites[k][1]].m_tNLSize;++ww)
			{
				double lfLoss = m_psmconf.m_Condition.m_vSelectedVarMod[m_pPeptide->m_AlphaPeptide.m_tModSites[k][1]].m_vlfMonoNeutralLoss_dif[ww];
				lfLoss *= MZMULTIPLIER;
				triple.clear();
				
				for(int j = m_pPeptide->m_AlphaPeptide.m_tModSites[k][0];j < (int)m_pPeptide->m_AlphaPeptide.m_tLength - 1;++j)
				{
					if(m_pPeptide->m_XLink.m_eXLinkType == 2)
					{	
						if(IonType.cType==0 && (j >= m_pPeptide->m_XLink.m_tAlphaSite && j < m_pPeptide->m_XLink.m_tBetaSite))
						{
							continue;
						}
						else if(IonType.cType==1 && (m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j > m_pPeptide->m_XLink.m_tAlphaSite && m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j <= m_pPeptide->m_XLink.m_tBetaSite))
						{
							continue;
						}
					}

					triple.nIonTypeOrder = i;
					triple.bAddToTag = false;
					triple.nPepPosOrder1 = j;
					triple.nPepPosOrder2 = -1;
					
					triple.bContainLinker = false;
					triple.bSameSite = false;
					triple.nAAnum = j + 1;
					if(IonType.cType == 0)
					{
						if(m_pPeptide->m_XLink.m_eXLinkType == 2 && j >= m_pPeptide->m_XLink.m_tBetaSite)
						{
							triple.bContainLinker = true;
							triple.nAAnum = j + 1;
						}
						else if(m_pPeptide->m_XLink.m_eXLinkType == 1 && j >= m_pPeptide->m_XLink.m_tAlphaSite)
						{
							triple.bContainLinker = true;
							triple.nAAnum = j + 1;
						}
						else if(m_pPeptide->m_XLink.m_eXLinkType == 3 && j >= m_pPeptide->m_XLink.m_tAlphaSite)
						{
							triple.bContainLinker = true;
							triple.nAAnum = j + 1 + nPepLen2;
						}
						else
						{
							triple.bContainLinker = false;
							triple.nAAnum = j + 1;
						}
					}
					else if(IonType.cType == 1)
					{
						if(m_pPeptide->m_XLink.m_eXLinkType == 2 && m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j <= m_pPeptide->m_XLink.m_tAlphaSite)
						{
							triple.bContainLinker = true;
							triple.nAAnum = j + 1;
						}
						else if(m_pPeptide->m_XLink.m_eXLinkType == 1 && m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j <= m_pPeptide->m_XLink.m_tAlphaSite)
						{
							triple.bContainLinker = true;
							triple.nAAnum = j + 1;
						}
						else if(m_pPeptide->m_XLink.m_eXLinkType == 3 && m_pPeptide->m_AlphaPeptide.m_tLength - 1 - j <= m_pPeptide->m_XLink.m_tAlphaSite)
						{
							triple.bContainLinker = true;
							triple.nAAnum = j + 1 + nPepLen2;
						}
						else
						{
							triple.bContainLinker = false;
							triple.nAAnum = j + 1;
						}
					}
					
					triple.bNTerm1 = false;
					triple.bNTerm2 = false;
					triple.nMz = IonType.cType == 0?bPepAAMass[j]:yPepAAMass[j];
					triple.nMz = int(triple.nMz-lfLoss);
					triple.nMz += (IonType.nCharge - 1) * npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					if(triple.nMz < 0)
						continue;
					triple.nMz /= IonType.nCharge;
					m_ionmz.push_back(triple);
					m_tionmzSize++;
					
					
				}
			}
		}
		
		if(m_pPeptide->m_bPair)
		{
			for(j = 0;j < (int)m_pPeptide->m_BetaPeptide.m_tLength - 1;++j)
			{
				triple.clear();
				triple.nIonTypeOrder = i;
				triple.bAddToTag = true;
				triple.nPepPosOrder1 = nPepLen1 + j;
				triple.nPepPosOrder2 = -1;
				
				triple.bContainLinker = false;
				triple.bSameSite = false;
				
				if(IonType.cType == 0)
				{
					if(m_pPeptide->m_XLink.m_eXLinkType == 3 && j >= m_pPeptide->m_XLink.m_tBetaSite)
					{
						triple.bContainLinker = true;
						triple.nAAnum = j + 1 + nPepLen1;
					}
						
					else
					{
						triple.bContainLinker = false;
						triple.nAAnum = j + 1 ;
					}
				}
				else if(IonType.cType == 1)
				{
					if(m_pPeptide->m_XLink.m_eXLinkType == 3 && m_pPeptide->m_BetaPeptide.m_tLength - 1 - j <= m_pPeptide->m_XLink.m_tBetaSite)
					{
						triple.bContainLinker = true;
						triple.nAAnum = j + 1 + nPepLen1;
					}
					else
					{
						triple.bContainLinker = false;
						triple.nAAnum = j + 1 ;
					}
				}
				
				triple.bNTerm1 = false;
				triple.bNTerm2 = false;
				triple.nMz = IonType.cType == 0?bPepAAMass[nPepLen1 + j]:yPepAAMass[nPepLen1 + j];
				triple.nMz += (IonType.nCharge - 1) * npmass_multi - IonType.nTotalLostVal;
				triple.nMz /= IonType.nCharge;
				++m_tionmzSize;
				
				m_ionmz.push_back(triple);
				
				m_SingleIonmz[m_tSingleIonMzSize].nIonTypeOrder = i;
				m_SingleIonmz[m_tSingleIonMzSize].bAddToTag = true;
				m_SingleIonmz[m_tSingleIonMzSize].nPepPosOrder1 = nPepLen1 + j;
				m_SingleIonmz[m_tSingleIonMzSize].nPepPosOrder2 = -1;
				m_SingleIonmz[m_tSingleIonMzSize].bNTerm1 = false;
				m_SingleIonmz[m_tSingleIonMzSize].bNTerm2 = false;
				m_SingleIonmz[m_tSingleIonMzSize].nMz = IonType.cType == 0?bSinglePepAAMass[nPepLen1 + j]:ySinglePepAAMass[nPepLen1 + j];
				m_SingleIonmz[m_tSingleIonMzSize].nMz += (IonType.nCharge - 1) * npmass_multi - IonType.nTotalLostVal;
				m_SingleIonmz[m_tSingleIonMzSize].nMz /= IonType.nCharge;
				++m_tSingleIonMzSize;
			}
			for(size_t k = 0;k < m_pPeptide->m_BetaPeptide.m_tModCnt;++k)
			{
				for(size_t ww = 0;ww <m_psmconf.m_Condition.m_vSelectedVarMod[m_pPeptide->m_BetaPeptide.m_tModSites[k][1]].m_tNLSize;++ww)
				{
					double lfLoss = m_psmconf.m_Condition.m_vSelectedVarMod[m_pPeptide->m_BetaPeptide.m_tModSites[k][1]].m_vlfMonoNeutralLoss_dif[ww];
					lfLoss *= MZMULTIPLIER;
					
					triple.clear();
					for(int j = m_pPeptide->m_BetaPeptide.m_tModSites[k][0];j < (int)m_pPeptide->m_BetaPeptide.m_tLength - 1;++j)
					{
						triple.nIonTypeOrder = i;
						triple.bAddToTag = false;
						triple.nPepPosOrder1 = nPepLen1 + j;
						triple.nPepPosOrder2 = -1;
						triple.bSameSite = false;
						triple.bContainLinker = false;
						if(IonType.cType == 0)
						{
							if(m_pPeptide->m_XLink.m_eXLinkType == 3 && j >= m_pPeptide->m_XLink.m_tBetaSite)
							{
								triple.bContainLinker = true;
								triple.nAAnum = j + 1 + nPepLen1;
							}
							else
							{
								triple.bContainLinker = false;
								triple.nAAnum = j + 1;
							}
						}
						else if(IonType.cType == 1)
						{
							if(m_pPeptide->m_XLink.m_eXLinkType == 3 && m_pPeptide->m_BetaPeptide.m_tLength - 1 - j <= m_pPeptide->m_XLink.m_tBetaSite)
							{
								triple.bContainLinker = true;
								triple.nAAnum = j + 1 + nPepLen1;
							}
							else
							{
								triple.bContainLinker = false;
								triple.nAAnum = j + 1;
							}
						}

						triple.bNTerm1 = false;
						triple.bNTerm2 = false;
						triple.nMz = IonType.cType == 0?bPepAAMass[nPepLen1 + j]:yPepAAMass[nPepLen1 + j];
						triple.nMz = int(triple.nMz-lfLoss);
						triple.nMz += (IonType.nCharge - 1) * npmass_multi;
						triple.nMz -= IonType.nTotalLostVal;
						if(triple.nMz < 0)
							continue;
						triple.nMz /= IonType.nCharge;
						m_ionmz.push_back(triple);
						m_tionmzSize++;
					}
				}
			}
			
		}
	
	}
	

	// for precursor (related) ions
	triple.clear();
	triple.nPepPosOrder1 = -1;
	triple.nPepPosOrder2 = -1;
	triple.bNTerm1 = false;
	triple.bNTerm2 = false;
	triple.bAddToTag = false;
	// todo : all deem false
	triple.bContainLinker = false;
	triple.bSameSite = false; 
	for(int i = 0;i < (int)m_psmconf.m_vIonTypes.size();++i)
	{
		const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[i];
		if(IonType.cType == 3)
		{
			triple.nIonTypeOrder = i;
			triple.nMz = nPepMass1 + (IonType.nCharge)*npmass_multi;
			triple.nMz -= IonType.nTotalLostVal;
			triple.nMz /= IonType.nCharge;
			triple.nAAnum = nPepLen1;
			m_tionmzSize++;
			m_ionmz.push_back(triple);
			if(nPepLen2 > 0)
			{
				triple.nIonTypeOrder = i;
				triple.nMz = nPepMass2 + (IonType.nCharge)*npmass_multi;
				triple.nMz -= IonType.nTotalLostVal;
				triple.nMz /= IonType.nCharge;
				triple.nAAnum = nPepLen2;
				m_tionmzSize++;
				m_ionmz.push_back(triple);
			}
		}
		else if(IonType.cType == 4 && nPepLen2 > 0)
		{
			triple.nIonTypeOrder = i;
			triple.nMz = nPepMass1 + nPepMass2 + nLinkerMass + (IonType.nCharge)*npmass_multi;
			triple.nMz -= IonType.nTotalLostVal;
			triple.nMz /= IonType.nCharge;
			triple.nAAnum = nPepLen1 + nPepLen2;
			m_tionmzSize++;
			m_ionmz.push_back(triple);
		}
	}
	
	// for precursor with activation
	// yP ions
	triple.clear();
	triple.nPepPosOrder1 = -1;
	triple.nPepPosOrder2 = -1;
	triple.bNTerm1 = false;
	triple.bNTerm2 = false;
	triple.bAddToTag = false;
	triple.bContainLinker = false;
	triple.bSameSite = false; 

	for(int m = 0;m < (int)m_psmconf.m_vIonTypes.size();++m)
	{
		const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[m];
		if(IonType.cType == 5)
		{
			triple.nIonTypeOrder = m;
			
			if(m_pPeptide->m_XLink.m_eXLinkType == 3)
			{
				for (int i = 1; i < m_pPeptide->m_XLink.m_tAlphaSite; i++)//add the yP ion on alpha pep
				{
					triple.nMz = 0;
					triple.nAAnum = 0;

					for (int j = i; j < nPepLen1; j++)
					{
						triple.nMz += nPepAAMass[j];
						triple.nAAnum++;
					}
					triple.nPepPosOrder1 = nPepLen1 - i -1;
					triple.nPepPosOrder2 = nPepLen1;
					triple.nMz += 2*nhmass_mono_multi + nomass_mono_multi; 

					triple.nMz += npmass_multi;
					triple.bContainLinker = true;
					
					triple.nMz += (IonType.nCharge - 1)*npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					triple.nMz /= IonType.nCharge;
					m_tionmzSize++;
					m_ionmz.push_back(triple);
				}
				for (int i = 1; i < m_pPeptide->m_XLink.m_tBetaSite; i++)//add the yP ion on beta pep
				{
					triple.nMz = 0;
					triple.nAAnum = 0;

					for (int j = i; j < nPepLen2; j++)
					{
						triple.nMz += nPepAAMass[nPepLen1+j];
						triple.nAAnum++;
					}
					triple.nPepPosOrder1 = nPepLen1 + nPepLen2 - i - 1;
					triple.nPepPosOrder2 = nPepLen1 + nPepLen2;
					triple.nMz += 2*nhmass_mono_multi + nomass_mono_multi; 

					triple.nMz += npmass_multi;
					triple.bContainLinker = true;

					triple.nMz += (IonType.nCharge - 1)*npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					triple.nMz /= IonType.nCharge;
					m_tionmzSize++;
					m_ionmz.push_back(triple);
				}
			}
		}
	}	

	// bP ions
	triple.clear();
	triple.nPepPosOrder1 = -1;
	triple.nPepPosOrder2 = -1;
	triple.bNTerm1 = false;
	triple.bNTerm2 = false;
	triple.bAddToTag = false;
	triple.bContainLinker = false;
	triple.bSameSite = false; 

	for(int m = 0;m < (int)m_psmconf.m_vIonTypes.size();++m)
	{
		const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[m];
		if(IonType.cType == 6)
		{
			triple.nIonTypeOrder = m;
			
			if(m_pPeptide->m_XLink.m_eXLinkType == 3)
			{
				for (int i = m_pPeptide->m_XLink.m_tAlphaSite - 1; i < nPepLen1-1; i++)//add the bP ion on alpha pep
				{
					triple.nMz = 0;
					triple.nAAnum = 0;

					for (int j = 0; j <= i; j++)
					{
						triple.nMz += nPepAAMass[j];
						triple.nAAnum++;
					}
					triple.nPepPosOrder1 = -1;
					triple.nPepPosOrder2 = i;

					triple.nMz += npmass_multi;
					triple.bContainLinker = true;
					
					triple.nMz += (IonType.nCharge - 1)*npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					triple.nMz /= IonType.nCharge;
					m_tionmzSize++;
					m_ionmz.push_back(triple);
/*					cout<<"triple.nMz: "<<triple.nMz<<endl;
					cout<<"IonType.nCharge: "<<IonType.nCharge<<endl;
					cout<<"nPepAAMass: "<<endl;
					for (int j = 0; j < i ; j++)
					{
						cout<<j<<"\t"<<nPepAAMass[j]<<endl;
					}
					system("pause");
*/				}
				for (int i = m_pPeptide->m_XLink.m_tBetaSite - 1; i < nPepLen2 - 1; i++)//add the bP ion on beta pep
				{
					triple.nMz = 0;
					triple.nAAnum = 0;
					
					for (int j = 0; j <= i; j++)
					{
						triple.nMz += nPepAAMass[nPepLen1+j];
						triple.nAAnum++;
					}
					triple.nPepPosOrder1 = nPepLen1;
					triple.nPepPosOrder2 = i + nPepLen1;

					triple.nMz += npmass_multi;
					triple.bContainLinker = true;

					triple.nMz += (IonType.nCharge - 1)*npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					triple.nMz /= IonType.nCharge;
					m_tionmzSize++;
					m_ionmz.push_back(triple);
				}
			}
		}
	}	

	//cout<<"start internal"<<endl;
	// for internal ions;
	vector<int> vInternalIonOrder;
	vInternalIonOrder.clear();
	for(int i = 0;i < (int)m_psmconf.m_vIonTypes.size();++i)
	{
		const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[i];
		if(IonType.cType == 2)
		{
			vInternalIonOrder.push_back(i);
		}
	}
	if(vInternalIonOrder.size() <= 0)
		return ;
	
	int nSum = 0;
	int nPepLen = nPepLen1 + nPepLen2;
	for(int i=0; i<nPepLen; ++i)
	{
		for(int j=i; j<nPepLen; ++j)
		{
			// get rid of boundary conditions
			/*
			if(i == 0 || i== nPepLen1-1 || i == nPepLen1 || i == nPepLen-1)
				continue;
			if(j == 0 || j== nPepLen1-1 || j == nPepLen1 || j == nPepLen-1)
				continue;
			*/
			
			//nPepAAMass
			// get the sum of residues
			nSum = 0;
			bool bNTerm1 = false;
			bool bNTerm2 = false;
			
			triple.clear();
			if(i == j)
				triple.bSameSite = true;
			else
				triple.bSameSite = false;
			
			 
			if(i < nPepLen1)
			{
				if( j < nPepLen1)
				{
					bNTerm1 = false;
					bNTerm2 = true;
					triple.nAAnum = 0;
					if(i == 0 || j == nPepLen1 - 1 )
						continue;
					
					// Ay_Ab
					nSum = 0;
					for(int k = i ; k<=j ; ++k)
					{
						nSum += nPepAAMass[k];
						triple.nAAnum ++;
					}
					
					nSum += npmass_multi;
					if(j == nPepLen1-1)
					{
						nSum += 2*nhmass_mono_multi + nomass_mono_multi;
					}
					
					triple.bContainLinker = false;
					if(m_pPeptide->m_XLink.m_eXLinkType == 2)
					{
						if(i <= m_pPeptide->m_XLink.m_tAlphaSite && j >= m_pPeptide->m_XLink.m_tBetaSite)
						{
							nSum += nLinkerMass;
							triple.bContainLinker = true;
						}
					}
					else if(m_pPeptide->m_XLink.m_eXLinkType == 3)
					{
						if(i <= m_pPeptide->m_XLink.m_tAlphaSite && j >= m_pPeptide->m_XLink.m_tAlphaSite)
						{
							nSum += (nLinkerMass + nPepMass2);
							triple.nAAnum += nPepLen2;
							triple.bContainLinker = true;
						}
					}
					
					triple.nPepPosOrder1 = nPepLen1 - 1 - i;
					triple.nPepPosOrder2 = j;
					
					triple.bNTerm1 = bNTerm1;
					triple.bNTerm2 = bNTerm2;
					triple.bAddToTag = false;
					for(int w = 0 ; w < vInternalIonOrder.size(); ++ w)
					{
						const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[vInternalIonOrder[w]];
						triple.nIonTypeOrder = vInternalIonOrder[w];
						triple.nMz = nSum + (IonType.nCharge - 1)*npmass_multi;
						triple.nMz -= IonType.nTotalLostVal;
						triple.nMz /= IonType.nCharge;
						m_tionmzSize++;
						m_ionmz.push_back(triple);
					}
				}
				else
				{
					
					if(i <= m_pPeptide->m_XLink.m_tAlphaSite)
					{
						if( j <= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite)
						{
							bNTerm1 = false;
							bNTerm2 = false;
							triple.nAAnum = 0;
							if(i != 0 && j != nPepLen1)
							{	
							
								// Ay_By
								nSum = 0;
								for(int k = i ; k<nPepLen1 ; ++k)
								{
									nSum += nPepAAMass[k];
									triple.nAAnum ++;
								}
								for(int k = j ; k<nPepLen ; ++k)
								{
									nSum += nPepAAMass[k];
									triple.nAAnum ++;
								}
							
								nSum += 2*(2*nhmass_mono_multi + nomass_mono_multi) + npmass_multi;
								nSum += nLinkerMass;
								triple.bContainLinker = true;
								
								triple.nPepPosOrder1 = nPepLen1 - 1 - i;
								triple.nPepPosOrder2 = nPepLen2 - 1 - (j - nPepLen1) + nPepLen1;
	
								triple.bNTerm1 = bNTerm1;
								triple.bNTerm2 = bNTerm2;
								triple.bAddToTag = false;
								/* todo :
								 * don't know how to calculate the ion mass
								 * 
								for(int w = 0 ; w < vInternalIonOrder.size(); ++ w)
								{
									const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[vInternalIonOrder[w]];
									triple.nIonTypeOrder = vInternalIonOrder[w];
									triple.nMz = nSum + (IonType.nCharge - 1)*npmass_multi;
									triple.nMz -= IonType.nTotalLostVal;
									triple.nMz /= IonType.nCharge;
									m_tionmzSize++;
									m_ionmz.push_back(triple);
								}
								*/
							}
						}
						if(j >= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite)
						{
							bNTerm1 = false;
							bNTerm2 = true;
							triple.nAAnum = 0;
							if(i != 0 && j != nPepLen - 1)
							{
								// Ay_Bb
								nSum = 0;
								for(int k = i ; k<= j  ; ++k)
								{
									nSum += nPepAAMass[k];
									triple.nAAnum ++;
								}
	
								nSum += 2*nhmass_mono_multi + nomass_mono_multi + npmass_multi;
							
								if(j == nPepLen - 1)
								{
									nSum += 2*nhmass_mono_multi + nomass_mono_multi;
								}
								
								nSum += nLinkerMass;
								triple.bContainLinker = true;
								
								triple.nPepPosOrder1 = nPepLen1 - 1 - i;
								triple.nPepPosOrder2 = (j - nPepLen1) + nPepLen1;
	
								triple.bNTerm1 = bNTerm1;
								triple.bNTerm2 = bNTerm2;
								triple.bAddToTag = false;
								for(int w = 0 ; w < vInternalIonOrder.size(); ++ w)
								{
									const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[vInternalIonOrder[w]];
									triple.nIonTypeOrder = vInternalIonOrder[w];
									triple.nMz = nSum + (IonType.nCharge - 1)*npmass_multi;
									triple.nMz -= IonType.nTotalLostVal;
									triple.nMz /= IonType.nCharge;
									m_tionmzSize++;
									m_ionmz.push_back(triple);
								}
	
							}
						}
					}
					if(i >= m_pPeptide->m_XLink.m_tAlphaSite)
					{
						if(j <= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite)
						{
							bNTerm1 = true;
							bNTerm2 = false;
							triple.nAAnum = 0;
							if(i != nPepLen1 - 1 && j != nPepLen1 )
							{
							
								// Ab_By
								nSum = 0;
								for(int k = 0 ; k<=i ; ++k)
								{
									nSum += nPepAAMass[k];
									triple.nAAnum++;
								}
								for(int k = j ; k<nPepLen ; ++k)
								{
									nSum += nPepAAMass[k];
									triple.nAAnum++;
								}
								
								nSum += (2*nhmass_mono_multi + nomass_mono_multi) + npmass_multi;
								
								if(i == nPepLen1 - 1)
								{
									nSum += 2*nhmass_mono_multi + nomass_mono_multi;
								}
								
								nSum += nLinkerMass;
								triple.bContainLinker = true;
								
								triple.nPepPosOrder1 = i;
								triple.nPepPosOrder2 = nPepLen2 - 1 - (j - nPepLen1) + nPepLen1;						
								triple.bNTerm1 = bNTerm1;
								triple.bNTerm2 = bNTerm2;
								triple.bAddToTag = false;
								for(int w = 0 ; w < vInternalIonOrder.size(); ++ w)
								{
									const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[vInternalIonOrder[w]];
									triple.nIonTypeOrder = vInternalIonOrder[w];
									triple.nMz = nSum + (IonType.nCharge - 1)*npmass_multi;
									triple.nMz -= IonType.nTotalLostVal;
									triple.nMz /= IonType.nCharge;
									m_tionmzSize++;
									m_ionmz.push_back(triple);
								}
							}


						}
						if(j >= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite)
						{
							bNTerm1 = true;
							bNTerm2 = true;
							triple.nAAnum = 0;
							if( i != nPepLen1 - 1 && j != nPepLen -1 )
							{
								// Ab_Bb
								nSum = 0;
								for(int k = 0 ; k<=i ; ++k)
								{
									nSum += nPepAAMass[k];
									triple.nAAnum++;
								}
								for(int k = nPepLen1 ; k<=j ; ++k)
								{
									nSum += nPepAAMass[k];
									triple.nAAnum++;
								}
								
								nSum += npmass_multi;
								
								if(i == nPepLen1 - 1)
								{
									nSum += 2*nhmass_mono_multi + nomass_mono_multi;
								}
								if(j == nPepLen - 1)
								{
									nSum += 2*nhmass_mono_multi + nomass_mono_multi;
								}
								
								nSum += nLinkerMass;
								triple.bContainLinker = true;
								/* todo :
								 * don't know how to calculate the ion mass
								 * 
								triple.nPepPosOrder1 = i;
								triple.nPepPosOrder2 = (j - nPepLen1) + nPepLen1;
								triple.bNTerm1 = bNTerm1;
								triple.bNTerm2 = bNTerm2;
								triple.bAddToTag = false;
								for(int w = 0 ; w < vInternalIonOrder.size(); ++ w)
								{
									const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[vInternalIonOrder[w]];
									triple.nIonTypeOrder = vInternalIonOrder[w];
									triple.nMz = nSum + (IonType.nCharge - 1)*npmass_multi;
									triple.nMz -= IonType.nTotalLostVal;
									triple.nMz /= IonType.nCharge;
									m_tionmzSize++;
									m_ionmz.push_back(triple);
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
				triple.nAAnum = 0;
				if( i == nPepLen1 || j == nPepLen - 1)
					continue;
				
				// By_Bb
				nSum = 0;
				for(int k = i ; k<=j ; ++k)
				{
					nSum += nPepAAMass[k];
					triple.nAAnum++;
				}
				
				nSum += npmass_multi;
				
				if(j == nPepLen-1)
				{
					nSum += 2*nhmass_mono_multi + nomass_mono_multi;
				}
				
				triple.bContainLinker = false;
				if((i <= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite) && (j >= nPepLen1 + m_pPeptide->m_XLink.m_tBetaSite))
				{
					nSum += (nLinkerMass+nPepMass1);
					triple.nAAnum += nPepLen1;
					triple.bContainLinker = true;
				}

				triple.nPepPosOrder1 = nPepLen2 - 1 - (i - nPepLen1) + nPepLen1;
				triple.nPepPosOrder2 = (j - nPepLen1) + nPepLen1;
				
				triple.bNTerm1 = bNTerm1;
				triple.bNTerm2 = bNTerm2;
				triple.bAddToTag = false;
				for(int w = 0 ; w < vInternalIonOrder.size(); ++ w)
				{
					const CIonTypeEx & IonType = m_psmconf.m_vIonTypes[vInternalIonOrder[w]];
					triple.nIonTypeOrder = vInternalIonOrder[w];
					triple.nMz = nSum + (IonType.nCharge - 1)*npmass_multi;
					triple.nMz -= IonType.nTotalLostVal;
					triple.nMz /= IonType.nCharge;
					m_tionmzSize++;
					m_ionmz.push_back(triple);
				}
			}
		}
	}
}
