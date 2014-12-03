#ifndef CALCBIOVARIABLE_H_
#define CALCBIOVARIABLE_H_

#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>
using namespace std;


#define PH_MIN 0		/* minimum pH value */
#define PH_MAX 14		/* maximum pH value */
#define MAXLOOP 2000	/* maximum number of iterations */
#define EPSI 0.0001		/* desired precision */
#define EXP10(x) (exp(2.30258509299405E0*(x)))

#define H_MONO 1.0078
#define H_AVRG 1.0079
#define O_MONO 15.9949
#define O_AVRG 15.9994
#define H2O_MONO (H_MONO * 2 + O_MONO)
#define H2O_AVRG (H_AVRG * 2 + O_AVRG)

#define NUM_AA 25

const double lfAAMass[NUM_AA][2] =
{
{ 71.03711, 71.0788 },
{ 156.10111, 156.1876 },
{ 114.04293, 114.1039 },
{ 115.02694, 115.0886 },
{ 103.00919, 103.1448 },
{ 129.04259, 129.1155 },
{ 128.05858, 128.1308 },
{ 57.02146, 57.0520 },
{ 137.05891, 137.1412 },
{ 113.08406, 113.1595 },
{ 113.08406, 113.1595 },
{ 128.09496, 128.1742 },
{ 131.04049, 131.1986 },
{ 147.06841, 147.1766 },
{ 97.05276, 97.1167 },
{ 87.03203, 87.0782 },
{ 101.04768, 101.1051 },
{ 186.07931, 186.2133 },
{ 163.06333, 163.1760 },
{ 99.06841, 99.1326 },
{ 173.03242, 173.1253 },
{ 113.04768, 113.1161 },
{ 113.08406, 113.1595 },
{ 111.03203, 111.1002 },
{ 103.00919, 103.1448 } };

const char szAAName[NUM_AA + 1] = "ARNDCEQGHILKMFPSTWYVJXUBZ";

const double lfAAGravy[NUM_AA] =
{ 1.800, -4.500, -3.500, -3.500, 2.500, -3.500, -3.500, -0.400, -3.200, 4.500, 3.800, -3.900,
		1.900, 2.800, -1.600, -0.800, -0.700, -0.900, -1.300, 4.200, 0, 0, 0, 0, 0 };

class CBioCalculator
{
public:

	CBioCalculator()
	{
		int i = 0;
		for (i = 0; i < NUM_AA; ++i)
		{
			m_mapMoleMonoMass[szAAName[i]] = lfAAMass[i][0];
			m_mapMoleAvrgMass[szAAName[i]] = lfAAMass[i][1];
		}
	}

	map<char, double> m_mapMoleMonoMass;
	map<char, double> m_mapMoleAvrgMass;

	double CalcMW(const string &s, bool bAvrg)
	{
		int len = (int) s.size();
		int i = 0;
		double res = 0;

		if (bAvrg)
		{
			for (i = 0; i < len; ++i)
				res += m_mapMoleAvrgMass[s[i]];
			res += H2O_AVRG;
		}
		else
		{
			for (i = 0; i < len; ++i)
				res += m_mapMoleMonoMass[s[i]];
			res += H2O_MONO;
		}
		return res;
	}

	double CalcGravy(const string & strSQ)
	{
		int i = 0;
		double nRes = 0;
		for (i = 0; i < (int) strSQ.size(); ++i)
		{
			nRes += lfAAGravy[(int) strSQ[i]];
		}
		return nRes;
	}

	double CalcPI(const string & strSQ)
	{
		static double cPk[26][5] =
		{
		{ 3.55, 7.59, 0.0, 0.0, 0.0 }, // A

				{ 3.55, 7.50, 0.0, 0.0, 0.0 }, // B

				{ 3.55, 7.50, 9.00, 9.00, 9.00 },// C

				{ 4.55, 7.50, 4.05, 4.05, 4.05 }, // D

				{ 4.75, 7.70, 4.45, 4.45, 4.45 }, // E

				{ 3.55, 7.50, 0.0, 0.0, 0.0 }, // F

				{ 3.55, 7.50, 0.0, 0.0, 0.0 }, // G

				{ 3.55, 7.50, 5.98, 5.98, 5.98 }, // H

				{ 3.55, 7.50, 0.0, 0.0, 0.0 }, // I

				{ 0.00, 0.00, 0.0, 0.0, 0.0 },// J

				{ 3.55, 7.50, 10.00, 10.00, 10.00 }, // K

				{ 3.55, 7.50, 0.0, 0.0, 0.0 }, // L

				{ 3.55, 7.00, 0.0, 0.0, 0.0 }, // M

				{ 3.55, 7.50, 0.0, 0.0, 0.0 }, // N

				{ 0.00, 0.00, 0.0, 0.0, 0.0 }, // O

				{ 3.55, 8.36, 0.0, 0.0, 0.0 }, // P

				{ 3.55, 7.50, 0.0, 0.0, 0.0 }, // Q

				{ 3.55, 7.50, 12.0, 12.0, 12.0 }, // R

				{ 3.55, 6.93, 0.0, 0.0, 0.0 },// S

				{ 3.55, 6.82, 0.0, 0.0, 0.0 }, // U

				{ 3.55, 7.44, 0.0, 0.0, 0.0 }, // V

				{ 3.55, 7.50, 0.0, 0.0, 0.0 }, // W

				{ 3.55, 7.50, 0.0, 0.0, 0.0 }, // X

				{ 3.55, 7.50, 10.00, 10.00, 10.00 }, // Y

				{ 3.55, 7.50, 0.0, 0.0, 0.0 } // Z
		};

		int i;
		int R, H, K, D, E, C, Y;
		int comp[26] =
		{ 0 };
		int nTermResidue;
		int cTermResidue;
		double phMin;
		double phMid;
		double phMax;
		double charge;
		double cter, nter, carg, chis, clys, casp, cglu, ccys, ctyr;

		R = 'R' - 'A';
		H = 'H' - 'A';
		K = 'K' - 'A';
		D = 'D' - 'A';
		E = 'E' - 'A';
		C = 'C' - 'A';
		Y = 'Y' - 'A';

		// Compute the amino-acid composition.

		int l = (int) strSQ.size();
		for (i = 0; i < l; i++)
			comp[strSQ[i] - 'A']++;

		// Look up N-terminal and C-terminal residue.

		nTermResidue = strSQ[0] - 'A';
		cTermResidue = strSQ[l - 1] - 'A';

		phMin = PH_MIN;
		phMax = PH_MAX;

		for (i = 0, charge = 1.0; i < MAXLOOP && (phMax - phMin) > EPSI; i++)
		{
			phMid = phMin + (phMax - phMin) / 2;

			// EXP10( x ) = exp( 2.30258509299405E0*x )
			cter = EXP10(-cPk[cTermResidue][0]) / (EXP10(-cPk[cTermResidue][0]) + EXP10(-phMid));
			nter = EXP10(-phMid) / (EXP10(-cPk[nTermResidue][1]) + EXP10(-phMid));

			carg = comp[R] * EXP10(-phMid) / (EXP10(-cPk[R][2]) + EXP10(-phMid));
			chis = comp[H] * EXP10(-phMid) / (EXP10(-cPk[H][2]) + EXP10(-phMid));
			clys = comp[K] * EXP10(-phMid) / (EXP10(-cPk[K][2]) + EXP10(-phMid));

			casp = comp[D] * EXP10(-cPk[D][2]) / (EXP10(-cPk[D][2]) + EXP10(-phMid));
			cglu = comp[E] * EXP10(-cPk[E][2]) / (EXP10(-cPk[E][2]) + EXP10(-phMid));

			ccys = comp[C] * EXP10(-cPk[C][2]) / (EXP10(-cPk[C][2]) + EXP10(-phMid));
			ctyr = comp[Y] * EXP10(-cPk[Y][2]) / (EXP10(-cPk[Y][2]) + EXP10(-phMid));

			charge = carg + clys + chis + nter - (casp + cglu + ctyr + ccys + cter);

			if (charge > 0.0)
				phMin = phMid;
			else
				phMax = phMid;
		}

		return phMid;
	}
};

#endif /* CALCBIOVARIABLE_H_ */
