#include "../include/sdk.h"
#include "Mass2PepIndex.h"
#include "PepFilter.h"


using namespace std;
using namespace proteomics_sdk;

namespace Mass2PepIndex
{
CPepFilter::CPepFilter(void)
:m_tMaxMissedClvg(DEFAULT_MAX_MISSED_SITES)
, m_tMinPepLength(DEFAULT_MIN_IDX_PEP_LEN)
, m_tMaxPepLength(DEFAULT_MAX_IDX_PEP_LEN)
, m_lfMinPepMass(DEFAULT_MIN_IDX_PEP_MASS)
, m_lfMaxPepMass(DEFAULT_MAX_IDX_PEP_MASS)
{
}

CPepFilter::~CPepFilter(void)
{
}

size_t CPepFilter::GetMaxMissedClvg( void )
{ 
	return m_tMaxMissedClvg;
}

size_t CPepFilter::GetMaxPepLength( void )
{
	return m_tMaxPepLength;
}

size_t CPepFilter::GetMinPepLength( void )
{
	return m_tMinPepLength;
}

double  CPepFilter::GetMaxPepMass( void )
{
	return m_lfMaxPepMass;
}

double  CPepFilter::GetMinPepMass( void )
{
	return m_lfMinPepMass;
}

void CPepFilter::SetMaxMissedClvg( size_t tMissedSites )
{
	m_tMaxMissedClvg = tMissedSites;
}

void CPepFilter::SetMaxPepLength( size_t tMax )
{
	if(tMax > DEFAULT_MAX_IDX_PEP_LEN || tMax<DEFAULT_MIN_IDX_PEP_LEN)
	{
		CErrInfo info("CPepFilter", "SetMaxPepLength", "Length invalid");
		char str[20];
//		itoa(tMax, str, 10);
		sprintf(str,"%d",tMax);
		info.Append("tMax = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

	m_tMaxPepLength = tMax;
}

void CPepFilter::SetMinPepLength( size_t tMin )
{	
	if(tMin > DEFAULT_MAX_IDX_PEP_LEN || tMin < DEFAULT_MIN_IDX_PEP_LEN)
	{
		CErrInfo info("CPepFilter", "SetMinPepLength", "Length invalid");
		char str[20];
//		itoa(tMin, str, 10);
		sprintf(str,"%d",tMin);
		info.Append("tMax = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
	
	m_tMinPepLength = tMin;
}

void CPepFilter::SetMaxPepMass( double lfMaxMass )
{
	if(lfMaxMass > DEFAULT_MAX_IDX_PEP_MASS || lfMaxMass <DEFAULT_MIN_IDX_PEP_MASS)
	{
		CErrInfo info("CPepFilter", "SetMaxPepMass", "Mass invalid");
		char str[20];
//		itoa((int)(lfMaxMass+0.5), str, 10);
		sprintf(str,"%d",(int)(lfMaxMass+0.5));
		info.Append("lfMaxMass = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

	m_lfMaxPepMass = lfMaxMass;
}

void CPepFilter::SetMinPepMass( double lfMinMass )
{
	if( lfMinMass > DEFAULT_MAX_IDX_PEP_MASS || lfMinMass < DEFAULT_MIN_IDX_PEP_MASS )
	{
		CErrInfo info("CPepFilter", "SetMinPepMass", "Mass invalid");
		char str[20];
//		itoa((int)(lfMinMass+0.5), str, 10);
		sprintf(str,"%d",(int)(lfMinMass+0.5));
		info.Append("lfMinMass = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

	m_lfMinPepMass = lfMinMass;
}

void CPepFilter::SetMaxPepMassAllowed(double lfMaxMass)
{
	if(lfMaxMass > DEFAULT_MAX_IDX_PEP_MASS || lfMaxMass <DEFAULT_MIN_IDX_PEP_MASS)
	{
		CErrInfo info("CPepFilter", "SetMaxPepMassAllowed", "Mass invalid");
		char str[20];
//		itoa((int)(lfMaxMass+0.5), str, 10);
		sprintf(str,"%d",(int)(lfMaxMass+0.5));
		info.Append("lfMaxMass = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

	m_lfMaxPepMassAllowed = lfMaxMass;
}
void CPepFilter::SetMinPepMassAllowed(double lfMinMass)
{
	if( lfMinMass > DEFAULT_MAX_IDX_PEP_MASS || lfMinMass < DEFAULT_MIN_IDX_PEP_MASS )
	{
		CErrInfo info("CPepFilter", "SetMinPepMassAllowed", "Mass invalid");
		char str[20];
//		itoa((int)(lfMinMass+0.5), str, 10);
		sprintf(str,"%d",(int)(lfMinMass+0.5));
		info.Append("lfMinMass = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
	
	m_lfMinPepMassAllowed = lfMinMass;
}
// Determining whether input Peptide'sequence is satisfied conditions.  true means passed. false means filtered.
bool CPepFilter::IsPassed(size_t tLen, double lfPepMass,size_t tMissedSite)
{
	bool bRlt = false;
	if (lfPepMass != 0) {
		bRlt = IsLowerMinPepMass(lfPepMass) || IsUpperMaxPepMass(lfPepMass);
	}
	if (tLen != 0) {
		bRlt = bRlt || IsUpperMaxPepLength(tLen) || IsLowerMinPepLength(tLen);
	}
	return !( bRlt || IsUpperMaxMissClvg(tMissedSite) );
}

bool CPepFilter::IsUpperMaxMissClvg(size_t tMissSite)
{
	return tMissSite > m_tMaxMissedClvg;
}

bool CPepFilter::IsUpperMaxPepLength(size_t tPepLen)
{
	return tPepLen > m_tMaxPepLength;
}

bool CPepFilter::IsUpperMaxPepMass(double lfMass)
{
	return lfMass > m_lfMaxPepMass;
}

bool CPepFilter::IsLowerMinPepLength(size_t tPepLen)
{
	return tPepLen < m_tMinPepLength;
}

bool CPepFilter::IsLowerMinPepMass(double lfMass)
{
	return lfMass < m_lfMinPepMass;
}
}
