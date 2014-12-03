#include "../include/sdk.h"
#include "../include/option.h"
#include "PepCalcFunc.h"

using namespace std; 
using namespace proteomics_sdk;

CPepCalcFunc::CPepCalcFunc(void)
:m_strPepSQ("")
,m_strAAPath("")
{
}

CPepCalcFunc::CPepCalcFunc(string strAAPath)
:m_strPepSQ("")
,m_strAAPath("")
{
	SetAAPath(strAAPath);
	SetMapAA();
}

CPepCalcFunc::~CPepCalcFunc(void)
{
}

void CPepCalcFunc::Init(std::string strAAPath)
{
	m_strPepSQ = "";
	m_strAAPath = "";
	SetAAPath(strAAPath);
	SetMapAA();
}

void CPepCalcFunc::SetAAPath(string strAAPath)
{
	m_strAAPath = strAAPath;
}

void CPepCalcFunc::SetMapAA()
{
	CAAConf aa_conf(m_strAAPath);
	m_mapAA = aa_conf.GetMapAAMass();
}

// Get the double-type mass of specified residue
double CPepCalcFunc::GetfAAMass(const char cAA, bool bMono)
{
	try
	{
		if ( cAA < 'A' || cAA > 'Z' ) return 0;
		if ( bMono ) 
			return ( m_mapAA.GetAA(cAA).m_lfMonoMass);
		else
			return ( m_mapAA.GetAA(cAA).m_lfAvrgMass);
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CPepCalcFunc", "GetfAAMass", "return double.");
		string str = " ";
		str[0] = cAA;		
		info.Append("cAA = " + str);		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPepCalcFunc", "GetfAAMass", "Caught an unkown exception,return double.!");
		string str = " ";
		str[0] = cAA;		
		info.Append("cAA = " + str);
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

// Get the int-type mass of specified residue
size_t CPepCalcFunc::GettAAMass(const char cAA, bool bMono)
{
	try
	{
		if ( cAA < 'A' || cAA > 'Z' ) return 0;
		if ( bMono ) 
			return ( m_mapAA.GetAA(cAA).m_tmonoMass);
		else
			return ( m_mapAA.GetAA(cAA).m_tavrgMass);
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CPepCalcFunc", "GetfAAMass", "return size_t.");
		string str = " ";
		str[0] = cAA;		
		info.Append("cAA = " + str);		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPepCalcFunc", "GetfAAMass", "Caught an unkown exception,return size_t.");
		string str = " ";
		str[0] = cAA;		
		info.Append("cAA = " + str);
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void CPepCalcFunc::SumSQResiduesMass(double&	fSum,const string& strSQ,bool bMono)
{
	fSum = 0 ;
	if ( strSQ.empty() )
	{
		CErrInfo info("CPepCalcFunc", "SumSQResiduesMass", "Strings should not be null!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

	try
	{
		string::const_iterator itAA;
		for( itAA = strSQ.begin(); itAA != strSQ.end(); ++itAA )
		{
			fSum += GetfAAMass( *itAA,bMono );
		}
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CPepCalcFunc", "SumSQResiduesMass", "return double.");
		info.Append("strSQ = " + strSQ);		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPepCalcFunc", "SumSQResiduesMass", "Caught an unkown exception,return double.");
		info.Append("strSQ = " + strSQ);
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}
void CPepCalcFunc::SumSQResiduesMass(size_t& tSum,const string& strSQ,bool bMono)
{
	tSum = 0 ;
	if ( strSQ.empty() )
	{
		CErrInfo info("CPepCalcFunc", "SumSQResiduesMass", "Strings should not be null!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}

	try
	{
		string::const_iterator itAA;
		for( itAA = strSQ.begin(); itAA != strSQ.end(); ++itAA )
		{
			tSum += GettAAMass( *itAA, bMono);
		}
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CPepCalcFunc", "SumSQResiduesMass", "return size_t.");
		info.Append("strSQ = " + strSQ);		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPepCalcFunc", "SumSQResiduesMass", "Caught an unkown exception,return size_t.");
		info.Append("strSQ = " + strSQ);
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}
// calc the precursor MH+ mass (add 18 +1)
void CPepCalcFunc::CalcPepMHMass(double&  fMH,const string& strSQ,bool bMono)
{
	try
	{
		SumSQResiduesMass(fMH,strSQ,bMono);
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CPepCalcFunc", "CalcPepMHMass", "return double.");
		info.Append("strSQ = " + strSQ);		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPepCalcFunc", "SumSQResiduesMass", "Caught an unkown exception,return double.");
		info.Append("strSQ = " + strSQ);
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}
void CPepCalcFunc::	CalcPepMHMass(size_t&  tMH,const string& strSQ,bool bMono)
{
	try
	{
		SumSQResiduesMass(tMH,strSQ,bMono);
	}
	catch(runtime_error &e)
	{
		CErrInfo info("CPepCalcFunc", "CalcPepMHMass", "return size_t.");
		info.Append("strSQ = " + strSQ);		
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CPepCalcFunc", "CalcPepMHMass", "Caught an unkown exception,return size_t.");
		info.Append("strSQ = " + strSQ);
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

