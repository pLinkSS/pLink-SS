
#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include "../include/sdk.h"
#include "../include/interface.h"
#include "XLinkPeptideEvaluater.h"

using namespace std;
using namespace proteomics_sdk;
#define yfu_new

//#define SCOREASEVALUE

//#define EMILY_DEBUG

CXLinkPeptideEvaluater::CXLinkPeptideEvaluater(void)
	:m_pResult(NULL)
{
	m_Scoefficient[0] = 0;
	m_Scoefficient[1] = 0;
	m_nEvalueNo = 0;

	for(int i = 0; i< 6; ++i)
		m_DT[i] = 0;
}

CXLinkPeptideEvaluater::~CXLinkPeptideEvaluater(void)
{
	Close();
}

void CXLinkPeptideEvaluater::Init(const CCondition &condition)
{
	if(!condition.ValidateAll())
	{
		CErrInfo info("CXLinkPeptideEvaluater", "Init", "The validation of inputing condition is failed.");
		throw invalid_argument(info.Get().c_str());
	}

	Close();
	m_Condition = condition;
}

void CXLinkPeptideEvaluater::Run(CSpectrum & spectrum, CXLinkMatchResult * pResult)
{
#ifdef SCOREASEVALUE
	for (int i = 0; i < pResult->m_vPeptideResults.size(); i++)
		pResult->m_vPeptideResults[i].m_lfEvalue = 1 - pResult->m_vPeptideResults[i].m_lfScore;
	return;
#endif

	if((NULL == pResult)||(pResult->m_vlfScores.empty()) || pResult->m_vPeptideResults.empty())
	{
		return;
	}
	
	sort(pResult->m_vPeptideResults.begin(), pResult->m_vPeptideResults.end(), CXLinkPepResult::Score_Greater);
	m_pResult = pResult;

	sort( pResult->m_vlfScores.begin() , pResult->m_vlfScores.end());

#ifdef EMILY_DEBUG
	FILE * fp;
	fp = fopen("evalue.txt","w");
	
	for(size_t i = 0;i < pResult->m_vlfScores.size() ; ++ i)
	{
		fprintf(fp,"%f\n",pResult->m_vlfScores[i]);
	}
	fclose(fp);
#endif	
	
	/*
	 * for refined - search load coef
	 * 
	 */
	if(spectrum.m_stEVCoef.lfCoef0 == 0 && spectrum.m_stEVCoef.lfCoef1 == 0)
	{
		_FitHighScores();
	}
	else
	{
		m_Scoefficient[0] = spectrum.m_stEVCoef.lfCoef0 ;
		m_Scoefficient[1] = spectrum.m_stEVCoef.lfCoef1 ;
	}

	/*
	 * for refined - search save coef
	 */
	
	spectrum.m_stEVCoef.lfCoef0 = m_Scoefficient[0];
	spectrum.m_stEVCoef.lfCoef1 = m_Scoefficient[1];
	spectrum.m_nEvalueNo = m_nEvalueNo;
	
	_ComputeExpectation(spectrum);

#ifdef EMILY_DEBUG
	FILE * fp;
	string strTmpFile("spec_score\\");
	strTmpFile += spectrum.m_strFilePath;
	strTmpFile += ".txt";
	fp = fopen(strTmpFile.c_str(),"w");
	if(fp!=NULL)
	{
		fprintf(fp,"%u\n",pResult->m_tRealCandidate);
		
		int i=0;
		for(i=0;i<pResult->m_vPeptideResults.size() && pResult->m_vPeptideResults[i].m_bEV==false ;++i);
	
		if(i < pResult->m_vPeptideResults.size())
		{
			fprintf(fp,"%.5f\n",pResult->m_vPeptideResults[i].m_lfScore);
		}
		else
		{
			fprintf(fp,"0.0\n");
		}
		
		for(i=0;i<pResult->m_vlfScores.size();++i)
		{
			fprintf(fp,"%.5f\n",pResult->m_vlfScores[i]);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}
#endif
	

	m_pResult = NULL;
	
}

void CXLinkPeptideEvaluater::_FitHighScores(void)
{
	m_Scoefficient[0] = 0;
	m_Scoefficient[1] = 0;

	for(int i = 0; i< 6; ++i)
		m_DT[i] = 0;
	
	int itmp=0;

	/*
	for (size_t i = m_pResult->m_vlfScores.size()-1; i>0; --i )
	{
		if ( m_pResult->m_vlfScores[i-1]== m_pResult->m_vlfScores[i])
			++itmp;
		else break;
	}
	*/
	
	// remove scores of true results
	
	int i ;
	for (i = int(m_pResult->m_vlfScores.size()-1); i>0 && (size_t)itmp < m_pResult->m_vPeptideResults.size(); --i )
	{
		if ( m_pResult->m_vlfScores[i] == m_pResult->m_vPeptideResults[itmp].m_lfScore)
			++itmp;
		else
			break;
	}
	 
	// remove redundant values
	int nCurId = 0;
	int nLen = int(m_pResult->m_vlfScores.size()) - itmp;
	if(nLen > 0)
		nCurId = 1;
	
	for( i = 1 ;i < nLen ; ++ i)
	{
		if(m_pResult->m_vlfScores[i] != m_pResult->m_vlfScores[i-1])
		{
			m_pResult->m_vlfScores[nCurId++] = m_pResult->m_vlfScores[i];  
		}
	}
	nLen = nCurId;
	
	/*
	for (; i>0; --i )
	{
		if ( m_pResult->m_vlfScores[i-1]== m_pResult->m_vlfScores[i])
			++itmp;
		else break;
	}
	*/
	
	size_t start = (size_t)( nLen*0.9 + 1);
	size_t n = nLen - start;
	
	m_nEvalueNo = n;
	//	<< "n = " << n << endl;
	
	if ( n > 10)//the number of scores cannot be too small # changed from 20 to 10 at 2013.10.25
	{
		double * logs=new double[n];
		double * logx=new double[n];
		double temp=log(0.1/n);
		for (size_t i=0;i<n;++i)
		{
			logx[i] = m_pResult->m_vlfScores[start+i];
			logs[i] = temp+log((double)(n-i));
		}
		
		_Isqt(logx, logs, n);

		delete[] logs;
		delete[] logx;
	} else {
	}
	
}


void CXLinkPeptideEvaluater::_ComputeExpectation(const CSpectrum & spectrum)
{
	/*
	if(m_pResult->m_vlfScores.empty())
	{
		CErrInfo info("CXLinkPeptideEvaluater", "_ComputeExpectation", "The best match list is empty.");
		throw invalid_argument(info.Get().c_str());
	}
	*/
	
	size_t i = 0;
	for(i = 0;i < m_pResult->m_vPeptideResults.size();++i)
	{
		double lfScore = m_pResult->m_vPeptideResults[i].m_lfScore;
		if ( lfScore < 0 )
			m_pResult->m_vPeptideResults[i].m_lfEvalue = 10000.0;
		else
		{
			double logs = 0;
			logs = m_Scoefficient[1] * lfScore + m_Scoefficient[0];

#ifdef yfu_new
			m_pResult->m_vPeptideResults[i].m_lfEvalue = exp(logs) * (double)m_pResult ->m_tRealCandidate;
//			cout<<"In evalue :"<<endl;
//			cout<<"m_Scoefficient[1]: "<<m_Scoefficient[1]<<" lfScore: "<<lfScore
//					<<" m_Scoefficient[0]: "<<m_Scoefficient[0]<<endl;
//			cout<<"logs: m_Scoefficient[1] * lfScore + m_Scoefficient[0]: "
//					<<m_Scoefficient[1] * lfScore + m_Scoefficient[0]<<endl;
//			cout<<"m_pResult ->m_tRealCandidate: "<<m_pResult ->m_tRealCandidate
//					<<" exp(logs) * (double)m_pResult ->m_tRealCandidate "<<exp(logs) * (double)m_pResult ->m_tRealCandidate
//					<<endl;
#else
			m_pResult->m_vPeptideResults[i].m_lfEvalue = exp(logs) * (double)m_pResult ->m_vlfScores.size();
#endif
			
			// modify by emily for refined - search
			if(m_Condition.m_bRSLoadInfo == false && m_pResult ->m_tScore < 100 )
			{
				m_pResult->m_vPeptideResults[i].m_lfEvalue = 10000.0;
			}

		}
		
		if(m_pResult->m_vPeptideResults[i].m_lfEvalue >= m_Condition.m_lfMaxEV)//E_Value大于等于1的肯定不是
			m_pResult->m_vPeptideResults[i].m_bEV = false;
		else
			m_pResult->m_vPeptideResults[i].m_bEV = true;
		
	}
	
}

void CXLinkPeptideEvaluater::_Isqt( double x[], double y[], size_t n)
{ 
	double xx=0.0;
	double yy=0.0;
	for (size_t i=0; i<=n-1; ++i)
	{
		xx+=x[i];
		yy+=y[i];
	}
	xx/=n;
	yy/=n;

	double e=0.0;
	double f=0.0;
	for (size_t i=0; i<=n-1; ++i)
	{
		double q=x[i]-xx;
		e=e+q*q;
		f=f+q*(y[i]-yy);
	}
	
	m_Scoefficient[1]=f/e;
	m_Scoefficient[0]=yy - m_Scoefficient[1]*xx;
	
	double q=0.0;
	double u=0.0;
	double p=0.0;
	double umax=0.0;
	double umin=1.0e+30;
	for (size_t i=0; i<=n-1; ++i)
	{
		double s=m_Scoefficient[1]*x[i]+m_Scoefficient[0];
		q=q+(y[i]-s)*(y[i]-s);
		p=p+(s-yy)*(s-yy);
		e=fabs(y[i]-s);
		if (e>umax) umax=e;
		if (e<umin) umin=e;
		u+=e;
	}
	u/=n;

	m_DT[0]=q;
	m_DT[1]=sqrt(q/n);
	m_DT[2]=p;
	m_DT[3]=umax;
	m_DT[4]=umin;
	m_DT[5]=u;

}

void CXLinkPeptideEvaluater::Close(void)
{
	m_Condition.clear();

	m_Scoefficient[0] = 0;
	m_Scoefficient[1] = 0;

	for(int i = 0; i< 6; ++i)
		m_DT[i] = 0;
}
