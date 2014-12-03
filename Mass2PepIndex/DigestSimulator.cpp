#include "DigestSimulator.h"

namespace Mass2PepIndex
{
double nPrime[100] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541};
CDigestSimulator::CDigestSimulator()
{
	try
	{
		m_vdfPrime.reserve(DEFAULT_MAX_IDX_PEP_LEN);
		double df = 0;
		for(size_t t=0; t<100; ++t)
		{
			df = log(nPrime[t]);
			m_vdfPrime.push_back(df);
		}
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CDigestSimulator", "CDigestSimulator");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDigestSimulator", "CDigestSimulator", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

CDigestSimulator::~CDigestSimulator()
{
}

void	CDigestSimulator::EmptyContent(void)
{
	m_vsDatBlock.clear();	
}

vector<PEP_INFO> &	CDigestSimulator::GetvSubPeptides() 
{
	return m_vsDatBlock;
}

void	CDigestSimulator::DigestPros2(string& strProteinSQ,const proteomics_sdk::CEnzyme& Enzyme, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<PEP_INFO_EX_TAG>	& m_vsDat, size_t tProPos, size_t tProID, size_t tDatNum)
{
	try
	{
		EmptyContent();
		const bool bNTerm = Enzyme.GetIsNTerm(); 
		const string strCleave = Enzyme.GetCleaveString();
		const string strNotCleave = Enzyme.GetNotCleave();
		const std::string::size_type size = strProteinSQ.length();
			
		double lfMasses = 0;
	
		vector<size_t> vtCleaveSite;
		vector<double> vlfMass;
	
		vtCleaveSite.push_back(0);
		vlfMass.push_back(0);
		
		if(!bNTerm)
		{
			for(std::string::size_type i = 0; i < size; ++i)
			{
				lfMasses += pepFunc.GetfAAMass(strProteinSQ[i], bMono);
				
				if(string::npos != strCleave.find(strProteinSQ[i]))
				{
					//vlfSiteMass.push_back(m_pepFunc.GetfAAMass(strProteinSQ[i], bMono));
					if(strNotCleave.find(strProteinSQ[i + 1]) == string::npos)		//C-term
					{
						vtCleaveSite.push_back(i+1);
						vlfMass.push_back(lfMasses);
	
						double lfTempMass = vlfMass[vlfMass.size()-1] - vlfMass[vlfMass.size()-2];
						size_t tTempLen = (vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-2]);
						
						if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
						{
							m_stDatExTagBlock.tPos = vtCleaveSite[vtCleaveSite.size()-2] + tProPos;
							m_stDatExTagBlock.cMiss = 0;
							m_stDatExTagBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER+0.5);
						
							m_stDatExTagBlock.cLen = (unsigned char)(tTempLen);
							m_stDatExTagBlock.tProID = tProID;
							m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
							
							//Godel code
							double dLogPrime = 0;
							//long long int lltLogPrime = 0;
							for(size_t t=0; t<m_stDatExTagBlock.cLen; ++t)
							{
								dLogPrime += (strProteinSQ[i-t] - 'A') * m_vdfPrime[t];
							}
							
							//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
							m_stDatExTagBlock.lfTag = dLogPrime;
	
							m_vsDat.push_back(m_stDatExTagBlock);				
						}
						
						//Miss cleave
						for (size_t j=1; j<=max_miss_cleavages&&j<(vtCleaveSite.size()-1); ++j)
						{				
							double lfTempMass = vlfMass[vlfMass.size()-1] - vlfMass[vlfMass.size()-2-j];
							size_t tTempLen = vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-2-j];
							
							if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
							{
								m_stDatExTagBlock.tPos = vtCleaveSite[vtCleaveSite.size()-2-j] + tProPos;
								m_stDatExTagBlock.cMiss = j;
								m_stDatExTagBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER+0.5);
								//m_vsDatPep2ProBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER);
								m_stDatExTagBlock.cLen = (unsigned char)(tTempLen);
								m_stDatExTagBlock.tProID = tProID;
								m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
								
								//Godel code
								double dLogPrime = 0;
								//long long int lltLogPrime = 0;
								for(size_t t=0; t<m_stDatExTagBlock.cLen; ++t)
								{
									dLogPrime += (strProteinSQ[i-t] - 'A') * m_vdfPrime[t];
								}
								
								//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
								m_stDatExTagBlock.lfTag = dLogPrime;
								
								m_vsDat.push_back(m_stDatExTagBlock);	
							}
						}
						
					}//end of if				
				}//end of if
				else if( i == (size-1) )
				{			
					double lfTempMass = lfMasses - vlfMass[vlfMass.size()-1];
					size_t tTempLen = i - vtCleaveSite[vtCleaveSite.size()-1] +1;
					
					if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
					{
						m_stDatExTagBlock.tPos = vtCleaveSite[vtCleaveSite.size()-1] + tProPos;
						m_stDatExTagBlock.cMiss = 0;
						m_stDatExTagBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER+0.5);
						
						m_stDatExTagBlock.cLen = (unsigned char)(tTempLen);
						m_stDatExTagBlock.tProID = tProID;
						m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
						
						//Godel code
						double dLogPrime = 0;
						//long long int lltLogPrime = 0;
						for(size_t t=0; t<m_stDatExTagBlock.cLen; ++t)
						{
							dLogPrime += (strProteinSQ[i-t] - 'A') * m_vdfPrime[t];
						}
						
						//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
						m_stDatExTagBlock.lfTag = dLogPrime;
						
						m_vsDat.push_back(m_stDatExTagBlock);		
					}
					
					//Miss cleave
					for (size_t j=1; j<=max_miss_cleavages&&j<(vtCleaveSite.size()); ++j)
					{			
						double lfTempMass = lfMasses - vlfMass[vlfMass.size()-1-j];
						size_t tTempLen = i - vtCleaveSite[vtCleaveSite.size()-1-j]+1;
							
						if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
						{
							m_stDatExTagBlock.tPos = vtCleaveSite[vtCleaveSite.size()-1-j] + tProPos;
							m_stDatExTagBlock.cMiss = j;
							m_stDatExTagBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER+0.5);
							//m_vsDatPep2ProBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER);
							m_stDatExTagBlock.cLen = (unsigned char)(tTempLen);
							m_stDatExTagBlock.tProID = tProID;
							m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
							
							//Godel code
							double dLogPrime = 0;
							//long long int lltLogPrime = 0;
							for(size_t t=0; t<m_stDatExTagBlock.cLen; ++t)
							{
								dLogPrime += (strProteinSQ[i-t] - 'A') * m_vdfPrime[t];
							}
							
							//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
							m_stDatExTagBlock.lfTag = dLogPrime;
							
							m_vsDat.push_back(m_stDatExTagBlock);		
						}
					}
	
					vtCleaveSite.push_back(i+1);
					vlfMass.push_back(lfMasses);				
				}//end of else if
			}//end of for
		}//end of if(!bNTerm)
		else
		{
			for(std::string::size_type i = 0; i < size; ++i)
			{
				lfMasses += pepFunc.GetfAAMass(strProteinSQ[i], bMono);

				if(string::npos != strCleave.find(strProteinSQ[i]))
				{
					//vlfSiteMass.push_back(m_pepFunc.GetfAAMass(strProteinSQ[i], bMono));
					if(strNotCleave.find(strProteinSQ[i + 1]) == string::npos)		//N-term
					{
						vtCleaveSite.push_back(i);
						vlfMass.push_back(lfMasses - pepFunc.GetfAAMass(strProteinSQ[i], bMono));
	
						double lfTempMass = vlfMass[vlfMass.size()-1] - vlfMass[vlfMass.size()-2];
						size_t tTempLen = vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-2];
	
						if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
						{
							m_stDatExTagBlock.tPos = vtCleaveSite[vtCleaveSite.size()-2] + tProPos;
							m_stDatExTagBlock.cMiss = 0;
							m_stDatExTagBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER+0.5);
						
							m_stDatExTagBlock.cLen = (unsigned char)(tTempLen);
							m_stDatExTagBlock.tProID = tProID;
							m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
							
							//Godel code
							double dLogPrime = 0;
							//long long int lltLogPrime = 0;
							for(size_t t=0; t<m_stDatExTagBlock.cLen; ++t)
							{
								dLogPrime += (strProteinSQ[i-t] - 'A') * m_vdfPrime[t];
							}
							
							//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
							m_stDatExTagBlock.lfTag = dLogPrime;
	
							m_vsDat.push_back(m_stDatExTagBlock);					
						}
	
						//Miss cleave
						for (size_t j=1; j<=max_miss_cleavages&&j<(vtCleaveSite.size()-1); ++j)
						{				
							double lfTempMass = vlfMass[vlfMass.size()-1] - vlfMass[vlfMass.size()-2-j];
							size_t tTempLen = vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-2-j];
	
							if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
							{
								m_stDatExTagBlock.tPos = vtCleaveSite[vtCleaveSite.size()-2-j] + tProPos;
								m_stDatExTagBlock.cMiss = j;
								m_stDatExTagBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER+0.5);
								//m_vsDatPep2ProBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER);
								m_stDatExTagBlock.cLen = (unsigned char)(tTempLen);
								m_stDatExTagBlock.tProID = tProID;
								m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
								
								//Godel code
								double dLogPrime = 0;
								//long long int lltLogPrime = 0;
								for(size_t t=0; t<m_stDatExTagBlock.cLen; ++t)
								{
									dLogPrime += (strProteinSQ[i-t] - 'A') * m_vdfPrime[t];
								}
								
								//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
								m_stDatExTagBlock.lfTag = dLogPrime;
								
								m_vsDat.push_back(m_stDatExTagBlock);				
							}
						}
					}
				}

				else if( i == (size-1) )
				{			
					double lfTempMass = lfMasses - vlfMass[vlfMass.size()-1];
					size_t tTempLen = i - vtCleaveSite[vtCleaveSite.size()-1] +1;
	
					if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
					{
						m_stDatExTagBlock.tPos = vtCleaveSite[vtCleaveSite.size()-1] + tProPos;
						m_stDatExTagBlock.cMiss = 0;
						m_stDatExTagBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER+0.5);
						
						m_stDatExTagBlock.cLen = (unsigned char)(tTempLen);
						m_stDatExTagBlock.tProID = tProID;
						m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
						
						//Godel code
						double dLogPrime = 0;
						//long long int lltLogPrime = 0;
						for(size_t t=0; t<m_stDatExTagBlock.cLen; ++t)
						{
							dLogPrime += (strProteinSQ[i-t] - 'A') * m_vdfPrime[t];
						}
						
						//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
						m_stDatExTagBlock.lfTag = dLogPrime;
						
						m_vsDat.push_back(m_stDatExTagBlock);					
					}
	
					//Miss cleave
					for (size_t j=1; j<=max_miss_cleavages&&j<(vtCleaveSite.size()); ++j)
					{				
						double lfTempMass = lfMasses - vlfMass[vlfMass.size()-1-j];
						size_t tTempLen = i - vtCleaveSite[vtCleaveSite.size()-1-j]+1;
	
						if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
						{
							m_stDatExTagBlock.tPos = vtCleaveSite[vtCleaveSite.size()-1-j] + tProPos;
							m_stDatExTagBlock.cMiss = j;
							m_stDatExTagBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER+0.5);
							//m_vsDatPep2ProBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER);
							m_stDatExTagBlock.cLen = (unsigned char)(tTempLen);
							m_stDatExTagBlock.tProID = tProID;
							m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
							
							//Godel code
							double dLogPrime = 0;
							//long long int lltLogPrime = 0;
							for(size_t t=0; t<m_stDatExTagBlock.cLen; ++t)
							{
								dLogPrime += (strProteinSQ[i-t] - 'A') * m_vdfPrime[t];
							}
							
							//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
							m_stDatExTagBlock.lfTag = dLogPrime;
							
							m_vsDat.push_back(m_stDatExTagBlock);			
						}
					}
					
					vtCleaveSite.push_back(i);
					vlfMass.push_back(lfMasses);				
				}
			}	
		}
		
		if('M' == strProteinSQ[0])
		{
			size_t tTempLen = 0;
			double lfTempMass = 0;
			if(vtCleaveSite.size() > 2)
			{
				tTempLen = vtCleaveSite[1]-1;
				lfTempMass = vlfMass[1] - pepFunc.GetfAAMass(strProteinSQ[0], bMono);
			}
			else
			{
				tTempLen = strProteinSQ.size() - 1;
				lfTempMass = lfMasses - pepFunc.GetfAAMass(strProteinSQ[0], bMono);
			}
	
			if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
			{
				m_stDatExTagBlock.tPos = 1 + tProPos;
				m_stDatExTagBlock.cMiss = 0;
				m_stDatExTagBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER+0.5);
				m_stDatExTagBlock.cLen = (unsigned char)(tTempLen);
				m_stDatExTagBlock.tProID = tProID;
				m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
	
				//Godel code
				double dLogPrime = 0;
				//long long int lltLogPrime = 0;
				for(size_t t=0; t<m_stDatExTagBlock.cLen; ++t)
				{
					dLogPrime += (strProteinSQ[tTempLen-t] - 'A') * m_vdfPrime[t];
				}
				
				//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
				m_stDatExTagBlock.lfTag = dLogPrime;
				
				m_vsDat.push_back(m_stDatExTagBlock);				
			}
	
			for (size_t j=1; j<=max_miss_cleavages&&j<(vtCleaveSite.size()-1); ++j)
			{	
//				if(strProteinSQ.size() == vtCleaveSite[j+1])
				double lfTempMass = vlfMass[j+1] - pepFunc.GetfAAMass(strProteinSQ[0], bMono);
				size_t tTempLen = vtCleaveSite[j+1] -1;
				
				if(strProteinSQ.size() - 1 == vtCleaveSite[j+1])//czhou
				{
					++tTempLen;
//					lfTempMass += pepFunc.GetfAAMass(strProteinSQ[vtCleaveSite[j+1]], bMono);
				}
				
				if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
				{
					m_stDatExTagBlock.tPos = 1 + tProPos;
					m_stDatExTagBlock.cMiss = j;
					m_stDatExTagBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER+0.5);
					//m_vsDatPep2ProBlock.tMass = (size_t)(lfTempMass*DEFAULT_MULTPLIER);
					m_stDatExTagBlock.cLen = (unsigned char)(tTempLen);
					m_stDatExTagBlock.tProID = tProID;
					m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
					
					//Godel code
					double dLogPrime = 0;
					//long long int lltLogPrime = 0;
					for(size_t t=0; t<m_stDatExTagBlock.cLen; ++t)
					{
						dLogPrime += (strProteinSQ[tTempLen-t] - 'A') * m_vdfPrime[t];
					}
					
					//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
					m_stDatExTagBlock.lfTag = dLogPrime;
					
					m_vsDat.push_back(m_stDatExTagBlock);	
				}
			}
		}
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CDigestSimulator", "DigestPros2");
		char str[20];
//		itoa(tProID, str, 10);
		sprintf(str,"%d",tProID);
		info.Append("tProID = " + (string)(str));
		info.Append("strProteinSQ = " + strProteinSQ);
//		itoa(m_vsDat.size(), str, 10);
		sprintf(str,"%d",m_vsDat.size());
		info.Append("m_vsDat.size() = " + (string)(str));
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDigestSimulator", "DigestPros2", "Caught an unkown exception!");
		char str[20];
//		itoa(tProID, str, 10);
		sprintf(str,"%d", tProID);
		info.Append("tProID = " + (string)(str));
		info.Append("strProteinSQ = " + strProteinSQ);
//		itoa(m_vsDat.size(), str, 10);
		sprintf(str,"%d",m_vsDat.size());
		info.Append("m_vsDat.size() = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}



void	CDigestSimulator::SemiDigestPros1(string& strProteinSQ,const proteomics_sdk::CEnzyme& Enzyme, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<PEP_INFO_EX_TAG>	& m_vsDat, size_t tProPos, size_t tProID, size_t tDatNum)
{
	try
	{
		EmptyContent();
	
		const bool bNTerm = Enzyme.GetIsNTerm(); 
		const string strCleave = Enzyme.GetCleaveString();
		const string strNotCleave = Enzyme.GetNotCleave();
		const std::string::size_type size = strProteinSQ.length();
	
		double lfMasses = 0;
	
		vector<size_t> vtCleaveSite;
		vector<double> vlfSumMass;
	
		vtCleaveSite.push_back(0);
		vlfSumMass.push_back(0);
	
		if (!bNTerm)
		{
			for(std::string::size_type i = 0; i < size; ++i)
			{
				lfMasses += pepFunc.GetfAAMass(strProteinSQ[i], bMono);
				vlfSumMass.push_back(lfMasses);
	
				
				if(string::npos != strCleave.find(strProteinSQ[i]))
				{
					if(strNotCleave.find(strProteinSQ[i + 1]) == string::npos)		//C-term
					{
						vtCleaveSite.push_back(i+1);
	
						for (size_t j=0; j<=max_miss_cleavages&&j<(vtCleaveSite.size()-1); ++j)
						{
							
							size_t t=(vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-1-j] + 1);
							t = (m_pFilter->GetMinPepLength()> t)?m_pFilter->GetMinPepLength():t;
	
							for (; (t<=(vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-2-j])) && (t<= m_pFilter->GetMaxPepLength()); ++t)
							{
								double lfTempMass = vlfSumMass[i+1]-vlfSumMass[i-t+1];					
	
								if ( (!m_pFilter) || (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>lfTempMass ) )
								{
									m_stDatExTagBlock.tMass = (size_t)( lfTempMass * DEFAULT_MULTPLIER + 0.5);
									m_stDatExTagBlock.tPos = (i-t+1) + tProPos;
									m_stDatExTagBlock.cMiss = j;
									m_stDatExTagBlock.cLen = (unsigned char)(t);
									m_stDatExTagBlock.tProID = tProID;
									m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
									
									//Godel code
									double dLogPrime = 0;
									//long long int lltLogPrime = 0;
									for(size_t t1=0; t1<m_stDatExTagBlock.cLen; ++t1)
									{
										dLogPrime += (strProteinSQ[i-t1] - 'A') * m_vdfPrime[t1];
									}
									
									//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
									m_stDatExTagBlock.lfTag = dLogPrime;
									
									m_vsDat.push_back(m_stDatExTagBlock);
								}				
							}
						}
					}
				}
	
				else
				{			
					for (size_t j=0; j<=max_miss_cleavages && j<vtCleaveSite.size(); ++j)
					{
						if ( (i - vtCleaveSite[vtCleaveSite.size() -1 - j] + 1) >= m_pFilter->GetMinPepLength() &&  (i - vtCleaveSite[vtCleaveSite.size() -1 - j] + 1) <= m_pFilter->GetMaxPepLength())
						{
							//double d1 = vlfSumMass[i];
							//double d2 = vlfSumMass[vtCleaveSite[vtCleaveSite.size() - 1 - j]];
							//size_t t = vtCleaveSite[vtCleaveSite.size() - 1 - j];
	
							double lfTempMass = vlfSumMass[i+1]-vlfSumMass[vtCleaveSite[vtCleaveSite.size() - 1 - j]];					
	
							if ( (!m_pFilter) || (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>lfTempMass ) )
							{
								m_stDatExTagBlock.tMass = (size_t)( lfTempMass * DEFAULT_MULTPLIER + 0.5);
								m_stDatExTagBlock.tPos = vtCleaveSite[vtCleaveSite.size() -1 - j]  + tProPos;
								m_stDatExTagBlock.cMiss = j;
								m_stDatExTagBlock.cLen = (unsigned char)(i - vtCleaveSite[vtCleaveSite.size() -1 - j] +1);
								m_stDatExTagBlock.tProID = tProID;
								m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
								
								//Godel code
								double dLogPrime = 0;
								//long long int lltLogPrime = 0;
								for(size_t t1=0; t1<m_stDatExTagBlock.cLen; ++t1)
								{
									dLogPrime += (strProteinSQ[i-t1] - 'A') * m_vdfPrime[t1];
								}
								
								//m_stDatExTagBlock.tTag = (size_t)(dLogPrime*DEFAULT_MULTPLIER+0.5);
								m_stDatExTagBlock.lfTag = dLogPrime;
								
								m_vsDat.push_back(m_stDatExTagBlock);
							}	
						}
										
					}
				}//end of else
			}//end of for
		}//end of if(!bNterm)	
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CDigestSimulator", "SemiDigestPros1");
		char str[20];
//		itoa(tProID, str, 10);
		sprintf(str,"%d",tProID);
		info.Append("tProID = " + (string)(str));
		info.Append("strProteinSQ = " + strProteinSQ);
//		itoa(m_vsDat.size(), str, 10);
		sprintf(str,"%d",m_vsDat.size());
		info.Append("m_vsDat.size() = " + (string)(str));
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDigestSimulator", "SemiDigestPros1", "Caught an unkown exception!");
		char str[20];
//		itoa(tProID, str, 10);
		sprintf(str,"%d",tProID);
		info.Append("tProID = " + (string)(str));
		info.Append("strProteinSQ = " + strProteinSQ);
//		itoa(m_vsDat.size(), str, 10);
		sprintf(str,"%d",m_vsDat.size());
		info.Append("m_vsDat.size() = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void	CDigestSimulator::NoneEnzymeDigestPros1(string& strProteinSQ, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<PEP_INFO_EX_TAG>	& m_vsDat, size_t tProPos, size_t tProID, size_t tDatNum)
{
	try
	{
		EmptyContent();

		const std::string::size_type size = strProteinSQ.length();	
		double lfMasses = 0;

		vector<double> vlfSumMass;
		vlfSumMass.push_back(0);

		for(std::string::size_type i = 0; i < size; ++i)
		{
			lfMasses += pepFunc.GetfAAMass(strProteinSQ[i], bMono);
			vlfSumMass.push_back(lfMasses);			

			for (size_t t = m_pFilter->GetMinPepLength(); t<m_pFilter->GetMaxPepLength() && t<=i+1; ++t)
			{
				double lfTempMass = vlfSumMass[i+1]-vlfSumMass[i-t+1];					

				if ( (!m_pFilter) || (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>lfTempMass ) )
				{
					m_stDatExTagBlock.tMass = (size_t)( lfTempMass * DEFAULT_MULTPLIER + 0.5);
					m_stDatExTagBlock.tPos = (i-t+1) + tProPos;
					m_stDatExTagBlock.cMiss = 0;
					m_stDatExTagBlock.cLen = (unsigned char)(t);
					m_stDatExTagBlock.tProID = tProID;
					m_stDatExTagBlock.cDatNum = (unsigned char)(tDatNum);
					
					//Godel code
					double dLogPrime = 0;
					for(size_t t1=0; t1<m_stDatExTagBlock.cLen; ++t1)
					{
						dLogPrime += (strProteinSQ[i-t1] - 'A') * m_vdfPrime[t1];
					}
					m_stDatExTagBlock.lfTag = dLogPrime;
					
					m_vsDat.push_back(m_stDatExTagBlock);
				}				
			}					
		}//end of for			
	}
	
	catch(runtime_error & e)
	{
		CErrInfo info("CDigestSimulator", "NoneEnzymeDigestPros1");
		char str[20];
//		itoa(tProID, str, 10);
		sprintf(str,"%d",tProID);
		info.Append("tProID = " + (string)(str));
		info.Append("strProteinSQ = " + strProteinSQ);
//		itoa(m_vsDat.size(), str, 10);
		sprintf(str,"%d",m_vsDat.size());
		info.Append("m_vsDat.size() = " + (string)(str));
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDigestSimulator", "NoneEnzymeDigestPros1", "Caught an unkown exception!");
		char str[20];
//		itoa(tProID, str, 10);
		sprintf(str,"%d",tProID);
		info.Append("tProID = " + (string)(str));
		info.Append("strProteinSQ = " + strProteinSQ);
//		itoa(m_vsDat.size(), str, 10);
		sprintf(str,"%d",m_vsDat.size());
		info.Append("m_vsDat.size() = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void	CDigestSimulator::ComputePepNumAtEachMassBin(string& strProteinSQ,const proteomics_sdk::CEnzyme& Enzyme, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<size_t> & vtPepMassNum, size_t tBitNum)
{
	try
	{
		EmptyContent();
	
		const bool bNTerm = Enzyme.GetIsNTerm(); 
		const string strCleave = Enzyme.GetCleaveString();
		const string strNotCleave = Enzyme.GetNotCleave();
		const std::string::size_type size = strProteinSQ.length();
	
		double lfMasses = 0;
	
		vector<size_t> vtCleaveSite;
		vector<double> vlfMass;
		vector<double> vlfSiteMass;
	
		vtCleaveSite.push_back(0);
		vlfMass.push_back(0);
	
		if (!bNTerm)
		{
			for(std::string::size_type i = 0; i < size; ++i)
			{
				lfMasses += pepFunc.GetfAAMass(strProteinSQ[i], bMono);
				
				if(string::npos != strCleave.find(strProteinSQ[i]))
				{
					//vlfSiteMass.push_back(m_pepFunc.GetfAAMass(strProteinSQ[i], bMono));
					if(strNotCleave.find(strProteinSQ[i + 1]) == string::npos)		//C-term
					{
						vtCleaveSite.push_back(i+1);
						vlfMass.push_back(lfMasses);
	
						double lfTempMass = vlfMass[vlfMass.size()-1] - vlfMass[vlfMass.size()-2];
						size_t tTempLen = (vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-2]);
	
						if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
						{
							size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
							++vtPepMassNum[ tMass >> tBitNum];	
							
						}
	
						//Miss cleave
						for (size_t j=1; j<=max_miss_cleavages&&j<(vtCleaveSite.size()-1); ++j)
						{				
							double lfTempMass = vlfMass[vlfMass.size()-1] - vlfMass[vlfMass.size()-2-j];
							size_t tTempLen = vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-2-j];
	
							if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
							{
								size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
								++vtPepMassNum[ tMass >> tBitNum];	
								
							}
						}
	
					}
				}
				
				else if( i == (size-1) )
				{			
					double lfTempMass = lfMasses - vlfMass[vlfMass.size()-1];
					size_t tTempLen = i - vtCleaveSite[vtCleaveSite.size()-1] +1;
	
					if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
					{
						size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
						++vtPepMassNum[ tMass >> tBitNum];	
						
					}
	
					//Miss cleave
					for (size_t j=1; j<=max_miss_cleavages&&j<(vtCleaveSite.size()); ++j)
					{				
						double lfTempMass = lfMasses - vlfMass[vlfMass.size()-1-j];
						size_t tTempLen = i - vtCleaveSite[vtCleaveSite.size()-1-j]+1;
	
						if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
						{
							size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
							++vtPepMassNum[ tMass >> tBitNum];	
						
						}
					}
					
					vtCleaveSite.push_back(i+1);
					vlfMass.push_back(lfMasses);	
				}
			}
		}
	
		else
		{
			for(std::string::size_type i = 0; i < size; ++i)
			{
				lfMasses += pepFunc.GetfAAMass(strProteinSQ[i], bMono);
				
				if(string::npos != strCleave.find(strProteinSQ[i]))
				{
					//vlfSiteMass.push_back(m_pepFunc.GetfAAMass(strProteinSQ[i], bMono));
					if(strNotCleave.find(strProteinSQ[i + 1]) == string::npos)		//N-term
					{
						vtCleaveSite.push_back(i);
						vlfMass.push_back(lfMasses - pepFunc.GetfAAMass(strProteinSQ[i], bMono));
	
						double lfTempMass = vlfMass[vlfMass.size()-1] - vlfMass[vlfMass.size()-2];
						size_t tTempLen = vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-2];
	
						if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
						{
							size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
							++vtPepMassNum[ tMass >> tBitNum];			
						}
	
						//Miss cleave
						for (size_t j=1; j<=max_miss_cleavages&&j<(vtCleaveSite.size()-1); ++j)
						{				
							double lfTempMass = vlfMass[vlfMass.size()-1] - vlfMass[vlfMass.size()-2-j];
							size_t tTempLen = vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-2-j];
	
							if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
							{
								size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
								++vtPepMassNum[ tMass >> tBitNum];				
							}
						}
					}
				}
				
				else if( i == (size-1) )
				{			
					double lfTempMass = lfMasses - vlfMass[vlfMass.size()-1];
					size_t tTempLen = i - vtCleaveSite[vtCleaveSite.size()-1] +1;
	
					if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
					{
						size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
						++vtPepMassNum[ tMass >> tBitNum];			
					}
	
					//Miss cleave
					for (size_t j=1; j<=max_miss_cleavages&&j<(vtCleaveSite.size()); ++j)
					{				
						double lfTempMass = lfMasses - vlfMass[vlfMass.size()-1-j];
						size_t tTempLen = i - vtCleaveSite[vtCleaveSite.size()-1-j]+1;
	
						if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
						{
							size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
							++vtPepMassNum[ tMass >> tBitNum];				
						}
					}
	
				}
			}
		}
		
		
		if('M' == strProteinSQ[0])
		{
			size_t tTempLen = 0;
			double lfTempMass = 0;
			if(vtCleaveSite.size() > 2)
			{
				tTempLen = vtCleaveSite[1]-1;
				lfTempMass = vlfMass[1] - pepFunc.GetfAAMass(strProteinSQ[0], bMono);
			}
			else
			{
				tTempLen = strProteinSQ.size() - 1;
				lfTempMass = lfMasses - pepFunc.GetfAAMass(strProteinSQ[0], bMono);
			}
				
			if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
			{
				size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
				++vtPepMassNum[ tMass >> tBitNum];	
			}
	
			for (size_t j=1; j<=max_miss_cleavages&&j<(vtCleaveSite.size()-1); ++j)
			{				
				double lfTempMass = vlfMass[j+1] - pepFunc.GetfAAMass(strProteinSQ[0], bMono);
				size_t tTempLen = vtCleaveSite[j+1]-1;
				
				if ( (!m_pFilter) || (m_pFilter->GetMinPepLength()<= tTempLen && m_pFilter->GetMaxPepLength() >= tTempLen) && (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>=lfTempMass ) )
				{
					size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
					++vtPepMassNum[ tMass >> tBitNum];	
				}
			}
		}
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CDigestSimulator", "ComputePepNumAtEachMassBin");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDigestSimulator", "ComputePepNumAtEachMassBin", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}

void	CDigestSimulator::ComputeSemiPepNumAtEachMassBin(string& strProteinSQ,const proteomics_sdk::CEnzyme& Enzyme, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<size_t> & vtPepMassNum, size_t tBitNum)
{
	try
	{
		const bool bNTerm = Enzyme.GetIsNTerm(); 
		const string strCleave = Enzyme.GetCleaveString();
		const string strNotCleave = Enzyme.GetNotCleave();
		const std::string::size_type size = strProteinSQ.length();
	
		double lfMasses = 0;
	
		vector<size_t> vtCleaveSite;
		vector<double> vlfSumMass;
	
		vtCleaveSite.push_back(0);
		vlfSumMass.push_back(0);
	
		if (!bNTerm)
		{
			for(std::string::size_type i = 0; i < size; ++i)
			{
				lfMasses += pepFunc.GetfAAMass(strProteinSQ[i], bMono);
				vlfSumMass.push_back(lfMasses);
	
				
				if(string::npos != strCleave.find(strProteinSQ[i]))
				{
					if(strNotCleave.find(strProteinSQ[i + 1]) == string::npos)		//C-term
					{
						vtCleaveSite.push_back(i+1);
	
						for (size_t j=0; j<=max_miss_cleavages&&j<(vtCleaveSite.size()-1); ++j)
						{
							
							size_t t=(vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-1-j] + 1);
							t = (m_pFilter->GetMinPepLength()> t)?m_pFilter->GetMinPepLength():t;
	
							for (; (t<=(vtCleaveSite[vtCleaveSite.size()-1] - vtCleaveSite[vtCleaveSite.size()-2-j])) && (t<= m_pFilter->GetMaxPepLength()); ++t)
							{
								double lfTempMass = vlfSumMass[i+1]-vlfSumMass[i-t+1];					
	
								if ( (!m_pFilter) || (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>lfTempMass ) )
								{
									size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
									//cout << tMass << ' ' << tBitNum << ' ' << (tMass >> tBitNum) << endl;
									++vtPepMassNum[ tMass >> tBitNum];
								}				
							}
						}
					}
				}
	
				else
				{			
					for (size_t j=0; j<=max_miss_cleavages && j<vtCleaveSite.size(); ++j)
					{
						if ( (i - vtCleaveSite[vtCleaveSite.size() -1 - j] + 1) >= m_pFilter->GetMinPepLength() &&  (i - vtCleaveSite[vtCleaveSite.size() -1 - j] + 1) <= m_pFilter->GetMaxPepLength())
						{
							double lfTempMass = vlfSumMass[i+1]-vlfSumMass[vtCleaveSite[vtCleaveSite.size() - 1 - j]];					
	
							if ( (!m_pFilter) || (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>lfTempMass ) )
							{
								size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);
								
								//cout << tMass << ' ' << tBitNum << ' ' << (tMass >> tBitNum) << endl;
								++vtPepMassNum[ tMass >> tBitNum];
								
							}	
						}
	
					}
				}//end of else
			}//end of for
		}//end of if(!bNterm)
	}
	catch(runtime_error & e)
	{
		CErrInfo info("CDigestSimulator", "ComputeSemiPepNumAtEachMassBin");
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDigestSimulator", "ComputeSemiPepNumAtEachMassBin", "Caught an unkown exception!");
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}
void	CDigestSimulator::ComputeNoneEnzymeAtEachMassBin(string& strProteinSQ, CPepCalcFunc &pepFunc, const bool bMono, const size_t max_miss_cleavages, CPepFilter* m_pFilter, vector<size_t> & vtPepMassNum, size_t tBitNum)
{
	try
	{
		EmptyContent();

		const std::string::size_type size = strProteinSQ.length();	
		double lfMasses = 0;

		vector<double> vlfSumMass;
		vlfSumMass.push_back(0);

		for(std::string::size_type i = 0; i < size; ++i)
		{
			lfMasses += pepFunc.GetfAAMass(strProteinSQ[i], bMono);
			vlfSumMass.push_back(lfMasses);			

			for (size_t t = m_pFilter->GetMinPepLength(); t<m_pFilter->GetMaxPepLength() && t<=i+1; ++t)
			{
				double lfTempMass = vlfSumMass[i+1]-vlfSumMass[i-t+1];					

				if ( (!m_pFilter) || (m_pFilter->GetMinPepMass()<= lfTempMass && m_pFilter->GetMaxPepMass()>lfTempMass ) )
				{
					size_t tMass = (size_t)( (lfTempMass - m_pFilter->GetMinPepMass() ) * DEFAULT_MULTPLIER + 0.5);					
					++vtPepMassNum[ tMass >> tBitNum];
				}				
			}					
		}//end of for			
	}
	
	catch(runtime_error & e)
	{
		CErrInfo info("CDigestSimulator", "ComputeNoneEnzymeAtEachMassBin");
		char str[20];
		info.Append("tProID = " + (string)(str));
		info.Append("strProteinSQ = " + strProteinSQ);
		info.Append("m_vsDat.size() = " + (string)(str));
		cerr << info.Get(e) << endl;
		throw runtime_error(info.Get(e).c_str());
	}
	catch(...)
	{
		CErrInfo info("CDigestSimulator", "ComputeNoneEnzymeAtEachMassBin", "Caught an unkown exception!");
		char str[20];
		info.Append("tProID = " + (string)(str));
		info.Append("strProteinSQ = " + strProteinSQ);
		info.Append("m_vsDat.size() = " + (string)(str));
		cerr << info.Get() << endl;
		throw runtime_error(info.Get().c_str());
	}
}
}
