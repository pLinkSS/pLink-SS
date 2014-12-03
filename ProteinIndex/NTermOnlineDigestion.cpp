#include "MSufSort.h"

#include "../include/sdk.h"
#include "../include/option.h"

#include "NTermOnlineDigestion.h"


namespace ProteinIndex
{

CNTermOnlineDigestion::CNTermOnlineDigestion()
{
	m_chDat = 0;LCP = 0;SA = 0;nClvy = 0;m_tMass =0; beg = 0;
}

CNTermOnlineDigestion::~CNTermOnlineDigestion()
{
	Close();
}

inline bool CNTermOnlineDigestion::_IsDiv(size_t p)
{
	return DIVPEP == m_chDat[p];
}

inline bool CNTermOnlineDigestion::_IsEnzyme(size_t p)
{
	return (string::npos != m_strCleave.find(m_chDat[p]) );
}

inline bool CNTermOnlineDigestion::_IsLeftTerm(size_t p)//because of reverse
{
	return DIVPEP == m_chDat[p] ;
}

inline bool CNTermOnlineDigestion::_IsRightTerm(size_t p)
{
	return (DIVPEP == m_chDat[p + 1]|| (('M' == m_chDat[p + 1] && DIVPEP == m_chDat[p + 2]) ));
}

inline bool CNTermOnlineDigestion::_IsLeft(size_t p) 
{
	return ( _IsLeftTerm(p) || _IsEnzyme(p));
}

inline bool CNTermOnlineDigestion::_IsRight(size_t p)
{
	return (_IsRightTerm(p) || _IsEnzyme(p));
}	

inline bool CNTermOnlineDigestion::_IsPep(size_t s, size_t e)
{
	size_t l = e - s + 1, Mass = (m_tMass[e] + (m_tMassMod - m_tMass[s - 1])) & m_tMassMod ;
	return (l >= m_pIndexItem.tMinPepLength && l <= m_pIndexItem.tMaxPepLength  && Mass >= m_pIndexItem.tMinPepMass && Mass <= m_pIndexItem.tMaxPepMass );
}


bool CNTermOnlineDigestion::InitPara(PINDEX_HEAD &pIndexHead,PINDEX_ITEM &pIndexItem, CMetaHandler &pMetaHandler)
{
	cout << "NTerm Init" << endl;
	
	m_pIndexHead = pIndexHead;
	m_pIndexItem = pIndexItem;
	m_pMetaHandler = pMetaHandler;
	if(m_pIndexItem.tMinPepMass < DEFAULT_MULTPLIER )
	{
		m_pIndexItem.tMinPepMass *= DEFAULT_MULTPLIER;
		m_pIndexItem.tMaxPepMass *= DEFAULT_MULTPLIER;
	}
	
	CEnzymeConf enzyme_conf(m_pIndexHead.stEnzymeList);
	if(!enzyme_conf.GetEnzyme(m_pIndexItem.vtEnzymeNames[0], enzyme))
	{
		proteomics_sdk::CErrInfo err_info("CNTermOnlineDigestion", "InitPara","Can't get enzyme!");
		throw runtime_error(err_info.Get().c_str());
	}
	m_strCleave = enzyme.GetCleaveString();
	m_strNotCleave = enzyme.GetNotCleave();

	
	CPepCalcFunc PepFunc;
	PepFunc.Init(m_pIndexHead.strAAList);
	m_tAAMass = new size_t[200];
	m_tAAMass[(int)DIVPEP] = 0;
	for(size_t t = 0; t < 26; ++t)
	{
		m_tAAMass[t + 'A'] = size_t( PepFunc.GetfAAMass((char)(t + 'A'), m_pIndexItem.bMono) * DEFAULT_MULTPLIER + 0.5);
	}
	
	m_uCurFileID = -1;
}

bool CNTermOnlineDigestion::Close()
{
	_DeleteAll();
}

void CNTermOnlineDigestion::_DeleteAll()
{
	try
	{
		if(NULL != m_chDat)
		{
			delete[] m_chDat;
			m_chDat = 0; 
		}
		if(NULL != LCP)
		{
			delete[] LCP; 
			LCP = 0;
		}
		if(NULL != SA)
		{
			delete[] SA; 
			SA = 0;
		}
		if(NULL != nClvy)
		{
			delete[] nClvy; 
			nClvy = 0; 
		}
		if(NULL != m_tMass)
		{
			delete[] m_tMass;
			m_tMass = 0;
		}
		if(NULL != beg)
		{
			delete[] beg; 
			beg = 0;
		}
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CNTermOnlineDigestion", "_DeleteAll");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CNTermOnlineDigestion", "_DeleteAll","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}



bool CNTermOnlineDigestion::GetNextAllPro(size_t tSufType)
{
	size_t clt = clock();
	cout << "GetNextAllPro is Begin..." << endl;
	
	if(++m_uCurFileID >= m_pMetaHandler.m_uFastaNum) return false;
	FILE *f = fopen(m_pMetaHandler.GenerFileName((uchar)m_uCurFileID,PRO_SQ_EXT).c_str(), "rb");
	if(NULL == f) return false;
	fseek(f, 0, SEEK_END);
	m_tLen = ftell(f);
	
	switch(tSufType)
	{
		case 0:
			m_chDat = new uchar[m_tLen + 10];
			break;
		case 3:
			m_chDat = new uchar[(m_tLen + 4) * 6 + 10];
			break;
		default:
			m_chDat = new uchar[m_tLen + 10];
			break;
	}
	
	fseek(f, 0, SEEK_SET);
	fread(m_chDat, sizeof(char), m_tLen, f);
	m_chDat[m_tLen] = 0;
	reverse(m_chDat, m_chDat+m_tLen);//for N Terminal digestion
	fclose(f);
	
	cout << "GetNextAllPro is end, the time is:" << clock() - clt << endl << endl;	
	
	return true;
}

bool CNTermOnlineDigestion::GetNextBeg()
{
	m_tNum = m_pMetaHandler.m_vmItems[m_uCurFileID].tEnd - m_pMetaHandler.m_vmItems[m_uCurFileID].tBeg; 
	size_t tn = 3 * ( m_tNum + 1);
	FILE *f = fopen(m_pMetaHandler.GenerFileName((uchar)m_uCurFileID,PRO_IDX_EXT).c_str(), "rb");
	if(NULL == f) return false;
	long *tmp = new long[tn];
	beg = new size_t[tn];
	fread(tmp,sizeof(long), tn , f);
	size_t t = 0;
	for(t =2; t < tn; t += 3)
	{
		beg[t/3] = tmp[t];
	}
	tn/=3;
	for(t = 0; t < tn - 2; ++t)//for N Terminal digestion,for the reverse of the proteins
	{
		beg[t] = (beg[tn - 1] - 1) - (beg[t + 1] - 2);//when reverse the proteins, the begin position of every protein is changed. 
	}
	beg[t] = (beg[tn - 1] - 1) - (beg[t + 1] - 1);//In the design of protein index, the first character of proteins is '?'
	reverse(beg,beg+ (tn-1));
	
//	for(t=0; t < tn; ++t)
//	{
//		cout << "beg:" << beg[t] << " ";
//		if(beg[t]) cout <<  m_chDat[beg[t] - 1] << " ";
//		cout << m_chDat[beg[t]] << " " << m_chDat[beg[t] + 1] << endl;
//	}
//	system("pause");
	
	delete[] tmp;tmp = 0;
	fclose(f);
	return true;
}

bool CNTermOnlineDigestion::GetNextLCP()
{
	FILE *fl = NULL;
	switch(m_pIndexItem.nCleaveWay)
	{
		case 0:
			fl = fopen(m_pMetaHandler.GenerFileName((uchar)m_uCurFileID,PRO_SQ_0_LCP).c_str(), "rb");
			if(NULL == fl) return false;
			LCP = new uchar[m_tLen + 2];
			fread(LCP, sizeof(uchar), m_tLen, fl);
			break;
		case 1:
			fl = fopen(m_pMetaHandler.GenerFileName((uchar)m_uCurFileID,PRO_SQ_2_LCP).c_str(), "rb");
			if(NULL == fl) return false;
			LCP = new uchar[m_tLen + 2];
			fread(LCP, sizeof(uchar), m_tLen, fl);
			break;
		case 2:
			fl = fopen(m_pMetaHandler.GenerFileName((uchar)m_uCurFileID,PRO_SQ_2_LCP).c_str(), "rb");
			if(NULL == fl) return false;
			LCP = new uchar[m_tLen + 2];
			fread(LCP, sizeof(uchar), m_tLen, fl);
			break;
		default:
			return false;
	}
	fclose(fl);
	return true;
}

bool CNTermOnlineDigestion::GetNextSA()
{
	FILE *fl = fopen(m_pMetaHandler.GenerFileName((uchar)m_uCurFileID,PRO_SQ_SA).c_str(), "rb");
	if(NULL == fl) return false;
	SA = new size_t[m_tLen + 2]; 
	if(fread((char*)SA, sizeof(size_t), m_tLen, fl) < m_tLen) return false;
	fclose(fl);
	return true;
}

void CNTermOnlineDigestion::_InitDigest()
{
	size_t clt = clock();
	cout << "_InitDigest is begin ..." << endl;
	
	m_tMass = new size_t[m_tLen + 3];
	nClvy = new size_t[m_tLen/4 + 3];
	ncp = 0;
	m_tMass[0] = 0,nClvy[ncp++] = 0;
	for(size_t t = 1; t < m_tLen - 1; ++t)//todo 
	{
		m_tMass[t] = (m_tMass[t - 1] + m_tAAMass[m_chDat[t]]) & m_tMassMod;
		if(_IsLeft(t) || _IsRight(t))
		{
			nClvy[ncp++] = t;
		}	
	}
	m_tMass[m_tLen] = m_tMass[m_tLen - 1] = m_tMass[m_tLen - 2];
	m_b = m_e = tbg = 0;
	tMaxClvy = m_pIndexItem.nMaxMissSite;
	cout << "_InitDigest is end , The time is:" << clock() - clt << " and ncp is:" << ncp << endl << endl;
}

void CNTermOnlineDigestion::_Print(size_t s, size_t e, uchar t,PEP_SQ &strPepSQ)
{
	char tmp = m_chDat[e +1];
	m_chDat[e + 1]	= 0;
	strPepSQ.strSQ = (char*)&m_chDat[s];
	reverse(strPepSQ.strSQ.begin(), strPepSQ.strSQ.end());
	m_chDat[e + 1] = tmp;	
	strPepSQ.dfMass = (((m_tMass[e] + (m_tMassMod + 1 - m_tMass[s - 1])) & m_tMassMod) * 1.0)/ DEFAULT_MULTPLIER;
	
	//the miss is wrong, however, in suffix array, we shouold don't care this
	strPepSQ.cMiss = t;
	
	//for protein infer
	strPepSQ.tInvertedFilesPos = s;
	
	strPepSQ.cEnd = 0;
	if(_IsLeftTerm(s - 1)) strPepSQ.cEnd |=2;//for Ntermal digestion, the protein has been reverse
	if(_IsRightTerm(e)) strPepSQ.cEnd |=1;
}

size_t CNTermOnlineDigestion::GetProID(size_t tp)
{
	size_t s = 0,e = m_tNum, mid = 0, rp = 0;
	while(s <= e)
	{
		mid =( s + e) >> 1;
		if(beg[mid] <= tp)
		{
			rp=mid,s=mid+1;
		}else e=mid-1;
	}	
	if(tp<=beg[rp+1]-2) return rp;
	return m_tNum;	
}

bool CNTermOnlineDigestion::GetPep2Pro(size_t tp, size_t tLen, vector<size_t> & vtProID)
{
	return false;
	vtProID.clear();
	do
	{	
		switch(m_pIndexItem.nCleaveWay)
		{
			case 0:
				if(_IsLeft(tp - 1))
				{
					vtProID.push_back(GetProID(tp));
				}
				break;
				
			case 1:
				if(_IsLeft(tp - 1) || _IsRight(tp + tLen -1)) vtProID.push_back(GetProID(tp));
				break;
				
			case 2:
				vtProID.push_back(GetProID(tp));
				break;
		}
		tp = SA[tp];
		if(tp >= m_tLen || LCP[tp] < tLen) break;
	}while(true);
	return true;
}

bool CNTermOnlineDigestion::GetSpecificPep(PEP_SQ &strPepSQ)
{
	bool fl = false;
	for(; m_b < ncp; ++m_b)
	{
		while(nClvy[m_b] + 1 >= beg[tbg + 1]) ++tbg;

		if(true == fl)
		{
			tMaxClvy = m_pIndexItem.nMaxMissSite;
			m_e = 0;
		}
		else fl = true;
		for(; m_e <= tMaxClvy && m_b + m_e + 1 < ncp; ++m_e)
		{
			if(nClvy[m_b + m_e + 1] >= beg[tbg + 1] - 1) break;// if not minus one,there will be rudundant of two, because the left of DIVPEP is also the  nClvy place				
			if(nClvy[m_b + m_e + 1] - nClvy[m_b] <= LCP[nClvy[m_b] + 1])
			{ 
				if(!_IsEnzyme(nClvy[m_b + m_e + 1])) ++tMaxClvy;
				continue;
			}
			if(_IsPep(nClvy[m_b] + 1, nClvy[m_b + m_e + 1]))
			{
				
				_Print(nClvy[m_b] + 1, nClvy[m_b + m_e + 1], (uchar)m_e, strPepSQ);
				
				++m_e;
				return true;
			}	
		}	
	}
						
	return false;
}

bool CNTermOnlineDigestion::_GetSemiSpecificPep1(PEP_SQ &strPepSQ)
{
	/*
	 *  In the suffix system, the first SQ of database begin with '?'
	 */ 
//	if(!m_b) ++m_b;
	
	bool fl = false;
	size_t t = 0;
	for(; m_b < ncp; ++m_b)
	{
		if(!_IsLeft(nClvy[m_b])) continue;
		while(nClvy[m_b] + 1 >= beg[tbg + 1]) ++tbg;
		tMaxClvy = m_pIndexItem.nMaxMissSite;
		for(t = 0; t <= tMaxClvy && m_b + t + 1 < ncp; ++t)
		{
			if(nClvy[m_b + t + 1] >= beg[tbg + 1] - 1) break;
			if(nClvy[m_b + t + 1] - nClvy[m_b] <= LCP[nClvy[m_b] + 1])
			{ 
				if(!_IsEnzyme(nClvy[m_b + t + 1])) ++tMaxClvy;
				continue;
			}			
		}
		--t;
		
		if(true == fl)m_e = 0;
		else fl = true;
		for( ; nClvy[m_b] + 1 + m_e  <  nClvy[m_b + t + 1]&& m_e < m_pIndexItem.tMaxPepLength; ++m_e)
		{
			if(nClvy[m_b] + 1 + m_e >= beg[tbg + 1] - 1) break;	
			if(_IsPep(nClvy[m_b] + 1, nClvy[m_b] + 1 +m_e) && !_IsRight(nClvy[m_b] + 1 + m_e))
			{
				_Print(nClvy[m_b] + 1, nClvy[m_b] + 1 + m_e , (uchar)0, strPepSQ);	
				
				++m_e;
				return true;
			}
		}
	}
	return false;
}

bool CNTermOnlineDigestion::_GetSemiSpecificPep2(PEP_SQ &strPepSQ)
{
	/*
	 *  In the suffix system, the first SQ of database begin with '?'
	 */ 
	if(!m_b) ++m_b;
	
	bool fl = false;
	for(; m_b < m_tLen; ++m_b)	
	{
		while(m_b >= beg[tbg + 1]) ++tbg;
		if(true == fl)
		{
			m_e = LCP[m_b];
			tMaxClvy = m_pIndexItem.nMaxMissSite;
			for(size_t t = 0; t < m_e; ++t)
			{
				if(_IsEnzyme(m_b + t)) --tMaxClvy;
			}
		}
		else fl = true;
		for( ; m_b + m_e  < m_tLen && m_e < m_pIndexItem.tMaxPepLength && tMaxClvy>=0; ++m_e)
		{
			if(m_b + m_e >= beg[tbg + 1] - 1) break;	
			if(_IsEnzyme(m_b + m_e)) --tMaxClvy;
			if(_IsPep(m_b, m_b +m_e) && _IsRight(m_b + m_e))
			{
				_Print(m_b, m_b + m_e , (uchar)0, strPepSQ);	
				
				++m_e;
				return true;
			}
		}
		
	}
	return false;
}

bool CNTermOnlineDigestion::GetSemiSpecificPep(PEP_SQ &strPepSQ)
{
	static bool flag = true;
	if(flag && _GetSemiSpecificPep2(strPepSQ)) return true;
	if(flag)
	{
		FILE *fl = fopen(m_pMetaHandler.GenerFileName((uchar)m_uCurFileID,PRO_SQ_0_LCP).c_str(), "rb");
		fread(LCP, sizeof(uchar), m_tLen, fl);
		fclose(fl);	
		m_b = m_e = tbg = 0;
		tMaxClvy = m_pIndexItem.nMaxMissSite;
		flag = false;
	}

	return _GetSemiSpecificPep1(strPepSQ);
	
}

bool CNTermOnlineDigestion::GetNonSpecificPep(PEP_SQ &strPepSQ)
{
	/*
	 *  In the suffix system, the first SQ of database begin with '?'
	 */ 
	if(!m_b) ++m_b;
	
	bool fl = false;
	for(; m_b < m_tLen; ++m_b)	
	{
		while(m_b >= beg[tbg + 1]) ++tbg;
		if(true == fl)m_e = LCP[m_b];
		else fl = true;
		for( ; m_b + m_e  < m_tLen && m_e < m_pIndexItem.tMaxPepLength; ++m_e)
		{
			if(m_b + m_e >= beg[tbg + 1] - 1) break;	
			if(_IsPep(m_b, m_b +m_e))
			{
				_Print(m_b, m_b + m_e , (uchar)0, strPepSQ);	
				
				++m_e;
				return true;
			}
		}
		
	}
	return false;
}

long long CNTermOnlineDigestion::StaticNonSpecificPep()
{
	long long fn = 0;
	size_t tbg = 0;
	for(m_b = 0; m_b < m_tLen; ++m_b)	
	{
		while(m_b >= beg[tbg + 1]) ++tbg;
		for(m_e = m_pIndexItem.tMinPepLength - 1 ; m_b + m_e  < m_tLen && m_e < m_pIndexItem.tMaxPepLength; ++m_e)
		{
			if(m_b + m_e >= beg[tbg + 1] - 1) break;	
			if(_IsPep(m_b, m_b +m_e))++fn;
		}
	}
	return fn;	
}

bool CNTermOnlineDigestion::GetNextBlock()
{
	_DeleteAll();
	if(false == GetNextAllPro(0)) return false;
	if(false == GetNextBeg())return false;
	if(false == GetNextLCP())return false;
	GetNextSA();
	_InitDigest();
	return true;
}

bool CNTermOnlineDigestion::GetNextPep(PEP_SQ &strPepSQ)
{
	switch(m_pIndexItem.nCleaveWay)
	{
		case 0:	
			return GetSpecificPep(strPepSQ);
		case 1:
			return GetSemiSpecificPep(strPepSQ);
		case 2:
			return GetNonSpecificPep(strPepSQ);

		default:
			break;
	}
	return false;
}	

bool CNTermOnlineDigestion::SuffixSort(size_t tSufType)
{
	size_t clt = 0;
	FILE *fl = NULL;
	switch(tSufType)
	{
		case 0:return false;
		case 1:break;
		case 2:break;
			
		case 3:
		{
			MSufSort * msufsort = new MSufSort();
			while(GetNextAllPro(tSufType))
			{				
				clt = clock();	
				m_chDat[--m_tLen] = 0;
				msufsort->Sort(m_chDat, m_tLen);
				
				SA = (size_t *)(m_chDat + ((m_tLen + 3) & ~3));
				m_chDat[m_tLen] = DIVPEP;//todo
				m_tLen = m_tLen + 1;					
				
				clt = clock();
				_GetLCP();
				
				clt = clock();
				fl = fopen(m_pMetaHandler.GenerFileName((uchar)m_uCurFileID,PRO_SQ_2_LCP).c_str(), "wb");
				fwrite((char*)LCP, sizeof(uchar), m_tLen, fl);
				fclose(fl);				
				
				clt = clock();
				_DealSpecial();
				fl = fopen(m_pMetaHandler.GenerFileName((uchar)m_uCurFileID,PRO_SQ_0_LCP).c_str(), "wb");
				fwrite((char*)LCP, sizeof(uchar), m_tLen, fl);					
				fclose(fl);
				
//				fl = fopen(m_pMetaHandler.GenerFileName((uchar)m_uCurFileID,PRO_SQ_SA).c_str(), "wb");
//				fwrite((char*)SA, sizeof(size_t), m_tLen, fl);
//				fclose(fl);
				
//				os.close();
				
				_DeleteAll();
			}
			delete msufsort;
			return true;
		}
		default:
			return false;
	}	
	return true;
}

void CNTermOnlineDigestion::_GetLCP()
{
	size_t clt = clock();
	cout << "_GetLCP is Begin...m_tLen is:" << m_tLen << endl;
	size_t *Rank = new size_t[m_tLen + 5];
	size_t t = 0, j = 0, tmp = 0;
	for(t = 0; t < m_tLen; ++t)
	{
		Rank[SA[t]] = t;
	}
	
	LCP = new unsigned char[m_tLen + 5];
	for(t = 0;t < m_tLen; ++t)
	{
		if(m_tLen - 1 == Rank[t])LCP[t]=0;
		else
		{
			if(t==0||LCP[t-1]<=1) j=0;
			else j=LCP[t-1]-1;
			for(tmp=SA[Rank[t]+1];m_chDat[t+j]==m_chDat[tmp+j] && j <= MAX_PEP_LENGTH;++j);//czhou, the max length of peptide
			LCP[t]=j;
		}
	}
	delete[] Rank;Rank = 0;

	cout << "_GetLCP is end, the time is:" << clock() - clt << endl << endl;	
}

void CNTermOnlineDigestion::_DealSpecial()
{
	size_t clt = clock();
	cout << "DealSpecial is Begin..." << endl;

	size_t *Rank = new size_t[m_tLen + 2];
	size_t t = 0;
	for(t = 0; t < m_tLen; ++t)
	{
		Rank[SA[t]] = t;
	}
	
	for(size_t i =0; i < m_tLen; ++i)
	{
		t = Rank[i];
		if(MIN_PEP_LENGTH - 1 > LCP[i] ) continue;
		if(_IsLeft(i - 1))//todo to judge i -1 < 0 is left
		{
			size_t p = 1;
			while(t + p < m_tLen)
			{						  
				int q = SA[t + p];							
				if(_IsLeft(q - 1)) break;
				if(LCP[q] < LCP[i]) LCP[i] = LCP[q];	
				if(LCP[q] == 0) break;
				++p;						  
			}	
			
		}
		if(_IsDiv(i + LCP[i]) && !_IsEnzyme(i + LCP[i] - 1)) --LCP[i];//蛋白质最后一个肽
	}	
	
	Rank[SA[0]] = m_tLen;//todo for not use AC algorithm
	for(t =m_tLen - 1; t >0 ; --t)
	{
		Rank[SA[t]] = SA[t - 1];
	}
	memcpy(SA, Rank, sizeof(size_t) * (m_tLen + 1));
	delete[] Rank;Rank = 0;
	cout << "DealSpecial is end, the time is:" << clock() - clt << endl << endl;		
}


}
