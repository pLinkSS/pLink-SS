#include "MSufSort.h"

#include "../include/sdk.h"
#include "../include/option.h"
#include "OnlineDigestionFactory.h"

#include "ProteinHandler.h"
namespace ProteinIndex
{
class CLineSuffix
{
public:
inline bool leq(int a1, int a2,   int b1, int b2) { // lexic. order for pairs
  return(a1 < b1 || a1 == b1 && a2 <= b2); 
}                                                   // and triples
inline bool leq(int a1, int a2, int a3,   int b1, int b2, int b3) {
  return(a1 < b1 || a1 == b1 && leq(a2,a3, b2,b3)); 
}

// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
static void radixPass(int* a, int* b, int* r, int n, int K) 
{ // count occurrences
  int* c = new int[K + 1];                          // counter array
  for (int i = 0;  i <= K;  i++) c[i] = 0;         // reset counters
  for (int i = 0;  i < n;  i++) c[r[a[i]]]++;    // count occurences
  for (int i = 0, sum = 0;  i <= K;  i++) { // exclusive prefix sums
     int t = c[i];  c[i] = sum;  sum += t;
  }
  for (int i = 0;  i < n;  i++) b[c[r[a[i]]]++] = a[i];      // sort
  delete [] c;
}

// find the suffix array SA of s[0..n-1] in {1..K}^n
// require s[n]=s[n+1]=s[n+2]=0, n>=2
void suffixArray(int* s, int* SA, int n, int K) {
	static int level = 0;
//	cout << level++ << endl;
  int n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2; 
  int* s12  = new int[n02 + 3];  s12[n02]= s12[n02+1]= s12[n02+2]=0; 
  int* SA12 = new int[n02 + 3]; SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;
  int* s0   = new int[n0];
  int* SA0  = new int[n0];
  
//  cout << "new is success..." << endl;
  // generate positions of mod 1 and mod  2 suffixes
  // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
  for (int i=0, j=0;  i < n+(n0-n1);  i++) if (i%3 != 0) s12[j++] = i;

  // lsb radix sort the mod 1 and mod 2 triples
  radixPass(s12 , SA12, s+2, n02, K);
  radixPass(SA12, s12 , s+1, n02, K);  
  radixPass(s12 , SA12, s  , n02, K);

  // find lexicographic names of triples
  int name = 0, c0 = -1, c1 = -1, c2 = -1;
  for (int i = 0;  i < n02;  i++) {
    if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) { 
      name++;  c0 = s[SA12[i]];  c1 = s[SA12[i]+1];  c2 = s[SA12[i]+2];
    }
    if (SA12[i] % 3 == 1) { s12[SA12[i]/3]      = name; } // left half
    else                  { s12[SA12[i]/3 + n0] = name; } // right half
  }

  // recurse if names are not yet unique
  if (name < n02) {
//  	cout << "DC3 is again..." << endl;
    suffixArray(s12, SA12, n02, name);
    // store unique names in s12 using the suffix array 
    for (int i = 0;  i < n02;  i++) s12[SA12[i]] = i + 1;
  } else // generate the suffix array of s12 directly
    for (int i = 0;  i < n02;  i++) SA12[s12[i] - 1] = i; 

  // stably sort the mod 0 suffixes from SA12 by their first character
  for (int i=0, j=0;  i < n02;  i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
  radixPass(s0, SA0, s, n0, K);

  // merge sorted SA0 suffixes and sorted SA12 suffixes
  for (int p=0,  t=n0-n1,  k=0;  k < n;  k++) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
    int i = GetI(); // pos of current offset 12 suffix
    int j = SA0[p]; // pos of current offset 0  suffix
    if (SA12[t] < n0 ? 
        leq(s[i],       s12[SA12[t] + n0], s[j],       s12[j/3]) :
        leq(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]))
    { // suffix from SA12 is smaller
      SA[k] = i;  t++;
      if (t == n02) { // done --- only SA0 suffixes left
        for (k++;  p < n0;  p++, k++) SA[k] = SA0[p];
      }
    } else { 
      SA[k] = j;  p++; 
      if (p == n0)  { // done --- only SA12 suffixes left
        for (k++;  t < n02;  t++, k++) SA[k] = GetI(); 
      }
    }  
  } 
  delete [] s12; delete [] SA12; delete [] SA0; delete [] s0; 
}

};

//static size_t val[] = {0,0,7103711,11453494,10300919,11502694,12904259,14706841,5702146,13705891,11308406,0,12809496,11308406,13104049,11404293,0,9705276,12805858,15610111,8703203,10104768,0,9906841,18607931,11308406,16306333,12855059};

CProteinHandler::CProteinHandler()
{
	m_fpProIDX = NULL;
	m_fpProDAT = NULL;
	m_chDat = 0;LCP = 0;SA = 0;nClvy = 0;m_tMass =0; beg = 0;
	memset(m_IDX, 0, sizeof(m_IDX));
	memset(m_DAT, 0, sizeof(m_DAT));
}

CProteinHandler::~CProteinHandler()
{
	Close();
}

void CProteinHandler::Init(const string & strWorkDir, const string &strDBName,const string &szOrgDBPath)
{
	 m_strWorkDir = strWorkDir;
	 m_strDBName = strDBName;
	 m_szOrgDBPath = szOrgDBPath;
	 m_vmItems.clear();
	 m_uCurFileID = -1;
	 m_pDigester = NULL;
}

void CProteinHandler::Close()
{	
	try
	{
//		if(MAP_FAILED != &m_MapIDX)
//		{
//			ACE_OS::close(m_MapIDX.handle());
//			m_MapIDX.unmap();
//			m_MapIDX.close();
//		}
//		
//		if(MAP_FAILED != &m_MapDAT)
//		{
//			ACE_OS::close(m_MapDAT.handle());
//			m_MapDAT.unmap();
//			m_MapDAT.close();
//		}

		for(uchar uFileID = 0; uFileID < m_vmItems.size(); ++uFileID)		
		{
			if(m_DAT[uFileID])
			{
				fclose(m_DAT[uFileID]);
			}
			m_DAT[uFileID] = NULL;
		}
		
		if(NULL != m_fpProIDX)
		{
			fclose(m_fpProIDX);
			m_fpProIDX = NULL;
		}
//		if(NULL != m_fpProDAT)
//		{
//			fclose(m_fpProDAT);
//			m_fpProDAT = NULL;
//		}
		m_uCurFileID = (uchar)-1;
		_DeleteAll();
		
		if(NULL != m_pDigester)
		{
			m_pDigester->Close();
			delete m_pDigester;
			m_pDigester = NULL;
		}
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "Close");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "Close","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CProteinHandler::_DeleteAll()
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
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "_DeleteAll");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "_DeleteAll","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CProteinHandler::ObtOStream(fstream &fIdx, fstream &fAC, fstream &fDE, fstream &fSQ, uchar uFileID)
{
	try
	{
		if(uFileID)
		{
			fIdx.close();
			fAC.close();
			fDE.close();
			fSQ.close();
		}
		
		fIdx.open(GenerFileName(uFileID,PRO_IDX_EXT).c_str(),ios::out | ios::trunc | ios::binary);
		fAC.open(GenerFileName(uFileID,PRO_AC_EXT).c_str(),ios::out | ios::trunc | ios::binary);
		fDE.open(GenerFileName(uFileID,PRO_DE_EXT).c_str(),ios::out | ios::trunc | ios::binary);
		fSQ.open(GenerFileName(uFileID,PRO_SQ_EXT).c_str(),ios::out | ios::trunc | ios::binary);
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "ObtOStream");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "ObtOStream","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CProteinHandler::InitPara(	PINDEX_HEAD pIndexHead,PINDEX_ITEM pIndexItem )//todo
{	
	CEnzymeConf enzyme_conf(pIndexHead.stEnzymeList);
	if(!enzyme_conf.GetEnzyme(pIndexItem.vtEnzymeNames[0], enzyme))
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitPara","Can't get enzyme!");
		throw runtime_error(err_info.Get().c_str());
	}
		
	CMetaHandler metaHandler;
	metaHandler.m_strWorkDir = m_strWorkDir;
	metaHandler.m_strDBName = m_strDBName;
	metaHandler.m_szOrgDBPath = m_szOrgDBPath;
	metaHandler.ReadMeta();
		
	COnlineDigestionFactory onlineDigestionFactory;
	m_pDigester = onlineDigestionFactory.GetOnlineDigestion(enzyme.GetIsNTerm());
	m_pDigester->InitPara(pIndexHead,pIndexItem,metaHandler);
	
	
}

inline bool CProteinHandler::IsDiv(size_t p)
{
	return DIVPEP == m_chDat[p];
}

inline bool CProteinHandler::IsEnzyme(size_t p)
{
	return (string::npos != enzyme.GetCleaveString().find(m_chDat[p]) );
}

inline bool CProteinHandler::IsLeftTerm(size_t p)
{
	return DIVPEP == m_chDat[p] || ('M' == m_chDat[p] && DIVPEP == m_chDat[p -1]) ;
}

inline bool CProteinHandler::IsRightTerm(size_t p)
{
	return (DIVPEP == m_chDat[p + 1]);
}

inline bool CProteinHandler::IsLeft(size_t p)//todo 
{
	return ( IsLeftTerm(p) || IsEnzyme(p));
}

inline bool CProteinHandler::IsRight(size_t p)
{
	return (IsRightTerm(p) || IsEnzyme(p));
}	

inline bool CProteinHandler::IsPep(size_t s, size_t e)
{
	size_t l = e - s + 1, Mass = (m_tMass[e] + (m_tMassMod - m_tMass[s - 1])) & m_tMassMod ;
	return (l >= m_pIndexItem.tMinPepLength && l <= m_pIndexItem.tMaxPepLength  && Mass >= m_pIndexItem.tMinPepMass && Mass <= m_pIndexItem.tMaxPepMass );
}

//inline bool CProteinHandler::IsPep(size_t s, size_t e)
//{
//	size_t l = e - s + 1, Mass = (m_tMass[e] + (m_tMassMod - m_tMass[s - 1])) & m_tMassMod ;
//	cout << "m_tMass[e]:" << m_tMass[e] << " m_tMass[s - 1]:" << m_tMass[s - 1] << endl;
//	system("pause");
//	cout << "l:" << l << " m_pIndexItem.tMinPepLength:" <<  m_pIndexItem.tMinPepLength << " m_pIndexItem.tMaxPepLength" << endl;
//	cout << "Mass:" << Mass << " m_pIndexItem.tMinPepMass" << m_pIndexItem.tMinPepMass << " m_pIndexItem.tMaxPepMass:" << m_pIndexItem.tMaxPepMass << endl;
//	
//	return (l >= m_pIndexItem.tMinPepLength && l <= m_pIndexItem.tMaxPepLength  && Mass >= m_pIndexItem.tMinPepMass && Mass <= m_pIndexItem.tMaxPepMass );
//}



bool CProteinHandler::GetNextAllPro(size_t tSufType)
{
	size_t clt = clock();
	cout << "GetNextAllPro is Begin..." << endl;
	
	if(++m_uCurFileID >= m_uFastaNum) return false;
	FILE *f = fopen(GenerFileName((uchar)m_uCurFileID,PRO_SQ_EXT).c_str(), "rb");
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
	
	fclose(f);
	
	cout << "GetNextAllPro is end, the time is:" << clock() - clt << endl << endl;	
	
	return true;
}

bool CProteinHandler::GetNextBeg()
{
	m_tNum = m_vmItems[m_uCurFileID].tEnd - m_vmItems[m_uCurFileID].tBeg; 
	size_t tn = 3 * ( m_tNum + 1);
	FILE *f = fopen(GenerFileName((uchar)m_uCurFileID,PRO_IDX_EXT).c_str(), "rb");
	long *tmp = new long[tn];
	beg = new size_t[tn];
	fread(tmp,sizeof(long), tn , f);
	for(size_t t =2; t < tn; t += 3)
	{
		beg[t/3] = tmp[t];
	}
	
	delete[] tmp;tmp = 0;
	fclose(f);
	return true;
}

bool CProteinHandler::GetNextLCP()
{
	FILE *fl = NULL;
	switch(m_pIndexItem.nCleaveWay)
	{
		case 0:
			fl = fopen(GenerFileName((uchar)m_uCurFileID,PRO_SQ_0_LCP).c_str(), "rb");
			LCP = new uchar[m_tLen + 2];
			fread(LCP, sizeof(uchar), m_tLen, fl);
			break;
		case 1:
			fl = fopen(GenerFileName((uchar)m_uCurFileID,PRO_SQ_2_LCP).c_str(), "rb");
			LCP = new uchar[m_tLen + 2];
			fread(LCP, sizeof(uchar), m_tLen, fl);
			break;
		case 2:
			fl = fopen(GenerFileName((uchar)m_uCurFileID,PRO_SQ_2_LCP).c_str(), "rb");
			LCP = new uchar[m_tLen + 2];
			fread(LCP, sizeof(uchar), m_tLen, fl);
			break;
		default:
			return false;
	}
	fclose(fl);
	return true;
}

bool CProteinHandler::GetNextSA()
{
	FILE *fl = fopen(GenerFileName((uchar)m_uCurFileID,PRO_SQ_SA).c_str(), "rb");
	SA = new size_t[m_tLen + 2]; 
	if(fread((char*)SA, sizeof(size_t), m_tLen, fl) < m_tLen) return false;
	fclose(fl);
	return true;
}

void CProteinHandler::InitDigest()
{
	size_t clt = clock();
	cout << "InitDigest is begin ..." << endl;
	
	m_tMass = new size_t[m_tLen + 3];
	nClvy = new size_t[m_tLen/4 + 3];
	ncp = 0;
	m_tMass[0] = 0,nClvy[ncp++] = 0;
	for(size_t t = 1; t < m_tLen - 1; ++t)//todo 
	{
		m_tMass[t] = (m_tMass[t - 1] + m_tAAMass[m_chDat[t]]) & m_tMassMod;
		if(IsLeft(t) || IsRight(t))
		{
				nClvy[ncp++] = t;
		}	
	}
	m_tMass[m_tLen] = m_tMass[m_tLen - 1] = m_tMass[m_tLen - 2];
	m_b = m_e = tbg = 0;
	tMaxClvy = m_pIndexItem.nMaxMissSite;
	cout << "InitDigest is end , The time is:" << clock() - clt << " and ncp is:" << ncp << endl << endl;
}

void CProteinHandler::Print(size_t s, size_t e, uchar t,PEP_SQ &strPepSQ)
{
	char tmp = m_chDat[e +1];
	m_chDat[e + 1]	= 0;
	strPepSQ.strSQ = (char*)&m_chDat[s];
	m_chDat[e + 1] = tmp;	
	
//	strPepSQ.dfMass = (m_tMass[e] + (m_tMassMod - m_tMass[s - 1])) & m_tMassMod ;
	strPepSQ.dfMass = (((m_tMass[e] + (m_tMassMod + 1 - m_tMass[s - 1])) & m_tMassMod) * 1.0)/ DEFAULT_MULTPLIER;
	
	
	//the miss is wrong, however, in suffix array, we shouold don't care this
	strPepSQ.cMiss = t;
	
	//for protein infer
	strPepSQ.tInvertedFilesPos = s;
	
	strPepSQ.cEnd = 0;
	if(IsLeftTerm(s - 1)) strPepSQ.cEnd |=1;
	if(IsRightTerm(e)) strPepSQ.cEnd |=2;
}

size_t CProteinHandler::GetProID(size_t tp)
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


bool CProteinHandler::GetSpecificPep(PEP_SQ &strPepSQ)
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
				if(!IsRight(nClvy[m_b + m_e + 1])) ++tMaxClvy;
				continue;
			}
			if(IsPep(nClvy[m_b] + 1, nClvy[m_b + m_e + 1]))
			{
				
				Print(nClvy[m_b] + 1, nClvy[m_b + m_e + 1], (uchar)m_e, strPepSQ);
				
				++m_e;
				return true;
			}	
		}	
	}
						
	return false;
}



bool CProteinHandler::_GetSemiSpecificPep1(PEP_SQ &strPepSQ)
{
	bool fl = false;
	size_t t = 0;
	for(; m_b < ncp; ++m_b)
	{
		if(!IsLeft(nClvy[m_b])) continue;
		while(nClvy[m_b] + 1 >= beg[tbg + 1]) ++tbg;
		tMaxClvy = m_pIndexItem.nMaxMissSite;
		for(t = 0; t < tMaxClvy && m_b + t + 1 < ncp; ++t)
		{
			if(nClvy[m_b + t + 1] >= beg[tbg + 1] - 1) break;
			if(nClvy[m_b + t + 1] - nClvy[m_b] <= LCP[nClvy[m_b] + 1])
			{ 
				if(!IsRight(nClvy[m_b + t + 1])) ++tMaxClvy;
				continue;
			}			
		}
		--t;
		
		if(true == fl)m_e = 0;
		else fl = true;
		for( ; nClvy[m_b] + 1 + m_e  <  nClvy[m_b + t + 1]&& m_e < m_pIndexItem.tMaxPepLength; ++m_e)
		{
			if(nClvy[m_b] + 1 + m_e >= beg[tbg + 1] - 1) break;	
			if(IsPep(nClvy[m_b] + 1, nClvy[m_b] + 1 +m_e) && !IsRight(nClvy[m_b] + 1 + m_e))
			{
				Print(nClvy[m_b] + 1, nClvy[m_b] + 1 + m_e , (uchar)0, strPepSQ);	
				
				++m_e;
				return true;
			}
		}
	}
	return false;
}

bool CProteinHandler::_GetSemiSpecificPep2(PEP_SQ &strPepSQ)
{
	bool fl = false;
	for(; m_b < m_tLen; ++m_b)	
	{
		while(m_b >= beg[tbg + 1]) ++tbg;
		if(true == fl)m_e = LCP[m_b];
		else fl = true;
		for( ; m_b + m_e  < m_tLen && m_e < m_pIndexItem.tMaxPepLength; ++m_e)
		{
			if(m_b + m_e >= beg[tbg + 1] - 1) break;	
			if(IsPep(m_b, m_b +m_e) && IsRight(m_b + m_e))
			{
				Print(m_b, m_b + m_e , (uchar)0, strPepSQ);	
				
				++m_e;
				return true;
			}
		}
		
	}
	return false;
}

bool CProteinHandler::GetSemiSpecificPep(PEP_SQ &strPepSQ)
{
	static bool flag = true;
	if(flag && _GetSemiSpecificPep2(strPepSQ)) return true;
	if(flag)
	{
		FILE *fl = fopen(GenerFileName((uchar)m_uCurFileID,PRO_SQ_0_LCP).c_str(), "rb");
		fread(LCP, sizeof(uchar), m_tLen, fl);
		fclose(fl);
		flag = false;
	}
	return _GetSemiSpecificPep1(strPepSQ);
	
}

bool CProteinHandler::GetNonSpecificPep(PEP_SQ &strPepSQ)
{
	bool fl = false;
	for(; m_b < m_tLen; ++m_b)	
	{
		while(m_b >= beg[tbg + 1]) ++tbg;
		if(true == fl)m_e = LCP[m_b];
		else fl = true;
		for( ; m_b + m_e  < m_tLen && m_e < m_pIndexItem.tMaxPepLength; ++m_e)
		{
			if(m_b + m_e >= beg[tbg + 1] - 1) break;	
			if(IsPep(m_b, m_b +m_e))
			{
				Print(m_b, m_b + m_e , (uchar)0, strPepSQ);	
				
				++m_e;
				return true;
			}
		}
		
	}
	return false;
}

long long CProteinHandler::StaticNonSpecificPep()
{
	long long fn = 0;
	size_t tbg = 0;
	for(m_b = 0; m_b < m_tLen; ++m_b)	
	{
		while(m_b >= beg[tbg + 1]) ++tbg;
		for(m_e = m_pIndexItem.tMinPepLength - 1 ; m_b + m_e  < m_tLen && m_e < m_pIndexItem.tMaxPepLength; ++m_e)
		{
			if(m_b + m_e >= beg[tbg + 1] - 1) break;	
			if(IsPep(m_b, m_b +m_e))++fn;
		}
	}
	return fn;	
}

bool CProteinHandler::GetPep2Pro(size_t tp, size_t tLen, vector<size_t> & vtProID)
{
//	vtProID.clear();
//	do
//	{	
//		switch(m_pIndexItem.nCleaveWay)
//		{
//			case 0:
//				if(IsLeft(tp - 1))
//				{
//					vtProID.push_back(GetProID(tp));
//				}
//				break;
//				
//			case 1:
//				if(IsLeft(tp - 1) || IsRight(tp + tLen -1)) vtProID.push_back(GetProID(tp));
//				break;
//				
//			case 2:
//				vtProID.push_back(GetProID(tp));
//				break;
//		}
//		tp = SA[tp];
//		if(tp >= m_tLen || LCP[tp] < tLen) break;
//	}while(true);
//	return true;
	
	return m_pDigester->GetPep2Pro(tp, tLen, vtProID);
}


bool CProteinHandler::GetNextBlock()
{
//	_DeleteAll();
//	if(false == GetNextAllPro(0)) return false;
//	GetNextBeg();
//	GetNextLCP();
//	GetNextSA();
//	InitDigest();
//	
//	return true;
	
	return m_pDigester->GetNextBlock();
}

bool CProteinHandler::GetNextPep(PEP_SQ &strPepSQ)
{
//	switch(m_pIndexItem.nCleaveWay)
//	{
//		case 0:	
//			return GetSpecificPep(strPepSQ);
//		case 1:
//			return GetSemiSpecificPep(strPepSQ);
//		case 2:
//			return GetNonSpecificPep(strPepSQ);
//
//		default:
//			break;
//	}
//	return false;
	
	return m_pDigester->GetNextPep(strPepSQ);
}	

bool CProteinHandler::SuffixSort(size_t tSufType)
{
	return m_pDigester->SuffixSort(tSufType);
//	size_t clt = 0;
//	FILE *fl = NULL;
//	switch(tSufType)
//	{
//		case 0:return false;
//		case 1:break;
//		case 2:break;
//			
//		case 3:
//		{
//			MSufSort * msufsort = new MSufSort();
//			while(GetNextAllPro(tSufType))
//			{
//				ofstream os("time.txt",ios::app);
//				
//				clt = clock();	
//				m_chDat[--m_tLen] = 0;
//				msufsort->Sort(m_chDat, m_tLen);
//				os << "The time of Create Suffix Array with Read File is:" << clock() - clt << endl << endl;
//				
//				SA = (size_t *)(m_chDat + ((m_tLen + 3) & ~3));
//				m_chDat[m_tLen] = DIVPEP;//todo
//				m_tLen = m_tLen + 1;					
//				
//				clt = clock();
//				_GetLCP();
//				os << "The time of GetLCP is:" << clock() - clt << endl << endl;
//				
//				clt = clock();
//				fl = fopen(GenerFileName((uchar)m_uCurFileID,PRO_SQ_2_LCP).c_str(), "wb");
//				fwrite((char*)LCP, sizeof(uchar), m_tLen, fl);
//				fclose(fl);
//				os << "The time of WriteLCP is:" << clock() - clt << endl << endl;
//				
//				
//				clt = clock();
//				_DealSpecial();
//				os << "The time of DealSpecial is:" << clock() - clt << endl << endl;
//
//				fl = fopen(GenerFileName((uchar)m_uCurFileID,PRO_SQ_0_LCP).c_str(), "wb");
//				fwrite((char*)LCP, sizeof(uchar), m_tLen, fl);					
//				fclose(fl);
//				
//				fl = fopen(GenerFileName((uchar)m_uCurFileID,PRO_SQ_SA).c_str(), "wb");
//				fwrite((char*)SA, sizeof(size_t), m_tLen, fl);
//				fclose(fl);
//				
//				os.close();
//			}
//			delete msufsort;
//			return true;
//		}
//		default:
//			return false;
//	}	
//	return true;
}

void CProteinHandler::_GetLCP()
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
//	LCP = new uchar[m_tLen + 5];
//	for(t = 0;t < m_tLen; ++t)
//	{
//		LCP[t] = (uchar)(h[t] >100?100:h[t]);
//	}
//	delete[] h;
	
	cout << "_GetLCP is end, the time is:" << clock() - clt << endl << endl;	
}

void CProteinHandler::_DealSpecial()
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
		if(IsLeft(i - 1))//todo to judge i -1 < 0 is left
		{
			size_t p = 1;
			while(t + p < m_tLen)
			{						  
					int q = SA[t + p];							
					if(IsLeft(q - 1)) break;
					if(LCP[q] < LCP[i]) LCP[i] = LCP[q];	
					if(LCP[q] == 0) break;
					++p;						  
			}	
			
		}
		if(IsDiv(i + LCP[i]) && !IsEnzyme(i + LCP[i] - 1)) --LCP[i];//蛋白质最后一个肽
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


//added by chihao
//version after 2012.4.4

bool CProteinHandler::InitDiskDatFileByProID(size_t tProID, size_t tStartPos, string &strIDX, string &strDAT/**/)
{
	try
	{
		for(uchar  uFileID = 0; uFileID < m_vmItems.size(); ++uFileID)	
		{
			if(tProID < m_vmItems[uFileID].tEnd)
			{
				if(m_uCurFileID != uFileID)
				{	
					if(-1 != m_uCurFileID)
					{
						fclose(m_fpProIDX);
						m_fpProIDX = NULL;
					}

					m_fpProIDX = fopen(GenerFileName(uFileID,strIDX).c_str(), "rb");
					m_fpProDAT = m_DAT[uFileID];
					
					
					if(!m_fpProIDX || !m_fpProDAT)
					{
						proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitDiskDatFileByProID","The file don't exist!");
						throw runtime_error(err_info.Get().c_str());
					}
					
					m_uCurFileID = uFileID;
				}				
				long offset = GetPosInFile(tProID);
				if(PRO_DE_EXT == strDAT) offset += sizeof(long) * 1;
				else if(PRO_SQ_EXT == strDAT) offset += sizeof(long) * 2;
				
				fseek(m_fpProIDX, offset,SEEK_SET);
				fread((char *)&offset,sizeof(long),1,m_fpProIDX);

				fseek(m_fpProDAT,offset + tStartPos,SEEK_SET);				
				
				return true;
			}
		}
		return false;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitDiskDatFileByProID");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitDiskDatFileByProID","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
	return true;
}
/*
//version before 2012.4.4

bool CProteinHandler::InitDiskDatFileByProID(size_t tProID, size_t tStartPos, string &strIDX, string &strDAT)
{
	try
	{
		for(uchar  uFileID = 0; uFileID < m_vmItems.size(); ++uFileID)	
		{
			if(tProID < m_vmItems[uFileID].tEnd)
			{
				if(m_uCurFileID != uFileID)
				{	
					if(-1 != m_uCurFileID)
					{
						fclose(m_fpProIDX);
						fclose(m_fpProDAT);
					}

					m_fpProIDX = fopen(GenerFileName(uFileID,strIDX).c_str(), "rb");
					m_fpProDAT = fopen(GenerFileName(uFileID,strDAT).c_str(), "rb");
					
					cout << "pro id openXXXXX " << GenerFileName(uFileID,strIDX) << endl;
					
					if(!m_fpProIDX || !m_fpProDAT)
					{
						proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitDiskDatFileByProID","The file don't exist!");
						throw runtime_error(err_info.Get().c_str());
					}
					
					m_uCurFileID = uFileID;
				}				
				long offset = GetPosInFile(tProID);
				if(PRO_DE_EXT == strDAT) offset += sizeof(long) * 1;
				else if(PRO_SQ_EXT == strDAT) offset += sizeof(long) * 2;
				
				fseek(m_fpProIDX, offset,SEEK_SET);
				fread((char *)&offset,sizeof(long),1,m_fpProIDX);

				fseek(m_fpProDAT,offset + tStartPos,SEEK_SET);				
				
				return true;
			}
		}
		return false;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitDiskDatFileByProID");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitDiskDatFileByProID","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
	return true;
}
*/

/*
void CProteinHandler::InitDiskDatFileByFileID(uchar uFileID, size_t tStartPos, string &strIDX, string &strDAT)
{
	try
	{
		if(m_uCurFileID != uFileID)
		{
			if(-1 != m_uCurFileID)
			{
				fclose(m_fpProIDX);
				fclose(m_fpProDAT);
			}
			m_fpProIDX = fopen(GenerFileName(uFileID,strIDX).c_str(), "rb");
			fseek(m_fpProIDX, 0, SEEK_SET);

			m_fpProDAT = fopen(GenerFileName(uFileID,strDAT).c_str(), "rb");
			fseek(m_fpProDAT, 0, SEEK_SET);
		}
		fseek(m_fpProDAT,tStartPos,SEEK_SET);

		m_uCurFileID = uFileID;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitDiskDatFileByProID");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitDiskDatFileByProID","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}
*/
//bool CProteinHandler::InitDiskFileByPepPos(size_t tProID, size_t tStartPos, string &strIDX, string &strDAT)
//{
//	try
//	{
//		for(uchar  uFileID = 0; uFileID < m_vmItems.size(); ++uFileID)
//		{
//			if( tProID < m_vmItems[uFileID].tEnd)
//			{
//				if(m_uCurFileID != uFileID)
//				{
//					if(-1 != m_uCurFileID)
//					{
//						fclose(m_fpProIDX);
//						m_fpProIDX = NULL;
//					}
//
//
//					m_fpProIDX = fopen(GenerFileName(uFileID,strIDX).c_str(), "rb");
//					m_fpProDAT = m_DAT[uFileID];
//
//					if(!m_fpProIDX || !m_fpProDAT)
//					{
//						proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitDiskDatFileByProID","The file don't exist!");
//						throw runtime_error(err_info.Get().c_str());
//					}
//
//					m_uCurFileID = uFileID;
//				}
//
//				fseek(m_fpProDAT, tStartPos,SEEK_SET);
//
//				return true;
//			}
//		}
//		return false;
//	}
//	catch(runtime_error &e)
//	{
//		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitMMapDatFileByProID");
//		throw runtime_error(err_info.Get(e).c_str());
//	}
//	catch(...)
//	{
//		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitMMapDatFileByProID","Unknown Error!");
//		throw runtime_error(err_info.Get().c_str());
//	}
//	return false;
//
//}


//bool CProteinHandler::InitMMapFileByPepPos(size_t tProID, size_t tStartPos, string &strIDX, string &strDAT)
//{
//	try
//	{
//		for(uchar  uFileID = 0; uFileID < m_vmItems.size(); ++uFileID)		
//		{
//			if( tProID < m_vmItems[uFileID].tEnd)
//			{
//				if(m_uCurFileID != uFileID)
//				{	
//					if(-1 != m_uCurFileID)
//					{
//						ACE_OS::close(m_MapIDX.handle());
//						m_MapIDX.unmap();
//						m_MapIDX.close();
//						
//						ACE_OS::close(m_MapDAT.handle());
//						m_MapDAT.unmap();
//						m_MapDAT.close();						
//					}
//
//
//					ACE_HANDLE handle = ACE_OS::open (GenerFileName(uFileID,strIDX).c_str(), O_RDONLY);
//					m_MapIDX.map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//					
//
//					handle = ACE_OS::open (GenerFileName(uFileID,strDAT).c_str(), O_RDONLY);
//					m_MapDAT.map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//					
//					
//					if(handle == ACE_INVALID_HANDLE || MAP_FAILED == &m_MapIDX || MAP_FAILED == &m_MapDAT )
//					{
//						proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitMMapDatFileByProID","The file don't exist!");
//						throw runtime_error(err_info.Get().c_str());
//					}
//					m_uCurFileID = uFileID;
//				}		
//				
//				ACE_OS::lseek(m_MapDAT.handle(), tStartPos,SEEK_SET);
//
//				return true;
//			}
//		}
//		return false;
//	}
//	catch(runtime_error &e)
//	{
//		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitMMapDatFileByProID");
//		throw runtime_error(err_info.Get(e).c_str());
//	}
//	catch(...)
//	{
//		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitMMapDatFileByProID","Unknown Error!");
//		throw runtime_error(err_info.Get().c_str());
//	}
//	return false;
//	
//}
//
//bool CProteinHandler::InitMMapDatFileByProID(size_t tProID, size_t tStartPos, string &strIDX, string &strDAT)
//{
//	try
//	{
//
//		for(uchar  uFileID = 0; uFileID < m_vmItems.size(); ++uFileID)		
//		{
//			if( tProID < m_vmItems[uFileID].tEnd)
//			{
//				if(m_uCurFileID != uFileID)
//				{	
//					if(-1 != m_uCurFileID)
//					{
//						ACE_OS::close(m_MapIDX.handle());
//						m_MapIDX.unmap();
//						m_MapIDX.close();
//						
//						ACE_OS::close(m_MapDAT.handle());
//						m_MapDAT.unmap();
//						m_MapDAT.close();						
//					}
//
//					ACE_HANDLE handle = ACE_OS::open (GenerFileName(uFileID,strIDX).c_str(), O_RDONLY);
//					m_MapIDX.map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//					
//					handle = ACE_OS::open (GenerFileName(uFileID,strDAT).c_str(), O_RDONLY);
//					m_MapDAT.map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//					
//					
//					if(handle == ACE_INVALID_HANDLE || MAP_FAILED == &m_MapIDX || MAP_FAILED == &m_MapDAT )
//					{
//						proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitMMapDatFileByProID","The file don't exist!");
//						throw runtime_error(err_info.Get().c_str());
//					}
//					m_uCurFileID = uFileID;
//				}		
//				
//
//				long offset = GetPosInFile(tProID);
//				if(PRO_DE_EXT == strDAT) offset += sizeof(long) * 1;
//				else if(PRO_SQ_EXT == strDAT) offset += sizeof(long) * 2;
//				
//				ACE_OS::lseek(m_MapIDX.handle(), offset, SEEK_SET);
//				ACE_OS::read(m_MapIDX.handle(), (char *)&offset, sizeof(long));
//				
//				ACE_OS::lseek(m_MapDAT.handle(), offset + tStartPos,SEEK_SET);
//				return true;
//			}
//		}
//		return false;
//	}
//	catch(runtime_error &e)
//	{
//		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitMMapDatFileByProID");
//		throw runtime_error(err_info.Get(e).c_str());
//	}
//	catch(...)
//	{
//		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitMMapDatFileByProID","Unknown Error!");
//		throw runtime_error(err_info.Get().c_str());
//	}
//	return false;
//}
//
//void CProteinHandler::InitMMapDatFileByFileID(uchar uFileID, size_t tStartPos, string &strIDX, string &strDAT)
//{
//	try
//	{
//		if(m_uCurFileID != uFileID)
//		{	
//			if(-1 != m_uCurFileID)
//			{
//				m_MapIDX.unmap();
//				m_MapIDX.close();
//				
//				m_MapDAT.unmap();
//				m_MapDAT.close();
//			}
//
//			ACE_HANDLE handle = ACE_OS::open (GenerFileName(uFileID,strIDX).c_str(), O_RDONLY);
//			m_MapIDX.map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//
//			handle = ACE_OS::open (GenerFileName(uFileID,strDAT).c_str(), O_RDONLY);
//			m_MapDAT.map(handle, (size_t)(-1),PROT_READ, ACE_MAP_SHARED);
//		}				
//	
//		ACE_OS::lseek(m_MapDAT.handle(), tStartPos, SEEK_SET);
//		
//		m_uCurFileID = uFileID;
//	}
//	catch(runtime_error &e)
//	{
//		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitMMapDatFileByFileID");
//		throw runtime_error(err_info.Get(e).c_str());
//	}
//	catch(...)
//	{
//		proteomics_sdk::CErrInfo err_info("CProteinHandler", "InitMMapDatFileByFileID","Unknown Error!");
//		throw runtime_error(err_info.Get().c_str());
//	}
//}

size_t CProteinHandler::GetPosInFile(size_t tProID)
{
	try
	{
		size_t t = m_uCurFileID;
		if(!(tProID >= m_vmItems[t].tBeg && tProID < m_vmItems[t].tEnd))
		{
			for(t = 0; t < m_vmItems.size(); ++t)//***		
			{
				if(/*tProID >= m_vmItems[uFileID].tBeg &&*/ tProID < m_vmItems[t].tEnd) break;
			}
		}
		return (tProID - m_vmItems[t].tBeg) * sizeof(long) * 3;;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "GetPosInFile");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "GetPosInFile","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

size_t CProteinHandler::GetEndPosInFile(size_t tProID)
{
	try
	{
		size_t t = m_uCurFileID;
		if(!(tProID >= m_vmItems[t].tBeg && tProID <= m_vmItems[t].tEnd))
		{
			for(t = 0; t < m_vmItems.size(); ++t)//***		
			{
				if(/*tProID >= m_vmItems[uFileID].tBeg &&*/ tProID <= m_vmItems[t].tEnd) break;
			}
		}
		return (tProID - m_vmItems[t].tBeg) * sizeof(long) * 3;;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "GetPosInFile");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "GetPosInFile","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CProteinHandler::InsertItem(BLOCK_ITEM vmItem)
{
	m_vmItems.push_back(vmItem);
}

void CProteinHandler::SetProInfor(size_t &tProNum, uchar & uFastaNum)
{
	m_tProNum = tProNum;
	m_uFastaNum = uFastaNum;
}

size_t CProteinHandler::GetProNum() const
{
	return m_tProNum;
}

string CProteinHandler::GenerFileName(uchar uFileID,string strApp)
{
	try
	{
		string str = _GenerFileName();
		char sz[10];
		sprintf(sz,".%u", uFileID);
		return str + string(sz) + strApp;
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "GenerFileName");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "GenerFileName","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

string CProteinHandler::_GenerFileName()
{
	return m_strWorkDir + m_strDBName;
}



void CProteinHandler::WriteMeta()
{
	_WriteHead();
	_WriteItem();
}

void CProteinHandler::ReadMeta()
{
	_ReadHead();
	_ReadItem();
}

void CProteinHandler::_WriteHead()
{
	try
	{
		string strMeta = _GenerFileName() + PRO_META_EXE;
		string strTxt = strMeta + ".txt";
		
		ofstream fMeta(strMeta.c_str(), ios::out | ios::trunc | ios::binary );
		ofstream fTxt(strTxt.c_str(), ios::out | ios::trunc );
		
		fMeta.write((char *)&m_HeadSize, sizeof(size_t));
		fMeta.write((char *)&m_uFastaNum, sizeof(uchar));
		fMeta.write((char *)&m_tProNum, sizeof(size_t));		
		
		fMeta.write(m_szOrgDBPath.c_str() , m_szOrgDBPath.size());
		char abc[PATH_MAX + 1];
		memset(abc,0, PATH_MAX + 1);
		fMeta.write(abc , PATH_MAX - m_szOrgDBPath.size());

		m_HeadSize = fMeta.tellp();
		fMeta.seekp(0, ios::beg);
		fMeta.write((char *)&m_HeadSize, sizeof(size_t));
		fMeta.close();
		
		fTxt << "[Head]" << endl;
		fTxt << "FastaNum = " << (int)(m_uFastaNum )<< endl;
		fTxt << "ProNum = " << m_tProNum << endl;
		fTxt << "OrgDBPath = " << m_szOrgDBPath << endl;
		fTxt.close();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "_WriteHead");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "_WriteHead","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CProteinHandler::_ReadHead()
{
	try
	{
		string strMeta = _GenerFileName() + PRO_META_EXE;
	
		ifstream fMeta(strMeta.c_str(), ios::binary );
		fMeta.seekg(0, ios::beg);
		fMeta.read((char *)&m_HeadSize, sizeof(size_t));
		fMeta.read((char *)&m_uFastaNum, sizeof(uchar));
		fMeta.read((char *)&m_tProNum, sizeof(size_t));
		char *szOrgDBPath = new char[PATH_MAX + 1];
		fMeta.read(szOrgDBPath, sizeof(char) * PATH_MAX);//****
		m_szOrgDBPath = szOrgDBPath;
		fMeta.close();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "_ReadHead");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "_ReadHead","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CProteinHandler::_WriteItem()
{
	try
	{
		string strMeta = _GenerFileName() + PRO_META_EXE;
		string strTxt = strMeta + ".txt";
		
		ofstream fMeta(strMeta.c_str(), ios::out  | ios::app |ios::binary );
		ofstream fTxt(strTxt.c_str(), ios::out | ios::app );
		
		fMeta.seekp(0, ios::end);
		fTxt.seekp(0, ios::end);
		
		size_t tDatSize = m_vmItems.size();
		for(size_t t = 0; t < tDatSize; ++t)
		{
			fMeta.write((char*)&m_vmItems[t].tBeg, sizeof(size_t));
			fMeta.write((char*)&m_vmItems[t].tEnd, sizeof(size_t));
			
			char sz[10];
			sprintf(sz,"%d",t);
			fTxt << "[dat" << string(sz) << "]" << endl;
			
			fTxt << "First Pro = " << m_vmItems[t].tBeg << endl;
			fTxt << "Last Pro = " << m_vmItems[t].tEnd << endl;
		}
		
		fMeta.close();
		fTxt.close();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "_WriteItem");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "_WriteItem","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

void CProteinHandler::_ReadItem()
{
	try
	{
		string strMeta = _GenerFileName() + PRO_META_EXE;
		
		ifstream fMeta(strMeta.c_str(), ios::binary );	
		fMeta.seekg(m_HeadSize, ios::beg);
		m_vmItems.clear();
		BLOCK_ITEM mItem;
		for(size_t t = 0; t < m_uFastaNum; ++t)
		{		
			fMeta.read((char*)&mItem.tBeg, sizeof(size_t));
			fMeta.read((char*)&mItem.tEnd, sizeof(size_t));
			m_vmItems.push_back(mItem);
		}
		fMeta.close();
	}
	catch(runtime_error &e)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "_ReadItem");
		throw runtime_error(err_info.Get(e).c_str());
	}
	catch(...)
	{
		proteomics_sdk::CErrInfo err_info("CProteinHandler", "_ReadItem","Unknown Error!");
		throw runtime_error(err_info.Get().c_str());
	}
}

}

