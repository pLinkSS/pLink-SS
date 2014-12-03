#ifndef _PROTEININDEX_H_INCLUDED_
#define _PROTEININDEX_H_INCLUDED_

#include <string>
#include <vector>
#include "../include/predefine.h"

using namespace std;

	//	File extenstion of every  index file type
#define PRO_META_EXE	".pro.meta"
#define PRO_IDX_EXT		".pro.IDX"
#define PRO_AC_EXT		".pro.AC"
#define PRO_DE_EXT		".pro.DE"
#define PRO_SQ_EXT		".pro.SQ"
#define PRO_SQ_0_LCP	".lcp.0.SQ"
#define PRO_SQ_2_LCP	".lcp.2.SQ"
#define PRO_SQ_SA		".sa.SQ"
#define FILE_FASTA		".fasta"

#define DEFAULT_MULTPLIER    100000
#define DATABASE_FILE_BUF_SIZE  40960

// The maximal DESCRIPTION length of protein entry
#define MAX_DE_LENGTH_PRO 100
#define MAX_AC_LENGTH_PRO 100

#define MAX_PEP_LENGTH	150
#define MIN_PEP_LENGTH	4

#ifndef MAX_SQ_LENGTH_PRO
#define MAX_SQ_LENGTH_PRO 3000
#endif

const size_t m_tMassMod = 2147483647L;
struct BLOCK_ITEM
{
	size_t tBeg, tEnd;//the range is [)
};

// Structure for protein index file
// Head structure of .pro.IDX
typedef struct _PRO_IDX_HEAD
{
	long        lProIdxHeadOffset;		//a pointer to the first POSITION POINTER BLOCK
	size_t		tProteinNum;			// number of protein entries
	char		szOrgDBPath[PATH_MAX];		// original FASTA path

} PRO_IDX_HEAD;

// POSITION POINTER BLOCK structure of .pro.IDX
typedef struct _PRO_IDX_POS_POINTER		// one item of the POSITION BLOCK in .pro.IDX file.
{
	long		lPosAC;			// starting postion of AC
	long		lPosDE;			// starting postion of DE
	long		lPosSQ;			// starting postion of SQ

} PRO_IDX_POS_POINTER;
//
//typedef struct _PEP_SQ
//{
//	double			dfMass;
//	string			strSQ;
//	unsigned char	cMiss;
//	size_t			tInvertedFilesPos;
//	unsigned char	cDatNum;
//	unsigned char	cEnd;//0 means nothing, 1 means left, 2 means right, 3 means left&right
//} PEP_SQ;


typedef struct _PINDEX_HEAD
{
	string  stEnzymeList;
	string  strAAList;
	size_t  nIndexNum; // The num of index existing
	string  strOutpath;
} PINDEX_HEAD;

typedef struct _PINDEX_ITEM
{
	string  strFastaFile;
	string  strDBName;   	    // The name of database, type string?
	int nEnzymeNum;
	vector<std::string> vtEnzymeNames;// The name of enzyme, type string?
	int nMaxMissSite;
	size_t tMinPepLength;
	size_t tMaxPepLength;
	size_t tMinPepMass;
	size_t tMaxPepMass;
	size_t tMaxMemSize;
	size_t tMassRange;
	size_t tMaxFileSize;
	size_t tOutTxt;
	size_t tDatFlow;
	size_t tDatsFlow;
	int nCleaveWay;
	int	nAutoReverse;
	bool	bMono;
	bool	bAvrg;
	bool	bPro;
	bool	bPep;
	bool	bMass2Pep;
	bool	bPep2Pro;
	size_t	bProIndexMap;
	bool	bSaveDebugFile;
	bool	bConsole;
	
}PINDEX_ITEM;

#endif

