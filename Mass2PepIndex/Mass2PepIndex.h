#ifndef MASS2PEPINDEX_H_
#define MASS2PEPINDEX_H_

#include "../ProteinIndex/ProteinIndex.h"

#ifndef PATH_MAX
#define PATH_MAX 260
#endif

#ifndef MAX_PATH
#define MAX_PATH 260
#endif

using namespace std;
 
namespace Mass2PepIndex
{
#define PEP_META   		".meta"
#define PEP_INFO_POS	".POS"
#define PEP_INFO_PID	".PID"
#define PEP_INFO_DAT	".PEP"
#define PEP_INFO_IDX	".IDX"
#define PEP_INFO_TMP	".Tmp"
#define PEP2PRO_DAT		".pep2pro.dat"

#define DEFAULT_MAX_MISSED_SITES 3
#define DEFAULT_MIN_IDX_PEP_LEN 4
#define DEFAULT_MAX_IDX_PEP_LEN 100
#define DEFAULT_MIN_IDX_PEP_MASS 300.0f
#define DEFAULT_MAX_IDX_PEP_MASS 10000.0f
#define DEFAULT_REDUNDANCY	0x80
#define PEP_NUM_EACH_PRO	50
#define DEFAULT_BIT_NUM		20
#define DEFAULT_PEP_RANGE	1000.0
#define DEFAULT_MEMORY_SIZE	536870912//512m
#define DEFAULT_MAX_TMP_PEP_FILE_NUM 1024
#define DEFAULT_MAX_SQ2PEP_TIME	100
#define CONSOLE
#define PEP2PROFILE

typedef struct _META_HEAD
{
	size_t	tOffset;	
	size_t	nCleaveWay;//czhou warning:can't change the place of nCleaveWay in this struct
	size_t	tMultiplier;
	size_t	tPepSQNum;
	size_t	tUniquePepSQNum;
	size_t 	tUniqueMassNum;
	size_t  tMinMass;
	size_t  tMaxMass;
	size_t	tMinLength;
	size_t	tMaxLength;
	char	cDBName[PATH_MAX];
	char	cEnzyme[PATH_MAX];
	unsigned char tIdxNum;
	unsigned char tMissCleaveSites;
	bool	bMono;
	bool	bPep2Pro;
}META_HEAD;

typedef struct _META_ITEM
{
	size_t tPepSQNum;
	size_t tUniquePepSQNum;
	size_t tUniqueMassNum;
	size_t tMinMass;
	size_t tMaxMass;
}META_ITEM;

typedef struct _PEP_INFO		
{
	size_t			tMass;
	size_t			tPos;		
	unsigned char	cLen;		
	unsigned char	cMiss;
} PEP_INFO;

typedef struct _PEP_INFO_EX		
{
	size_t			tMass;
	size_t			tPos;		
	size_t			tProID;
	size_t			tInvertedFilesPos;
	unsigned char	cLen;		
	unsigned char	cMiss;
	unsigned char	cDatNum;
} PEP_INFO_EX;

typedef struct _PEP_INFO_EX_TAG		
{
	size_t			tMass;
	size_t			tPos;		
	size_t			tProID;
	size_t			tInvertedFilesPos;
	double		lfTag;
	//size_t			tTag;
	unsigned char	cLen;		
	unsigned char	cMiss;
	unsigned char	cDatNum;
	unsigned char	cEnd;//0 means nothing, 1 means left, 2 means right, 3 means left&right
} PEP_INFO_EX_TAG;

typedef struct _PEP2PRO	
{
	size_t			tProNum;
	vector<size_t>	vtProID;
} PEP2PRO;



//typedef struct _PEP_SQ_EX		
//{
//	string			strSQ; 
//	double			dfMass;
//	size_t			tInvertedFilesPos;
//	unsigned char	cMiss;
//	unsigned char	cDatNum;
//} PEP_SQ_EX	;





typedef struct _HEAP_ITEM
{
	size_t tSer;
	PEP_INFO_EX_TAG sDatExTagBlock;
}HEAP_ITEM;

}
#endif /*MASS2PEPINDEX_H_*/

