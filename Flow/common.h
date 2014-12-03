#ifndef COMMON_H_
#define COMMON_H_

#define RANDOM_TIME_SEED 10000

struct SPEC_WND_INFO
{
	int nIndex;
	double lfMassMin;
	double lfMassMax;
};

bool Score_Less(const CPeptideResult & pep1, const CPeptideResult & pep2);
bool XLINK_Score_Less(const CXLinkPepResult & pep1, const CXLinkPepResult & pep2);
bool XLINK_Score_More(const CXLinkPepResult & pep1, const CXLinkPepResult & pep2);
bool MIN_MASS_LESS(const SPEC_SORTER_INDEX_INFO & a, const SPEC_SORTER_INDEX_INFO & b);
bool XLINK_MIN_MASS_LESS(const SPEC_WND_INFO & a, const SPEC_WND_INFO & b);
bool XLINK_OPENFLOW_Score_Less(const CXLinkOpenPepResult & pep1, const CXLinkOpenPepResult & pep2);
bool XLINK_OPENFLOW_Score_More(const CXLinkOpenPepResult & pep1, const CXLinkOpenPepResult & pep2);
#endif /*COMMON_H_*/
