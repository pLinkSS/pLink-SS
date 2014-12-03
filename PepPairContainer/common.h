#ifndef PEPPAIRCONTAINER_COMMON_H_
#define PEPPAIRCONTAINER_COMMON_H_

namespace proteomics_search{

#define WRITE_BUFFER_SIZE 10000
#define READ_BUFFER_SIZE 10000

#define MAX_SORT_UNIT 500000

#define LOW_LIMIT_FOR_QSORT 3

struct PEP_PAIR_ITEM
{
	bool bPair;
	int tPepId1;
	int tPepId2;
	double	lfMass;
	PEP_PAIR_ITEM()
	{
		bPair = true;
		tPepId1 = 0;
		tPepId2 = 0;
		lfMass = 0.0;
	}

	bool operator < (const PEP_PAIR_ITEM &item) const {
		return (lfMass < item.lfMass);
	}

	bool operator > (const PEP_PAIR_ITEM &item) const {
		return (lfMass > item.lfMass);
	}
};

}

#endif /*PEPPAIRCONTAINER_COMMON_H_*/
