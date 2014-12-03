#ifndef MASS2PEPINDEXCREATOR_H_
#define MASS2PEPINDEXCREATOR_H_

#include "Mass2PepIndex.h"
//using namespace ProteinIndex;
namespace Mass2PepIndex{

class CMass2PepIndexCreator
{
public:
	virtual ~CMass2PepIndexCreator(){};
	
	virtual void Init(const PINDEX_HEAD & , const PINDEX_ITEM &, bool) = 0;//init member
	virtual void Close()= 0;
	
	virtual void WriteIndex()= 0;
	
	
};

}
#endif /*MASS2PEPINDEXCREATOR_H_*/
