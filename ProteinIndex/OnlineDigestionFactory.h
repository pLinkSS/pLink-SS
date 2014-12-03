#ifndef ONLINEDIGESTIONFACTORY_H_
#define ONLINEDIGESTIONFACTORY_H_

#include "OnlineDigestion.h"

namespace ProteinIndex
{

class COnlineDigestionFactory
{
public:
	COnlineDigestionFactory();
	virtual ~COnlineDigestionFactory();
	COnlineDigestion * GetOnlineDigestion(bool);
	
};

}

#endif /*ONLINEDIGESTIONFACTORY_H_*/
