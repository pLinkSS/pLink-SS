#include "CTermOnlineDigestion.h"
#include "NTermOnlineDigestion.h"
#include "OnlineDigestionFactory.h"

namespace ProteinIndex
{

COnlineDigestionFactory::COnlineDigestionFactory()
{
}

COnlineDigestionFactory::~COnlineDigestionFactory()
{
}

COnlineDigestion * COnlineDigestionFactory::GetOnlineDigestion(bool bType)
{
	if(bType) return new CNTermOnlineDigestion();
	else return new CCTermOnlineDigestion();
//	return new CCTermOnlineDigestion();	
}

}
