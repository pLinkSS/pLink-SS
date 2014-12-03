
#include "Inferor.h"
#include "XLinkInferor.h"
#include "InferorFactory.h"


CInferorFactory::CInferorFactory()
{
}

CInferorFactory::~CInferorFactory()
{
}

CInferor * CInferorFactory::GetInferor(INFEROR_TYPE eType)
{
	if(eType == XLINK_INFEROR)
	{
		return new CXLinkInferor();
	}
	else
	{
		return new CXLinkInferor();
	}
}
