#include <string>
#include <iostream>
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "SimpleAppender.h"

using namespace std;

namespace proteomics_sdk
{

CSimpleAppender::CSimpleAppender()
{
}

CSimpleAppender::~CSimpleAppender()
{
}
void CSimpleAppender::Init(const CCondition &condition)
{
}
void CSimpleAppender::Out(const string &strOut)
{
	cout<<strOut<<endl;
}
void CSimpleAppender::Close(void)
{
}
}
