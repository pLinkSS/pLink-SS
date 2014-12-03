#include <iostream>
#include <string>
 
#include "../include/predefine.h"
#include "../include/sdk.h"
#include "../include/interface.h"
#include "../include/option.h"

#include "../XLinkPepResultFilter/XLinkPepResultFilter.h"
#include "FilterConf.h"

#include "XLinkResultFilterInterface.h"
#include "XLinkResultFilterFactory.h"

using namespace std;
using namespace proteomics_sdk;
int main(int argc, char * argv[])
{
	//cout << argc << endl;

	if(argc != 2)
		return 0;
	
	COptionTool * pOption(NULL);
	pOption = new COptionTool("filter", argv[1]);

	if(!pOption)
    {
    	cout << "can't open pOption file : " << argv[1] << " on [filter]" << endl;
    }
	
	int nFilterType = pOption->GetInteger("FilterMethod", 0);
	delete pOption;
	
	CXLinkResultFilterInterface * filter = NULL;
	CXLinkResultFilterFactory filterFact;
	filter = filterFact.GetFilter(nFilterType);
	if(filter == NULL)
	{
		cout << "can't get filter of type " << nFilterType << endl;
		return 1;
	}
	
	//cout << "init" << endl;
	filter->Init(argv[1]);
	//cout << "Run" << endl;
	filter->Run();
	//cout << "Close" << endl;
	filter->Close();
	//cout << "Over" << endl;
	delete filter;
	
}
