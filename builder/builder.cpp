#include "../ResultParse/SampleInfor.h"

#include <iostream>

using namespace std;

int main(int argc, char * argv[])
{
	try
	{
		//cout <<  "argc = " << argc << endl;
		CSampleInfor Sample;
		if (argc == 2)
			Sample.Run(argv[1], "");
		else if (argc > 2)
			Sample.Run(argv[1], argv[2]);
	} catch (runtime_error & e)
	{
		cout << e.what() << endl;
	}
	return 0;
}
