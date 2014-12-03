#ifndef FILE_PTR_H_
#define FILE_PTR_H_
#include <cstdio>

using namespace std;

class File_ptr
{
	FILE * p;
public:
	File_ptr(const char* filename, const char* mode)
	{
		p = 0;
		p = fopen(filename, mode);
	}
	~File_ptr()
	{
		if(p)
			fclose(p);
	}
	operator FILE*()
	{
		return p;
	}
	bool isOpen()
	{
		return (p!=0);
	}
};

#endif /*FILE_PTR_H_*/
