#ifndef ACECONFIGTOOL_H_
#define ACECONFIGTOOL_H_

#include "../include/option.h"
using namespace proteomics_sdk;

class CACEConfigTool
{
public:
	CACEConfigTool();
	void Open(std::string strFileName);
	std::string GetString(const std::string strSection, const std::string strKey, const std::string strDefault= "");
	int GetInteger(const std::string strSection, const std::string strKey, const int nDefault=0);
	size_t GetSizeT(const std::string strSection, const std::string strKey, const size_t tDefault = 0);
	bool GetBool(const std::string strSection, const std::string strKey, const bool bDefault = false);
	void Close();
	virtual ~CACEConfigTool();

private:
	COptionTool * m_pOption;
};

#endif /*ACECONFIGTOOL_H_*/
