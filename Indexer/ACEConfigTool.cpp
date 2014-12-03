#include <stdexcept>
#include <string>
#include <iostream>
#include "../include/sdk.h"
#include "ACEConfigTool.h"

using namespace std;

CACEConfigTool::CACEConfigTool()
{
//	impExp = NULL;
}

CACEConfigTool::~CACEConfigTool()
{
	Close();
}

void CACEConfigTool::Open(string strFileName)
{
	FILE * pDBConf = NULL;
	pDBConf = fopen(strFileName.c_str(), "r");

	if( NULL == pDBConf)
	{
		CErrInfo info("CACEConfigTool", "Load", "fopen returns NULL");
		info.Append("strFileName=" + strFileName);
		throw runtime_error(info.Get().c_str());
	}

	fclose(pDBConf);

	m_pOption = new COptionTool("total", strFileName.c_str());

}

void CACEConfigTool::Close()
{
	if(NULL != m_pOption)
	{
		delete m_pOption;
		m_pOption = NULL;
	}
}

std::string CACEConfigTool::GetString(const std::string strSection, const std::string strKey, const std::string strDefault)
{
	return m_pOption->GetString(strSection.c_str(), strKey.c_str(), strDefault.c_str());
}

int CACEConfigTool::GetInteger(const std::string strSection, const std::string strKey, const int nDefault)
{
	return m_pOption->GetInteger(strSection.c_str(), strKey.c_str(), nDefault);
}

size_t CACEConfigTool::GetSizeT(const std::string strSection, const std::string strKey, const size_t tDefault)
{
	return m_pOption->GetSizeT(strSection.c_str(), strKey.c_str(),tDefault);
}

bool CACEConfigTool::GetBool(const std::string strSection, const std::string strKey, const bool bDefault)
{
	return m_pOption->GetBool(strSection.c_str(), strKey.c_str(), bDefault);
}
