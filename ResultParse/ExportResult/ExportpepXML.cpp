#include "ExportpepXML.h"
#include "../BasicFunction/ReadWrite.h"

using namespace bio_analysis;
//extern ofstream flog;
extern ostringstream osspBuildLog;
CExportpepXML::CExportpepXML()
{

}

CExportpepXML::~CExportpepXML()
{
}

void CExportpepXML::Export(const IntegratedInfo & HaveToInterData, const CConf & conf)
{
	//cout << "ExportpepXML" << endl;
	osspBuildLog << "ExportpepXML" << endl;
	//pBuildLog("ExportpepXML");
}

