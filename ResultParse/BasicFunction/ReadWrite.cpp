#include "ReadWrite.h"

using namespace bio_analysis;

char SearchEngineName[10][20] = { "pFind", "Mascot", "SEQUEST", "SQT" };

extern ostringstream osspBuildLog;

void pBuildLog(const string & strOut) {
	ifstream fin("pBuild.log");

	ostringstream strTemp;

	if (fin.good()) {
		string strVal;

		while (GetLine(fin, strVal)) {
			strTemp << strVal << endl;
		}
	}

	fin.close();

	ofstream fout("pBuild.log");
	fout << strOut;
	fout << strTemp.str();
	fout.close();
}

bool ISFolder(const string & strPath) {
	struct stat info;
	stat(strPath.c_str(), &info);

	if (S_ISDIR(info.st_mode)) {
		//printf("This is a directory");
		return true;
	}

	return false;
}

void RemoveFile(const string & strFile) {
	if (access(strFile.c_str(), 0) != 0) {
		//printf("error");
		return;
	}

	if (remove(strFile.c_str()) != 0)
		perror("remove");
}

void RemoveFolder(const string & strFolder) {
	DIR * dir = opendir(strFolder.c_str());

	if (NULL == dir) {
		CErrInfo info("ReadWrite.cpp", "RemoveFolder",
				"Cannot open the folder: " + strFolder);
		throw runtime_error(info.Get());
		printf("error");
		return;
	}
	struct dirent * de;
	for (; NULL != (de = readdir(dir));) {
		if (de->d_name[0] == '.')
			continue;

		string FileName;
		if (strFolder[strFolder.size() - 1] == '/')
			FileName = strFolder + de->d_name;
		else
			FileName = strFolder + "/" + de->d_name;

		if (ISFolder(FileName))
			RemoveFolder(FileName);
		else
			RemoveFile(FileName);
	}

	rmdir(strFolder.c_str());
}

void RemoveFolderOrFile(string strPath) {
	if (strPath[strPath.size() - 1] == '/') {
		strPath = strPath.substr(0, strPath.size() - 1);
	}

	if (ISFolder(strPath))
		RemoveFolder(strPath);
	else
		RemoveFile(strPath);
}

void PrintProteinDatabase(FILE * fout, const string & strPro) {
	size_t t = 0;
	for (t = 0; t < strPro.size(); t++) {
		if (t != 0 && t % 100 == 0) {
			fprintf(fout, "\n");
		}
		fprintf(fout, "%c", strPro[t]);
	}
	//if (t % 100 != 0)
	fprintf(fout, "\n");
}

void CheckIndexFile(const string & strIndexFile) {
	ifstream fin(strIndexFile.c_str());
	if (!fin.good()) {
		CErrInfo info("ReadWrite.cpp", "CheckIndexFile",
				"Cannot open the file: " + strIndexFile);
		throw runtime_error(info.Get());
	}

	string strVal;
	GetLine(fin, strVal);

	if ("DATASIZE" != strVal.substr(0, 8)) {
		CErrInfo info("ReadWrite.cpp", "CheckIndexFile",
				"This file is not IndexFile" + strIndexFile);
		throw runtime_error(info.Get());
	}

	fin.close();
}

bool CheckFileExist(const string & strFilePath) {
	ifstream fin(strFilePath.c_str());
	if (!fin.good()) {
		fin.close();
		return false;
	}

	fin.close();
	return true;
}

string GetTimePre(const CConf & conf) {
#ifdef DEBUG
	string strPre;

	char pre[PATH_MAX] =
	{	'\0'};

	struct tm * timeinfo;
	timeinfo = localtime(&conf.time_now);

	sprintf(pre, "%04d-%02d-%02d-%02d-%02d-%02d_", 1900 + timeinfo->tm_year, 1 + timeinfo->tm_mon,
			timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

	strPre = pre;
	return strPre;
#endif
	return "";
}
bool CheckSimpleSame(const CConf & conf, const string & strParserSimple) {
	char chrInput[PATH_MAX] = { 0 };

	string strInput = "";

	int Total_Samples = conf.m_vEngineInputInfo.size();

	sprintf(chrInput, "Total_Samples=%d", Total_Samples);
	strInput += chrInput;

	for (int i = 0; i < Total_Samples; i++) {
		sprintf(chrInput, "[Sample%d]", i + 1);
		strInput += chrInput;

		sprintf(chrInput, "EngineType=%s",
				SearchEngineName[conf.m_vEngineInputInfo[i].m_Type]);
		strInput += chrInput;

		size_t tsize = conf.m_vEngineInputInfo[i].m_vResultPath.size();
		sprintf(chrInput, "SubItems="PRI_SIZE_T, tsize);
		strInput += chrInput;

		for (size_t j = 0; j < tsize; j++) {
			sprintf(chrInput, "SubItem"PRI_SIZE_T"=%s,%d", j + 1,
					conf.m_vEngineInputInfo[i].m_vResultPath[j].c_str(),
					conf.m_vEngineInputInfo[i].m_Valid[j]);
			//cout << chrInput << endl;
			strInput += chrInput;
		}
	}

	sprintf(chrInput, "[END]");
	strInput += chrInput;

	ifstream fin(strParserSimple.c_str());
	if (!fin.good()) {
		return false;
		//		CErrInfo info("ReadWrite.cpp", "CheckSimpleSame", "Cannot open the file: "
		//				+ strParserSimple);

		//throw runtime_error(info.Get());
	}

	string strVal;
	GetLine(fin, strVal);

	fin.close();

	if (strVal == strInput)
		return true;

	return false;
}

bool CheckFolderExist(const string & strForderPath) {
	DIR * dir = opendir(strForderPath.c_str());
	if (NULL == dir) {
		closedir(dir);
		return false;
	}
	closedir(dir);

	//return false;//todo 永远返回没有索引
	return true;
}

void SetDefault(CConf & conf) {
	conf.m_SampleID = 1;
	conf.m_vEngineInputInfo.clear();
	conf.m_vProDBPath.clear();
	conf.m_vDecoyTags.clear();
	conf.m_vDecoyTags.push_back("REVERSE");

	conf.m_ReParse = 0;

	conf.m_OutPutForder = "C:\\";
	conf.m_OutPutFile = "Task1";
	conf.m_OutPutName = "C:\\Task1";
	conf.ExportFormat.clear();
	conf.ExportFormat.push_back(ST_HTMLFormat);
	conf.ExportFileType.clear();

	conf.strExportPath = "C:\\Task1";

	conf.m_Filter.m_UseFilter = 1;
	conf.m_Filter.m_lfPepMassLowerBound = 800.0;
	conf.m_Filter.m_lfPepMassUpperBound = 100000000.0;
	conf.m_Filter.m_tLengthLowerBound = 1;
	conf.m_Filter.m_tLengthUpperBound = 1000000;
	conf.m_Filter.m_Sepatrate = false;
	conf.m_Filter.m_lfPepTolBase.clear();
	conf.m_Filter.m_lfPepTolBase.push_back(0.0);
	conf.m_Filter.m_lfPepTolLowerBound.clear();
	conf.m_Filter.m_lfPepTolLowerBound.push_back(-10.0);
	conf.m_Filter.m_lfPepTolUpperBound.clear();
	conf.m_Filter.m_lfPepTolUpperBound.push_back(10);
	conf.m_Filter.m_strPepTolType = "ppm";

	for (int i = 1; i <= 16; i++) {
		conf.m_Filter.m_vChargeState.push_back(i);
	}

	conf.m_Filter.m_tRankLimit = 1;
	conf.m_Filter.m_bRedundant = 1;
	conf.m_Filter.m_bFixedDeltCn = true;
	conf.m_Filter.m_lfFixedDeltCn = 0.1;
	conf.m_Filter.m_bUseFDR = true;
	conf.m_Filter.m_lfFDR = 0.01;
	conf.m_Filter.m_FDRFormula = 0;
	conf.m_Filter.m_ScoreMax = 100000000;
	conf.m_Filter.m_ScoreMin = 0.000000;

	conf.m_Filter.m_nDistinctPepLimit = 1;
	conf.m_Filter.m_nDistinctSpecLimit = 1;

	conf.m_Filter.m_strModSites.clear();
	conf.m_Filter.m_strCTerminal.clear();
	conf.m_Filter.m_strNTerminal.clear();

	conf.m_ReParse = 0;
	conf.m_strPNovoParam.clear();
	//conf.m_bReserveDecoy = 1;
	conf.m_bReserveDecoy = 0;
	conf.m_bCrossLink = 0;
	conf.m_vLinkID.clear();
	conf.m_vLinkType.clear();
}

void ReadFiltration_old(CFiltration & m_Filter, ifstream & fin) {
	//	SetDefaultCFiltration(m_Filter);
	string strTemp;
	ReadUntilString(fin, "[filtration]");

	ThrowFirstOfRemain(fin, "UseFilter", strTemp);
	m_Filter.m_UseFilter = atoi(strTemp.c_str());

	ThrowFirstOfRemain(fin, "PepMassLower", strTemp);
	sscanf(strTemp.c_str(), "%lf", &m_Filter.m_lfPepMassLowerBound);
	ThrowFirstOfRemain(fin, "PepMassUpper", strTemp);
	sscanf(strTemp.c_str(), "%lf", &m_Filter.m_lfPepMassUpperBound);

	ThrowFirstOfRemain(fin, "LengthLower", strTemp);
	sscanf(strTemp.c_str(), PRI_SIZE_T, &m_Filter.m_tLengthLowerBound);
	ThrowFirstOfRemain(fin, "LengthUpper", strTemp);
	sscanf(strTemp.c_str(), PRI_SIZE_T, &m_Filter.m_tLengthUpperBound);

	vector<string> vStrTemp;
	ThrowFirstOfRemain(fin, "PepTolBase", strTemp);
	SplitStringByComma(strTemp, vStrTemp);

	for (size_t t = 0; t < vStrTemp.size(); t++) {
		double lfTmp;
		sscanf(vStrTemp[t].c_str(), "%lf", &lfTmp);
		m_Filter.m_lfPepTolBase.push_back(lfTmp);
	}
	ThrowFirstOfRemain(fin, "PepTolLower", strTemp);

	SplitStringByComma(strTemp, vStrTemp);

	for (size_t t = 0; t < vStrTemp.size(); t++) {
		double lfTmp;
		sscanf(vStrTemp[t].c_str(), "%lf", &lfTmp);
		m_Filter.m_lfPepTolLowerBound.push_back(lfTmp);
	}

	ThrowFirstOfRemain(fin, "PepTolUpper", strTemp);
	SplitStringByComma(strTemp, vStrTemp);

	for (size_t t = 0; t < vStrTemp.size(); t++) {
		double lfTmp;
		sscanf(vStrTemp[t].c_str(), "%lf", &lfTmp);
		m_Filter.m_lfPepTolUpperBound.push_back(lfTmp);
	}

	ThrowFirstOfRemain(fin, "PepTolType", m_Filter.m_strPepTolType);

	ThrowFirstOfRemain(fin, "ChargeState", strTemp);

	vector<string> vstrTemp;
	vstrTemp.clear();

	SplitStringByComma(strTemp, vstrTemp);

	m_Filter.m_vChargeState.clear();
	for (size_t t = 0; t < vstrTemp.size(); t++)
		m_Filter.m_vChargeState.push_back(atoi(vstrTemp[t].c_str()));
	vstrTemp.clear();

	ThrowFirstOfRemain(fin, "RankLimit", strTemp);
	sscanf(strTemp.c_str(), PRI_SIZE_T, &m_Filter.m_tRankLimit);

	ThrowFirstOfRemain(fin, "Redundant", strTemp);
	m_Filter.m_bRedundant = atoi(strTemp.c_str());

	ThrowFirstOfRemain(fin, "FixedDeltCn", strTemp);
	m_Filter.m_bFixedDeltCn = atoi(strTemp.c_str());

	ThrowFirstOfRemain(fin, "DeltCnLower", strTemp);
	sscanf(strTemp.c_str(), "%lf", &m_Filter.m_lfFixedDeltCn);
	ThrowFirstOfRemain(fin, "UseFDR", strTemp);
	m_Filter.m_bUseFDR = atoi(strTemp.c_str());
	ThrowFirstOfRemain(fin, "FDR", strTemp);
	sscanf(strTemp.c_str(), "%lf", &m_Filter.m_lfFDR);

	ThrowFirstOfRemain(fin, "FDRFormula", strTemp);
	sscanf(strTemp.c_str(), "%d", &m_Filter.m_FDRFormula);

	ThrowFirstOfRemain(fin, "ScoreMin", strTemp);
	sscanf(strTemp.c_str(), "%lf", &m_Filter.m_ScoreMin);

	ThrowFirstOfRemain(fin, "ScoreMax", strTemp);
	sscanf(strTemp.c_str(), "%lf", &m_Filter.m_ScoreMax);

	ThrowFirstOfRemain(fin, "DistinctPep", strTemp);
	sscanf(strTemp.c_str(), PRI_SIZE_T, &m_Filter.m_nDistinctPepLimit);

	ThrowFirstOfRemain(fin, "DistinctSpec", strTemp);
	sscanf(strTemp.c_str(), PRI_SIZE_T, &m_Filter.m_nDistinctSpecLimit);

	ThrowFirstOfRemain(fin, "ModSites", strTemp);
	m_Filter.m_strModSites.clear();
	for (size_t t = 0; t < strTemp.length(); t++) {
		if (strTemp[t] >= 'A' && strTemp[t] <= 'Z')
			m_Filter.m_strModSites += strTemp[t];
	}
	/////////////////////////////////////////////////////////
	//todo
	ThrowFirstOfRemain(fin, "CTerminal", strTemp);
	m_Filter.m_strCTerminal.clear();

	for (size_t t = 0; t < strTemp.length(); t++) {
		if (strTemp[t] >= 'A' && strTemp[t] <= 'Z')
			m_Filter.m_strCTerminal += strTemp[t];
	}

	ThrowFirstOfRemain(fin, "NTerminal", strTemp);
	m_Filter.m_strNTerminal.clear();
	for (size_t t = 0; t < strTemp.length(); t++) {
		if (strTemp[t] >= 'A' && strTemp[t] <= 'Z')
			m_Filter.m_strNTerminal += strTemp[t];
	}

	/////////////////////////////////////////////////////////
}

void WriteFiltration(const CFiltration & m_Filter, FILE * fp) {
	fprintf(fp, "[filtration]\n");
	fprintf(fp, "UseFilter=%d\n", m_Filter.m_UseFilter);
	fprintf(fp, "PepMassLower=%lf\n", m_Filter.m_lfPepMassLowerBound);
	fprintf(fp, "PepMassUpper=%lf\n", m_Filter.m_lfPepMassUpperBound);
	fprintf(fp, "LengthLower="PRI_SIZE_T"\n", m_Filter.m_tLengthLowerBound);
	fprintf(fp, "LengthUpper="PRI_SIZE_T"\n", m_Filter.m_tLengthUpperBound);
	fprintf(fp, "Separate=%d\n", m_Filter.m_Sepatrate);

	fprintf(fp, "PepTolBase=");
	fprintf(fp, "%lf", m_Filter.m_lfPepTolBase[0]);
	for (size_t t = 1; t < m_Filter.m_lfPepTolBase.size(); t++) {
		fprintf(fp, ",%lf", m_Filter.m_lfPepTolBase[t]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "PepTolLower=");
	fprintf(fp, "%lf", m_Filter.m_lfPepTolLowerBound[0]);
	for (size_t t = 1; t < m_Filter.m_lfPepTolLowerBound.size(); t++) {
		fprintf(fp, ",%lf", m_Filter.m_lfPepTolLowerBound[t]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "PepTolUpper=");
	fprintf(fp, "%lf", m_Filter.m_lfPepTolUpperBound[0]);
	for (size_t t = 1; t < m_Filter.m_lfPepTolUpperBound.size(); t++) {
		fprintf(fp, ",%lf", m_Filter.m_lfPepTolUpperBound[t]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "PepTolType=%s\n", m_Filter.m_strPepTolType.c_str());

	fprintf(fp, "ChargeState=");
	for (size_t t = 0; t < m_Filter.m_vChargeState.size() - 1; t++)
		fprintf(fp, "%d,", m_Filter.m_vChargeState[t]);
	fprintf(fp, "%d\n",
			m_Filter.m_vChargeState[m_Filter.m_vChargeState.size() - 1]);

	fprintf(fp, "RankLimit="PRI_SIZE_T"\n", m_Filter.m_tRankLimit);

	fprintf(fp, "Redundant=%d\n", m_Filter.m_bRedundant);

	fprintf(fp, "FixedDeltCn=%d\n", m_Filter.m_bFixedDeltCn);

	fprintf(fp, "DeltCnLower=%lf\n", m_Filter.m_lfFixedDeltCn);

	fprintf(fp, "UseFDR=%d\n", m_Filter.m_bUseFDR);

	fprintf(fp, "FDR=%lf\n", m_Filter.m_lfFDR);

	fprintf(fp, "FDRFormula=%d\n", m_Filter.m_FDRFormula);

	fprintf(fp, "ScoreMin=%lf\n", m_Filter.m_ScoreMin);
	fprintf(fp, "ScoreMax=%lf\n", m_Filter.m_ScoreMax);

	fprintf(fp, "DistinctPep="PRI_SIZE_T"\n", m_Filter.m_nDistinctPepLimit);
	fprintf(fp, "DistinctSpec="PRI_SIZE_T"\n", m_Filter.m_nDistinctSpecLimit);

	fprintf(fp, "ModSites=%s\n", m_Filter.m_strModSites.c_str());
	////////////////////////////////////////////////////////////////
	fprintf(fp, "CTerminal=%s\n", m_Filter.m_strCTerminal.c_str());
	fprintf(fp, "NTerminal=%s\n", m_Filter.m_strNTerminal.c_str());
	////////////////////////////////////////////////////////////////
	fprintf(fp, "[END]\n");
}

void WriteFiltrationHTML(const CFiltration & m_Filter, FILE * fp) {
	fprintf(fp, "<br>[filtration]\n");
	fprintf(fp, "<BR>PepMassLower=%lf\n", m_Filter.m_lfPepMassLowerBound);
	fprintf(fp, "<BR>PepMassUpper=%lf\n", m_Filter.m_lfPepMassUpperBound);
	fprintf(fp, "<BR>LengthLower="PRI_SIZE_T"\n", m_Filter.m_tLengthLowerBound);
	fprintf(fp, "<BR>LengthUpper="PRI_SIZE_T"\n", m_Filter.m_tLengthUpperBound);
	fprintf(fp, "<br>Separate=%d\n", m_Filter.m_Sepatrate);
	fprintf(fp, "<BR>PepTolBase=");
	fprintf(fp, "%lf", m_Filter.m_lfPepTolBase[0]);
	for (size_t t = 1; t < m_Filter.m_lfPepTolBase.size(); t++) {
		fprintf(fp, ",%lf", m_Filter.m_lfPepTolBase[t]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "<BR>PepTolLower=");
	fprintf(fp, "%lf", m_Filter.m_lfPepTolLowerBound[0]);
	for (size_t t = 1; t < m_Filter.m_lfPepTolLowerBound.size(); t++) {
		fprintf(fp, ",%lf", m_Filter.m_lfPepTolLowerBound[t]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "<BR>PepTolUpper=");
	fprintf(fp, "%lf", m_Filter.m_lfPepTolUpperBound[0]);
	for (size_t t = 1; t < m_Filter.m_lfPepTolUpperBound.size(); t++) {
		fprintf(fp, ",%lf", m_Filter.m_lfPepTolUpperBound[t]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "<BR>PepTolType=%s\n", m_Filter.m_strPepTolType.c_str());
	fprintf(fp, "<BR>ChargeState=");
	for (size_t t = 0; t < m_Filter.m_vChargeState.size() - 1; t++)
		fprintf(fp, "%d,", m_Filter.m_vChargeState[t]);
	fprintf(fp, "%d\n",
			m_Filter.m_vChargeState[m_Filter.m_vChargeState.size() - 1]);
	fprintf(fp, "<BR>Redundant=%d\n", m_Filter.m_bRedundant);
	fprintf(fp, "<BR>FixedDeltCn=%d\n", m_Filter.m_bFixedDeltCn);
	fprintf(fp, "<BR>DeltCnLower=%lf\n", m_Filter.m_lfFixedDeltCn);
	fprintf(fp, "<BR>FDR=%lf\n", m_Filter.m_lfFDR);
	if (m_Filter.m_FDRFormula == 2)
		fprintf(fp, "<BR>FDRFormula=2*F/(T+F)");
	else
		fprintf(fp, "FDRFormula=F/T");
	fprintf(fp, "<BR>ScoreMin=%lf\n", m_Filter.m_ScoreMin);
	fprintf(fp, "<BR>ScoreMax=%lf\n", m_Filter.m_ScoreMax);
	fprintf(fp, "<BR>DistinctPep="PRI_SIZE_T"\n", m_Filter.m_nDistinctPepLimit);
	fprintf(fp, "<BR>ModSites=%s\n", m_Filter.m_strModSites.c_str());
}

void ReadCConf_old(CConf & conf, const string & m_strParamFilePath) {
	ifstream fin(m_strParamFilePath.c_str());
	if (!fin.good()) {
		CErrInfo info("ReadWrite.cpp", "ReadCConf",
				"Cannot open the file: " + m_strParamFilePath);
		throw runtime_error(info.Get());
	}

	DataSetInputPathInfo SearchEngineTemp;
	string strTemp;

	int Total_Samples, SubItems, Valid;

	ThrowFirstOfRemain(fin, "SampleID", conf.m_SampleID);
	ThrowFirstOfRemain(fin, "Total_Samples", strTemp);

	sscanf(strTemp.c_str(), "%d", &Total_Samples);
	for (int i = 1; i <= Total_Samples; i++) {
		SearchEngineTemp.m_nDataSetID = i;
		ThrowFirstOfRemain(fin, "EngineType", strTemp);
		SearchEngineTemp.m_Type = GetEnginType(strTemp);
		ThrowFirstOfRemain(fin, "SubItems", strTemp);
		sscanf(strTemp.c_str(), "%d", &SubItems);
		for (int j = 1; j <= SubItems; j++) {
			GetLine(fin, strTemp);
			size_t pos1 = strTemp.find_first_of('=');
			size_t pos2 = strTemp.find_first_of(',');

			SearchEngineTemp.m_vResultPath.push_back(
					strTemp.substr(pos1 + 1, pos2 - pos1 - 1));

			strTemp = strTemp.substr(pos2 + 1);
			sscanf(strTemp.c_str(), "%d", &Valid);
			SearchEngineTemp.m_Valid.push_back(Valid);
		}
		conf.m_vEngineInputInfo.push_back(SearchEngineTemp);
		SearchEngineTemp.m_Valid.clear();
		SearchEngineTemp.m_vResultPath.clear();
	}

	ReadUntilString(fin, "[DATABASE INFORMATION]");
	ThrowFirstOfRemain(fin, "Fasta_File_Path", strTemp);

	GetNSplitStringByComma(strTemp, conf.m_vProDBPath);
	ThrowFirstOfRemain(fin, "Decoy_Tags", strTemp);

	GetNSplitStringByComma(strTemp, conf.m_vDecoyTags);
	ReadFiltration_old(conf.m_Filter, fin);

	ThrowFirstOfRemain(fin, "OutPutPath", conf.m_OutPutForder);

	DIR * dir = opendir(conf.m_OutPutForder.c_str());
	if (NULL == dir) {
		CErrInfo info("ReadWrite.cpp", "ReadCConf",
				"Cannot open the forder: " + conf.m_OutPutForder);
		throw runtime_error(info.Get());
	}

	ThrowFirstOfRemain(fin, "OutPutName", conf.m_OutPutFile);

	if (conf.m_OutPutForder[conf.m_OutPutForder.length() - 1] != '/')
		conf.m_OutPutForder += '/';
	conf.m_OutPutName = conf.m_OutPutForder + conf.m_OutPutFile;

	ThrowFirstOfRemain(fin, "ExportFormat", strTemp);

	vector<string> vstrTemp;
	SplitStringByComma(strTemp, vstrTemp);
	for (size_t t = 0; t < vstrTemp.size(); t++) {
		if (0 == atoi(vstrTemp[t].c_str()))
			conf.ExportFormat.push_back(ST_HTMLFormat);
		else if (1 == atoi(vstrTemp[t].c_str()))
			conf.ExportFormat.push_back(ST_PBULIDPREFormat);
		else if (2 == atoi(vstrTemp[t].c_str()))
			conf.ExportFormat.push_back(ST_PEPXMLFormat);
		else
			conf.ExportFormat.push_back(ST_NEXT_Format);
	}

	fin.close();
}

void ReadCConf(CConf & conf, const string & m_strParamFilePath) {
	SetDefault(conf);
	ifstream fin(m_strParamFilePath.c_str());
	if (!fin.good()) {
		CErrInfo info("ReadWrite.cpp", "ReadCConf",
				"Cannot open the file: " + m_strParamFilePath);
		throw runtime_error(info.Get());
	}
	DataSetInputPathInfo SearchEngineTemp;
	string strTemp, strResult;
	vector<string> vStrTemp;
	int Total_Samples, SubItems, Valid;

	while (GetLine(fin, strTemp)) {
		if (!GetBeforeEqual(strTemp, strResult))
			continue;

		if (strResult == "SampleID") {
			conf.m_SampleID = strTemp;
		} else if (strResult == "Total_Samples") {
			sscanf(strTemp.c_str(), "%d", &Total_Samples);
			conf.m_vEngineInputInfo.clear();

			for (int i = 1; i <= Total_Samples; i++) {
				SearchEngineTemp.m_nDataSetID = i;
				ThrowFirstOfRemain(fin, "EngineType", strTemp);

				SearchEngineTemp.m_Type = GetEnginType(strTemp);
				ThrowFirstOfRemain(fin, "SubItems", strTemp);
				sscanf(strTemp.c_str(), "%d", &SubItems);

				for (int j = 1; j <= SubItems; j++) {
					GetLine(fin, strTemp);
					size_t pos1 = strTemp.find_first_of('=');
					size_t pos2 = strTemp.find_first_of(',');

					SearchEngineTemp.m_vResultPath.push_back(
							strTemp.substr(pos1 + 1, pos2 - pos1 - 1));

					strTemp = strTemp.substr(pos2 + 1);
					sscanf(strTemp.c_str(), "%d", &Valid);
					SearchEngineTemp.m_Valid.push_back(Valid);
				}

				conf.m_vEngineInputInfo.push_back(SearchEngineTemp);
				SearchEngineTemp.m_Valid.clear();
				SearchEngineTemp.m_vResultPath.clear();
			}
		} else if (strResult == "Fasta_File_Path") {
			conf.m_vProDBPath.clear();
			GetNSplitStringByComma(strTemp, conf.m_vProDBPath);
		}

		else if (strResult == "Decoy_Tags") {
			conf.m_vDecoyTags.clear();
			GetNSplitStringByComma(strTemp, conf.m_vDecoyTags);
		}

		else if (strResult == "UseFilter") {
			conf.m_Filter.m_UseFilter = atoi(strTemp.c_str());
		}

		else if (strResult == "PepMassLower") {
			sscanf(strTemp.c_str(), "%lf",
					&conf.m_Filter.m_lfPepMassLowerBound);
		}

		else if (strResult == "PepMassUpper") {
			sscanf(strTemp.c_str(), "%lf",
					&conf.m_Filter.m_lfPepMassUpperBound);
		}

		else if (strResult == "LengthLower") {
			sscanf(strTemp.c_str(), PRI_SIZE_T, &conf.m_Filter.m_tLengthLowerBound);
		}

		else if (strResult == "LengthUpper") {
			sscanf(strTemp.c_str(), PRI_SIZE_T, &conf.m_Filter.m_tLengthUpperBound);
		} else if (strResult == "Separate") {
			conf.m_Filter.m_Sepatrate = atoi(strTemp.c_str());
		}

		else if (strResult == "PepTolBase") {
			SplitStringByComma(strTemp, vStrTemp);
			conf.m_Filter.m_lfPepTolBase.clear();

			for (size_t t = 0; t < vStrTemp.size(); t++) {
				double lfTmp;
				sscanf(vStrTemp[t].c_str(), "%lf", &lfTmp);
				conf.m_Filter.m_lfPepTolBase.push_back(lfTmp);
			}

		}

		else if (strResult == "PepTolLower") {
			SplitStringByComma(strTemp, vStrTemp);
			conf.m_Filter.m_lfPepTolLowerBound.clear();

			for (size_t t = 0; t < vStrTemp.size(); t++) {
				double lfTmp;
				sscanf(vStrTemp[t].c_str(), "%lf", &lfTmp);
				conf.m_Filter.m_lfPepTolLowerBound.push_back(lfTmp);
			}

		}

		else if (strResult == "PepTolUpper") {
			SplitStringByComma(strTemp, vStrTemp);
			conf.m_Filter.m_lfPepTolUpperBound.clear();

			for (size_t t = 0; t < vStrTemp.size(); t++) {
				double lfTmp;
				sscanf(vStrTemp[t].c_str(), "%lf", &lfTmp);
				conf.m_Filter.m_lfPepTolUpperBound.push_back(lfTmp);
			}

		}

		else if (strResult == "PepTolType") {
			conf.m_Filter.m_strPepTolType = strTemp;
		}

		else if (strResult == "ChargeState") {
			SplitStringByComma(strTemp, vStrTemp);
			conf.m_Filter.m_vChargeState.clear();

			for (size_t t = 0; t < vStrTemp.size(); t++)
				conf.m_Filter.m_vChargeState.push_back(
						atoi(vStrTemp[t].c_str()));

			vStrTemp.clear();
		}

		else if (strResult == "RankLimit") {
			sscanf(strTemp.c_str(), PRI_SIZE_T, &conf.m_Filter.m_tRankLimit);
		}

		else if (strResult == "Redundant") {
			conf.m_Filter.m_bRedundant = atoi(strTemp.c_str());
		}

		else if (strResult == "FixedDeltCn") {
			conf.m_Filter.m_bFixedDeltCn = atoi(strTemp.c_str());
		}

		else if (strResult == "DeltCnLower") {
			sscanf(strTemp.c_str(), "%lf", &conf.m_Filter.m_lfFixedDeltCn);
		}

		else if (strResult == "UseFDR") {
			conf.m_Filter.m_bUseFDR = atoi(strTemp.c_str());
		}

		else if (strResult == "FDR") {
			sscanf(strTemp.c_str(), "%lf", &conf.m_Filter.m_lfFDR);
		}

		else if (strResult == "FDRFormula") {
			sscanf(strTemp.c_str(), "%d", &conf.m_Filter.m_FDRFormula);
		}

		else if (strResult == "ScoreMin") {
			sscanf(strTemp.c_str(), "%lf", &conf.m_Filter.m_ScoreMin);
		}

		else if (strResult == "ScoreMax") {
			sscanf(strTemp.c_str(), "%lf", &conf.m_Filter.m_ScoreMax);
		}

		else if (strResult == "DistinctPep") {
			sscanf(strTemp.c_str(), PRI_SIZE_T, &conf.m_Filter.m_nDistinctPepLimit);
		}

		else if (strResult == "DistinctSpec") {
			sscanf(strTemp.c_str(), PRI_SIZE_T, &conf.m_Filter.m_nDistinctSpecLimit);
		} else if (strResult == "ModSites") {
			conf.m_Filter.m_strModSites.clear();
			for (size_t t = 0; t < strTemp.length(); t++) {
				if (strTemp[t] >= 'A' && strTemp[t] <= 'Z')
					conf.m_Filter.m_strModSites += strTemp[t];
			}

		} else if (strResult == "CTerminal") {
			conf.m_Filter.m_strCTerminal.clear();
			for (size_t t = 0; t < strTemp.length(); t++) {
				if (strTemp[t] >= 'A' && strTemp[t] <= 'Z')
					conf.m_Filter.m_strCTerminal += strTemp[t];
			}
		} else if (strResult == "NTerminal") {
			conf.m_Filter.m_strNTerminal.clear();
			for (size_t t = 0; t < strTemp.length(); t++) {
				if (strTemp[t] >= 'A' && strTemp[t] <= 'Z')
					conf.m_Filter.m_strNTerminal += strTemp[t];
			}
		} else if (strResult == "OutPutPath") {
			conf.m_OutPutForder = strTemp;
			DIR * dir = opendir(conf.m_OutPutForder.c_str());
			if (NULL == dir) {
				CErrInfo info("ReadWrite.cpp", "ReadCConf",
						"Cannot open the folder: " + conf.m_OutPutForder);
				throw runtime_error(info.Get());
			}
		}

		else if (strResult == "OutPutName") {
			conf.m_OutPutFile = strTemp;
			if (conf.m_OutPutForder[conf.m_OutPutForder.length() - 1] != cSlash)
				conf.m_OutPutForder += cSlash;
			conf.m_OutPutName = conf.m_OutPutForder + conf.m_OutPutFile;
		}

		else if (strResult == "ExportPath") {
			conf.strExportPath = strTemp;
		} else if (strResult == "ExportFormat") {
			vector<string> vstrTemp;
			SplitStringByComma(strTemp, vstrTemp);
			conf.ExportFormat.clear();

			for (size_t t = 0; t < vstrTemp.size(); t++) {
				if (0 == atoi(vstrTemp[t].c_str()))
					conf.ExportFormat.push_back(ST_INTERFACEEXPORT);
				else if (1 == atoi(vstrTemp[t].c_str()))
					conf.ExportFormat.push_back(ST_HTMLFormat);
				else if (2 == atoi(vstrTemp[t].c_str()))
					conf.ExportFormat.push_back(ST_PBULIDPREFormat);
				else if (3 == atoi(vstrTemp[t].c_str()))
					conf.ExportFormat.push_back(ST_PFINDFormat);
				else if (4 == atoi(vstrTemp[t].c_str()))
					conf.ExportFormat.push_back(ST_PEPXMLFormat);
				else
					conf.ExportFormat.push_back(ST_NEXT_Format);
			}

		}

		else if (strResult == "ExportFile") {
			conf.ExportFileType = strTemp;
		}

		else if (strResult == "ReParse") {
			conf.m_ReParse = atoi(strTemp.c_str());
		}

		else if (strResult == "PNovoParamPath") {
			conf.m_strPNovoParam = strTemp;
		} else if (strResult == "ReserveDecoy") {
			conf.m_bReserveDecoy = atoi(strTemp.c_str());
		}

		else if (strResult == "CrossLink") {
			conf.m_bCrossLink = atoi(strTemp.c_str());
		} else if (strResult == "LinkerID") {
			SplitStringByComma(strTemp, vStrTemp);
			conf.m_vLinkID.clear();

			for (size_t t = 0; t < vStrTemp.size(); t++)
				conf.m_vLinkID.push_back(atoi(vStrTemp[t].c_str()));
			vStrTemp.clear();
		} else if (strResult == "LinkerType") {
			SplitStringByComma(strTemp, vStrTemp);
			conf.m_vLinkType.clear();

			for (size_t t = 0; t < vStrTemp.size(); t++)
				conf.m_vLinkType.push_back(atoi(vStrTemp[t].c_str()));
			vStrTemp.clear();
		}

	}

	fin.close();
}

void WriteCConf(const CConf & conf, const string & strFilePath) {
	FILE * fp = fopen(strFilePath.c_str(), "w");
	if (!fp) {
		CErrInfo info("ReadWrite.cpp", "Run",
				"Cannot open the file: " + strFilePath);
		throw runtime_error(info.Get());
	}
	fprintf(fp, "SampleID=Task1\n");

	int Total_Samples = conf.m_vEngineInputInfo.size();
	fprintf(fp, "Total_Samples=%d\n", Total_Samples);

	for (int i = 0; i < Total_Samples; i++) {
		fprintf(fp, "[Sample%d]\n", i + 1);
		fprintf(fp, "EngineType=%s\n",
				SearchEngineName[conf.m_vEngineInputInfo[i].m_Type]);

		size_t tsize = conf.m_vEngineInputInfo[i].m_vResultPath.size();
		fprintf(fp, "SubItems="PRI_SIZE_T"\n", tsize);
		for (size_t j = 0; j < tsize; j++)
			fprintf(fp, "SubItem"PRI_SIZE_T"=%s,%d\n", j + 1,
					conf.m_vEngineInputInfo[i].m_vResultPath[j].c_str(),
					conf.m_vEngineInputInfo[i].m_Valid[j]);

	}

	fprintf(fp, "[END]\n");

	fprintf(fp, "[DATABASE INFORMATION]\n");
	size_t tsize = conf.m_vProDBPath.size();
	fprintf(fp, "Fasta_File_Path="PRI_SIZE_T, tsize);
	for (size_t j = 0; j < tsize; j++)
		fprintf(fp, ",%s", conf.m_vProDBPath[j].c_str());
	fprintf(fp, "\n");

	tsize = conf.m_vDecoyTags.size();
	fprintf(fp, "Decoy_Tags="PRI_SIZE_T, tsize);
	for (size_t j = 0; j < tsize; j++)
		fprintf(fp, ",%s", conf.m_vDecoyTags[j].c_str());
	fprintf(fp, "\n");

	fprintf(fp, "[END]\n");

	WriteFiltration(conf.m_Filter, fp);

	fprintf(fp, "OutPutPath=%s\n", conf.m_OutPutForder.c_str());
	fprintf(fp, "OutPutName=%s\n", conf.m_OutPutFile.c_str());
	fprintf(fp, "ExportPath=%s\n", conf.strExportPath.c_str());

	fprintf(fp, "ExportFormat=");
	for (size_t j = 0; j < conf.ExportFormat.size() - 1; j++)
		fprintf(fp, "%d,", conf.ExportFormat[j]);

	fprintf(fp, "%d\n", conf.ExportFormat[conf.ExportFormat.size() - 1]);

	fprintf(fp, "ExportFile=%s\n", conf.ExportFileType.c_str());

	fprintf(fp, "[END]\n");

	fclose(fp);
}
void WriteCConfHTML(FILE * fp, const CConf & conf) {
	fprintf(fp, "<BR>SampleID=Task1\n");
	int Total_Samples = conf.m_vEngineInputInfo.size();
	fprintf(fp, "<BR>Total_Samples=%d\n", Total_Samples);
	for (int i = 0; i < Total_Samples; i++) {
		fprintf(fp, "<BR>[Sample%d]\n", i + 1);
		fprintf(fp, "<BR>EngineType=%s\n",
				SearchEngineName[conf.m_vEngineInputInfo[i].m_Type]);
		size_t tsize = conf.m_vEngineInputInfo[i].m_vResultPath.size();
		for (size_t j = 0; j < tsize; j++)
			fprintf(fp, "<BR>SubItem"PRI_SIZE_T"=%s,%d\n", j + 1,
					conf.m_vEngineInputInfo[i].m_vResultPath[j].c_str(),
					conf.m_vEngineInputInfo[i].m_Valid[j]);
	}
	fprintf(fp, "<BR>[END]\n");
	fprintf(fp, "<BR>[DATABASE INFORMATION]\n");
	size_t tsize = conf.m_vProDBPath.size();
	fprintf(fp, "<BR>Fasta_File_Path="PRI_SIZE_T, tsize);
	for (size_t j = 0; j < tsize; j++)
		fprintf(fp, ",%s", conf.m_vProDBPath[j].c_str());
	fprintf(fp, "\n");
	tsize = conf.m_vDecoyTags.size();
	fprintf(fp, "<BR>Decoy_Tags="PRI_SIZE_T, tsize);
	for (size_t j = 0; j < tsize; j++)
		fprintf(fp, ",%s", conf.m_vDecoyTags[j].c_str());
	fprintf(fp, "\n");
	fprintf(fp, "<BR>[END]\n");
	WriteFiltrationHTML(conf.m_Filter, fp);
}

void ReadSingeSpectra(const SPECTRA_PATH & pairPath,
		CMatchSpectraInfo & SpectraTemp) {
	ifstream fin(pairPath.first.c_str());
	if (!fin.good()) {
		CErrInfo info("ReadWrite.cpp", "ReadSingeSpectra",
				"Cannot open the file: " + pairPath.first);
		throw runtime_error(info.Get());
	}
	fin.seekg(pairPath.second, ios::beg);
	SpectraTemp.m_vPeptides.clear();
	ReadSpectra(fin, SpectraTemp);
	fin.close();
}

void ReadSingeSpectraFisrtPep(const SPECTRA_PATH & pairPath,
		CMatchSpectraInfo & SpectraTemp) {
	ifstream fin(pairPath.first.c_str());
	if (!fin.good()) {
		CErrInfo info("ReadWrite.cpp", "ReadSingeSpectraFisrtPep",
				"Cannot open the file: " + pairPath.first);
		cout << info.Get() << endl;
		throw runtime_error(info.Get());
	}

	fin.seekg(pairPath.second, ios::beg);

	SpectraTemp.m_vPeptides.clear();
	for (size_t t = 0; t < SpectraTemp.m_vPeptides.size(); t++) {
		SpectraTemp.m_vPeptides[t].m_vMod.clear();
		SpectraTemp.m_vPeptides[t].m_vProteinAC.clear();
		SpectraTemp.m_vPeptides[t].m_vlfScores.clear();
	}

	SpectraTemp.m_vPeptides.clear();

	string strTemp;

	GetLine(fin, strTemp);
	sscanf(strTemp.c_str(), "[Spectrum%d]", &SpectraTemp.m_nSpectraNum);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.m_nDataSetID = atoi(strTemp.c_str());

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.m_nFileID = atoi(strTemp.c_str());

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.EnginType = GetEnginType(strTemp);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.m_strFileName = strTemp;

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	sscanf(strTemp.c_str(), "%d", &SpectraTemp.m_nCharge);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	sscanf(strTemp.c_str(), "%lf", &SpectraTemp.m_lfMH);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);

	size_t ValidCandidate;
	sscanf(strTemp.c_str(), PRI_SIZE_T, &ValidCandidate);
	CMatchPeptideInfo mp;
	if (ValidCandidate < 1) {
		return;
	}

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	size_t ScoreNum;
	sscanf(strTemp.c_str(), PRI_SIZE_T, &ScoreNum);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	double score;
	sscanf(strTemp.c_str(), "%lf", &score);
	mp.m_vlfScores.push_back(score);

	if (ScoreNum == 2) {
		GetLine(fin, strTemp);
		ReturnAfternEqual(strTemp);
		sscanf(strTemp.c_str(), "%lf", &score);
		mp.m_vlfScores.push_back(score);
	}

	if (SpectraTemp.EnginType == ST_PFIND) {
		double lfTemp = mp.m_vlfScores[0];
		mp.m_vlfScores[0] = mp.m_vlfScores[1];
		mp.m_vlfScores[1] = lfTemp;
	}

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	sscanf(strTemp.c_str(), "%lf", &mp.m_lfCalc_MH);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	sscanf(strTemp.c_str(), "%lf", &mp.m_lfDelta);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	sscanf(strTemp.c_str(), "%lf", &mp.m_lfPPM);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	mp.m_cPrev = strTemp[0];
	mp.m_cNext = strTemp[1];

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	mp.m_strSQ = strTemp;

	GetLine(fin, strTemp);
	size_t pos = strTemp.find_first_of('=');
	strTemp = strTemp.substr(pos + 1);
	GetNSplitStringByComma(strTemp, mp.m_vProteinAC);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	size_t ModificationNum;
	sscanf(strTemp.c_str(), PRI_SIZE_T, &ModificationNum);

	CModificationSiteInfo ModificationTemp;
	for (size_t i = 0; i < ModificationNum; i++) {
		GetLine(fin, strTemp);
		char * ModificationTempchr = new char[strTemp.length() + 1];
		sscanf(strTemp.c_str(), "%s "PRI_SIZE_T, ModificationTempchr,
				&ModificationTemp.m_tPos);
		ModificationTemp.m_strModName = ModificationTempchr;
		delete[] ModificationTempchr;
		mp.m_vMod.push_back(ModificationTemp);
	}

	SpectraTemp.m_vPeptides.push_back(mp);

	mp.m_vMod.clear();
	mp.m_vProteinAC.clear();
	mp.m_vlfScores.clear();
	mp.clear();

	fin.close();

	return;
}

int ReadSpectra(ifstream & fin, string & strBuf,
		CMatchSpectraInfo & SpectraTemp) {
	for (size_t t = 0; t < SpectraTemp.m_vPeptides.size(); t++) {
		SpectraTemp.m_vPeptides[t].m_vMod.clear();
		SpectraTemp.m_vPeptides[t].m_vProteinAC.clear();
		SpectraTemp.m_vPeptides[t].m_vlfScores.clear();
		SpectraTemp.m_vPeptides[t].clear();
	}

	SpectraTemp.m_vPeptides.clear();
	SpectraTemp.clear();

	string strTemp;
	bool eof = GetLine(fin, strBuf, strTemp);
	if (eof == false)
		return 0;
	sscanf(strTemp.c_str(), "[Spectrum%d]", &SpectraTemp.m_nSpectraNum);

	GetLine(fin, strBuf, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.m_nDataSetID = atoi(strTemp.c_str());

	GetLine(fin, strBuf, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.m_nFileID = atoi(strTemp.c_str());

	GetLine(fin, strBuf, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.EnginType = GetEnginType(strTemp);

	GetLine(fin, strBuf, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.m_strFileName = strTemp;

	GetLine(fin, strBuf, strTemp);
	ReturnAfternEqual(strTemp);
	sscanf(strTemp.c_str(), "%d", &SpectraTemp.m_nCharge);

	GetLine(fin, strBuf, strTemp);
	ReturnAfternEqual(strTemp);
	sscanf(strTemp.c_str(), "%lf", &SpectraTemp.m_lfMH);

	GetLine(fin, strBuf, strTemp);
	ReturnAfternEqual(strTemp);
	size_t ValidCandidate;
	sscanf(strTemp.c_str(), PRI_SIZE_T, &ValidCandidate);

	CMatchPeptideInfo mp;
	for (size_t t = 1; t <= ValidCandidate; t++) {
		GetLine(fin, strBuf, strTemp);
		ReturnAfternEqual(strTemp);
		size_t ScoreNum;
		sscanf(strTemp.c_str(), PRI_SIZE_T, &ScoreNum);

		GetLine(fin, strBuf, strTemp);
		ReturnAfternEqual(strTemp);
		double score;
		sscanf(strTemp.c_str(), "%lf", &score);
		mp.m_vlfScores.push_back(score);
		if (ScoreNum == 2) {
			GetLine(fin, strBuf, strTemp);
			ReturnAfternEqual(strTemp);
			sscanf(strTemp.c_str(), "%lf", &score);
			mp.m_vlfScores.push_back(score);
		}

		if (SpectraTemp.EnginType == ST_PFIND) {
			double lfTemp = mp.m_vlfScores[0];
			mp.m_vlfScores[0] = mp.m_vlfScores[1];
			mp.m_vlfScores[1] = lfTemp;
		}

		GetLine(fin, strBuf, strTemp);
		ReturnAfternEqual(strTemp);
		sscanf(strTemp.c_str(), "%lf", &mp.m_lfCalc_MH);

		GetLine(fin, strBuf, strTemp);
		ReturnAfternEqual(strTemp);
		sscanf(strTemp.c_str(), "%lf", &mp.m_lfDelta);

		GetLine(fin, strBuf, strTemp);
		ReturnAfternEqual(strTemp);
		sscanf(strTemp.c_str(), "%lf", &mp.m_lfPPM);

		GetLine(fin, strBuf, strTemp);
		ReturnAfternEqual(strTemp);
		mp.m_cPrev = strTemp[0];
		mp.m_cNext = strTemp[1];

		GetLine(fin, strBuf, strTemp);
		ReturnAfternEqual(strTemp);
		mp.m_strSQ = strTemp;

		GetLine(fin, strBuf, strTemp);
		size_t pos = strTemp.find_first_of('=');
		strTemp = strTemp.substr(pos + 1);
		GetNSplitStringByComma(strTemp, mp.m_vProteinAC);

		GetLine(fin, strBuf, strTemp);
		ReturnAfternEqual(strTemp);
		size_t ModificationNum;
		sscanf(strTemp.c_str(), PRI_SIZE_T, &ModificationNum);

		CModificationSiteInfo ModificationTemp;
		for (size_t i = 0; i < ModificationNum; i++) {
			GetLine(fin, strBuf, strTemp);
			char * ModificationTempchr = new char[strTemp.length() + 1];
			sscanf(strTemp.c_str(), "%s "PRI_SIZE_T, ModificationTempchr,
					&ModificationTemp.m_tPos);
			ModificationTemp.m_strModName = ModificationTempchr;
			delete[] ModificationTempchr;
			mp.m_vMod.push_back(ModificationTemp);
		}

		SpectraTemp.m_vPeptides.push_back(mp);

		mp.m_vMod.clear();
		mp.m_vProteinAC.clear();
		mp.m_vlfScores.clear();
		mp.clear();

	}

	return 1;

}

int ReadSpectra(ifstream & fin, CMatchSpectraInfo & SpectraTemp) {
	for (size_t t = 0; t < SpectraTemp.m_vPeptides.size(); t++) {
		SpectraTemp.m_vPeptides[t].m_vMod.clear();
		SpectraTemp.m_vPeptides[t].m_vProteinAC.clear();
		SpectraTemp.m_vPeptides[t].m_vlfScores.clear();
	}

	SpectraTemp.m_vPeptides.clear();

	string strTemp;
	bool eof = GetLine(fin, strTemp);
	if (eof == false)
		return 0;

	sscanf(strTemp.c_str(), "[Spectrum%d]", &SpectraTemp.m_nSpectraNum);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.m_nDataSetID = atoi(strTemp.c_str());

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.m_nFileID = atoi(strTemp.c_str());

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.EnginType = GetEnginType(strTemp);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	SpectraTemp.m_strFileName = strTemp;

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	sscanf(strTemp.c_str(), "%d", &SpectraTemp.m_nCharge);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	sscanf(strTemp.c_str(), "%lf", &SpectraTemp.m_lfMH);

	GetLine(fin, strTemp);
	ReturnAfternEqual(strTemp);
	size_t ValidCandidate;
	sscanf(strTemp.c_str(), PRI_SIZE_T, &ValidCandidate);

	CMatchPeptideInfo mp;
	for (size_t t = 1; t <= ValidCandidate; t++) {
		GetLine(fin, strTemp);
		ReturnAfternEqual(strTemp);
		size_t ScoreNum;
		sscanf(strTemp.c_str(), PRI_SIZE_T, &ScoreNum);

		GetLine(fin, strTemp);
		ReturnAfternEqual(strTemp);
		double score;
		sscanf(strTemp.c_str(), "%lf", &score);
		mp.m_vlfScores.push_back(score);
		if (ScoreNum == 2) {
			GetLine(fin, strTemp);
			ReturnAfternEqual(strTemp);
			sscanf(strTemp.c_str(), "%lf", &score);
			mp.m_vlfScores.push_back(score);
		}

		if (SpectraTemp.EnginType == ST_PFIND) {
			double lfTemp = mp.m_vlfScores[0];
			mp.m_vlfScores[0] = mp.m_vlfScores[1];
			mp.m_vlfScores[1] = lfTemp;
		}

		GetLine(fin, strTemp);
		ReturnAfternEqual(strTemp);
		sscanf(strTemp.c_str(), "%lf", &mp.m_lfCalc_MH);

		GetLine(fin, strTemp);
		ReturnAfternEqual(strTemp);
		sscanf(strTemp.c_str(), "%lf", &mp.m_lfDelta);

		GetLine(fin, strTemp);
		ReturnAfternEqual(strTemp);
		sscanf(strTemp.c_str(), "%lf", &mp.m_lfPPM);

		GetLine(fin, strTemp);
		ReturnAfternEqual(strTemp);
		mp.m_cPrev = strTemp[0];
		mp.m_cNext = strTemp[1];

		GetLine(fin, strTemp);
		ReturnAfternEqual(strTemp);
		mp.m_strSQ = strTemp;

		GetLine(fin, strTemp);
		size_t pos = strTemp.find_first_of('=');
		strTemp = strTemp.substr(pos + 1);
		GetNSplitStringByComma(strTemp, mp.m_vProteinAC);

		GetLine(fin, strTemp);
		ReturnAfternEqual(strTemp);
		size_t ModificationNum;
		sscanf(strTemp.c_str(), PRI_SIZE_T, &ModificationNum);

		CModificationSiteInfo ModificationTemp;
		for (size_t i = 0; i < ModificationNum; i++) {
			GetLine(fin, strTemp);
			char * ModificationTempchr = new char[strTemp.length() + 1];
			sscanf(strTemp.c_str(), "%s "PRI_SIZE_T, ModificationTempchr,
					&ModificationTemp.m_tPos);
			ModificationTemp.m_strModName = ModificationTempchr;
			delete[] ModificationTempchr;
			mp.m_vMod.push_back(ModificationTemp);
		}

		SpectraTemp.m_vPeptides.push_back(mp);

		mp.m_vMod.clear();
		mp.m_vProteinAC.clear();
		mp.m_vlfScores.clear();
		mp.clear();

	}

	return 1;

}

void OutPutSpectraFirstPep(FILE * fout, const CMatchSpectraInfo & SpectraTemp) {
	fprintf(fout, "[Spectrum%d]\n", SpectraTemp.m_nSpectraNum);
	fprintf(fout, "DataSetID=%d\n", SpectraTemp.m_nDataSetID);
	fprintf(fout, "FileID=%d\n", SpectraTemp.m_nFileID);

	fprintf(fout, "SearchEngine=%s\n", SearchEngineName[SpectraTemp.EnginType]);
	//fprintf(fout, "Input=%s\n", GetScanNum(SpectraTemp.m_strFileName).c_str());//todo Get ScanNum
	fprintf(fout, "Input=%s\n", SpectraTemp.m_strFileName.c_str()); //todo Get ScanNum
	fprintf(fout, "Charge=%d\n", SpectraTemp.m_nCharge);
	fprintf(fout, "MH=%lf\n", SpectraTemp.m_lfMH);

	if (SpectraTemp.m_vPeptides.size() < 1) {
		fprintf(fout, "ValidCandidate=0\n");
		return;
	}

	fprintf(fout, "ValidCandidate=1\n");
	size_t k = 0;
	fprintf(fout, "NO"PRI_SIZE_T"_ScoreNum="PRI_SIZE_T"\n", k + 1,
			SpectraTemp.m_vPeptides[k].m_vlfScores.size());
	if (SpectraTemp.EnginType == ST_PFIND)
		fprintf(fout, "NO"PRI_SIZE_T"_Score=%lf\n", k + 1,
				SpectraTemp.m_vPeptides[k].m_vlfScores[1]);
	else
		fprintf(fout, "NO"PRI_SIZE_T"_Score=%lf\n", k + 1,
				SpectraTemp.m_vPeptides[k].m_vlfScores[0]);

	if (SpectraTemp.m_vPeptides[k].m_vlfScores.size() == 2) {
		if (SpectraTemp.EnginType == ST_PFIND)
			fprintf(fout, "NO"PRI_SIZE_T"_EValue=%lf\n", k + 1,
					SpectraTemp.m_vPeptides[k].m_vlfScores[0]); //only pFind
		else
			fprintf(fout, "NO"PRI_SIZE_T"_EValue=%lf\n", k + 1,
					SpectraTemp.m_vPeptides[k].m_vlfScores[1]); //only pFind
	}

	fprintf(fout, "NO"PRI_SIZE_T"_Mass=%lf\n", k + 1,
			SpectraTemp.m_vPeptides[k].m_lfCalc_MH);
	fprintf(fout, "NO"PRI_SIZE_T"_Delta=%lf\n", k + 1,
			SpectraTemp.m_vPeptides[k].m_lfDelta);
	fprintf(fout, "NO"PRI_SIZE_T"_PPM=%lf\n", k + 1, SpectraTemp.m_vPeptides[k].m_lfPPM);

	fprintf(fout, "NO"PRI_SIZE_T"_PreviousAndNext=%c%c\n", k + 1,
			SpectraTemp.m_vPeptides[k].m_cPrev,
			SpectraTemp.m_vPeptides[k].m_cNext);

	fprintf(fout, "NO"PRI_SIZE_T"_SQ=%s\n", k + 1,
			SpectraTemp.m_vPeptides[k].m_strSQ.c_str());

	fprintf(fout, "NO"PRI_SIZE_T"_Proteins="PRI_SIZE_T, k + 1,
			SpectraTemp.m_vPeptides[k].m_vProteinAC.size());
	for (size_t j = 0; j < SpectraTemp.m_vPeptides[k].m_vProteinAC.size();
			j++) {
		fprintf(fout, ",%s",
				SpectraTemp.m_vPeptides[k].m_vProteinAC[j].c_str());
	}
	fprintf(fout, "\n");

	fprintf(fout, "NO"PRI_SIZE_T"_ModificationNum="PRI_SIZE_T"\n", k + 1,
			SpectraTemp.m_vPeptides[k].m_vMod.size());

	for (size_t j = 0; j < SpectraTemp.m_vPeptides[k].m_vMod.size(); j++) {
		fprintf(fout, "%s "PRI_SIZE_T"\n",
				SpectraTemp.m_vPeptides[k].m_vMod[j].m_strModName.c_str(),
				SpectraTemp.m_vPeptides[k].m_vMod[j].m_tPos);
	}

}

void OutPutSpectra(FILE * fout, const CMatchSpectraInfo & SpectraTemp,
		FILE * ferrorfigure_Da, FILE * ferrorfigure_ppm, const CConf & conf) {
	static int spectrN = 1;
	fprintf(fout, "[Spectrum%d]\n", SpectraTemp.m_nSpectraNum);
	fprintf(fout, "DataSetID=%d\n", SpectraTemp.m_nDataSetID);

	fprintf(fout, "FileID=%d\n", SpectraTemp.m_nFileID);

	fprintf(fout, "SearchEngine=%s\n", SearchEngineName[SpectraTemp.EnginType]);

	fprintf(fout, "Input=%s\n", SpectraTemp.m_strFileName.c_str());
	fprintf(fout, "Charge=%d\n", SpectraTemp.m_nCharge);
	fprintf(fout, "MH=%lf\n", SpectraTemp.m_lfMH);

	fprintf(fout, "ValidCandidate="PRI_SIZE_T"\n", SpectraTemp.m_vPeptides.size());
	if (SpectraTemp.m_vPeptides.size() > 0) {
		fprintf(ferrorfigure_Da, "%d ", spectrN++);
		fprintf(ferrorfigure_ppm, "%d ", spectrN++);
		fprintf(ferrorfigure_Da, "%lf ",
				SpectraTemp.m_vPeptides[0].m_vlfScores[0]);
		fprintf(ferrorfigure_ppm, "%lf ",
				SpectraTemp.m_vPeptides[0].m_vlfScores[0]);
		size_t nIsRev = IsReverse(SpectraTemp, conf);
		if (nIsRev == 2 || nIsRev == 3 || nIsRev == 4) {
			fprintf(ferrorfigure_Da, "%lf 0\n",
					SpectraTemp.m_vPeptides[0].m_lfDelta);
			fprintf(ferrorfigure_ppm, "%lf 0\n",
					SpectraTemp.m_vPeptides[0].m_lfPPM);
		} else {
			fprintf(ferrorfigure_Da, "%lf 1\n",
					SpectraTemp.m_vPeptides[0].m_lfDelta);
			fprintf(ferrorfigure_ppm, "%lf 1\n",
					SpectraTemp.m_vPeptides[0].m_lfPPM);
		}
	}

	for (size_t k = 0; k < SpectraTemp.m_vPeptides.size(); k++) {
		fprintf(fout, "NO"PRI_SIZE_T"_ScoreNum="PRI_SIZE_T"\n", k + 1,
				SpectraTemp.m_vPeptides[k].m_vlfScores.size());

		if (SpectraTemp.EnginType == ST_PFIND) {
			fprintf(fout, "NO"PRI_SIZE_T"_Score=%lf\n", k + 1,
					SpectraTemp.m_vPeptides[k].m_vlfScores[1]);
		} else {
			fprintf(fout, "NO"PRI_SIZE_T"_Score=%lf\n", k + 1,
					SpectraTemp.m_vPeptides[k].m_vlfScores[0]);
		}

		if (SpectraTemp.m_vPeptides[k].m_vlfScores.size() == 2) {
			if (SpectraTemp.EnginType == ST_PFIND) {
				fprintf(fout, "NO"PRI_SIZE_T"_EValue=%lf\n", k + 1,
						SpectraTemp.m_vPeptides[k].m_vlfScores[0]); //only pFind
			} else
				fprintf(fout, "NO"PRI_SIZE_T"_EValue=%lf\n", k + 1,
						SpectraTemp.m_vPeptides[k].m_vlfScores[1]); //only pFind
		}

		fprintf(fout, "NO"PRI_SIZE_T"_Mass=%lf\n", k + 1,
				SpectraTemp.m_vPeptides[k].m_lfCalc_MH);

		fprintf(fout, "NO"PRI_SIZE_T"_Delta=%lf\n", k + 1,
				SpectraTemp.m_vPeptides[k].m_lfDelta);
		//		if (conf.m_Filter.m_strPepTolType == "Da")
		//		{
		//			if (IsReverse(SpectraTemp, conf))
		//				fprintf(ferrorfigure, "%lf 0\n", SpectraTemp.m_vPeptides[k].m_lfDelta);
		//			else
		//				fprintf(ferrorfigure, "%lf 1\n", SpectraTemp.m_vPeptides[k].m_lfDelta);
		//		}
		//		else
		//		{
		//			if (IsReverse(SpectraTemp, conf))
		//				fprintf(ferrorfigure, "%lf 0\n", SpectraTemp.m_vPeptides[k].m_lfPPM);
		//			else
		//				fprintf(ferrorfigure, "%lf 1\n", SpectraTemp.m_vPeptides[k].m_lfPPM);
		//		}

		fprintf(fout, "NO"PRI_SIZE_T"_PPM=%lf\n", k + 1,
				SpectraTemp.m_vPeptides[k].m_lfPPM);

		fprintf(fout, "NO"PRI_SIZE_T"_PreviousAndNext=%c%c\n", k + 1,
				SpectraTemp.m_vPeptides[k].m_cPrev,
				SpectraTemp.m_vPeptides[k].m_cNext);

		fprintf(fout, "NO"PRI_SIZE_T"_SQ=%s\n", k + 1,
				SpectraTemp.m_vPeptides[k].m_strSQ.c_str());

		fprintf(fout, "NO"PRI_SIZE_T"_Proteins="PRI_SIZE_T"", k + 1,
				SpectraTemp.m_vPeptides[k].m_vProteinAC.size());

		for (size_t j = 0; j < SpectraTemp.m_vPeptides[k].m_vProteinAC.size();
				j++) {
			fprintf(fout, ",%s",
					SpectraTemp.m_vPeptides[k].m_vProteinAC[j].c_str());
		}
		fprintf(fout, "\n");

		fprintf(fout, "NO"PRI_SIZE_T"_ModificationNum="PRI_SIZE_T"\n", k + 1,
				SpectraTemp.m_vPeptides[k].m_vMod.size());

		for (size_t j = 0; j < SpectraTemp.m_vPeptides[k].m_vMod.size(); j++) {
			fprintf(fout, "%s "PRI_SIZE_T"\n",
					SpectraTemp.m_vPeptides[k].m_vMod[j].m_strModName.c_str(),
					SpectraTemp.m_vPeptides[k].m_vMod[j].m_tPos);
		}
	}
}

void OutPutCondition(FILE * fp, const CConditionInfo & cond,
		const int & TotalSpectra) {
	//todo
	fprintf(fp, "[Meta]\n");
	fprintf(fp, "[Search]\n");

	fprintf(fp, "InputPath=%s\n", cond.m_strInputPath.c_str());
	fprintf(fp, "Database="PRI_SIZE_T"", cond.m_vDatabase.size());
	for (size_t j = 0; j < cond.m_vDatabase.size(); j++)
		fprintf(fp, ",%s", cond.m_vDatabase[j].c_str());
	fprintf(fp, "\n");

	fprintf(fp, "Time=UNKNOWN\n");
	fprintf(fp, "Enzyme=%s\n", cond.m_Enzyme.c_str());
	fprintf(fp, "Fixed_modifications="PRI_SIZE_T"", cond.m_vFixMod.size());

	for (size_t j = 0; j < cond.m_vFixMod.size(); j++)
		fprintf(fp, ",%s", cond.m_vFixMod[j].c_str());
	fprintf(fp, "\n");

	fprintf(fp, "Variable_modifications="PRI_SIZE_T"", cond.m_vVarMod.size());
	for (size_t j = 0; j < cond.m_vVarMod.size(); j++)
		fprintf(fp, ",%s", cond.m_vVarMod[j].c_str());
	fprintf(fp, "\n");

	if (1 == cond.m_bPepMono)
		fprintf(fp, "Peptide_Mass=Monoisotopic\n");
	else
		fprintf(fp, "Peptide_Mass=Average\n");
	//cond.m_lfPepBase = 0;//todo 注意这个暂时赋予的一个值

	fprintf(fp, "Peptide_Mass_Tolerance_Base=%lf\n", 0.0);
	//cond.m_lfFragTol = 0;

	fprintf(fp, "Peptide_Mass_Tolerance=%lf\n", 0.0);
	//cond.m_strPepTolType = "UNKONW";

	fprintf(fp, "Peptide_Mass_Tolerance_Type=%s\n", "UNKONW");

	if (1 == cond.m_bFragMono)
		fprintf(fp, "Fragment_Mass=Monoisotopic\n");
	else
		fprintf(fp, "Fragment_Mass=Average\n");
	//cond.m_lfFragBase = 0;

	fprintf(fp, "Fragment_Mass_Tolerance_Base=%lf\n", 0.0);
	//cond.m_lfFragTol = 0;

	fprintf(fp, "Fragment_Mass_Tolerance=%lf\n", 0.0);
	//cond.m_strFragTolType = "UNKNOW";

	fprintf(fp, "Fragment_Mass_Tolerance_Type=%s\n", "UNKNOW");
	//cond.m_nMaxMiss = 0;

	fprintf(fp, "Max_Missed_Cleavages=%d\n", 0);
	//cond.m_strInstrument = "UNKNOW";

	fprintf(fp, "Instrument_type=%s\n", "UNKNOW");
	fprintf(fp, "[Total]\n");
	fprintf(fp, "Spectra=%d\n", TotalSpectra);
	fprintf(fp, "Proteins=0\n");
	fprintf(fp, "FPR=0.00\n");
	fprintf(fp, "Threshold=0.00\n");
}

void ReadIndexFile(const string & strIndexFile,
		vector<OneDATASET> & AllDataSetParser) {
	ifstream fin(strIndexFile.c_str());
	if (!fin.good()) {
		CErrInfo info("ReadWrite.cpp", "ReadIndexFile",
				"Cannot open the file: " + strIndexFile);
		throw runtime_error(info.Get());
	}

	string strVal;
	ThrowFirstOfRemain(fin, "DATASIZE", strVal);

	size_t Rt;
	sscanf(strVal.c_str(), PRI_SIZE_T, &Rt);
	for (size_t datasize = 0; datasize < Rt; datasize++) {
		OneDATASET DataSetOne;
		map<int, vector<SPECTRAINFO> > DataSetSpectra;

		ThrowFirstOfRemain(fin, "DATANUM", strVal);
		ThrowFirstOfRemain(fin, "SearchEngine", strVal);

		DataSetOne.first.EngineType = GetEnginType(strVal);
		//conf.m_vEngineInputInfo

		ThrowFirstOfRemain(fin, "CHARGESIZE", strVal);

		size_t ChargeSize;
		sscanf(strVal.c_str(), PRI_SIZE_T, &ChargeSize);

		for (size_t chsi = 0; chsi < ChargeSize; chsi++) {
			ThrowFirstOfRemain(fin, "CHARGE", strVal);
			ThrowFirstOfRemain(fin, "SPECTRANUM", strVal);

			size_t spectranum;
			sscanf(strVal.c_str(), PRI_SIZE_T, &spectranum);

			for (size_t spnum = 0; spnum < spectranum; spnum++) {
				GetLine(fin, strVal);
				SPECTRA_PATH pathTemp;

				long int move;

				size_t pos = strVal.find_last_of(' ');
				pathTemp.first = strVal.substr(0, pos);
				sscanf(strVal.substr(pos + 1).c_str(), "%ld", &move);

				pathTemp.second = move;

				CMatchSpectraInfo SpectraTemp;
				ReadSingeSpectraFisrtPep(pathTemp, SpectraTemp);

				GetIndex(DataSetSpectra, SpectraTemp, pathTemp);
			}
		}

		DataSetOne.second = DataSetSpectra;
		AllDataSetParser.push_back(DataSetOne);

	}
	fin.close();
}

void WriteParserSimple(const string & strSimpleFile, const CConf & conf) {
	char chrInput[1000] = { 0 };

	string strInput = "";

	FILE * fp = fopen(strSimpleFile.c_str(), "w");
	int Total_Samples = conf.m_vEngineInputInfo.size();

	sprintf(chrInput, "Total_Samples=%d", Total_Samples);
	strInput += chrInput;

	for (int i = 0; i < Total_Samples; i++) {
		sprintf(chrInput, "[Sample%d]", i + 1);
		strInput += chrInput;

		sprintf(chrInput, "EngineType=%s",
				SearchEngineName[conf.m_vEngineInputInfo[i].m_Type]);
		strInput += chrInput;

		size_t tsize = conf.m_vEngineInputInfo[i].m_vResultPath.size();
		sprintf(chrInput, "SubItems="PRI_SIZE_T"", tsize);
		strInput += chrInput;

		for (size_t j = 0; j < tsize; j++) {
			sprintf(chrInput, "SubItem"PRI_SIZE_T"=%s,%d", j + 1,
					conf.m_vEngineInputInfo[i].m_vResultPath[j].c_str(),
					conf.m_vEngineInputInfo[i].m_Valid[j]);
			strInput += chrInput;
		}
	}

	sprintf(chrInput, "[END]");
	strInput += chrInput;

	fprintf(fp, "%s\n", strInput.c_str());

	fclose(fp);
}

void WriteIndexFile(const string & strIndexFile,
		const vector<OneDATASET> & AllDataSetParser) {
	FILE * fout = fopen(strIndexFile.c_str(), "w");

	fprintf(fout, "DATASIZE="PRI_SIZE_T"\n", AllDataSetParser.size());

	for (size_t Rt = 0; Rt < AllDataSetParser.size(); Rt++) {
		fprintf(fout, "DATANUM="PRI_SIZE_T"\n", Rt + 1);
		fprintf(fout, "SearchEngine=%s\n",
				SearchEngineName[AllDataSetParser[Rt].first.EngineType]);
		fprintf(fout, "CHARGESIZE="PRI_SIZE_T"\n", AllDataSetParser[Rt].second.size());

		for (map<int, vector<SPECTRAINFO> >::const_iterator it =
				AllDataSetParser[Rt].second.begin();
				it != AllDataSetParser[Rt].second.end(); it++) {
			fprintf(fout, "CHARGE=%d\n", it->first);
			fprintf(fout, "SPECTRANUM="PRI_SIZE_T"\n", it->second.size());

			for (size_t t = 0; t < it->second.size(); t++) {
				fprintf(fout, "%s %ld\n", it->second[t].first.first.c_str(),
						it->second[t].first.second);
			}
		}
	}

	fprintf(fout, "[END]\n");

	fclose(fout);
}
void WriteSingePeptide(FILE * fout, const CMatchPeptideInfo & PepTemp,
		const CMatchSpectraInfo & SpectraTemp) {
	if (SpectraTemp.EnginType == ST_PFIND)
		fprintf(fout, "%10.10e\t", EValueToSmall(PepTemp.m_vlfScores[0]));

	else
		fprintf(fout, "%lf\t", PepTemp.m_vlfScores[0]);

	fprintf(fout, "%lf\t", PepTemp.m_lfCalc_MH);
	fprintf(fout, "%lf\t", PepTemp.m_lfDelta);

	for (size_t t = 0; t < PepTemp.m_vMod.size(); t++) {
		fprintf(fout, ""PRI_SIZE_T" %c(%s) ; ", PepTemp.m_vMod[t].m_tPos + 1,
				PepTemp.m_strSQ[PepTemp.m_vMod[t].m_tPos],
				PepTemp.m_vMod[t].m_strModName.c_str());
	}

	if (PepTemp.m_vMod.size() == 0)
		fprintf(fout, "null");

	fprintf(fout, "\t");
	fprintf(fout, "%d\t", SpectraTemp.m_nDataSetID);

	fprintf(fout, "%d.%s\t", SpectraTemp.m_nFileID,
			SearchEngineName[SpectraTemp.EnginType]);
	fprintf(fout, "1\t");

	fprintf(fout, "%s", PepTemp.m_vProteinAC[0].c_str());
	for (size_t p = 1; p < PepTemp.m_vProteinAC.size(); p++)
		fprintf(fout, "/%s", PepTemp.m_vProteinAC[p].c_str());
	fprintf(fout, "\n");
}

void WritePeptide_Spec(FILE * fout, const CMatchSpectraInfo & SpectraTemp,
		const size_t & t) {
	fprintf(fout, "\t"PRI_SIZE_T"\t", t + 1);

	fprintf(fout, "%s\t%d\t", SpectraTemp.m_strFileName.c_str(),
			SpectraTemp.m_nCharge);

	CMatchPeptideInfo PepTemp = SpectraTemp.m_vPeptides[0];

	WriteSingePeptide(fout, PepTemp, SpectraTemp);
}

void WriteSingePeptideforCSV(FILE * fout, const CMatchPeptideInfo & PepTemp,
		const CMatchSpectraInfo & SpectraTemp) {
	if (SpectraTemp.EnginType == ST_PFIND)
		fprintf(fout, "%10.10e,", EValueToSmall(PepTemp.m_vlfScores[0]));

	else
		fprintf(fout, "%lf,", PepTemp.m_vlfScores[0]);

	fprintf(fout, "%lf,", PepTemp.m_lfCalc_MH);
	fprintf(fout, "%lf,", PepTemp.m_lfDelta);

	for (size_t t = 0; t < PepTemp.m_vMod.size(); t++) {
		fprintf(fout, ""PRI_SIZE_T" %c(%s) ; ", PepTemp.m_vMod[t].m_tPos + 1,
				PepTemp.m_strSQ[PepTemp.m_vMod[t].m_tPos],
				PepTemp.m_vMod[t].m_strModName.c_str());
	}

	if (PepTemp.m_vMod.size() == 0)
		fprintf(fout, "null");
	fprintf(fout, ",");
	fprintf(fout, "%d,", SpectraTemp.m_nDataSetID);
	fprintf(fout, "%d.%s,", SpectraTemp.m_nFileID,
			SearchEngineName[SpectraTemp.EnginType]);

	fprintf(fout, "0,");
	fprintf(fout, "-1,");
	fprintf(fout, "1,");

	fprintf(fout, "%s", PepTemp.m_vProteinAC[0].c_str());
	for (size_t p = 1; p < PepTemp.m_vProteinAC.size(); p++)
		fprintf(fout, "/%s", PepTemp.m_vProteinAC[p].c_str());

	fprintf(fout, "\n");
}

void GetProtein_Peptide(FILE * fout, const CMatchSpectraInfo & SpectraTemp,
		const int & i) {
	fprintf(fout, "\t%d\t", i + 1);

	fprintf(fout, "%s\t%d\t", SpectraTemp.m_strFileName.c_str(),
			SpectraTemp.m_nCharge);

	fprintf(fout, "%c.%s.%c\t", SpectraTemp.m_vPeptides[0].m_cPrev,
			SpectraTemp.m_vPeptides[0].m_strSQ.c_str(),
			SpectraTemp.m_vPeptides[0].m_cNext);

	CMatchPeptideInfo PepTemp = SpectraTemp.m_vPeptides[0];

	WriteSingePeptide(fout, PepTemp, SpectraTemp);

}
void GetProteinHead(FILE * fout, const string & strPro,
		const vector<SPECTRAINFO> & SpecPath) {
	set<string> setUniPep;

	for (size_t t = 0; t < SpecPath.size(); t++)
		setUniPep.insert(SpecPath[t].second.m_strFirstPeptide);

	fprintf(fout, "\t%s\t%lf\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t\n", strPro.c_str(), 0.0,
			SpecPath.size(), setUniPep.size());
}

double CalcCoverage(FILE * fout, const CMatchProteinInfo & ProteinTemp,
		map<string, vector<SPECTRAINFO> > & protein_peptide) {
	double lfRet = 0.0;
	const string & strSQ = ProteinTemp.m_strSQ;

	if (protein_peptide[ProteinTemp.m_strAC].size() == 0) {
		fprintf(fout, "null\n");
		return lfRet;
	}

	int * pCov = new int[strSQ.size()];
	memset(pCov, 0, sizeof(int) * strSQ.size());

	for (size_t i = 0; i < protein_peptide[ProteinTemp.m_strAC].size(); ++i) {
		const string & strPep =
				protein_peptide[ProteinTemp.m_strAC][i].second.m_strFirstPeptide;

		int size = strPep.size();
		size_t nPos1 = 0, nPos2 = 0;

		while (1) {
			if ((nPos2 = strSQ.find(strPep, nPos1)) != string::npos) {
				fill(pCov + nPos2, pCov + nPos2 + size, true);
				nPos1 = nPos2 + 1;
			} else {
				break;
			}
		}
	}

	lfRet = count(pCov, pCov + strSQ.size(), 1) / (double) strSQ.size();

	int tag = 0;
	for (size_t i = 0; i < strSQ.size(); i++) {
		if (pCov[i] == 1) {
			fprintf(fout, ""PRI_SIZE_T";", i);
			tag = 1;
		}
	}

	if (tag == 1)
		fprintf(fout, "\n");
	else
		fprintf(fout, "null\n");

	delete[] pCov;

	return lfRet;
}

void ReadProteinANDCalculation(const CConf & conf,
		map<string, CProteinInfo> & ProteinInfo,
		const map<string, vector<SPECTRAINFO> > & m_protein_peptide) {
	//cout << "Read Protein Database..." << endl;
	osspBuildLog << "Read Protein Database..." << endl;

	map<string, vector<SPECTRAINFO> > protein_peptide = m_protein_peptide;
	CBioCalculator BioCalculator;

	string strProDatabase = conf.m_OutPutName;
	strProDatabase = conf.m_outPutForder_Java
			+ "ProteinDatabase_pBuild.protein";
	FILE * fout = fopen(strProDatabase.c_str(), "w");

	string strFASTA = conf.m_outPutForder_Index + conf.m_OutPutFile + ".fasta";
	//cout << strFASTA << endl;
	FILE * fFASTA = fopen(strFASTA.c_str(), "w");

	int ProteinNum = 0;

	for (size_t t = 0; t < conf.m_vProDBPath.size(); t++) {
		ifstream fin(conf.m_vProDBPath[t].c_str());
		if (!fin.good()) {
			CErrInfo info("ReadWrite.cpp", "ReadProteinDATABASE",
					"Cannot open the file: " + conf.m_vProDBPath[t]);
			throw runtime_error(info.Get());
		}

		string ProteinSQ = "";
		string ProteinAC = "";
		string ProteinDE = "";
		string strRead = "";

		GetLine(fin, strRead);

		size_t pos = strRead.find_first_of(' ');
		ProteinAC = strRead.substr(1, pos - 1);
		ProteinDE = strRead.substr(pos + 1);

		while (GetLine(fin, strRead)) {
			if (strRead[0] == '>') {
				CMatchProteinInfo ProteinTemp;
				ProteinTemp.m_strAC = ProteinAC;
				ProteinTemp.m_strDE = ProteinDE;
				ProteinTemp.m_strSQ = ProteinSQ;

				fprintf(fout, "%s\n", ProteinAC.c_str());
				fprintf(fout, "%s\n", ProteinDE.c_str());
				PrintProteinDatabase(fout, ProteinSQ);

				CProteinInfo ProteinInfoTemp;
				ProteinInfoTemp.m_Coverage = CalcCoverage(fout, ProteinTemp,
						protein_peptide);
				ProteinInfoTemp.m_MW = BioCalculator.CalcMW(ProteinTemp.m_strSQ,
						false);
				ProteinInfoTemp.m_pI = BioCalculator.CalcPI(
						ProteinTemp.m_strSQ);
				ProteinInfoTemp.m_strDE = ProteinDE;

				ProteinInfo[ProteinAC] = ProteinInfoTemp;

				fprintf(fout, "%.2lf %lf %lf\n",
						ProteinInfoTemp.m_Coverage * 100, ProteinInfoTemp.m_MW,
						ProteinInfoTemp.m_pI);
				if (ProteinInfoTemp.m_Coverage > 0) {
					fprintf(fFASTA, ">%s %s\n", ProteinAC.c_str(),
							ProteinDE.c_str());
					//	fprintf(fFASTA, "%s\n", ProteinSQ.c_str());
					PrintProteinDatabase(fFASTA, ProteinSQ);
				}

				++ProteinNum;
				fprintf(fout, "[END%d]\n", ProteinNum);
				ProteinSQ = "";
				pos = strRead.find_first_of(' ');

				ProteinAC = strRead.substr(1, pos - 1);
				ProteinDE = strRead.substr(pos + 1);
			} else {
				for (size_t strReadt = 0; strReadt < strRead.size();
						strReadt++) {
					if (strRead[strReadt] >= 'A' && strRead[strReadt] <= 'Z')
						ProteinSQ += strRead[strReadt];
				}
				//ProteinSQ += strRead.substr(0, strRead.size());
			}
		}

		CMatchProteinInfo ProteinTemp;
		ProteinTemp.m_strAC = ProteinAC;
		ProteinTemp.m_strDE = ProteinDE;
		ProteinTemp.m_strSQ = ProteinSQ;

		fprintf(fout, "%s\n", ProteinAC.c_str());
		fprintf(fout, "%s\n", ProteinDE.c_str());
		PrintProteinDatabase(fout, ProteinSQ);
		//////////////////////////////////////////////////////////////////////////////
		CProteinInfo ProteinInfoTemp;
		ProteinInfoTemp.m_Coverage = CalcCoverage(fout, ProteinTemp,
				protein_peptide);
		ProteinInfoTemp.m_MW = BioCalculator.CalcMW(ProteinTemp.m_strSQ, false);
		ProteinInfoTemp.m_pI = BioCalculator.CalcPI(ProteinTemp.m_strSQ);
		ProteinInfoTemp.m_strDE = ProteinDE;

		if (ProteinInfoTemp.m_Coverage > 0) {
			fprintf(fFASTA, ">%s %s\n", ProteinAC.c_str(), ProteinDE.c_str());
			//fprintf(fFASTA, "%s\n", ProteinSQ.c_str());
			PrintProteinDatabase(fFASTA, ProteinSQ);
		}

		ProteinInfo[ProteinAC] = ProteinInfoTemp;
		fprintf(fout, "%.2lf %lf %lf\n", ProteinInfoTemp.m_Coverage * 100,
				ProteinInfoTemp.m_MW, ProteinInfoTemp.m_pI);
		++ProteinNum;
		fprintf(fout, "[END%d]\n", ProteinNum);
		fin.close();
	}

	fclose(fout);

	fclose(fFASTA);
}

//void ReadFromIntermediateFromWriteSinglePeptide(istringstream & iss, string & Score,
//		string & Calc_M, string & Delta_M, string & ppm, string & Mod_Sites, string & SampleID,
//		string & Engine, string & MatchedIons, string & MissCleaveNum, string & Rank,
//		string & Proteins)
//{
//	iss >> Score >> Calc_M >> Delta_M >> ppm >> Mod_Sites >> SampleID >> Engine;
//
//	iss >> MatchedIons >> MissCleaveNum >> Rank >> Proteins;
//}
void ReadFromIntermediateFromWriteSinglePeptide(istringstream & iss,
		string & Score, string & Calc_M, string & Delta_M, string & ppm,
		string & Mod_Sites, string & SampleID, string & Engine,
		string & MatchedIons, string & MissCleaveNum, string & Rank,
		string & Proteins) {
	//cout << iss.str() << endl;
	getline(iss, Score, '\t');
	//cout << Score << endl;
	getline(iss, Calc_M, '\t');
	//cout << Calc_M << endl;
	getline(iss, Delta_M, '\t');
	//cout << Delta_M << endl;
	getline(iss, ppm, '\t');
	//cout << ppm << endl;
	getline(iss, Mod_Sites, '\t');
	//cout << Mod_Sites << endl;
	getline(iss, SampleID, '\t');
	//cout << SampleID << endl;
	getline(iss, Engine, '\t');
	//cout << Engine << endl;
	getline(iss, MatchedIons, '\t');
	//cout << MatchedIons << endl;
	getline(iss, MissCleaveNum, '\t');
	//cout << MissCleaveNum << endl;
	getline(iss, Rank, '\t');
	//cout << Rank << endl;
	getline(iss, Proteins, '\n');
	//cout << Proteins << endl;

	//cout << "=======================================" << endl;

	//iss >> Score >> Calc_M >> Delta_M >> ppm >> Mod_Sites >> SampleID >> Engine;

	//iss >> MatchedIons >> MissCleaveNum >> Rank >> Proteins;
}

void ReadFromIntermediateFileSpectra1(istringstream & iss, string & order,
		string & Spectrum, string & PeptideNum, string & UniquePepNum,
		string & Samples, string & Score, string & Condition) {
	//iss >> order >> Spectrum >> PeptideNum >> UniquePepNum >> Samples >> Score >> Condition;
	getline(iss, order, '\t');
	getline(iss, Spectrum, '\t');
	getline(iss, PeptideNum, '\t');
	getline(iss, UniquePepNum, '\t');
	getline(iss, Samples, '\t');
	getline(iss, Score, '\t');
	getline(iss, Condition, '\n');
}

void ReadFromIntermediateFileSpectra2(istringstream & iss, string & strOrder,
		string & Peptide, string & Score, string & Calc_M, string & Delta_M,
		string & ppm, string & Mod_Sites, string & SampleID, string & Engine,
		string & MatchedIons, string & MissCleaveNum, string & Rank,
		string & Proteins) {
	string strTmp;
	//iss >> strTmp >> strOrder >> Peptide;
	getline(iss, strTmp, '\t');
	getline(iss, strOrder, '\t');
	getline(iss, Peptide, '\t');

	ReadFromIntermediateFromWriteSinglePeptide(iss, Score, Calc_M, Delta_M, ppm,
			Mod_Sites, SampleID, Engine, MatchedIons, MissCleaveNum, Rank,
			Proteins);
}

void ReadFromIntermediateFilePeptide1(istringstream & iss, string & order,
		string & SQ, string & Spectra, string & Samples, string & Score,
		string & Proteins) {
	//iss >> order >> SQ >> Spectra >> Samples >> Score >> Proteins;
	getline(iss, order, '\t');
	getline(iss, SQ, '\t');
	getline(iss, Spectra, '\t');
	getline(iss, Samples, '\t');
	getline(iss, Score, '\t');
	getline(iss, Proteins, '\n');
}

void ReadFromIntermediateFilePeptide2(istringstream & iss, string & strOrder,
		string & Spectra, string & Score, string & Calc_M, string & Delta_M,
		string & ppm, string & Mod_Sites, string & SampleID, string & Engine,
		string & MatchedIons, string & MissCleaveNum, string & Rank,
		string & Proteins) {
	string strTmp;
	//iss >> strTmp >> strOrder >> Spectra;
	getline(iss, strTmp, '\t');
	getline(iss, strOrder, '\t');
	getline(iss, Spectra, '\t');

	ReadFromIntermediateFromWriteSinglePeptide(iss, Score, Calc_M, Delta_M, ppm,
			Mod_Sites, SampleID, Engine, MatchedIons, MissCleaveNum, Rank,
			Proteins);
}

void ReadFromIntermediateFileProtein1(istringstream & iss, string & order,
		string & ProteinAC, string & MW, string & pI, string & Coverage,
		string & UniquePepNum, string & SpecNum, string & NonModifiedSpecNum,
		string & ModifiedSpecNum, string & UniqueModifiedPepNum,
		string& Samples, string & Description) {

	//todo 改成遇到/t 或者 /n 停止读入，而不是空格或者其他的。
	//cout << iss.str() << endl;
	//	iss >> order >> ProteinAC >> MW >> pI >> Coverage >> UniquePepNum >> SpecNum
	//			>> NonModifiedSpecNum >> ModifiedSpecNum >> UniqueModifiedPepNum >> Samples;
	getline(iss, order, '\t');
	getline(iss, ProteinAC, '\t');
	getline(iss, MW, '\t');
	getline(iss, pI, '\t');
	getline(iss, Coverage, '\t');
	getline(iss, UniquePepNum, '\t');
	getline(iss, SpecNum, '\t');
	getline(iss, NonModifiedSpecNum, '\t');
	getline(iss, ModifiedSpecNum, '\t');
	getline(iss, UniqueModifiedPepNum, '\t');
	getline(iss, Samples, '\t');
	getline(iss, Description, '\n');
	//	size_t t = order.size() + 1;
	//	t += ProteinAC.size() + 1;
	//	t += MW.size() + 1;
	//	t += pI.size() + 1;
	//	t += Coverage.size() + 1;
	//	t += UniquePepNum.size() + 1;
	//	t += SpecNum.size() + 1;
	//	t += NonModifiedSpecNum.size() + 1;
	//	t += ModifiedSpecNum.size() + 1;
	//	t += UniqueModifiedPepNum.size() + 1;
	//	t += Samples.size() + 1;
	//	Description = iss.str().substr(t);
	//	//cout << "----------------------------------" << endl;
	//	//cout << Description << endl;
}

void ReadFromIntermediateFileProtein2(istringstream & iss, string & strOrder,
		string & Spectrum, string & Sequence, string & Score, string & Calc_M,
		string & Delta_M, string & ppm, string & Mod_Sites, string & SampleID,
		string & Engine, string & MatchedIons, string & MissCleaveNum,
		string & Rank, string & Proteins) {

	string strTmp;
	//iss >> strTmp >> strOrder >> Spectrum >> Sequence;
	getline(iss, strTmp, '\t');
	getline(iss, strOrder, '\t');
	getline(iss, Spectrum, '\t');
	getline(iss, Sequence, '\t');

	ReadFromIntermediateFromWriteSinglePeptide(iss, Score, Calc_M, Delta_M, ppm,
			Mod_Sites, SampleID, Engine, MatchedIons, MissCleaveNum, Rank,
			Proteins);
}

string TILELFormat(const string & titlex) {
	string title = titlex;

	size_t pos1 = title.find_first_of('.');
	string stra1 = title.substr(0, pos1); //scan

	title = title.substr(pos1 + 1);
	size_t pos2 = title.find_first_of('.');
	string stra2 = title.substr(0, pos2); //scan1

	title = title.substr(pos2 + 1);
	size_t pos3 = title.find_first_of('.');
	string stra3 = title.substr(0, pos3); //scan2

	title = title.substr(pos3 + 1);

	char chr[1000];
	sprintf(chr, "%s%06d.%06d.%s", stra1.c_str(), atoi(stra2.c_str()),
			atoi(stra3.c_str()), title.c_str());

	string ans = chr;
	return ans;

}

void ReadDataFORTILTE(const string & strPath) {
	map<string, string> outTmp;
	ifstream fin(strPath.c_str());
	ofstream fout("out.mgf");
	string strIn;
	char ss[1000];
	pair<string, string> pairTemp;
	while (fin.getline(ss, 1000, '\n')) {
		strIn = ss;
		if (strIn[0] == 'T') {
			string t = TILELFormat(strIn);
			//fout << t << endl;
			pairTemp.first = t;

			pairTemp.second += t + '\n';
		} else if (strIn[0] == 'E') {
			pairTemp.second += strIn + '\n';
			outTmp.insert(pairTemp);
			pairTemp.second.clear();
		} else {
			pairTemp.second += strIn + '\n';
		}
	}
	for (map<string, string>::iterator it = outTmp.begin(); it != outTmp.end();
			it++) {
		fout << it->second;
	}
	fin.close();
	fout.close();
}
