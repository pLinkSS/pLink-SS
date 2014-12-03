/************************************************************
 * 1.肽段，蛋白个数的重复计算
 * 2.数据的读写 ACE memoryMAP
 * 3.disciple>2时蛋白数可能有问题，可以参考测试报告。
 * 4.//todo 2009-9-13 今天测试了一个大的数据，发现
 * 后来越来越慢，然后在浩哥的指导下，测试了一个小的数据，发现内
 * 存使用增加的很快，一方面可能是内存泄漏，也可能是建索引花费了内存。
 * 5.由于spectra文件中map的第一个关键字存放的是去掉.charge.dta,
 * 导致后来在dispeptide>1的时候，比较时采用的不是一种形式，而且
 * 在spiltTitle函数中，取scan号也会改变。
 * 6.先用java写一个参数输入的界面，结果用网页的形式来显示。这样可以先让用户
 * 使用起来同时，网页的形式很好实现，也可以做的很漂亮。实现起来比较简单，以后有
 * 时间可以再做对话框的形式。
 * 7.todo 给内存和外存的形式(让用户自己选择)。
 * 8.todo 注意INDEXFILTER文件大于2G
 ********************************************************/
/*
 * 每次第一次运行新的数据时，建立索引文件，并将索引文件写到硬盘里，以便下次解析时直接调用这些索引
 * 1.重新过滤默认读写的参数文件名 Filter.REFITER
 * 2.改参数之后，新建一个参数文件，去读该参数文件，同时也可以选择自己需要的参数文件
 * 3.读已经建立索引的文件(已经解析完成的索引文件)
 * 4.重新过滤
 *********************************************************/

/*
 * 2010年6月23日
 * 这几天删除了输出的工厂模式，直接需要输出那种格式，就调用哪个函数
 * 首先，将以后可能使用到的所有结果都输出到到一个中间文件中，之后需要哪种格式的
 * 文件，直接去读取这个文件就可以，可以避免很多不必要的计算，同时，也实现了程序运
 * 行完成之后也可以导出结果的功能，之后都是要导出哪个文件，必须要在程序运行开始时
 * 就指定这样不符合人们的一般思维
 * 同时，pQuant格式文件的输出嵌入到了Protein文件的输出函数中，避免重复计算
 *********************************************************/

/*
 * 2010年6月23日
 * 每次读谱是都打开一次文件，是不是很耗时？可以把文件全部打开，以后直接去读，不用每次都打开文件
 *
 *
 *********************************************************/

/*
 Correlation between MH and HZ

 MH = M + P;
 MZ = (M + C * P) / C;

 neutral Ion
 b = [M1] + P;
 y = [M2] + H2O + P;
 MH = [M3] + H2O + P;
 [M3] = [M1] + [M2]; //todo


 * y = residue_y + H2 + OH;
 * b = residue_b + H;
 * precursor = residue_pre + H2O + H;
 * residue_pre = residue_y + residue_b;


 residue mass = sum(AA);
 molecule mass = residue mass + H20;
 MH = molecule mass + p;
 MZ = (molecule mass + p * c) / c;
 */
//#define DEBUG
#include "SampleInfor.h"

using namespace bio_analysis;

const long Maxfseek = 99999999;
const long MaxLongInt = 999999999;

extern char SearchEngineName[10][20];

ostringstream osspBuildLog;

void CSampleInfor::_OutPutCConf(CConf & conf) {
	string ParmerTemp = conf.m_OutPutName + ".pbuild";
//	cout<<"Jeremy: "<<ParmerTemp<<endl;
	WriteCConf(conf, ParmerTemp);
}

void CSampleInfor::_ReadparametersAndConditions(CConf & conf) {
	string parametersAndConditions = conf.m_outPutForder_Index
			+ "parametersAndConditions.txt";
	ifstream fin(parametersAndConditions.c_str());
	if (!fin.good()) {
		CErrInfo info("SampleInfor.cpp", "_ReadparametersAndConditions",
				"Cannot open the file: " + parametersAndConditions);
		throw runtime_error(info.Get());
	}

	//cout << conf.m_vEngineInputInfo.size() << endl;
	for (size_t t = 0; t < conf.m_vEngineInputInfo.size(); t++) {
		vector<CConditionInfo> m_vConditions;
		//这里之前写错了，只考虑了sample的个数，没有考虑sample里subItem的个数
		for (size_t j = 0; j < conf.m_vEngineInputInfo[t].m_vResultPath.size();
				j++) {
			string strTemp;
			ThrowFirstOfRemain(fin, "RawInputPath", strTemp);
			CConditionInfo cond;
			cond.m_strInputPath = strTemp;
			//			if (conf.m_vEngineInputInfo.size() <= t)
			//				break;
			m_vConditions.push_back(cond);
		}
		conf.m_vEngineInputInfo[t].m_vConditions = m_vConditions;
	}

	fin.close();
}

void CSampleInfor::_Makdir(CConf & conf) {
	conf.m_OutPutForder += GetTimePre(conf) + conf.m_OutPutFile.c_str() + cSlash;

#ifdef WIN32
	mkdir(conf.m_OutPutForder.c_str());
#else			
	mkdir(conf.m_OutPutForder.c_str(),0777);
#endif	

	conf.m_OutPutName = conf.m_OutPutForder;
	conf.m_OutPutName += conf.m_OutPutFile;
	osspBuildLog << conf.m_OutPutForder << endl;
}

void CSampleInfor::Run(const string & strParamFilePath,
		const string & strIndexFile) {
	if (strIndexFile == "-denovo") {
		//m_CDeNovo.DenovoRun(strParamFilePath);
		return;
	}
	vector<OneDATASET> AllDataSetParser;
	CConf conf;
	start = clock();
	ReadCConf(conf, strParamFilePath);
	conf.m_bDataAnalysis.clear();

	if (strIndexFile == "dataanalysisforhfchen" || strIndexFile == "-denovo") {
		conf.m_bDataAnalysis = "dataanalysisforhfchen";
		conf.m_ReParse = 1;
	}

	osspBuildLog
			<< "------------------------------------------------------------"
			<< endl;
	osspBuildLog << "pBuild 2.0.0.5" << endl;
	osspBuildLog << "Release date: April 2, 2011" << endl;
	osspBuildLog << endl;
	conf.time_now = time(NULL);
	osspBuildLog << asctime(localtime(&conf.time_now));

	cout << endl;
	cout << "##############################" << endl;
	cout << "#" << endl;
	cout << "#  pBuild 2.0.0.5" << endl;
	cout << "#  Release date: April 2, 2011" << endl;
	cout << "#" << endl;
	cout << "##############################" << endl;

	cout << endl;
	cout << "Reading parameters ..." << endl;
	cout << "	Path: " << strParamFilePath << endl;

	_Makdir(conf);
	_OutPutCConf(conf);

	conf.m_outPutForder_Index = conf.m_OutPutName + "_pBuildINDEX" + cSlash;
	conf.m_outPutForder_pLabel = conf.m_OutPutName + "_pLabel" + cSlash;
	conf.m_outPutForder_Java = conf.m_OutPutName + "_Java" + cSlash;
	conf.m_outPutForder_Charge = conf.m_outPutForder_Java + "Charge" + cSlash;
	conf.m_OutPutForder_ErrorFigure = conf.m_outPutForder_Java + "ErrorFigure" + cSlash;

	if (strIndexFile == "export") {
		/*export the result to the folder you select*/
		CExportpBuildPre ExportTemp;
		string ExportPathAndName = conf.strExportPath;
		size_t pos = ExportPathAndName.find_last_of(cSlash);

		osspBuildLog << "Export Result to " << ExportPathAndName.substr(0, pos)
				<< endl;

		if (!CheckFolderExist(ExportPathAndName.substr(0, pos))) {
			CErrInfo info("SampleInfor.cpp", "Run",
					"Cannot open the folder: " + ExportPathAndName);
			throw runtime_error(info.Get());
		}

		ExportTemp.Export(conf, ExportPathAndName, conf.ExportFileType);

		end = clock();
		osspBuildLog << "Time = " << (double) (end - start) / CLOCKS_PER_SEC
				<< "s" << endl;
		osspBuildLog << "Run Over" << endl;
		pBuildLog(osspBuildLog.str());

		return;
	}

	if (strIndexFile == "delete") {
		/*delete some jurk file*/
		string strIndexFile = conf.m_outPutForder_Index;
		osspBuildLog << "Delete Folder : " << conf.m_outPutForder_Index << endl;
		RemoveFolderOrFile(strIndexFile);
		return;
	}

	RemoveFolderOrFile(conf.m_outPutForder_Charge);

#ifdef WIN32
	mkdir(conf.m_outPutForder_Java.c_str());
	mkdir(conf.m_outPutForder_Charge.c_str());
	mkdir(conf.m_OutPutForder_ErrorFigure.c_str());
#else			
	mkdir(conf.m_outPutForder_Java.c_str(),0777);
	mkdir(conf.m_outPutForder_Charge.c_str(),0777);
	mkdir(conf.m_OutPutForder_ErrorFigure.c_str(),0777);
#endif	

	bool IndexExist = false;
	string strParserSimple = conf.m_outPutForder_Index + conf.m_OutPutFile
			+ ".samples";
	if (conf.m_ReParse == 0 && CheckFolderExist(conf.m_outPutForder_Index)
			&& CheckSimpleSame(conf, strParserSimple)) {

		string strIndexFile2 = conf.m_outPutForder_Index + "INDEX.REFILTER";

		if (CheckFileExist(strIndexFile2)) {
#ifdef DEBUG
			cout << "IndexExist" << endl;
#endif
			CheckIndexFile(strIndexFile2);
			osspBuildLog << "Read Index File..." << endl;
			ReadIndexFile(strIndexFile2, AllDataSetParser);
			_ReadparametersAndConditions(conf);
			IndexExist = true;
		}
	}

	CDataAnalysis m_DataAnalysis(conf);
	if (IndexExist == false) {
#ifdef DEBUG
		cout << "Index Not Exist" << endl;
#endif

#ifdef WIN32
		mkdir(conf.m_outPutForder_Index.c_str());
#else			
		mkdir(conf.m_outPutForder_Index.c_str(),0775);
#endif	

		_Parser(AllDataSetParser, conf, m_DataAnalysis);
		string strIndexFileOther = conf.m_outPutForder_Index + "INDEX.REFILTER"; //todo 注意这里大于2G
		string strSimpleFile = conf.m_outPutForder_Index + conf.m_OutPutFile
				+ ".samples";

		WriteParserSimple(strSimpleFile, conf);
		WriteIndexFile(strIndexFileOther, AllDataSetParser);
	}

	if (conf.m_Filter.m_UseFilter == false) {
		/*if you do not want to filter the result, then merge the file to one*/
		CMergeTopFind CMergeTopFindT;
		CMergeTopFindT.pFindMerge_Export(AllDataSetParser, conf);
		AllDataSetParser.clear();
	} else {
		/*filter and export the result*/
		vector<OneDATASET> AllDataSetFilter;
		_Filter(AllDataSetParser, AllDataSetFilter, conf);
#ifdef DEBUG
		cout << "Filter Over" << endl;
#endif
		AllDataSetParser.clear();
		IntegratedInfo HaveToInterData;
		m_DataIntegrate.DataTOIntegrate(AllDataSetFilter, conf,
				HaveToInterData);
		//AllDataSetFilter.clear();
#ifdef DEBUG
		cout << "DataIntegrate Over" << endl;
#endif
		CExportToIntermediateFile m_ExportToIntermediateFile;
		m_ExportToIntermediateFile.ExportToIntermediateFile(HaveToInterData,
				conf);
		cout << "Writing ..." << endl;
#ifdef DEBUG
		cout << "ExportToIntermediateFile Over" << endl;
#endif
#ifdef  dataanalysisforhfchen
		if (conf.m_bDataAnalysis == "dataanalysisforhfchen")
		{
			m_DataAnalysis.OutPutEachResult(HaveToInterData, conf); //todo only for hfchen
			//m_DataAnalysis.getwhy(HaveToInterData, conf, 0);
			//m_DataAnalysis.getwhy(HaveToInterData, conf, 1);
			//m_DataAnalysis.getwhy(HaveToInterData, conf, 2);
			//m_DataAnalysis.getwhy(HaveToInterData, conf, 3);
			//m_DataAnalysis.getwhy(HaveToInterData, conf, 4);
			//m_DataAnalysis.getVenn();
			//m_DataAnalysis.DataAnalysis(HaveToInterData, conf)
		}
#endif
		if (conf.m_bCrossLink == 1 || conf.m_bCrossLink == 2) {
			CMergeTopFind CMergeTopFindT;
			CMergeTopFindT.pFindMerge_Export(AllDataSetFilter, conf);
		}
		for (size_t t = 0; t < conf.ExportFormat.size(); t++) {
			//cout << conf.ExportFormat[t] << endl;
			if (conf.ExportFormat[t] == ST_HTMLFormat) {
				CExportHTML ExportTemp;
				string ExportPathAndName = conf.m_OutPutForder
						+ conf.m_OutPutFile;
				ExportTemp.Export(HaveToInterData, conf);
			}
			if (conf.ExportFormat[t] == ST_PFINDFormat
					&& conf.m_bCrossLink == 0) {
				CMergeTopFind CMergeTopFindT;
				CMergeTopFindT.pFindMerge_Export(AllDataSetFilter, conf);
			}
		}

		if (strIndexFile == "-denovo") {
			map<string, vector<string> > Spec_Pep;
			for (map<string, vector<SPECTRAINFO> >::iterator it =
					HaveToInterData.spec_peptide.begin();
					it != HaveToInterData.spec_peptide.end(); it++) {
				for (size_t t = 0; t < it->second.size(); t++) {
					Spec_Pep[it->first].push_back(
							it->second[t].second.m_strFirstPeptide);
				}
			}
			//m_CDeNovo.DenovoRun(conf.m_strPNovoParam);
		}
		AllDataSetFilter.clear();
		HaveToInterData.clear();
		_ExportResult(conf);
#ifdef DEBUG
		cout << "ExportResult Over" << endl;
#endif
	}
	_Close(conf);

#ifdef DEBUG
	cout << "Run Over" << endl;
#endif

#ifdef  dataanalysisforhfchen
	if (conf.m_bDataAnalysis == "dataanalysisforhfchen")
	m_DataAnalysis.DataAnalysis(conf); // TODO onlyforhfchen
#endif

	osspBuildLog << "Run Over" << endl;
	pBuildLog(osspBuildLog.str());
#ifdef DEBUG
	cout << "END" << endl;
#endif
	cout << "done!" << endl;
}

void CSampleInfor::_Close(const CConf & conf) {
	string PathTemp = conf.m_OutPutName;
	string time = PathTemp + ".time.txt";
	FILE * ftime = fopen(time.c_str(), "w");

	end = clock();
	fprintf(ftime, "Time = %lfs\n", (double) (end - start) / CLOCKS_PER_SEC);
	osspBuildLog << "Time = " << (double) (end - start) / CLOCKS_PER_SEC << "s"
			<< endl;
	fclose(ftime);
}

void CSampleInfor::_Parser(vector<OneDATASET> & AllDataSetParser, CConf & conf,
		CDataAnalysis & m_DataAnalysis) {
	osspBuildLog << "Parse..." << endl;
	cout << "Parse..." << endl;

	string parametersAndConditions = conf.m_outPutForder_Index
			+ "parametersAndConditions.txt";
	FILE * fparametersAndConditionsout = fopen(parametersAndConditions.c_str(),
			"w");

	int index = 1;
	string IndexFile = conf.m_outPutForder_Index + "index1.INDEX";
	FILE * fout = fopen(IndexFile.c_str(), "w");
	SPECTRA_PATH SpecPathTemp;

	SpecPathTemp.first = IndexFile;
	fprintf(fparametersAndConditionsout, "Total_Samples="PRI_SIZE_T"\n",
			conf.m_vEngineInputInfo.size());
	for (size_t i = 0; i < conf.m_vEngineInputInfo.size(); i++) {
		fprintf(fparametersAndConditionsout, "[Sample"PRI_SIZE_T"]\n", i + 1);
		osspBuildLog << "   Sample " << i + 1 << endl;

		cout << "	Running Sample " << i + 1 << " ..." << endl;

		OneDATASET pairTemp;
		map<int, vector<SPECTRAINFO> > DataSetTemp;
		CResultParserFactory parsefac;
		CResultParser * p = parsefac.GetParser(
				conf.m_vEngineInputInfo[i].m_Type);
		char OutPathTemp_Da[PATH_MAX] = { 0 };
		sprintf(OutPathTemp_Da, "Sample"PRI_SIZE_T"_Da.errorfigure", i + 1);
		string strFerrorfigure_Da = conf.m_OutPutForder_ErrorFigure
				+ OutPathTemp_Da;
		FILE * ferrorfigure_Da = fopen(strFerrorfigure_Da.c_str(), "w");
		char OutPathTemp_ppm[PATH_MAX] = { 0 };
		sprintf(OutPathTemp_ppm, "Sample"PRI_SIZE_T"_ppm.errorfigure", i + 1);
		string strFerrorfigure_ppm = conf.m_OutPutForder_ErrorFigure
				+ OutPathTemp_ppm;
		FILE * ferrorfigure_ppm = fopen(strFerrorfigure_ppm.c_str(), "w");
		fprintf(fparametersAndConditionsout, "EngineType=%s\n",
				SearchEngineName[conf.m_vEngineInputInfo[i].m_Type]);
		fprintf(fparametersAndConditionsout, "SubItems="PRI_SIZE_T"\n",
				conf.m_vEngineInputInfo[i].m_vResultPath.size());
		for (size_t j = 0; j < conf.m_vEngineInputInfo[i].m_vResultPath.size();
				j++) {
			fprintf(fparametersAndConditionsout, "SubItem"PRI_SIZE_T"=%s\n", j + 1,
					conf.m_vEngineInputInfo[i].m_vResultPath[j].c_str());

			if (0 == conf.m_vEngineInputInfo[i].m_Valid[j]) {
				fprintf(fparametersAndConditionsout, "TotalSpec=0\n");
				fprintf(fparametersAndConditionsout, "RawInputPath=NULL\n");
				continue;
			}
			CParseredDataSetInfo pbres;
			osspBuildLog << "   " << j + 1 << "/"
					<< conf.m_vEngineInputInfo[i].m_vResultPath.size() << "   -"
					<< conf.m_vEngineInputInfo[i].m_vResultPath[j] << " ";
#ifdef DEBUG
			cout << "   " << j + 1 << "/" << conf.m_vEngineInputInfo[i].m_vResultPath.size()
			<< "   -" << conf.m_vEngineInputInfo[i].m_vResultPath[j] << " ";
#endif
			p->Parse(conf.m_vEngineInputInfo[i].m_vResultPath[j], pbres);
			conf.m_vEngineInputInfo[i].m_vConditions.push_back(
					pbres.m_vConditions[0]);

			for (map<int, vector<CMatchSpectraInfo> >::iterator it =
					pbres.m_mapResultData.begin();
					it != pbres.m_mapResultData.end(); it++) {
				for (size_t t = 0; t < pbres.m_mapResultData[it->first].size();
						t++) {
					pbres.m_mapResultData[it->first][t].m_nFileID = j + 1;
					pbres.m_mapResultData[it->first][t].m_nDataSetID = i + 1;
					pbres.m_mapResultData[it->first][t].EnginType =
							conf.m_vEngineInputInfo[i].m_Type;

					SpecPathTemp.second = ftell(fout);
					long move = SpecPathTemp.second;

					OutPutSpectra(fout, pbres.m_mapResultData[it->first][t],
							ferrorfigure_Da, ferrorfigure_ppm, conf);
#ifdef  dataanalysisforhfchen
					if (conf.m_bDataAnalysis == "dataanalysisforhfchen")
					m_DataAnalysis.CompleteSpectraInfo.push_back(
							pbres.m_mapResultData[it->first][t]);
#endif

					GetIndex(DataSetTemp, pbres.m_mapResultData[it->first][t],
							SpecPathTemp);

					if (move >= Maxfseek) {
						char chrTemp[PATH_MAX] = { '\0' };
						sprintf(chrTemp, "index%d.INDEX", ++index);
						string IndexFile = conf.m_outPutForder_Index + chrTemp;

						fprintf(fout, "[END]\n");
						fclose(fout);
						fout = fopen(IndexFile.c_str(), "w");
						SpecPathTemp.first = IndexFile;
					}
				}
			}

			pbres.m_vConditions[0].EngineType =
					conf.m_vEngineInputInfo[i].m_Type;
			pairTemp.first = pbres.m_vConditions[0];

			fprintf(fparametersAndConditionsout, "TotalSpec=%d\n",
					pbres.m_vConditions[0].m_nTotalSpec);
			fprintf(fparametersAndConditionsout, "RawInputPath=%s\n",
					pbres.m_vConditions[0].m_strInputPath.c_str());
		}

		fclose(ferrorfigure_Da);
		fclose(ferrorfigure_ppm);

		pairTemp.second = DataSetTemp;
		AllDataSetParser.push_back(pairTemp);

		delete p;
	}
	fclose(fparametersAndConditionsout);
	fprintf(fout, "[END]\n");
	fclose(fout);
}

void CSampleInfor::_ModSiteFilter(OneDATASET & ResultDataSet,
		const CConf & conf) {
	OneDATASET NewResultDataSet;

	for (map<int, vector<SPECTRAINFO> >::const_iterator it =
			ResultDataSet.second.begin(); it != ResultDataSet.second.end();
			it++) {
		for (size_t t = 0; t < it->second.size(); t++) {
			CMatchSpectraInfo SpectraTemp;
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);

			if (SpectraTemp.m_vPeptides.size() == 0)
				continue;

			for (size_t k = 0; k < SpectraTemp.m_vPeptides[0].m_vMod.size();
					k++) {
				char chrtmp;
				if (SpectraTemp.m_vPeptides[0].m_vMod[k].m_tPos == 0)
					chrtmp =
							SpectraTemp.m_vPeptides[0].m_strSQ[SpectraTemp.m_vPeptides[0].m_vMod[k].m_tPos];
				else
					chrtmp =
							SpectraTemp.m_vPeptides[0].m_strSQ[SpectraTemp.m_vPeptides[0].m_vMod[k].m_tPos
									- 1];

				if (conf.m_Filter.m_strModSites.find(chrtmp) != string::npos) {
					NewResultDataSet.second[it->first].push_back(it->second[t]);
					break;
				}
			}
		}
	}
	ResultDataSet.second.clear();
	ResultDataSet.second = NewResultDataSet.second;
}

void CSampleInfor::_TerminalFilter(OneDATASET & ResultDataSet,
		const CConf & conf) {
	OneDATASET NewResultDataSet;

	for (map<int, vector<SPECTRAINFO> >::const_iterator it =
			ResultDataSet.second.begin(); it != ResultDataSet.second.end();
			it++) {
		for (size_t t = 0; t < it->second.size(); t++) {
			CMatchSpectraInfo SpectraTemp;
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);

			if (SpectraTemp.m_vPeptides.size() == 0)
				continue;

			size_t pos1, pos2;
			pos1 = 0;
			pos2 = 0;
			if (conf.m_Filter.m_strCTerminal.size() > 0)
				pos1 = conf.m_Filter.m_strCTerminal.find(
						SpectraTemp.m_vPeptides[0].m_cPrev);

			if (conf.m_Filter.m_strNTerminal.size() > 0)
				pos2 = conf.m_Filter.m_strNTerminal.find(
						SpectraTemp.m_vPeptides[0].m_cNext);

			if (pos1 != string::npos && pos2 != string::npos)
				NewResultDataSet.second[it->first].push_back(it->second[t]);

		}
	}
	ResultDataSet.second.clear();
	ResultDataSet.second = NewResultDataSet.second;
}

void CSampleInfor::_CrossLinkFilter(vector<OneDATASET> & AllDataSetParser,
		const CConf & conf) {
	for (size_t t = 0; t < AllDataSetParser.size(); t++) {

		OneDATASET DataSetTemp;
		map<int, vector<SPECTRAINFO> > mapCrossLink;
		for (map<int, vector<SPECTRAINFO> >::const_iterator it =
				AllDataSetParser[t].second.begin();
				it != AllDataSetParser[t].second.end(); it++) {
			for (size_t i = 0; i < it->second.size(); i++) {
				if (it->second[i].second.m_strFirstPeptide == "")
					continue;

				int nReagentType = GetReagentType(
						it->second[i].second.m_strFirstPeptide);
				int nCrossType = GetCrossType(
						it->second[i].second.m_strFirstPeptide);
				if (ISinVector(it->first, conf.m_Filter.m_vChargeState) == false
						|| ISinVector(nCrossType, conf.m_vLinkType) == false
						|| ISinVector(nReagentType, conf.m_vLinkID) == false)
					continue;
				int type = GetCrossLinkType(it->first, nCrossType,
						nReagentType);

				mapCrossLink[type].push_back(it->second[i]);
			}

		}
		AllDataSetParser[t].second = mapCrossLink;
#ifdef DEBUG
		for (map<int, vector<SPECTRAINFO> >::const_iterator it = AllDataSetParser[t].second.begin(); it
				!= AllDataSetParser[t].second.end(); it++)
		{
			cout << it->first << " " << it->second.size() << endl;
		}
#endif
	}
}

void CSampleInfor::_DeCrossLinkFilter(vector<OneDATASET> & AllDataSetParser,
		const CConf & conf) {
	for (size_t t = 0; t < AllDataSetParser.size(); t++) {
		OneDATASET DataSetTemp;
		map<int, vector<SPECTRAINFO> > mapCrossLink;
		for (map<int, vector<SPECTRAINFO> >::const_iterator it =
				AllDataSetParser[t].second.begin();
				it != AllDataSetParser[t].second.end(); it++) {
			for (size_t i = 0; i < it->second.size(); i++) {
				if (it->second[i].second.m_strFirstPeptide == "")
					continue;
				int type = DeCrossLinkType(it->first);
				mapCrossLink[type].push_back(it->second[i]);
			}

		}
		AllDataSetParser[t].second = mapCrossLink;
#ifdef DEBUG
		for (map<int, vector<SPECTRAINFO> >::const_iterator it = AllDataSetParser[t].second.begin(); it
				!= AllDataSetParser[t].second.end(); it++)
		{
			cout << it->first << " " << it->second.size() << endl;
		}
#endif
	}
}

void CSampleInfor::_Filter(vector<OneDATASET> & AllDataSetParser,
		vector<OneDATASET> & AllDataSetFilter, const CConf & conf) {
	osspBuildLog << "Filter..." << endl;

	cout << "Filter ..." << endl;

	if (conf.m_Filter.m_UseFilter == false) {

		CMergeTopFind CMergeTopFindT;
		CMergeTopFindT.pFindMerge_Export(AllDataSetParser, conf);
		return;
	}

	if (conf.m_Filter.m_strModSites.size() > 0) //这里大于0之前写成了大于1
			{

		for (size_t t = 0; t < AllDataSetParser.size(); t++) {
			_ModSiteFilter(AllDataSetParser[t], conf);
		}
	}

	if (conf.m_Filter.m_strCTerminal.size() > 0
			|| conf.m_Filter.m_strNTerminal.size() > 0) {

		for (size_t t = 0; t < AllDataSetParser.size(); t++) {
			_TerminalFilter(AllDataSetParser[t], conf);
		}
	}

	if (conf.m_bCrossLink == 1) {
		_CrossLinkFilter(AllDataSetParser, conf);
	}

	for (size_t t = 0; t < AllDataSetParser.size(); t++) {
		cout << "	Running Sample " << t + 1 << " ..." << endl;
		OneDATASET pairTemp;
		map<int, vector<SPECTRAINFO> > DataSetTemp;

		for (map<int, vector<SPECTRAINFO> >::const_iterator it =
				AllDataSetParser[t].second.begin();
				it != AllDataSetParser[t].second.end(); it++) {
			vector<SPECTRAINFO> SpectraInfoTemp;
			pair<int, int> pairSimpleCharge;

			pairSimpleCharge.first = t + 1;
			pairSimpleCharge.second = it->first;
//			cout << "it->second" << it->second.size() << endl;
			m_FilterUtility.Filter(it->second, SpectraInfoTemp, conf,
					AllDataSetParser[t].first.EngineType, pairSimpleCharge);
			//todo 这里用到了condition里的变量
			DataSetTemp[it->first] = SpectraInfoTemp;
		}

		pairTemp.first = AllDataSetParser[t].first;
		pairTemp.second = DataSetTemp;

		AllDataSetFilter.push_back(pairTemp);
	}

	if (conf.m_bCrossLink == 1) {
		_DeCrossLinkFilter(AllDataSetFilter, conf);
	}
}

void CSampleInfor::_ExportResult(const CConf & conf) {
	for (size_t t = 0; t < conf.ExportFormat.size(); t++) {
		if (conf.ExportFormat[t] == ST_INTERFACEEXPORT) {
			CExportForInterface ExportTemp;
			ExportTemp.Export(conf);
		} else if (conf.ExportFormat[t] == ST_PBULIDPREFormat) {
			CExportpBuildPre ExportTemp;
			string ExportPathAndName = conf.m_OutPutForder + conf.m_OutPutFile;
			ExportTemp.Export(conf, ExportPathAndName, conf.ExportFileType);
		}
	}
}
