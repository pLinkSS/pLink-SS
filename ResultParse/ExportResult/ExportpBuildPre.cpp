#include "ExportpBuildPre.h"

using namespace bio_analysis;

extern char SearchEngineName[4][20];

extern ostringstream osspBuildLog;

CExportpBuildPre::CExportpBuildPre()
{
}

CExportpBuildPre::~CExportpBuildPre()
{
}

void CExportpBuildPre::_Spectra_Export(const CConf & conf, const string & ExportPathAndName)
{
	string ExportSpectra_Export = ExportPathAndName + ".spectra.xls";
	FILE * fout = fopen(ExportSpectra_Export.c_str(), "w");

	fprintf(fout, "#\tSpectrum\tPeptideNum\tUniquePepNum\tSamples\tScore\tCondition\n");
	fprintf(fout,
			"*\t#,Spec\tPeptide\tSampleID\tMod_Sites\tScore\tCalc_M\tDelta_M\tppm\tProteins\n");

	string ExportSpectra = conf.m_outPutForder_Index + "IntermediateFile.spectra";
	ifstream fin(ExportSpectra.c_str());
	if (!fin.good())
	{
		CErrInfo info("ExportpBuildPre.cpp", "_Spectra_Export", "Cannot open the file: "
				+ ExportSpectra);
		throw runtime_error(info.Get());
	}

	string Spectrum, Condition, strOrder, Score, Peptide, SampleID, Mod_Sites, Proteins, strTmp,
			Engine, order, PeptideNum, UniquePepNum, Samples, Rank, Calc_M, Delta_M, ppm, strVal,
			MatchedIons, MissCleaveNum;

	GetLine(fin, strVal);
	GetLine(fin, strVal);

	while (GetLine(fin, strVal))
	{
		if (strVal[0] == '*')
		{
			istringstream iss(strVal);

			ReadFromIntermediateFileSpectra2(iss, strOrder, Peptide, Score, Calc_M, Delta_M, ppm,
					Mod_Sites, SampleID, Engine, MatchedIons, MissCleaveNum, Rank, Proteins);

			size_t pos1 = Peptide.find_first_of('.');
			size_t pos2 = Peptide.find_last_of('.');
			Peptide = Peptide.substr(pos1 + 1, pos2 - pos1 - 1);

			fprintf(fout, "*\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", strOrder.c_str(),
					Peptide.c_str(), SampleID.c_str(), Mod_Sites.c_str(), Score.c_str(),
					Calc_M.c_str(), Delta_M.c_str(), ppm.c_str(), Proteins.c_str());
		}

		else
		{
			istringstream iss(strVal);

			ReadFromIntermediateFileSpectra1(iss, order, Spectrum, PeptideNum, UniquePepNum,
					Samples, Score, Condition);

			fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", order.c_str(), Spectrum.c_str(),
					PeptideNum.c_str(), UniquePepNum.c_str(), Samples.c_str(), Score.c_str(),
					Condition.c_str());
		}
	}

	fin.close();
	fclose(fout);
}

void CExportpBuildPre::_Peptide_Export(const CConf & conf, const string & ExportPathAndName)
{
	string ExportPeptide_Export = ExportPathAndName + ".peptide.xls";
	FILE * fout = fopen(ExportPeptide_Export.c_str(), "w");

	fprintf(fout, "Order\tSequence\n");
	fprintf(fout,
			"\tOrder\tSpectrum\tScore\tCalc_M\tDelta\tppm\tModification\tSampleID\tEngine\tRank\tProteins\n");

	string ExportPeptide = conf.m_outPutForder_Index + "IntermediateFile.peptide";
	ifstream fin(ExportPeptide.c_str());
	if (!fin.good())
	{
		CErrInfo info("ExportpBuildPre.cpp", "_Peptide_Export", "Cannot open the file: "
				+ ExportPeptide);
		throw runtime_error(info.Get());
	}

	string Spectrum, SQ, Spectra, Condition, strOrder, Score, Peptide, SampleID, Mod_Sites,
			MatchedIons, MissCleaveNum, Proteins, strTmp, Engine, order, PeptideNum, UniquePepNum,
			Samples, Rank, Calc_M, Delta_M, ppm, strVal;

	GetLine(fin, strVal);
	GetLine(fin, strVal);

	while (GetLine(fin, strVal))
	{
		if (strVal[0] == '*')
		{
			istringstream iss(strVal);

			ReadFromIntermediateFilePeptide2(iss, strOrder, Spectra, Score, Calc_M, Delta_M, ppm,
					Mod_Sites, SampleID, Engine, MatchedIons, MissCleaveNum, Rank, Proteins);

			fprintf(fout, "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", strOrder.c_str(),
					Spectra.c_str(), Score.c_str(), Calc_M.c_str(), Delta_M.c_str(), ppm.c_str(),
					Mod_Sites.c_str(), SampleID.c_str(), Engine.c_str(), Rank.c_str(),
					Proteins.c_str());
		}

		else
		{
			istringstream iss(strVal);

			ReadFromIntermediateFilePeptide1(iss, order, SQ, Spectra, Samples, Score, Proteins);

			fprintf(fout, "%s\t%s\n", order.c_str(), SQ.c_str());
		}
	}

	fin.close();
	fclose(fout);
}

//void CExportpBuildPre::_Protein_Export(const CConf & conf, const string & ExportPathAndName)
//{
//	string ExportProtein_Export = ExportPathAndName + ".protein.xls";
//
//	FILE * fout = fopen(ExportProtein_Export.c_str(), "w");
//	string ExportProtein = conf.m_outPutForder_Index +"IntermediateFile.protein";
//
//	ifstream fin(ExportProtein.c_str());
//	if (!fin.good())
//	{
//		CErrInfo info("ExportForInterface.cpp", "_Protein_Export", "Cannot open the file: "
//				+ ExportProtein);
//		throw runtime_error(info.Get());
//	}
//
//	string Spectrum, SQ, Spectra, Condition, strOrder, Score, Peptide, SampleID, Mod_Sites,
//			MatchedIons, MissCleaveNum, Proteins, strTmp, Engine, order, PeptideNum, UniquePepNum,
//			Samples, Rank, Calc_M, Delta_M, ppm, strVal, Sequence, ProteinAC, MW, pI, Coverage,
//			SpecNum, NonModifiedSpecNum, ModifiedSpecNum, UniqueModifiedPepNum, Description;
//
//	GetLine(fin, strVal);
//	GetLine(fin, strVal);
//
//	bool SameSubSet = false;
//
//	while (GetLine(fin, strVal))
//	{
//		if (strVal.substr(0, 5) == "[END]")
//			break;
//
//		else if (strVal[0] == '*')
//		{
//			if (SameSubSet == true)
//				continue;
//
//			istringstream iss(strVal);
//
//			ReadFromIntermediateFileProtein2(iss, strOrder, Spectrum, Sequence, Score, Calc_M,
//					Delta_M, ppm, Mod_Sites, SampleID, Engine, MatchedIons, MissCleaveNum, Rank,
//					Proteins);
//
//			size_t pos1 = Sequence.find_first_of('.');
//			size_t pos2 = Sequence.find_last_of('.');
//			Sequence = Sequence.substr(pos1 + 1, pos2 - pos1 - 1);
//
//			fprintf(fout, "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
//					strOrder.c_str(), Spectrum.c_str(), Sequence.c_str(), Score.c_str(),
//					Calc_M.c_str(), Delta_M.c_str(), ppm.c_str(), Mod_Sites.c_str(),
//					SampleID.c_str(), Engine.c_str(), MatchedIons.c_str(), MissCleaveNum.c_str(),
//					Rank.c_str(), Proteins.c_str());
//		}
//
//		else if (strVal.substr(0, 8) == "[1Title]")
//		{
//
//			fprintf(
//					fout,
//					"Order\tProteinAC\tMW\tpI\tUniquePepNum\tCoverage(%%)\tSpecNum\tNon-ModifiedSpecNum\tModifiedSpecNum\tUniqueModifiedPepNum\tDescription\n");
//
//		}
//		else if (strVal.substr(0, 8) == "[2Title]")
//		{
//
//			fprintf(
//					fout,
//					"\tOrder\tSpectrum\tSequence\tScore\tCalc_M\tDelta_M\tppm\tModification\tSample\tEngine\tMatchedIons\tMissCleaveNum\tRank\tProteins\n");
//
//			SameSubSet = false;
//		}
//		else if (strVal.substr(0, 8) == "[3Title]")
//		{
//
//			fprintf(
//					fout,
//					"\tProteinAC\tMW\tpI\tUniquePepNum\tCoverage(%%)\tSpecNum\tNon-ModifiedSpecNum\tModifiedSpecNum\tUniqueModifiedPepNum\tDescription\n");
//
//			SameSubSet = true;
//		}
//		else
//		{
//			istringstream iss(strVal);
//
//			ReadFromIntermediateFileProtein1(iss, order, ProteinAC, MW, pI, Coverage, UniquePepNum,
//					SpecNum, NonModifiedSpecNum, ModifiedSpecNum, UniqueModifiedPepNum, Samples,
//					Description);
//
//			fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", order.c_str(),
//					ProteinAC.c_str(), MW.c_str(), pI.c_str(), UniquePepNum.c_str(),
//					Coverage.c_str(), SpecNum.c_str(), NonModifiedSpecNum.c_str(),
//					ModifiedSpecNum.c_str(), UniqueModifiedPepNum.c_str(), Description.c_str());
//		}
//	}
//	fprintf(fout, "\n\n");
//
//	while (GetLine(fin, strVal))
//	{
//		if (strVal.substr(0, 14) == "Total_proteins")
//			fprintf(fout, "\n");
//
//		fprintf(fout, "%s\n", strVal.c_str());
//
//		if (strVal.substr(0, 13) == "---summary---")
//			fprintf(fout, "\n");
//	}
//	fin.close();
//	fclose(fout);
//
//	string ExportProteinSimple_Export = ExportPathAndName + ".proteins_simple.xls";
//	fout = fopen(ExportProteinSimple_Export.c_str(), "w");
//
//	string ExportProteinSimple = conf.m_outPutForder_Index + conf.m_OutPutFile + ".proteins_simple";
//	ifstream fin2(ExportProteinSimple.c_str());
//	if (!fin2.good())
//	{
//		CErrInfo info("ExportForInterface.cpp", "_Protein_java", "Cannot open the file: "
//				+ ExportProteinSimple);
//		throw runtime_error(info.Get());
//	}
//
//	while (GetLine(fin2, strVal))
//	{
//		fprintf(fout, "%s\n", strVal.c_str());
//	}
//
//	fin2.close();
//	fclose(fout);
//}

bool CMP_Pair(const pair<string, int> & a, const pair<string, int> & b)
{
	if (a.second > b.second)
		return true;
	return false;
}

void CExportpBuildPre::_Protein_Export(const CConf & conf, const string & ExportPathAndName)
{
	string ExportProtein = conf.m_outPutForder_Index + "IntermediateFile.protein";

	map<string, int> map_path_Protein;
	vector<pair<string, int> > ProteinVector;
	ifstream fin(ExportProtein.c_str(), ios_base::binary);
	if (!fin.good())
	{
		CErrInfo info("ExportpBuildPre.cpp", "_Protein_Export", "Cannot open the file: "
				+ ExportProtein);
		throw runtime_error(info.Get());
	}

	string Spectrum, SQ, Spectra, Condition, strOrder, Score, Peptide, SampleID, Mod_Sites,
			MatchedIons, MissCleaveNum, Proteins, strTmp, Engine, order, PeptideNum, UniquePepNum,
			Samples, Rank, Calc_M, Delta_M, ppm, strVal, Sequence, ProteinAC, MW, pI, Coverage,
			SpecNum, NonModifiedSpecNum, ModifiedSpecNum, UniqueModifiedPepNum, Description;

	//ofstream out2("ou.out");
	while (GetLine(fin, strVal))
	{
		if (strVal.substr(0, 5) == "[END]")
			break;
		if (strVal.substr(0, 8) == "[1Title]")
		{
			int move = fin.tellg();
			GetLine(fin, strVal);
			istringstream iss(strVal);
			ReadFromIntermediateFileProtein1(iss, order, ProteinAC, MW, pI, Coverage, UniquePepNum,
					SpecNum, NonModifiedSpecNum, ModifiedSpecNum, UniqueModifiedPepNum, Samples,
					Description);
			map_path_Protein[ProteinAC] = move;
			pair<string, int> pairTemp;
			pairTemp.first = ProteinAC;
			pairTemp.second = atoi(UniquePepNum.c_str());
			ProteinVector.push_back(pairTemp);
			//out2 << map_path_Protein[ProteinAC] << endl;
			//out2 << ProteinAC<< "" <<  UniquePepNum << endl;
		}
	}
	//out2.close();

	fin.close();
	sort(ProteinVector.begin(), ProteinVector.end(), CMP_Pair);
	//	ofstream out("t.out");
	//	for (size_t t = 0; t < ProteinVector.size(); t++)
	//	{
	//		cout << ProteinVector[t].first << " " << ProteinVector[t].second << endl;
	//		out << ProteinVector[t].first << " " << ProteinVector[t].second << endl;
	//	}
	//	out.close();

	ifstream fin2(ExportProtein.c_str(), ios_base::binary);
	if (!fin2.good())
	{
		CErrInfo info("ExportpBuildPre.cpp", "_Protein_Export", "Cannot open the file: "
				+ ExportProtein);
		throw runtime_error(info.Get());
	}

	string ExportProtein_Export = ExportPathAndName + ".protein.xls";
	FILE * fout = fopen(ExportProtein_Export.c_str(), "w");

	int nPepNum = 0;
	for (size_t t = 0; t < ProteinVector.size(); t++)
	{
		fin2.seekg(map_path_Protein[ProteinVector[t].first], ios_base::beg);
		fprintf(
				fout,
				"Order\tProteinAC\tMW\tpI\tUniquePepNum\tCoverage(%%)\tSpecNum\tNon-ModifiedSpecNum\tModifiedSpecNum\tUniqueModifiedPepNum\tDescription\n");
		bool SameSubSet = false;
		while (GetLine(fin2, strVal))
		{
			if (strVal.substr(0, 5) == "[END]")
				break;
			else if (strVal.substr(0, 8) == "[1Title]")
				break;
			else if (strVal[0] == '*')
			{
				if (SameSubSet == true)
					continue;

				istringstream iss(strVal);

				ReadFromIntermediateFileProtein2(iss, strOrder, Spectrum, Sequence, Score, Calc_M,
						Delta_M, ppm, Mod_Sites, SampleID, Engine, MatchedIons, MissCleaveNum,
						Rank, Proteins);

				size_t pos1 = Sequence.find_first_of('.');
				size_t pos2 = Sequence.find_last_of('.');
				Sequence = Sequence.substr(pos1 + 1, pos2 - pos1 - 1);

				fprintf(fout, "\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
						nPepNum++, Spectrum.c_str(), Sequence.c_str(), Score.c_str(),
						Calc_M.c_str(), Delta_M.c_str(), ppm.c_str(), Mod_Sites.c_str(),
						SampleID.c_str(), Engine.c_str(), MatchedIons.c_str(),
						MissCleaveNum.c_str(), Rank.c_str(), Proteins.c_str());
			}
			else if (strVal.substr(0, 8) == "[2Title]")
			{

				fprintf(
						fout,
						"\tOrder\tSpectrum\tSequence\tScore\tCalc_M\tDelta_M\tppm\tModification\tSample\tEngine\tMatchedIons\tMissCleaveNum\tRank\tProteins\n");

				SameSubSet = false;
				nPepNum = 1;
			}
			else if (strVal.substr(0, 8) == "[3Title]")
			{

				fprintf(
						fout,
						"\tProteinAC\tMW\tpI\tUniquePepNum\tCoverage(%%)\tSpecNum\tNon-ModifiedSpecNum\tModifiedSpecNum\tUniqueModifiedPepNum\tDescription\n");

				SameSubSet = true;
			}
			else
			{
				istringstream iss(strVal);

				ReadFromIntermediateFileProtein1(iss, order, ProteinAC, MW, pI, Coverage,
						UniquePepNum, SpecNum, NonModifiedSpecNum, ModifiedSpecNum,
						UniqueModifiedPepNum, Samples, Description);

				if (SameSubSet == true)
					fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", order.c_str(),
							ProteinAC.c_str(), MW.c_str(), pI.c_str(), UniquePepNum.c_str(),
							Coverage.c_str(), SpecNum.c_str(), NonModifiedSpecNum.c_str(),
							ModifiedSpecNum.c_str(), UniqueModifiedPepNum.c_str(),
							Description.c_str());
				else
					fprintf(fout, ""PRI_SIZE_T"\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", t + 1,
							ProteinAC.c_str(), MW.c_str(), pI.c_str(), UniquePepNum.c_str(),
							Coverage.c_str(), SpecNum.c_str(), NonModifiedSpecNum.c_str(),
							ModifiedSpecNum.c_str(), UniqueModifiedPepNum.c_str(),
							Description.c_str());
			}
		}
	}
	/*****************************************************************************************
	 *************************************************************************************/
	fprintf(fout, "\n\n");
	while (GetLine(fin2, strVal))
	{
		if (strVal.substr(0, 14) == "---summary---")
		{

			fprintf(fout, "%s\n", strVal.c_str());
			fprintf(fout, "\n");
			break;
		}
	}
	while (GetLine(fin2, strVal))
	{
		if (strVal.substr(0, 14) == "Total_proteins")
			fprintf(fout, "\n");
		fprintf(fout, "%s\n", strVal.c_str());
	}
	fin2.close();
	fclose(fout);

	/***************************************************************************************/
	string ExportProteinSimple_Export = ExportPathAndName + ".proteins_sample.xls";
	fout = fopen(ExportProteinSimple_Export.c_str(), "w");

	string ExportProteinSimple = conf.m_outPutForder_Index + conf.m_OutPutFile + ".proteins_sample";
	ifstream fin3(ExportProteinSimple.c_str());
	if (!fin3.good())
	{
		CErrInfo info("ExportpBuildPre.cpp", "_Protein_java", "Cannot open the file: "
				+ ExportProteinSimple);
		throw runtime_error(info.Get());
	}

	while (GetLine(fin3, strVal))
	{
		fprintf(fout, "%s\n", strVal.c_str());
	}

	fin3.close();
	fclose(fout);
}

void CExportpBuildPre::_Protein_Export_NoUniPepSort(const CConf & conf,
		const string & ExportPathAndName)
{
	string ExportProtein_Export = ExportPathAndName + ".protein.xls";

	FILE * fout = fopen(ExportProtein_Export.c_str(), "w");
	string ExportProtein = conf.m_outPutForder_Index + "IntermediateFile.protein";

	ifstream fin(ExportProtein.c_str());
	if (!fin.good())
	{
		CErrInfo info("ExportpBuildPre.cpp", "_Protein_Export", "Cannot open the file: "
				+ ExportProtein);
		throw runtime_error(info.Get());
	}

	string Spectrum, SQ, Spectra, Condition, strOrder, Score, Peptide, SampleID, Mod_Sites,
			MatchedIons, MissCleaveNum, Proteins, strTmp, Engine, order, PeptideNum, UniquePepNum,
			Samples, Rank, Calc_M, Delta_M, ppm, strVal, Sequence, ProteinAC, MW, pI, Coverage,
			SpecNum, NonModifiedSpecNum, ModifiedSpecNum, UniqueModifiedPepNum, Description;

	GetLine(fin, strVal);
	GetLine(fin, strVal);

	bool SameSubSet = false;

	while (GetLine(fin, strVal))
	{
		if (strVal.substr(0, 5) == "[END]")
			break;

		else if (strVal[0] == '*')
		{
			if (SameSubSet == true)
				continue;

			istringstream iss(strVal);

			ReadFromIntermediateFileProtein2(iss, strOrder, Spectrum, Sequence, Score, Calc_M,
					Delta_M, ppm, Mod_Sites, SampleID, Engine, MatchedIons, MissCleaveNum, Rank,
					Proteins);

			size_t pos1 = Sequence.find_first_of('.');
			size_t pos2 = Sequence.find_last_of('.');
			Sequence = Sequence.substr(pos1 + 1, pos2 - pos1 - 1);

			fprintf(fout, "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					strOrder.c_str(), Spectrum.c_str(), Sequence.c_str(), Score.c_str(),
					Calc_M.c_str(), Delta_M.c_str(), ppm.c_str(), Mod_Sites.c_str(),
					SampleID.c_str(), Engine.c_str(), MatchedIons.c_str(), MissCleaveNum.c_str(),
					Rank.c_str(), Proteins.c_str());
		}

		else if (strVal.substr(0, 8) == "[1Title]")
		{

			fprintf(
					fout,
					"Order\tProteinAC\tMW\tpI\tUniquePepNum\tCoverage(%%)\tSpecNum\tNon-ModifiedSpecNum\tModifiedSpecNum\tUniqueModifiedPepNum\tDescription\n");

		}
		else if (strVal.substr(0, 8) == "[2Title]")
		{

			fprintf(
					fout,
					"\tOrder\tSpectrum\tSequence\tScore\tCalc_M\tDelta_M\tppm\tModification\tSample\tEngine\tMatchedIons\tMissCleaveNum\tRank\tProteins\n");

			SameSubSet = false;
		}
		else if (strVal.substr(0, 8) == "[3Title]")
		{

			fprintf(
					fout,
					"\tProteinAC\tMW\tpI\tUniquePepNum\tCoverage(%%)\tSpecNum\tNon-ModifiedSpecNum\tModifiedSpecNum\tUniqueModifiedPepNum\tDescription\n");

			SameSubSet = true;
		}
		else
		{
			istringstream iss(strVal);

			ReadFromIntermediateFileProtein1(iss, order, ProteinAC, MW, pI, Coverage, UniquePepNum,
					SpecNum, NonModifiedSpecNum, ModifiedSpecNum, UniqueModifiedPepNum, Samples,
					Description);

			fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", order.c_str(),
					ProteinAC.c_str(), MW.c_str(), pI.c_str(), UniquePepNum.c_str(),
					Coverage.c_str(), SpecNum.c_str(), NonModifiedSpecNum.c_str(),
					ModifiedSpecNum.c_str(), UniqueModifiedPepNum.c_str(), Description.c_str());
		}
	}
	fprintf(fout, "\n\n");

	while (GetLine(fin, strVal))
	{
		if (strVal.substr(0, 14) == "Total_proteins")
			fprintf(fout, "\n");

		fprintf(fout, "%s\n", strVal.c_str());

		if (strVal.substr(0, 13) == "---summary---")
			fprintf(fout, "\n");
	}
	fin.close();
	fclose(fout);
	/***************************************************************************************/
	string ExportProteinSimple_Export = ExportPathAndName + ".proteins_sample.xls";
	fout = fopen(ExportProteinSimple_Export.c_str(), "w");

	string ExportProteinSimple = conf.m_outPutForder_Index + conf.m_OutPutFile + ".proteins_sample";
	ifstream fin2(ExportProteinSimple.c_str());
	if (!fin2.good())
	{
		CErrInfo info("ExportpBuildPre.cpp", "_Protein_java", "Cannot open the file: "
				+ ExportProteinSimple);
		throw runtime_error(info.Get());
	}

	while (GetLine(fin2, strVal))
	{
		fprintf(fout, "%s\n", strVal.c_str());
	}

	fin2.close();
	fclose(fout);
}

void CExportpBuildPre::_plabel_Export(const CConf & conf, const string & ExportPathAndName)
{
	string ExportpLabel_Name = conf.m_outPutForder_pLabel + "plabel_Name.txt";
	ifstream fin(ExportpLabel_Name.c_str());
	if (!fin.good())
	{
		CErrInfo info("ExportpBuildPre.cpp", "_plabel_Export", "Cannot open the file: "
				+ ExportpLabel_Name);
		throw runtime_error(info.Get());
	}

	size_t pos = ExportPathAndName.find_last_of('/');
	string pLabelForder = ExportPathAndName.substr(0, pos) + "/pLabel/";

#ifdef WIN32
	mkdir(pLabelForder.c_str());
#else			
	mkdir(pLabelForder.c_str(),0775);
#endif	

	string strpLabelName = pLabelForder + ExportPathAndName.substr(pos + 1);

	int num = 1;
	string strVal;

	while (GetLine(fin, strVal))
	{
		ifstream fin2(strVal.c_str());
		if (!fin2.good())
		{
			CErrInfo info("ExportpBuildPre.cpp", "_plabel_Export", "Cannot open the file: "
					+ strVal);
			throw runtime_error(info.Get());
		}

		string ExportppLabel_Export = strpLabelName;
		char chr[100];
		sprintf(chr, "_%d.plabel", num++);
		ExportppLabel_Export += chr;
		FILE * fout = fopen(ExportppLabel_Export.c_str(), "w");

		string strTemp;

		while (GetLine(fin2, strTemp))
		{
			fprintf(fout, "%s\n", strTemp.c_str());
		}

		fin2.close();
		fclose(fout);
	}
	fin.close();
}

void CExportpBuildPre::_pQuant_Export(const CConf & conf, const string & ExportPathAndName)
{
	string ExportpQuant_Export = ExportPathAndName + "_00.pquant";
	FILE * fout = fopen(ExportpQuant_Export.c_str(), "w");

	string ExportQuant = conf.m_outPutForder_Index + "IntermediateFile.pquant";
	ifstream fin(ExportQuant.c_str());
	if (!fin.good())
	{
		CErrInfo info("ExportpBuildPre.cpp", "_pQuant_Export", "Cannot open the file: "
				+ ExportQuant);
		throw runtime_error(info.Get());
	}

	string strVal;
	int cnt = 0;
	int fcnt = 1;
	while (GetLine(fin, strVal))
	{
		if (strVal[0] == 'P' && cnt > 2000)
		{
			fclose(fout);
			char chr[100];
			sprintf(chr, "_%02d", fcnt++);
			ExportpQuant_Export = ExportPathAndName + chr + ".pquant";
			fout = fopen(ExportpQuant_Export.c_str(), "w");
			cnt = 0;

			fprintf(fout, "Order\tProteinAC\tCoverage\tPeptideNum\tUniquePepNum\n");
			fprintf(
					fout,
					"\tOrder\tUniqueProNum\tSpectrum\tCharge\tSequence\tScore\tCalc_M\tDelta\tModification\tSampleID\tEngine\tRank\tProteins\n");
		}
		fprintf(fout, "%s\n", strVal.c_str());
		cnt++;
	}

	fin.close();
	fclose(fout);
}

void CExportpBuildPre::_FASTA_Export(const CConf & conf, const string & ExportPathAndName)
{
	string ExportFASTA_Export = ExportPathAndName + ".fasta";
	FILE * fout = fopen(ExportFASTA_Export.c_str(), "w");

	string ExportFASTA = conf.m_outPutForder_Index + conf.m_OutPutFile + ".fasta";
	ifstream fin(ExportFASTA.c_str());
	if (!fin.good())
	{
		CErrInfo info("ExportpBuildPre.cpp", "_pQuant_Export", "Cannot open the file: "
				+ ExportFASTA);
		throw runtime_error(info.Get());
	}

	string strVal;

	while (GetLine(fin, strVal))
	{
		fprintf(fout, "%s\n", strVal.c_str());
	}

	fin.close();
	fclose(fout);
}

void CExportpBuildPre::Export(const CConf & conf, const string & ExportPathAndName,
		const string & ExportFileType)
{
	//cout << "pBuild Export..." << endl;
	osspBuildLog << "pBuild Export..." << endl;

	if (ExportFileType.find('0') != string::npos)
	{
		_Spectra_Export(conf, ExportPathAndName);
	}

	if (ExportFileType.find('1') != string::npos)
	{
		_Peptide_Export(conf, ExportPathAndName);
	}

	if (ExportFileType.find('2') != string::npos)
	{
		_Protein_Export(conf, ExportPathAndName);
	}

	if (ExportFileType.find('3') != string::npos)
	{
		_pQuant_Export(conf, ExportPathAndName);
	}

	if (ExportFileType.find('4') != string::npos)
	{
		_plabel_Export(conf, ExportPathAndName);
	}
	if (ExportFileType.find('5') != string::npos)
	{
		_FASTA_Export(conf, ExportPathAndName);
	}

}
