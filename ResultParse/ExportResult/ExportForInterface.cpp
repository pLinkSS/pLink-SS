#include "ExportForInterface.h"

using namespace bio_analysis;

extern char SearchEngineName[4][20];

extern ostringstream osspBuildLog;

CExportForInterface::CExportForInterface()
{
}

CExportForInterface::~CExportForInterface()
{
}

void CExportForInterface::_Spectra_java(const CConf & conf)
{
	conf.m_OutPutFile + ".spectra.txt";
	string ExportSpectra_java = conf.m_outPutForder_Java + conf.m_OutPutFile + ".spectra.txt";

	FILE * fout = fopen(ExportSpectra_java.c_str(), "w");
	fprintf(fout, "#\tSpectrum\tUni_Pet\tSamples\tScore\tCondition\n");
	fprintf(fout,
			"*\t#,Spec\tPeptide\tSampleID\tMod_Sites\tScore\tCalc_M\tDelta_M\tppm\tProteins\n");

	string ExportSpectra = conf.m_outPutForder_Index + "IntermediateFile.spectra";
	ifstream fin(ExportSpectra.c_str());
	if (!fin.good())
	{
		CErrInfo info("ExportForInterface.cpp", "_Spectra_java", "Cannot open the file: "
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

			fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\n", order.c_str(), Spectrum.c_str(),
					UniquePepNum.c_str(), Samples.c_str(), Score.c_str(), Condition.c_str());
		}
	}
	fin.close();
	fclose(fout);
}

void CExportForInterface::_Peptide_java(const CConf & conf)
{
	conf.m_OutPutFile + ".peptide.txt";

	string ExportPeptide_java = conf.m_outPutForder_Java + conf.m_OutPutFile + ".peptide.txt";
	FILE * fout = fopen(ExportPeptide_java.c_str(), "w");

	fprintf(fout, "#\tSQ\tSpectra\tSamples\tScore\tProtein\n");
	fprintf(fout, "*\t#,Pep\tSampleID\tSpectra\tMod_Sites\tScore\tCalc_M\tDelta_M\tppm\tProteins\n");

	string ExportPeptide = conf.m_outPutForder_Index + "IntermediateFile.peptide";
	ifstream fin(ExportPeptide.c_str());
	if (!fin.good())
	{
		CErrInfo info("ExportForInterface.cpp", "_Peptide_java", "Cannot open the file: "
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

			fprintf(fout, "*\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", strOrder.c_str(),
					SampleID.c_str(), Spectra.c_str(), Mod_Sites.c_str(), Score.c_str(),
					Calc_M.c_str(), Delta_M.c_str(), ppm.c_str(), Proteins.c_str());
		}

		else
		{
			istringstream iss(strVal);

			ReadFromIntermediateFilePeptide1(iss, order, SQ, Spectra, Samples, Score, Proteins);

			fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\n", order.c_str(), SQ.c_str(), Spectra.c_str(),
					Samples.c_str(), Score.c_str(), Proteins.c_str());
		}
	}

	fin.close();
	fclose(fout);
}

void CExportForInterface::_Protein_java(const CConf & conf)
{
	string ExportProtein_java = conf.m_outPutForder_Java + conf.m_OutPutFile + ".protein.txt";
	FILE * fout = fopen(ExportProtein_java.c_str(), "w");

	fprintf(fout, "#\tProteinAC\tCoverage\tSpectra\tUni_pep\tSamples\n");
	fprintf(fout,
			"*\t#\tPeptide\tSpectrum\tSampleID\tModification\tScore\tCalc_M\tDelta_M\tppm\tProteins\n");

	string ExportProtein = conf.m_outPutForder_Index + "IntermediateFile.protein";
	ifstream fin(ExportProtein.c_str());
	if (!fin.good())
	{
		CErrInfo info("ExportForInterface.cpp", "_Protein_java", "Cannot open the file: "
				+ ExportProtein);
		throw runtime_error(info.Get());
	}

	string Spectrum, SQ, Spectra, Condition, strOrder, Score, Peptide, SampleID, Mod_Sites,
			MatchedIons, MissCleaveNum, Proteins, strTmp, Engine, order, PeptideNum, UniquePepNum,
			Samples, Rank, Calc_M, Delta_M, ppm, strVal, Sequence, ProteinAC, MW, pI, Coverage,
			SpecNum, NonModifiedSpecNum, ModifiedSpecNum, UniqueModifiedPepNum, Description;

	GetLine(fin, strVal);
	GetLine(fin, strVal);

	while (GetLine(fin, strVal))
	{
		if (strVal.substr(0, 5) == "[END]")
			break;
		if (strVal[0] == '*')
		{
			istringstream iss(strVal);

			ReadFromIntermediateFileProtein2(iss, strOrder, Spectrum, Sequence, Score, Calc_M,
					Delta_M, ppm, Mod_Sites, SampleID, Engine, MatchedIons, MissCleaveNum, Rank,
					Proteins);

			size_t pos1 = Sequence.find_first_of('.');
			size_t pos2 = Sequence.find_last_of('.');
			Sequence = Sequence.substr(pos1 + 1, pos2 - pos1 - 1);

			fprintf(fout, "*\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", strOrder.c_str(),
					Sequence.c_str(), Spectrum.c_str(), SampleID.c_str(), Mod_Sites.c_str(),
					Score.c_str(), Calc_M.c_str(), Delta_M.c_str(), ppm.c_str(), Proteins.c_str());
		}

		else if (strVal.substr(0, 8) == "[1Title]" || strVal.substr(0, 8) == "[2Title]"
				|| strVal.substr(0, 8) == "[3Title]")
			continue;

		else
		{
			istringstream iss(strVal);

			ReadFromIntermediateFileProtein1(iss, order, ProteinAC, MW, pI, Coverage, UniquePepNum,
					SpecNum, NonModifiedSpecNum, ModifiedSpecNum, UniqueModifiedPepNum, Samples,
					Description);

			fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\n", order.c_str(), ProteinAC.c_str(),
					Coverage.c_str(), SpecNum.c_str(), UniquePepNum.c_str(), Samples.c_str());
		}
	}

	fin.close();
	fclose(fout);
}

void CExportForInterface::Export(const CConf & conf)
{
	//cout << "Interface Data..." << endl;
	osspBuildLog << "Interface Data..." << endl;

	_Spectra_java(conf);
	_Peptide_java(conf);
	_Protein_java(conf);
}
