#include "ExportHTML.h"

using namespace bio_analysis;

//extern ofstream flog;
extern ostringstream osspBuildLog;
CExportHTML::CExportHTML()
{
}

CExportHTML::~CExportHTML()
{
}

void CExportHTML::_GetLink(FILE * fout, const CConf & conf)
{
	fprintf(fout, "<br>\n<a href=\"%s\" target_blank>%s</a>\n<br>\n", conf.m_OutPutForder.c_str(),
			conf.m_OutPutName.c_str());
	string ExportSpectraHTML = conf.m_OutPutFile + ".spectra.htm";
	string ExportSummaryHTML = conf.m_OutPutFile + ".summary.htm";
	string ExportPeptideHTML = conf.m_OutPutFile + ".peptide.htm";
	string ExportProteinHTML = conf.m_OutPutFile + ".protein.htm";
	fprintf(fout, "<br>\n");
	fprintf(fout, "<table border='1' align=\"center\">\n");
	fprintf(fout, "<tr>\n");
	fprintf(fout, "<td><a href=\"%s\">Summary Result</a>\n", ExportSummaryHTML.c_str());
	fprintf(fout, "<td><a href=\"%s\">Spectra Result</a>\n", ExportSpectraHTML.c_str());
	fprintf(fout, "<td><a href=\"%s\">Peptide Result</a>\n", ExportPeptideHTML.c_str());
	fprintf(fout, "<td><a href=\"%s\">Protein Result</a>\n", ExportProteinHTML.c_str());
	fprintf(fout, "\n</tr>\n");
	fprintf(fout, "</table>\n");
	fprintf(fout, "<br>\n");
}

//void CExportHTML::_GetSpectraLink(FILE * fout, const CConf & conf)
//{
//	string ExportSpectraHTML = conf.m_OutPutName;
//	ExportSpectraHTML += ".spectra.htm";
//	fprintf(fout, "<td><a href=\"%s\">Spectra WebSite</a>",
//			ExportSpectraHTML.c_str());
//}
//void CExportHTML::_GetPeptideLink(FILE * fout, const CConf & conf)
//{
//	string ExportPeptideHTML = conf.m_OutPutName;
//	ExportPeptideHTML += ".peptide.htm";
//	fprintf(fout, "<td><a href=\"%s\">Peptide WebSite</a>",
//			ExportPeptideHTML.c_str());
//}
//
//void CExportHTML::_GetProteinLink(FILE * fout, const CConf & conf)
//{
//	string ExportProteinHTML = conf.m_OutPutName;
//	ExportProteinHTML += ".protein.htm";
//	fprintf(fout, "<td><a href=\"%s\">Protein WebSite</a>",
//			ExportProteinHTML.c_str());
//}

void CExportHTML::_GetProteinHead(FILE * fout, const string & strPro,
		const vector<SPECTRAINFO> & SpecPath)
{
	set<string> setUniPep;
	for (size_t t = 0; t < SpecPath.size(); t++)
		setUniPep.insert(SpecPath[t].second.m_strFirstPeptide);
	fprintf(fout, "\t<FONT color=red>%s</FONT>\t%lf\t"PRI_SIZE_T"\t"PRI_SIZE_T"\t\n", strPro.c_str(), 0.0,
			SpecPath.size(), setUniPep.size());
}

void CExportHTML::_Protein(const CConf & conf, const IntegratedInfo & HaveToInterData)
{
	string ExportProteinHTML = conf.m_OutPutName;
	ExportProteinHTML += ".protein.htm";
	FILE * fout = fopen(ExportProteinHTML.c_str(), "w");
	CMatchSpectraInfo SpectraTemp;
	int order = 1;
	fprintf(fout, "<HTML><HEAD><TITLE>%s-pBuild2.0</TITLE>", conf.m_OutPutFile.c_str());
	fprintf(fout, "<META http-equiv=Content-Type content=\"text/html; charset=GB2312\"></HEAD>\n");
	fprintf(fout, "<BODY bgcolor=#F0F8FF>\n");
	_GetLink(fout, conf);
	map<string, set<string> > m_mapSameSet = HaveToInterData.m_mapSameSet;
	map<string, set<string> > m_mapSubSet = HaveToInterData.m_mapSubSet;
	map<string, vector<SPECTRAINFO> > protein_peptide = HaveToInterData.protein_peptide;
	for (map<string, vector<SPECTRAINFO> >::const_iterator it =
			HaveToInterData.protein_peptide.begin(); it != HaveToInterData.protein_peptide.end(); it++)
	{
		if (HaveToInterData.m_setInvalidProteins.find(it->first)
				!= HaveToInterData.m_setInvalidProteins.end())
			continue;
		if (HaveToInterData.m_deleteProtein.find(it->first)
				!= HaveToInterData.m_deleteProtein.end())
			continue;
		set<string> setUniPep;
		vector<CMatchSpectraInfo> vSpectra;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);
			setUniPep.insert(SpectraTemp.m_vPeptides[0].m_strSQ);
			vSpectra.push_back(SpectraTemp);
		}
		fprintf(fout, "<TABLE><TBODY>\n");
		fprintf(fout, "<TR>\n");
		fprintf(fout, "<TD>%d</TD>", order++);
		fprintf(fout, "<TD><FONT color=red>%s</FONT></TD>", it->first.c_str());
		fprintf(fout, "<TD>%lf</TD>\n<TD>"PRI_SIZE_T"</TD>\n<TD>"PRI_SIZE_T"</TD>\n", 0.0, it->second.size(),
				setUniPep.size());
		fprintf(fout, "</TR>\n");
		fprintf(fout, "</TBODY></TABLE>\n");
		fprintf(fout, "<PRE>\n");
		for (set<string>::const_iterator t = m_mapSameSet[it->first].begin(); t
				!= m_mapSameSet[it->first].end(); t++)
		{
			fprintf(fout, "SameSet");
			_GetProteinHead(fout, *t, protein_peptide[*t]);
		}
		for (set<string>::const_iterator t = m_mapSubSet[it->first].begin(); t
				!= m_mapSubSet[it->first].end(); t++)
		{
			fprintf(fout, "SubSet");
			_GetProteinHead(fout, *t, protein_peptide[*t]);
		}
		for (size_t t = 0; t < it->second.size(); t++)
		{
			fprintf(fout, "\t"PRI_SIZE_T"\t", t + 1);
			fprintf(fout, "%s\t%d\t", vSpectra[t].m_strFileName.c_str(), vSpectra[t].m_nCharge);
			fprintf(fout, "%c.%s.%c\t", vSpectra[t].m_vPeptides[0].m_cPrev,
					vSpectra[t].m_vPeptides[0].m_strSQ.c_str(), vSpectra[t].m_vPeptides[0].m_cNext);
			CMatchPeptideInfo PepTemp = vSpectra[t].m_vPeptides[0];
			WriteSingePeptide(fout, PepTemp, vSpectra[t]);
		}
		fprintf(fout, "</PRE>\n");
		fprintf(fout, "<HR>\n");
	}
	nProtein_Group_Num = order - 1;
	fprintf(fout, "</BODY></HTML>\n");
	fclose(fout);
}

void CExportHTML::_Spectra(const CConf & conf, const IntegratedInfo & HaveToInterData)
{
	string ExportSpectraHTML = conf.m_OutPutName;
	ExportSpectraHTML += ".spectra.htm";
	FILE * fout = fopen(ExportSpectraHTML.c_str(), "w");
	fprintf(fout, "<HTML><HEAD><TITLE>%s-pBuild2.0</TITLE>", conf.m_OutPutFile.c_str());
	fprintf(fout, "<META http-equiv=Content-Type content=\"text/html; charset=GB2312\"></HEAD>\n");
	fprintf(fout, "<BODY bgcolor=#F0F8FF>\n");
	_GetLink(fout, conf);
	//	fprintf(fout, "Order\tTitle\tPeptideNum\tUniquePepNum\tCharge\n");
	//	fprintf(
	//			fout,
	//			"\tOrder\tSequence\tScore\tCalc_M\tDelta\tModification\tSampleID\tEngine\tRank\tProteins\n");
	int order = 1;
	bool xinjing = true;
	CMatchSpectraInfo SpectraTemp;
	for (map<string, vector<SPECTRAINFO> >::const_iterator it =
			HaveToInterData.spec_peptide.begin(); it != HaveToInterData.spec_peptide.end(); it++)
	{
		if (HaveToInterData.m_deleteSpectra.find(it->first)
				!= HaveToInterData.m_deleteSpectra.end())
			continue;
		vector<CMatchSpectraInfo> vSpectra;
		set<string> setUniPep;
		set<int> setIDs;
		for (size_t t = 0; t < it->second.size(); t++)
		{
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);
			setUniPep.insert(SpectraTemp.m_vPeptides[0].m_strSQ);
			setIDs.insert(SpectraTemp.m_nDataSetID);
			vSpectra.push_back(SpectraTemp);
		}
		fprintf(fout, "<TABLE><TBODY>\n");
		fprintf(fout, "<TR>\n");
		fprintf(fout, "<TD>%d</TD>", order++);
		//		fprintf(fout, "<TD><FONT color=red>%s</FONT> %s</TD>", SplitTitle(it->first).c_str(),
		//				it->first.c_str());
		fprintf(fout, "<TD><FONT color=red>%s</FONT></TD>", it->first.c_str());
		fprintf(fout, "<TD>"PRI_SIZE_T"</TD>\n<TD>"PRI_SIZE_T"</TD>\n<TD>%d</TD>\n", it->second.size(),
				setUniPep.size(), vSpectra[0].m_nCharge);
		fprintf(fout, "</TR>\n");
		fprintf(fout, "</TBODY></TABLE>\n");
		fprintf(fout, "<PRE>\n");
		for (size_t t = 0; t < it->second.size(); t++)
		{
			if (xinjing == true)
			{
				fprintf(fout, "#\t");
				xinjing = false;
			}
			else
			{
				fprintf(fout, "*\t");
				xinjing = true;
			}
			fprintf(fout, ""PRI_SIZE_T"\t%c.%s.%c\t", t + 1, vSpectra[t].m_vPeptides[0].m_cPrev,
					vSpectra[t].m_vPeptides[0].m_strSQ.c_str(), vSpectra[t].m_vPeptides[0].m_cNext);
			WriteSingePeptide(fout, vSpectra[t].m_vPeptides[0], vSpectra[t]);
		}
		fprintf(fout, "</PRE>\n");
		fprintf(fout, "<HR>\n");
	}
	fprintf(fout, "</BODY></HTML>\n");
	fclose(fout);
}

void CExportHTML::_Peptide(const CConf & conf, const IntegratedInfo & HaveToInterData)
{
	string ExportPeptideHTML = conf.m_OutPutName;
	ExportPeptideHTML += ".peptide.htm";
	FILE * fout = fopen(ExportPeptideHTML.c_str(), "w");
	fprintf(fout, "<HTML><HEAD><TITLE>%s-pBuild2.0</TITLE>", conf.m_OutPutFile.c_str());
	fprintf(fout, "<META http-equiv=Content-Type content=\"text/html; charset=GB2312\"></HEAD>\n");
	fprintf(fout, "<BODY bgcolor=#F0F8FF>\n");
	_GetLink(fout, conf);
	//	fprintf(fout, "Order\tSequence\n");
	//	fprintf(
	//			fout,
	//			"\tOrder\tSpectrum\tCharge\tScore\tCalc_M\tDelta\tModification\tSampleID\tEngine\tRank\tProteins\n");
	CMatchSpectraInfo SpectraTemp;
	int order = 1;
	for (map<string, vector<SPECTRAINFO> >::const_iterator it =
			HaveToInterData.peptide_spec.begin(); it != HaveToInterData.peptide_spec.end(); it++)
	{
		if (HaveToInterData.m_deletePeptide.find(it->first)
				!= HaveToInterData.m_deletePeptide.end())
			continue;
		fprintf(fout, "<TABLE><TBODY>\n");
		fprintf(fout, "<TR>\n");
		fprintf(fout, "<TD>%d</TD>", order++);
		fprintf(fout, "<TD><FONT color=red>%s</FONT></TD>", it->first.c_str());
		fprintf(fout, "</TR>\n");
		fprintf(fout, "</TBODY></TABLE>\n");
		fprintf(fout, "<PRE>\n");
		for (size_t t = 0; t < it->second.size(); t++)
		{
			ReadSingeSpectraFisrtPep(it->second[t].first, SpectraTemp);
			WritePeptide_Spec(fout, SpectraTemp, t);
		}
		fprintf(fout, "</PRE>\n");
		fprintf(fout, "<HR>\n");
	}
	fprintf(fout, "</BODY></HTML>\n");
	fclose(fout);

}

void CExportHTML::_Summary(const CConf & conf, const IntegratedInfo & HaveToInterData)
{
	string ExportSummaryHTML = conf.m_OutPutName;
	ExportSummaryHTML += ".summary.htm";
	FILE * fout = fopen(ExportSummaryHTML.c_str(), "w");
	fprintf(fout, "<HTML><HEAD><TITLE>%s-pBuild2.0</TITLE>", conf.m_OutPutFile.c_str());
	fprintf(fout, "<META http-equiv=Content-Type content=\"text/html; charset=GB2312\"></HEAD>\n");
	fprintf(fout, "<BODY bgcolor=#F0F8FF>pBuild 2.0\n");
	_GetLink(fout, conf);
	fprintf(fout, "Spectra_Num = "PRI_SIZE_T"<br>\n", HaveToInterData.spec_peptide.size()
			- HaveToInterData.m_deleteSpectra.size());
	fprintf(fout, "Peptide_Num = "PRI_SIZE_T"<br>\n", HaveToInterData.peptide_spec.size()
			- HaveToInterData.m_deletePeptide.size());
	fprintf(fout, "Protein_Num = "PRI_SIZE_T"<br>\n", HaveToInterData.protein_peptide.size()
			- HaveToInterData.m_deleteProtein.size());
	fprintf(fout, "Protein_Group_Num = "PRI_SIZE_T"\n", nProtein_Group_Num);
	fprintf(fout, "<br><br><br>\n");
	//WriteFiltrationHTML(conf.m_Filter, fout);
	WriteCConfHTML(fout, conf);
	fprintf(fout, "</BODY></HTML><br>\n");
	fclose(fout);
}
void CExportHTML::Export(const IntegratedInfo & HaveToInterData, const CConf & conf)
{
	//cout << "HTML Export..." << endl;
	osspBuildLog << "HTML Export..." << endl;
	//pBuildLog("HTML Export...");
	_Spectra(conf, HaveToInterData);
	_Peptide(conf, HaveToInterData);
	_Protein(conf, HaveToInterData);
	_Summary(conf, HaveToInterData);
}

