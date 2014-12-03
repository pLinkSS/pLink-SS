#!/usr/bin/env python
#-*-coding:utf8 -*-
import os
import sys
import time
import string
import fnmatch

import Path

def _GetFastaInfo(file_path):
	fasta_info = []
	file = open(file_path, 'r')
	str = file.readline()
	while 1:
		str = file.readline()
		if not str:
			break
		fasta_info.append(str.split(','))
	file.close()

	return fasta_info

def _GetIndexInfo(file_path):
	Index_info = []
	file = open(file_path, 'r')
	str = file.readline()
	while 1:
		str = file.readline()
		if not str:
			break
		Index_info.append(str.split(','))
	file.close()

	return Index_info

def _GetSpecFDRInfo(file_path):
	FDR_info = []
	file = open(file_path, 'r')
	while 1:
		str = file.readline()
		if not str:
			break
		[k, v]=str.split('=')
		FDR_info.append(v)
	file.close()

	return FDR_info

def _GetCompareInfo(file_path):

	file = open(file_path, 'r')
	str = file.readline()
	spec_compare_info = str.split(',')[1:]
	str = file.readline()
	pep_compare_info = str.split(',')[1:]
	str = file.readline()
	pro_compare_info = str.split(',')[1:]

	return spec_compare_info, pep_compare_info, pro_compare_info

def _WriteHTMLHead(file, title):
	file.write('<html>\n')
	file.write('\n')
	file.write('<head>\n')
	file.write('<meta http-equiv=Content-Type content=\"text/html;>\n')
	file.write('<meta name=Generator content=\"pFind Studio 2.6\">\n')
	file.write('<title>pFind Studio MS/MS Analysis Report - %s</title>\n' % title)
	file.write('<style>\n')
	file.write('  body {font-family:Times New Roman,serif;font-size:16.0pt;}\n')
	file.write('  .title1 {text-align:center;font-size:32.0pt;font-weight:bold}\n')
	file.write('  .version1 {text-align:center;font-size:12.0pt;font-style:italic;}\n')
	file.write('  .section1 {text-align:left;font-size:24.0pt;font-weight:bold;}\n')
	file.write('  .figure_legend1 {text-align:center;font-size:16.0pt;}\n')
	file.write('  .table_legend1 {text-align:center;font-size:16.0pt;}\n')
	file.write('  .table1 {border-collapse:collapse;border:none;font-size:16.0pt;}\n')
	file.write('  .table_head1 {border:solid black 1.0pt;background:#8DB3E2;padding:0cm 5.4pt 0cm 5.4pt;text-align:center;color:white;}\n')
	file.write('  .table_data1 {border:solid black 1.0pt; padding:0cm 5.4pt 0cm 5.4pt;text-align:center;}\n')
	file.write('  .Reference1 {font-style:italic;}\n')
	file.write('</style>\n')
	file.write('</head>\n')
	file.write('\n')

def WriteGeneral(filepath, initparam, sample_list, report_path):
	sample_num = len(sample_list)
	xlinktype = 'inter'
	# get information
	fasta_info_file = 'fastainfo.csv'
	fasta_info = _GetFastaInfo(os.path.join(report_path, fasta_info_file))

	index_info_file = 'indexinfo.csv'
	index_info = _GetIndexInfo(os.path.join(report_path, index_info_file))

	compare_info_file='compare.csv'
	if sample_num > 1:
		[spec_compare_info, pep_compare_info, pro_compare_info] = _GetCompareInfo(os.path.join(report_path, compare_info_file))

	raw_num_list=[]
	fdr_info_list =[]
	id_num_list=[]
	for i in range(0, sample_num):
		raw_num_list.append(len(sample_list[i].mgf_path))
		spectra_title = initparam['sample%d.spectra.title' % (i+1)]
		spec_fdr_info_file = os.path.join(report_path, 'sample%d' % (i+1), 'spec_u1_'+spectra_title+'_'+xlinktype+'.statistics.txt')
		spec_fdr_Info = _GetSpecFDRInfo(spec_fdr_info_file)
		fdr_info_list.append(spec_fdr_Info)
		pfind_id_num=string.atoi(spec_fdr_Info[0])
		id_num_list.append(pfind_id_num)

	# draw figures and tables
	nTable = 0
	file = open(filepath, 'w')
	_WriteHTMLHead(file, 'General Information')
	file.write('<body>\n')
	file.write('<p id=\'top\' class=\'title1\'>pLink MS/MS Analysis Report</p>\n')
	file.write('\n')
	file.write('<p class=\'version1\'>%s automatically created by pLink 1.0</p>\n' % time.strftime('%Y-%m-%d',time.localtime(time.time())))
	file.write('<p>&nbsp;</p>\n')

	file.write('<p class=section1>1. Database</p>\n')
	file.write('\n')
	file.write('<p>&nbsp;</p>\n')
	file.write('\n')

	nTable = nTable+1
	file.write('<p class=\'table_legend1\'>Table %d Protein database information</p>\n' % nTable)
	file.write('\n')
	file.write('<div align=center>\n')
	file.write('\n')
	file.write('<table class=\'table1\'>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Fasta file</td>\n')
	file.write('  <td class=\'table_head1\'># Proteins</td>\n')
	file.write('  <td class=\'table_head1\'>File size (MB)</td>\n')
	file.write('  <td class=\'table_head1\'>Create time</td>\n')
	file.write('  <td class=\'table_head1\'>Last modified time</td>\n')
	file.write(' </tr>\n')
	for L in fasta_info:
		file.write(' <tr>\n')
		file.write('  <td class=\'table_data1\'>%s</td>\n' % L[0])
		file.write('  <td class=\'table_data1\'>%s</td>\n' % L[1])
		file.write('  <td class=\'table_data1\'>%s</td>\n' % L[2])
		file.write('  <td class=\'table_data1\'>%s</td>\n' % L[3])
		file.write('  <td class=\'table_data1\'>%s</td>\n' % L[4])
		file.write(' </tr>\n')
	file.write('</table>\n')
	file.write('\n')
	file.write('<a href=\'%s\'>Excel File ...</a>\n' % fasta_info_file)
	file.write('</div>\n')
	file.write('\n')
	file.write('<p>&nbsp;</p>\n')
	file.write('\n')

	nTable = nTable+1
	file.write('<p class=\'table_legend1\'>Table %d pLink index information (From pIndex)</p>\n' % nTable)
	file.write('\n')
	file.write('<div align=center>\n')
	file.write('\n')
	file.write('<table class=\'table1\'>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Enzyme</td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % index_info[0][0])
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Max missed site</td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % index_info[0][1])
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Cleavage way</td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % index_info[0][2])
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Min peptide length</td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % index_info[0][3])
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Max peptide length</td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % index_info[0][4])
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Min mass</td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % index_info[0][5])
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Max mass</td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % index_info[0][6])
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'># Redundant peptides</td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % index_info[0][7])
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'># Unique peptides</td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % index_info[0][8])
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Create time</td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % index_info[0][9])
	file.write(' </tr>\n')
	file.write('</table>\n')
	file.write('\n')
	file.write('<a href=\'%s\'>Excel File ...</a>\n' % index_info_file)
	file.write('</div>\n')
	file.write('<p align=center><a href=\'#top\'>return to top</a></p>\n')
	file.write('\n')
	
	file.write('<div align=center>\n')
	nTable = nTable+1

	file.write('<p class=section1>2. Identification Results (From pLink and pBuild) </p>\n')
	file.write('\n')
	file.write('<p>&nbsp;</p>\n')
	file.write('\n')
	
	file.write('<p class=\'table_legend1\'>Table %d Idenfication Results (From pFind and pBuild)</p>\n' % nTable)
	file.write('\n')
	file.write('<table class=\'table1\'>\n')
	file.write(' <tr>\n')
	file.write('   <td class=\'table_head1\' rowspan=2>Sample ID</td>\n')
	file.write('   <td class=\'table_head1\' rowspan=2>#Raw</td>\n')
	file.write('   <td class=\'table_head1\' rowspan=2>#Total Spectra</td>\n')
	file.write('   <td class=\'table_head1\' colspan=4>pFind results </td>\n')
	file.write('   <td class=\'table_head1\'  rowspan=2>Detail Information</td>\n')
	file.write('  </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\' >#Spectra</td>\n')
	file.write('  <td class=\'table_head1\'>#Peptides</td>\n')
	file.write('  <td class=\'table_head1\'>#Proteins</td>\n')
	file.write('  <td class=\'table_head1\'>#Groups</td>\n')
	file.write(' </tr>\n')
	for td in range(0, sample_num):
		file.write(' <tr>\n')
		file.write('  <td class=\'table_data1\'>%d</td>\n' % (td+1))
		file.write('  <td class=\'table_data1\'>%d</td>\n' % raw_num_list[td])
		file.write('  <td class=\'table_data1\'>%d</td>\n' % id_num_list[td])
		file.write('  <td class=\'table_data1\'>%s</td>\n' % fdr_info_list[td][0])
		file.write('  <td class=\'table_data1\'>%s</td>\n' % fdr_info_list[td][1])
		file.write('  <td class=\'table_data1\'>%s</td>\n' % fdr_info_list[td][2])
		file.write('  <td class=\'table_data1\'>%s</td>\n' % fdr_info_list[td][3])
		url = './sample%d/report.html' % (td+1)
		file.write('  <td class=\'table_data1\'><a href=\'%s\'>Click to See Section 3.%d...</a></td>\n' % (url, (td+1)))
		file.write(' </tr>\n')
	file.write('</table>\n')
	file.write('</div>\n')
	file.write('<p align=center><a href=\'#top\'>return to top</a></p>\n')

	if sample_num>1:
		file.write('\n')
		file.write('<p class=section1>3. Comparison</p>\n')
		file.write('\n')

		file.write('<div align=center>\n')
		nTable = nTable+1
		file.write('<p class=\'table_legend1\'>Table %d Compare the Union and Intersection (From pBuild)</p>\n' % nTable)
		file.write('\n')
		file.write('<table class=\'table1\'>\n')
		file.write(' <tr>\n')
		for td in range(0, sample_num):
			file.write('  <td class=\'table_head1\'>Sample %d</td>\n' % (td+1))
		file.write('  <td class=\'table_head1\'>#Spectra</td>\n')
		file.write('  <td class=\'table_head1\'>#peptide</td>\n')
		file.write('  <td class=\'table_head1\'>#protein</td>\n')
		file.write(' </tr>\n')
		line_num = 2**sample_num-1
		for tr in range(0, line_num):
			file.write(' <tr>\n')
			v = tr + 1
			for s in  range(0, sample_num):
				yes = v % 2
				v = v / 2
				if 1==yes:
					file.write('  <td class=\'table_data1\'>Yes</td>\n')
				else:
					file.write('  <td class=\'table_data1\'>No</td>\n')
			file.write('  <td class=\'table_data1\'>%s</td>\n' % spec_compare_info[tr])
			file.write('  <td class=\'table_data1\'>%s</td>\n' % pep_compare_info[tr])
			file.write('  <td class=\'table_data1\'>%s</td>\n' % pro_compare_info[tr])
			file.write(' </tr>\n')
		file.write('</table>\n')
		file.write('</div>\n')
		file.write('\n')
	file.write('\n')
	file.write('<p align=center><a href=\'#top\'>return to top</a></p>\n')
	file.write('\n')

	file.write('</body>\n')
	file.write('\n')
	file.write('</html>\n')
	file.close()
	return

def WriteSample(filepath, SampleNo, report_path, sample_title, spectra_list, RunningTimeImg, initparam):
	# get information
	pair_path = os.path.join(report_path, 'pair')
	img_path = os.path.join(report_path, 'img')
	xlinktype = 'inter'

	spec_fdr_info_file = 'spec_u1_' + sample_title+'_'+xlinktype+'.statistics.txt'
	spec_u1_fdr_Info = _GetSpecFDRInfo(os.path.join(report_path, spec_fdr_info_file))

	# draw figures and tables
	nFig = 0
	nTable = 0
	file = open(filepath, 'w')
	title = '2.%d Detail Information of Sample %d' % (SampleNo,SampleNo)
	_WriteHTMLHead(file, title)
	file.write('<body>\n')
	file.write('<p id=\'top\' class=\'title1\'>%s</p>\n' % title)

	file.write('<p class=section1>2.%d.1 Data</p>\n' % SampleNo)
	file.write('\n')
	file.write('<p>&nbsp;</p>\n')
	file.write('\n')

	nTable = nTable+1
	file.write('<p class=\'table_legend1\'>Table 2.%d.%d Database search results (Frome pFind and pBuild)</p>\n' % (SampleNo,nTable))
	file.write('\n')
	file.write('<div align=center>\n')
	file.write('\n')
	file.write('<table class=\'table1\'>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>&nbsp;</td>\n')
	file.write('  <td class=\'table_head1\'># Spectra</td>\n')
	file.write('  <td class=\'table_head1\'># Peptides</td>\n')
	file.write('  <td class=\'table_head1\'># U1 proteins (groups)</td>\n')
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'> Spectral level</br>(FDR&lt;=5%) </td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % spec_u1_fdr_Info[0])
	file.write('  <td class=\'table_data1\'>%s</td>\n' % spec_u1_fdr_Info[1])
	file.write('  <td class=\'table_data1\'>%s (%s) </td>\n' % (spec_u1_fdr_Info[2],spec_u1_fdr_Info[3]))
	file.write(' </tr>\n')
	file.write('</table>\n')
	file.write('\n')
	file.write('</div>\n')
	file.write('\n')

	file.write('<p>&nbsp;</p>\n')
	
	nTable = nTable+1
	file.write('<p class=\'table_legend1\'>Figure 2.%d.%d FDR</p>\n' % (SampleNo,nTable))
	file.write('\n')
	file.write('<div align=center>\n')
	file.write('\n')
	file.write('<p id=\'fdr\' class=\'figure_legend1\'> </p>\n')
	file.write('<p align=center><img src="img/fdr0.png"></p>\n')
	
	file.write('<p>&nbsp;</p>\n')
	file.write('\n')
	
	nTable = nTable+1
	file.write('<p class=\'table_legend1\'>Figure 2.%d.%d Tolerance of all precursors</p>\n' % (SampleNo,nTable))
	file.write('\n')
	file.write('<div align=center>\n')
	file.write('\n')
	file.write('<p id=\'alltol\' class=\'figure_legend1\'> </p>\n')
	file.write('<p align=center><img src="img/allprecur0.png"></p>\n')
	
	file.write('<p>&nbsp;</p>\n')
	file.write('\n')

	
	nTable = nTable+1
	file.write('<p class=\'table_legend1\'>Table 2.%d.%d Results of each raw</p>\n' % (SampleNo,nTable))
	file.write('\n')
	file.write('<div align=center>\n')
	file.write('\n')
	file.write('<table class=\'table1\'>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Raw file</td>\n')
	file.write('  <td class=\'table_head1\'>FDR curve</td>\n')
	file.write('  <td class=\'table_head1\'>Tolerance of precursors</td>\n')
	file.write('  <td class=\'table_head1\'>Search results</td>\n')
	file.write('  </tr>\n')
	raw_idx=0
	for L in spectra_list:
		raw_idx=raw_idx+1
		raw_title = L
		file.write(' <tr>\n')
		file.write('  <td class=\'table_data1\'>%s</td>\n' % raw_title)
		file.write('  <td class=\'table_data1\'><a href=\'raw%d_normal.html#fdr\'>Click to show...</a> </td>\n' % raw_idx)
		file.write('  <td class=\'table_data1\'><a href=\'raw%d_normal.html#alltol\'>Click to show...</a> </td>\n' % raw_idx)
		file.write('  <td class=\'table_data1\'><a href=\'raw%d_normal.html#result\'>Click to show...</a> </td>\n' % raw_idx)
		file.write(' </tr>\n')

		raw_normal_path = os.path.join(report_path, 'raw%d_normal.html' % raw_idx)
		WriteRawNormal(raw_normal_path, report_path, SampleNo, sample_title, raw_title, raw_idx, initparam)
	
	file.write('</table>\n')
	file.write('\n')
	file.write('</div>\n')

	file.write('<p>&nbsp;</p>\n')

	nTable = nTable+1
	table_legend = 'Table 2.%d.%d Combined Identified Results (From pBuild)\n' % (SampleNo,nTable)
	IdProTablepath = os.path.join(report_path, 'IdProTable_combine.html')
	protein_Info_file_name = '%s_%s_combine.protein.xls' % (sample_title, xlinktype)
	protein_Info_file_path = os.path.join(report_path,protein_Info_file_name)
	WriteIdResultTable(IdProTablepath, sample_title, -1, SampleNo, protein_Info_file_path, protein_Info_file_name, table_legend, report_path, initparam) #2013.5.3 modified

	IdPepTablepath = os.path.join(report_path, 'IdPepTable_combine.html')
	peptide_Info_file_name = '%s_%s_combine.peptide.xls' % (sample_title, xlinktype)
	peptide_Info_file_path = os.path.join(report_path,peptide_Info_file_name)
	WriteIdResultTable(IdPepTablepath, sample_title, -1, SampleNo, peptide_Info_file_path, peptide_Info_file_name, table_legend, report_path, initparam) #2013.5.3 modified

	IdSpecTablepath = os.path.join(report_path, 'IdSpecTable_combine.html')
	spec_Info_file_name = '%s_%s_combine.spectra.xls' % (sample_title, xlinktype)
	spec_Info_file_path = os.path.join(report_path,spec_Info_file_name)
	WriteIdResultTable(IdSpecTablepath, sample_title, -1, SampleNo, spec_Info_file_path, spec_Info_file_name, table_legend, report_path, initparam) #2013.5.3 modified
	
	spec_fdr_info_file = 'spec_u1_%s_%s_combine.statistics.txt' % (sample_title, xlinktype)
	spec_u1_fdr_Info = _GetSpecFDRInfo(os.path.join(report_path, spec_fdr_info_file))

	file.write('<p class=\'table_legend1\'>Table 2.%d.%d Combined Identified Result List (From pBuild)</p>\n' % (SampleNo,nTable))
	file.write('\n')
	file.write('<div align=center>\n')
	file.write('\n')
	file.write('<table class=\'table1\'>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>&nbsp;</td>\n')
	file.write('  <td class=\'table_head1\'># Spectra</td>\n')
	file.write('  <td class=\'table_head1\'># Peptides</td>\n')
	file.write('  <td class=\'table_head1\'># U1 proteins (groups)</td>\n')
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'> Spectral level</br>(FDR&lt;=5%) </td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % spec_u1_fdr_Info[0])
	file.write('  <td class=\'table_data1\'>%s</td>\n' % spec_u1_fdr_Info[1])
	file.write('  <td class=\'table_data1\'>%s (%s) </td>\n' % (spec_u1_fdr_Info[2],spec_u1_fdr_Info[3]))
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Result Details</td>\n')
	file.write('  <td class=\'table_data1\'><a href=\'./IdSpecTable_combine.html\'>Click here</a></td>\n')
	file.write('  <td class=\'table_data1\'><a href=\'./IdPepTable_combine.html\'>Click here</a></td>\n')
	file.write('  <td class=\'table_data1\'><a href=\'./IdProTable_combine.html\'>Click here</a></td>\n')
	file.write(' </tr>\n')
	file.write('</table>\n')
	file.write('\n')
	file.write('</div>\n')


	nTable = nTable+1
	table_legend = 'Table 2.%d.%d Identified protein list (From pBuild)\n' % (SampleNo,nTable)
	IdProTablepath = os.path.join(report_path, 'IdProTable.html')
	protein_Info_file_name = '%s_%s.protein.xls' % (sample_title, xlinktype)
	protein_Info_file_path = os.path.join(report_path,protein_Info_file_name)
	WriteIdResultTable(IdProTablepath, sample_title, -1, SampleNo, protein_Info_file_path, protein_Info_file_name, table_legend, report_path, initparam) #2013.5.3 modified
	file.write('<p class=\'table_legend1\'>Table 2.%d.%d Identified protein list (From pBuild)</p>\n' % (SampleNo,nTable))
	file.write('\n')
	file.write('<div align=center>\n')
	file.write('<a href=\'./IdProTable.html\'>Click to Show ...</a>\n')
	file.write('</div>\n')

	nTable = nTable+1
	table_legend = 'Table 2.%d.%d Identified peptide list (From pBuild)\n' % (SampleNo,nTable)
	IdPepTablepath = os.path.join(report_path, 'IdPepTable.html')
	peptide_Info_file_name = '%s_%s.peptide.xls' % (sample_title, xlinktype)
	peptide_Info_file_path = os.path.join(report_path,peptide_Info_file_name)
	WriteIdResultTable(IdPepTablepath, sample_title, -1, SampleNo, peptide_Info_file_path, peptide_Info_file_name, table_legend, report_path, initparam) #2013.5.3 modified

	file.write('<p class=\'table_legend1\'>Table 2.%d.%d Identified peptide list (From pBuild)</p>\n' % (SampleNo,nTable))
	file.write('\n')
	file.write('<div align=center>\n')
	file.write('<a href=\'./IdPepTable.html\'>Click to Show ...</a>\n')
	file.write('</div>\n')
	file.write('<p align=center><a href=\'#top\'>return to top</a>  |  <a href=\'../general.html\'>return to Introduction</a></p>\n')
	file.write('\n')
	
	file.write('</body>\n')
	file.write('\n')
	file.write('</html>\n')
	file.close()
	return

def WritePSMInfor(psm_html_file,sample_no,raw_no,spectra_title,psm_path,spectra_name_list):
	
	file = open(psm_html_file,'w');
	h = 'Peptide Spectra Matching Detail of Raw (' +spectra_title  +')'
	_WriteHTMLHead(file, h)
	file.write('<body>\n')
	file.write('<p id=\'top\' class=\'table_legend1\'>%s</p>\n' % h)
	file.write('\n')
	file.write('<p align=center><a href=\'./report.html\'>return to Sample %d</a>  |  <a href=\'../general.html\'>return to Introduction</a></p>\n' % sample_no)
	file.write('<div align=center>\n')
	file.write('\n')
	file.write('<table class=\'table1\'>\n')
	file.write('<tr>\n')
	file.write('<td class = \'table_head1\'>')
	file.write('psm')
	file.write('</td>\n')
	file.write('</tr>\n')
	for cur_spectrum_name in spectra_name_list:
		file.write('<tr>\n')
		file.write('<td class = \'table_data1\'>')
		file.write('<p id=\'%s\' class=\'figure_legend1\'> </p>\n' % cur_spectrum_name)
		file.write('<p align=center><img src="./psm/%s.png"></p>\n' % (cur_spectrum_name))
		file.write('</td>\n')
		file.write('</tr>\n')

	file.write('</table>\n')
	file.write('\n')
	file.write('</div>\n')
	file.write('<p align=center><a href=\'#top\'>return to top</a>  |  <a href=\'./report.html\'>return to Sample %d</a>  |  <a href=\'../general.html\'>return to Introduction</a></p>\n' % sample_no)
	file.write('\n')
	file.write('</body>\n')
	file.write('\n')
	file.write('</html>\n')

	file.close();
	return ;
	
def WriteRawNormal(filepath, report_path, SampleNo, sample_title, spectra_title, raw_no, initparam):
	xlinktype = 'inter'
	spec_fdr_info_file = 'spec_u1_%s_%d_%s.statistics.txt' % (sample_title, raw_no, xlinktype)
	spec_u1_fdr_Info = _GetSpecFDRInfo(os.path.join(report_path, spec_fdr_info_file))

	table_legend = 'Identified protein list of No.%d spectra from Sample %s (From pBuild)\n' % (raw_no,sample_title)
	IdProTablepath = os.path.join(report_path, 'IdProTable_%d.html' % raw_no)
	protein_Info_file_name = sample_title+'_%d_%s.protein.xls' % (raw_no, xlinktype)
	protein_Info_file_path = os.path.join(report_path,protein_Info_file_name)
	WriteIdResultTable(IdProTablepath, sample_title,raw_no,SampleNo, protein_Info_file_path, protein_Info_file_name, table_legend, report_path, initparam) #2013.5.3 modified
        	
	table_legend = 'Identified peptide list of No.%d spectra from Sample %s (From pBuild)\n' % (raw_no,sample_title)
	IdPepTablepath = os.path.join(report_path, 'IdPepTable_%d.html' % raw_no)
	peptide_Info_file_name = sample_title+'_%d_%s.peptide.xls' % (raw_no, xlinktype)
	peptide_Info_file_path = os.path.join(report_path,peptide_Info_file_name)
	WriteIdResultTable(IdPepTablepath, sample_title,raw_no,SampleNo, peptide_Info_file_path, peptide_Info_file_name, table_legend, report_path, initparam) #2013.5.3 modified

	table_legend = 'Identified spectra list of No.%d spectra from Sample %s (From pBuild)\n' % (raw_no,sample_title)
	IdSpecTablepath = os.path.join(report_path, 'IdSpecTable_%d.html' % raw_no)
	spectra_Info_file_name = sample_title+'_%d_%s.spectra.xls' % (raw_no, xlinktype)
	spectra_Info_file_path = os.path.join(report_path,spectra_Info_file_name)
	WriteIdResultTable(IdSpecTablepath, sample_title,raw_no,SampleNo, spectra_Info_file_path, spectra_Info_file_name, table_legend, report_path, initparam) #2013.5.3 modified

	file = open(filepath, 'w')
	h = 'Normal Search Detail of Raw (' +spectra_title  +')'
	_WriteHTMLHead(file, h)
	file.write('<body>\n')
	file.write('<p id=\'top\' class=\'title1\'>%s</p>\n' % h)

	file.write('<p id=\'fdr\' class=\'figure_legend1\'> </p>\n')
	file.write('<p align=center><img src="img/fdr%d.png"></p>\n' % raw_no)
	file.write('<p class=\'figure_legend1\'>Figure 1 FDR curve</p>\n')

	file.write('<p>&nbsp;</p>\n')
	file.write('\n')
	
	file.write('<p id=\'alltol\' class=\'figure_legend1\'> </p>\n')
	file.write('<p align=center><img src="img/allprecur%d.png"></p>\n' % raw_no)
	file.write('<p class=\'figure_legend1\'>Figure 2 Tolerance of all precursors</p>\n')

	file.write('<p>&nbsp;</p>\n')
	file.write('\n')
	
	file.write('<p id=\'result\' class=\'table_legend1\'>Table 1 Database search results (Frome pFind and pBuild)</p>\n' )
	file.write('\n')
	file.write('<div align=center>\n')
	file.write('\n')
	file.write('<table class=\'table1\'>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>&nbsp;</td>\n')
	file.write('  <td class=\'table_head1\'># Spectra</td>\n')
	file.write('  <td class=\'table_head1\'># Peptides</td>\n')
	file.write('  <td class=\'table_head1\'># U1 proteins (groups)</td>\n')
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'> Spectral level</br>(FDR&lt;=5%) </td>\n')
	file.write('  <td class=\'table_data1\'>%s</td>\n' % spec_u1_fdr_Info[0])
	file.write('  <td class=\'table_data1\'>%s</td>\n' % spec_u1_fdr_Info[1])
	file.write('  <td class=\'table_data1\'>%s (%s) </td>\n' % (spec_u1_fdr_Info[2],spec_u1_fdr_Info[3]))
	file.write(' </tr>\n')
	file.write(' <tr>\n')
	file.write('  <td class=\'table_head1\'>Result Details</td>\n')
	file.write('  <td class=\'table_data1\'><a href=\'./IdSpecTable_%d.html\'>Click here</a></td>\n' % raw_no)
	file.write('  <td class=\'table_data1\'><a href=\'./IdPepTable_%d.html\'>Click here</a></td>\n' % raw_no)
	file.write('  <td class=\'table_data1\'><a href=\'./IdProTable_%d.html\'>Click here</a></td>\n' % raw_no)
	file.write(' </tr>\n')
	file.write('</table>\n')
	file.write('\n')
	file.write('</div>\n')

	file.write('\n')
	file.write('<p align=center><a href=\'#top\'>return to top</a>  |  <a href=\'./report.html\'>return to Sample %d</a>  |  <a href=\'../general.html\'>return to Introduction</a></p>\n' % SampleNo)
	file.write('\n')
	file.write('</body>\n')
	file.write('\n')
	file.write('</html>\n')
	file.close()
	return

def RmIllegalCharsForFile(strTitle): # added at 2013.5.3 by fan, remove the illegal Chars for File name
        strTitle = strTitle.replace("\\", " ")
        strTitle = strTitle.replace("/", " ")
        strTitle = strTitle.replace(":", " ")
        strTitle = strTitle.replace("*", " ")
        strTitle = strTitle.replace("?", " ")
        strTitle = strTitle.replace("<", " ")
        strTitle = strTitle.replace(">", " ")
        strTitle = strTitle.replace("|", " ")
        strTitle = strTitle.replace("\"", " ")

        return strTitle

def WriteIdResultTable(filepath, sample_title, raw_no ,SampleNo, protein_Info_file_path, protein_Info_file_name, table_legend, report_path, initparam):#2013.5.3 modified
	inputfile = open(protein_Info_file_path, 'r')
	file = open(filepath, 'w')
	_WriteHTMLHead(file, table_legend)
	file.write('<body>\n')
	file.write('<p id=\'top\' class=\'table_legend1\'>%s</p>\n' % table_legend)
	file.write('\n')
	file.write('<p align=center><a href=\'./report.html\'>return to Sample %d</a>  |  <a href=\'../general.html\'>return to Introduction</a></p>\n' % SampleNo)
	file.write('<div align=center>\n')
	file.write('\n')
	file.write('<table class=\'table1\'>\n')

	while 1:
		str = inputfile.readline()
		if not str:
			break
		L = str.split('\t')
		cl = 'none'
		if len(L) >= 2:
			if L[0] == 'Order' or L[1] == 'Order':
				cl = 'table_head1'
			else:
				cl = 'table_data1'
	
		file.write('<tr>\n')


##                if L[1].isdigit():
##                        FileName = RmIllegalCharsForFile(L[2]) #modified at 2013.5.3
##                        TotalFileName = os.path.join(report_path, 'psm', FileName)
##                        TotalFileName = TotalFileName + ".png"
##
##                        file.write('<td class = \'%s\'>' % cl)
##                        if not os.path.isfile(TotalFileName):
##                                print ('%s doesn\'t exist. Skip.' % TotalFileName)
##                                file.write(FileName)                                
##                        else:
##                                file.write('<a href=\'./psm/%s.png\'>%s</a>\n' % (FileName, FileName));
##                        file.write('</td>\n')   
                                
		for item in L:                        
			dtaname = '';
			pos = item.find('.dta');
			pos2 = item.find('Scan')        # modified at 2013.5.7, add the support to another scan type
			pos3 = item.find('scan')		# modified at 2013.9.26, add the support to another scan type, bug fix
			
			if (pos >= 0 and pos == len(item) - 4) or (pos2 >= 0) or (pos3 >= 0): # modified at 2013.5.7
				dtaname = RmIllegalCharsForFile(item);
			else:
				dtaname = '';
			
			if (initparam['drawpsm'] != 'true'): #added at 2014.2.11
				dtaname = ''
		
			file.write('<td class = \'%s\'>' % cl)
			if dtaname != '':
				file.write('<a href=\'./psm/%s.png\'>%s</a>\n' % (dtaname,dtaname));
			else:
				file.write(item)
			file.write('</td>\n')
			
		file.write('</tr>\n')

	file.write('</table>\n')
	file.write('\n')
	file.write('<a href=\'%s\'>Excel File ...</a>\n' % protein_Info_file_name)
	
	if fnmatch.fnmatch(protein_Info_file_name,'*.spectra.xls'):
		plabel_file_name = '%s%d.plabel' % (sample_title,raw_no);
		file.write('<a href=\'%s\'>pLabel File ...</a>\n' % plabel_file_name)
		
	file.write('</div>\n')
	file.write('<p align=center><a href=\'#top\'>return to top</a>  |  <a href=\'./report.html\'>return to Sample %d</a>  |  <a href=\'../general.html\'>return to Introduction</a></p>\n' % SampleNo)
	file.write('\n')
	file.write('</body>\n')
	file.write('\n')
	file.write('</html>\n')
	file.close()
	inputfile.close()
	return

