#!/usr/bin/env python
#-*-coding:utf8 -*-

import os
import sys
import string
import platform
import time

import Path
import Filter
import pBuild # for the compare only, use the old py file

def _GetVennInfo(file_path):
	spec_venn_info = []
	pep_venn_info = []
	pro_venn_info = []

	file = open(file_path, 'r')

	while 1:
		str = file.readline()
		if not str:
			break
		if str == '[END]\n':
			break

	venn_info_now = []
	while 1:
		str = file.readline()
		if not str:
			break
		if str == '[Spectra]\n':
			venn_info_now = spec_venn_info
		else:
			if str == '[Peptide]\n':
				venn_info_now = pep_venn_info
			else:
				if str == '[Protein]\n':
					venn_info_now = pro_venn_info
				else:
					if str != '[END]\n':
						venn_info_now.append(str[:-1].split(' '))

	file.close()

	return spec_venn_info,pep_venn_info,pro_venn_info


def _WriteCompareResult(file_path, sample_num, spec_venn_info,pep_venn_info,pro_venn_info):
	total_line_number = 2** sample_num - 1
	sub_spec = {}
	sub_pep = {}
	sub_pro = {}
	for i in range(1, total_line_number+1):
		sub_spec[i] = '0'
		sub_pep[i] = '0'
		sub_pro[i] = '0'

	for L in spec_venn_info:
		idx =string.atoi(L[0])
		sub_spec[idx] = L[1]

	for L in pep_venn_info:
		idx =string.atoi(L[0])
		sub_pep[idx] = L[1]

	for L in pro_venn_info:
		idx =string.atoi(L[0])
		sub_pro[idx] = L[1]

	file = open(file_path, 'w')
	file.write('spectra')
	for i in range(1, total_line_number+1):
		file.write(',%s' % sub_spec[i])
	file.write('\n')

	file.write('peptide')
	for i in range(1, total_line_number+1):
		file.write(',%s' % sub_pep[i])
	file.write('\n')

	file.write('protein')
	for i in range(1, total_line_number+1):
		file.write(',%s' % sub_pro[i])
	file.write('\n')

	file.close()

def Run(initparam,database_fasta,sample_num,compare_path):
	bin_path = initparam['bin.path']
	Path.Check(compare_path)

	[total_output_path,tmp] = os.path.split(compare_path)
	pFind_list = []
	for i in range(1,sample_num+1):
		key = 'sample%d' % i
		title = initparam[key+'.spectra.title'] + "_inter" # modified at 2013.12.23
		pfind_result_file = os.path.join(total_output_path,'%d.sample' % i,'search',title,'%s_DataSet1_pfind.txt' % title)
		pFind_list.append([pfind_result_file])

	SampleID = 'compare'
	compare_result_path = os.path.join(compare_path,SampleID,SampleID+'_Java','Sample_VennDiagram.txt')

        if platform.system() == 'Linux':
                bat_file=os.path.join(compare_path,'normal.bash') # modified at 2013.12.23
                file=open(bat_file, 'w')
    
        if platform.system() == 'Windows':
                bat_file=os.path.join(compare_path,'normal.bat') # modified at 2013.12.23
                file=open(bat_file, 'w')
                file.write('%s\n' % bin_path[0:2])

	file.write('cd %s\n' % bin_path)

	pBuild_Param={}
	pBuild_Param['SampleID'] = SampleID
	pBuild_Param['Fasta_File_Path'] = database_fasta
	pBuild_Param['OutPutPath'] = compare_path
	pBuild_Param['OutPutName'] = SampleID
	pBuild_Param['ExportPath'] = compare_path
	pBuild_Param['ExportFormat'] = '0,2'
	pBuild_Param['Redundant']='1'
	pBuild_Param['DistinctPep']='1'
	pBuild_Param['PepTolBase']='0'
	pBuild_Param['PepTolLower']='-100'
	pBuild_Param['PepTolUpper']='100'
	pBuild_Param['PepTolType']='Da'

        pbuild_file_path = os.path.join(compare_path,'compare.pbuild')
	pBuild.Write(pbuild_file_path, pBuild_Param, pFind_list)

	if os.path.isfile(compare_result_path) == 0:
		print 'Step : Compare results in different samples by pBuild'
		file.write('%s %s\n' % (os.path.join(bin_path,'builder'), pbuild_file_path))
	file.close()
        if platform.system() == 'Linux':
                os.system('chmod 766 %s' % bat_file)
	os.system(bat_file)

	[spec_venn_info,pep_venn_info,pro_venn_info] = _GetVennInfo(compare_result_path)
	compare_file=os.path.join(compare_path, 'compare.csv')
	_WriteCompareResult(compare_file, sample_num, spec_venn_info,pep_venn_info,pro_venn_info)

	compare_time = os.stat(compare_result_path).st_mtime
	#compare_time = time.time();
	
	return compare_time
