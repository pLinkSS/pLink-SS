#!/usr/bin/env python
# -*-coding:utf8 -*-
import os
import platform
import time
import string
import ConfigParser
import Path
import ClusterSetting
	
def _GetProteinSeqNum(FastaInfo_file):
	Index_info = []
	ffile = open(FastaInfo_file, 'r')
	ss = ffile.readline()
	while 1:
		ss = ffile.readline()
		if not ss:
			break
		Index_info.append(ss.split(','))
	ffile.close()
	proteinSeqNum = string.atoi(Index_info[0][1]);
	return proteinSeqNum

def _WriteFastaInfo(FastaInfo_file, meta_file, database_fasta):
	cfg = ConfigParser.ConfigParser()
	cfg.read(meta_file)
	v = cfg.items('Head')
	for x in v:
		if 'pronum' == x[0]:
			pronum = string.atoi(x[1]) / 2

	filesize = os.stat(database_fasta).st_size / (1024.0 * 1024.0)
	ctime = os.stat(database_fasta).st_ctime
	mtime = os.stat(database_fasta).st_mtime

	ffile = open(FastaInfo_file, 'w')
	ffile.write('Fasta file,# Proteins,File size (MB),Create time,Last modified time\n')
	ffile.write('%s,' % database_fasta)
	ffile.write('%d,' % pronum)
	ffile.write('%.1f,' % filesize)
	ffile.write('%s,' % time.ctime(ctime))
	ffile.write('%s,' % time.ctime(mtime))
	ffile.write('\n')
	ffile.close()

	return

def _WriteIndexInfo(IndexInfo_file, index_param, pep_meta_file):
	Cleavage = {0:'full', 1:'semi-specific', 2:'non-specific'}

	fasta_file = index_param['fasta_file']
	enzyme1 = index_param['enzyme1']
	max_miss_site = index_param['max_miss_site']
	cleave_way = string.atoi(index_param['cleave_way'])
	min_pep_length = index_param['min_pep_length']
	max_pep_length = index_param['max_pep_length']
	
	min_mass = index_param['min_mass']
	max_mass = index_param['max_mass']

	cfg = ConfigParser.ConfigParser()
	cfg.read(pep_meta_file)
	v = cfg.items('Head')
	for x in v:
		if 'pep_sq_num' == x[0]:
			pep_sq_num = string.atoi(x[1])
		if 'unique_pep_sq_num' == x[0]:
			unique_pep_sq_num = string.atoi(x[1])

	ctime = os.stat(fasta_file).st_ctime
	
	ifile = open(IndexInfo_file, 'w')
	ifile.write('Enzyme,Max missed site,Cleavage way,Min peptide length,Max peptide length,Min mass,Max mass,# Redundant peptides,# Unique peptides,Create time\n')
	ifile.write('%s,' % enzyme1)
	ifile.write('%s,' % max_miss_site)
	ifile.write('%s,' % Cleavage[cleave_way])
	ifile.write('%s,' % min_pep_length)
	ifile.write('%s,' % max_pep_length)
	ifile.write('%s,' % min_mass)
	ifile.write('%s,' % max_mass)
	ifile.write('%d,' % pep_sq_num)
	ifile.write('%d,' % unique_pep_sq_num)
	ifile.write('%s,' % time.ctime(ctime))
	ifile.write('\n')
	ifile.close()

	return

def _ConstructpIndex(initparam):
	pindex_param = {}
	
	# total section
	total = []
	# out_path
	bin_path = initparam['bin.path']
	out_path = os.path.join(bin_path, 'index')
	Path.Check(out_path)
	if platform.system()=='Linux':
		out_path = out_path+'/'
	if platform.system()=='Windows':
		out_path = out_path+'\\'
	total.append(('out_path', out_path))
	# aa_list
	aa_list = os.path.join(bin_path, 'aa.ini')
	Path.FileCheck(aa_list)
	total.append(('aa_list', aa_list))
	# enzyme_list
	enzyme_list = os.path.join(bin_path, 'enzyme.ini')
	Path.FileCheck(enzyme_list)
	total.append(('enzyme_list', enzyme_list))
	# db_num
	total.append(('db_num', '1'))
	pindex_param['[total]'] = total
	
	# db1 section
	db1 = []
	# fasta_file
	fasta_file = initparam['database.path']
	db1.append(('fasta_file', fasta_file))
	# db_name
	db_name = initparam['database.name']
	db1.append(('db_name', db_name))
	# enzyme_number
	db1.append(('enzyme_number', '1'))
	# enzyme1
	enzyme1 = initparam['enzyme.name']
	db1.append(('enzyme1', enzyme1))
	# max_miss_site
	if not 'max_miss_site' in initparam:
		max_miss_site = '2'
	else:
		max_miss_site = initparam['max_miss_site']
	db1.append(('max_miss_site', max_miss_site))
	# min_pep_length
	if not 'min_pep_length' in initparam:
		min_pep_length = '4'
	else:
		min_pep_length = initparam['min_pep_length']
	db1.append(('min_pep_length', min_pep_length))
	# max_pep_length
	if not 'max_pep_length' in initparam:
		max_pep_length = '100'
	else:
		max_pep_length = initparam['max_pep_length']	
	db1.append(('max_pep_length', max_pep_length))
	# min_mass
	if not 'min_mass' in initparam:
		min_mass = '400'
	else:
		min_mass = initparam['min_mass']	
	db1.append(('min_mass', min_mass))
	# max_mass
	if not 'max_mass' in initparam:
		max_mass = '10000'
	else:
		max_mass = initparam['max_mass']
	db1.append(('max_mass', max_mass))
	# mono
	db1.append(('mono', '1'))
	# avrg
	db1.append(('avrg', '1'))
	# cleave_way
	if not 'cleave_way' in initparam:
		cleave_way = '0'
	else:
		cleave_way = initparam['cleave_way']
	db1.append(('cleave_way', cleave_way))
	# auto_reverse
	db1.append(('auto_reverse', '1'))
	# max_mem_size
        if not 'max_mem_size' in initparam: #added at 2014.2.20
            max_mem_size = '236870912'
        else:
            max_mem_size = initparam['max_mem_size']
	db1.append(('max_mem_size', max_mem_size))
	# Remark
	db1.append(('Remark', ''))
	pindex_param['[db1]'] = db1
	
	return pindex_param

def _WritepIndex(pindex_param, pindex_file_path):
	sections = ['[total]', '[db1]']
	pindex_file = open(pindex_file_path, 'w')
	for section in sections:
		pindex_file.write(section+'\n')
		for elem in pindex_param[section]:
			pindex_file.write(elem[0]+'='+elem[1]+'\n')
	return

def Run(initparam, index_dir_path):
	# 0.index
	Path.Check(index_dir_path)
	
	# release/index
	bin_dir_path = initparam['bin.path']
	index_save_path = os.path.join(bin_dir_path, 'index')
	Path.Check(index_save_path)
	
	# *.pindex
	database_name = initparam['database.name']
	pindex_file_path = os.path.join(index_dir_path, database_name + '.pindex')
	if os.path.isfile(pindex_file_path):
		index_param_time = os.stat(pindex_file_path).st_mtime
	else:
		index_param_time = time.time()
		
	enzyme_name = initparam['enzyme.name']  # attention
	pro_meta_file = os.path.join(index_save_path, database_name + '.pro.meta.txt')
	pep_meta_file = os.path.join(index_save_path, database_name + '.' + enzyme_name + '.mono.meta.txt')
	FastaInfo_file = os.path.join(index_dir_path, 'fastainfo.csv')
	IndexInfo_file = os.path.join(index_dir_path, 'indexinfo.csv')
	if Path.CheckTime(pep_meta_file, index_param_time):
		pindex_param = _ConstructpIndex(initparam)
		_WritepIndex(pindex_param, pindex_file_path)
		
		if platform.system() == 'Linux':
			bat_file_path = os.path.join(index_dir_path, 'index.bash')
			bat_file = open(bat_file_path, 'w')
			bat_file.write('export PATH=%s:$PATH\n' % ClusterSetting.MPIPath)  # modified 2012.6.11
			bat_file.write('export LD_LIBRARY_PATH=/%s:$LD_LIBRARY_PATH\n' % ClusterSetting.pLinkBinPath)
		
		if platform.system() == 'Windows':
			bat_file_path = os.path.join(index_dir_path, 'index.bat')
			bat_file = open(bat_file_path, 'w')
			bat_file.write('%s\n' % bin_dir_path[0:2])
			
		bat_file.write('cd "%s"\n' % bin_dir_path)  # modified 2012.6.2 ''	
		bat_file.write('"%s" "%s"' % (os.path.join(bin_dir_path, 'Indexer'), pindex_file_path))
		bat_file.close()
		
		if platform.system() == 'Linux':
			os.system('chmod 766 %s' % bat_file_path)
		
		# execution
		os.system('"%s"' % bat_file_path)
		
		fasta_file_path = initparam['database.path']
		_WriteFastaInfo(FastaInfo_file, pro_meta_file, fasta_file_path)
		db1 = pindex_param['[db1]'];
		db1_param = {}
		for elem in db1:
			db1_param[elem[0]] = elem[1]
		_WriteIndexInfo(IndexInfo_file, db1_param, pep_meta_file)
	
	index_time = os.stat(pep_meta_file).st_mtime
	pro_seq_num = _GetProteinSeqNum(FastaInfo_file)
	return pro_seq_num, index_time
