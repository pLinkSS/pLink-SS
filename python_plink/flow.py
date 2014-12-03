#!/usr/bin/env python
#-*-coding:utf8 -*-
import os
import sys
import time
import string
import platform
import Param
import Path
import Index
import Normal
import Report
import Sampleinfo
import Compare
import Disulfide
from multiprocessing import cpu_count

def main():
	if len(sys.argv) < 2:
		para_file_path = 'pLink.ini'
	else:
		para_file_path = sys.argv[1]
	
# time control
       	if platform.system() == 'Windows':
                if string.atoi(time.strftime('%Y', time.localtime(time.time()))) > 2014:
                        print 'Sorry, pLink expired after 2013. Please register again. Thanks for your support to pLink.'
                        return
# param
	try:
		initparam = Param.Read(para_file_path)
		ModifyInfo = Param.ReadModify(initparam)
		AAInfo = Param.ReadAA(initparam)
		LinkerInfo = Param.ReadLinker(initparam)
		output_dir_path = initparam['output.path']
		Path.Check(output_dir_path)
	except Exception,data:
		print 'Failed to read the param:', data
		return

# start log
	start_time = time.time()
	log_path = os.path.join(output_dir_path,'pLink.log.txt')
	log_file = open(log_path,'a')
	log_file.write('%s\n' % ('*'*40))
	log_file.write('start time: %s\n' % time.ctime(start_time))
	log_file.write('parameter file: %s\n' % para_file_path)

# index
	try:
		index_dir_path = os.path.join(output_dir_path,'0.index')
		log_file.write('\nstep 0: index\n')
		log_file.write('\tdirectory: %s\n' % index_dir_path)

		t1 = time.time()
		[pro_seq_num, index_time] = Index.Run(initparam, index_dir_path)
		t2 = time.time()
		del_time = t2-t1
		
		search_mode = 0
		# search mode according to protein seq number
		if pro_seq_num < 100:
			search_mode = 0
		else:
			search_mode = 1

		#added at 2013.5.20, SUMO will use all mode
		if 'sumo' in initparam:
			search_mode = 0
		
		#added at 2013.9.30, search_mode can be selected
		if 'search_mode' in initparam:
			if initparam['search_mode'] == '1':
				search_mode = 1
			elif initparam['search_mode'] == '0':
				if pro_seq_num > 100:
					print 'Notice: exhaustion mode has been selected. This may cost long time for large database.'
					print 'Program will continue in 2 seconds.'
					time.sleep(2)
				search_mode = 0
			elif initparam['search_mode'] == '2': #For sumo
				search_mode = 2
			elif initparam['search_mode'] == '3': #For tri-pep
				search_mode = 3
			elif initparam['search_mode'] == '4': #For tri-pep ion index
				search_mode = 4
				
		#added at 2013.12.23, check processor number for search_mode
		if 'processor_num' in initparam:
			if int(initparam['processor_num']) > 1:
				if search_mode == 0:
					print 'Sorry, exhaustion mode can not use multi-thread currently. Will use 1 thread.'
					print 'Program will continue in 2 seconds.'
					time.sleep(2)
					initparam['processor_num'] = '1'

		#added at 2013.12.23, check the max processor number the computer can use
		if 'processor_num' in initparam:
			if int(initparam['processor_num']) > cpu_count():
				print 'Notice: process_num is larger than the max processor number of this computer. Will use max processor number - 1.'
				print 'Program will continue in 2 seconds.'
				time.sleep(2)
				initparam['processor_num'] = cpu_count() - 1
		
		print '0.index: ',time.ctime(index_time)
		log_file.write('\tlast modified time: %s\n' % time.ctime(index_time))
		log_file.write('\trunning time: %.1f (s)\n' % del_time)
		
		time_file = os.path.join(index_dir_path,'running_time.txt')
		Path.WriteTime(time_file,index_time,del_time)
		
		print 'protein sequence number = ', pro_seq_num
		print 'search mode =', search_mode
		log_file.write('\tprotein sequence number : %d\n' % pro_seq_num);
		log_file.write('\tsearch mode : %d\n' % search_mode);
		initparam['search_mode'] = str(search_mode)
		
	except Exception, data:
		print 'Failed to run Index:',data
		log_file.write('\tFailed to run Index: %s\n' % data)
		log_file.close()
		return
	
	if 'SS' == initparam['linker.name1'] or 'SS_0' == initparam['linker.name1']:
        # Skip the tri-pep identification
		if not ( search_mode == 4 or search_mode == 3 ):
			# for disulfide bond search and filter
			if pro_seq_num <= 20:
				initparam['fdr_strategy'] = 1
			else:
				initparam['fdr_strategy'] = 0
			initparam['origin.output.path'] = output_dir_path
			if 'filter_peptide_tol_type' not in initparam:
				initparam['filter_peptide_tol_type'] = 'ppm'
			else:
				initparam['filter_peptide_tol_type'] = initparam['filter_peptide_tol_type'].lower()
			if 'filter_peptide_tol_ub' not in initparam:
				initparam['filter_peptide_tol_ub'] = 10
	
			Disulfide.main(log_file, initparam)
			end_log(log_file, log_path, start_time)
			return

	# Samples
	Sampleinfo_list = []
	sample_num = string.atoi(initparam['sample.num'])
	for sno in range(1, sample_num+1):
		log_file.write('%s\n' % ('-'*30))

		# param
		key = 'sample%d' % sno
		initparam['spectra.instrument'] = initparam[key+'.spectra.instrument']
		initparam['spectra.format'] = initparam[key+'.spectra.format']
		initparam['spectra.path'] = initparam[key+'.spectra.path']
		initparam['spectra.title'] = initparam[key+'.spectra.title']
		initparam['output.path'] = os.path.join(output_dir_path,'%d.sample' % sno)
		Path.Check(initparam['output.path'])

		# MGF list
		file_list = [];
		spectra_type = initparam['spectra.format'];
		if platform.system() == 'Linux':
			spec_format = ('*.' + spectra_type).upper();#modified at 2012.3.2
			Path.Walk(initparam['spectra.path'],spec_format,file_list);
		spec_format = ('*.' + spectra_type).lower();
		Path.Walk(initparam['spectra.path'],spec_format,file_list);
		
		# todo for DTA
		spectra_list = [];
		for i in file_list:
			spectra_list.append(os.path.join(i[0],i[1]));
		
		print 'spectra list in sample %d are :' % sno ;
		log_file.write('spectra list in sample %d are :' % sno);
		for i in spectra_list:
			print i;
			log_file.write('%s\n' % i);
		
		try:
			# Search
			search_dir_path=os.path.join(initparam['output.path'],'search')
			log_file.write('\nstep %d: search\n' % sno)
			log_file.write('\tdirectory: %s\n' % search_dir_path)
			
			t1 = time.time()
			search_time = Normal.Run(initparam, spectra_list, search_mode, search_dir_path)
			t2 = time.time()
			del_time = t2-t1

			print '1.search: ', time.ctime(search_time)
			log_file.write('\tlast modified time: %s\n' % time.ctime(search_time))
			log_file.write('\trunning time: %.1f (s)\n' % del_time)
			time_file = os.path.join(search_dir_path,'running_time.txt')
			Path.WriteTime(time_file,search_time,del_time)
		except Exception,data:
			print 'Failed to run normal search:', data
			log_file.write('\tFailed to run normal search: %s\n' % data)
			log_file.close()
			return
		
		spectra_title = initparam['spectra.title']
		SampleIDList = []
		for i in range(0, len(spectra_list)):
			SampleIDList.append('%s_%d_inter' % (spectra_title, (i+1)))
		SampleIDList.append('%s_inter' % spectra_title)
		SampleIDList.append('%s_inter_combine' % spectra_title)
		s = Sampleinfo.Sampleinfo(initparam['spectra.title'],initparam['spectra.path'],SampleIDList,spectra_list,search_dir_path)
		Sampleinfo_list.append(s)
	
	# compare annotated at 2013.12.23, compare not use
	sno = sample_num
	compare_path = ''
	if sample_num>1:
		try:
			sno = sno+1
			compare_path = os.path.join(output_dir_path,'%d.compare' % sno) # modified at 2013.12.23
			log_file.write('%s\n' % ('-'*30))
			log_file.write('\nstep %d: compare\n' % (sample_num+1))
			log_file.write('\tdirectory: %s\n' % compare_path)

			t1 = time.time()
			database_fasta = initparam['database.path'] # modified at 2013.12.23
			compare_time = Compare.Run(initparam,database_fasta,sample_num,compare_path) 
			t2 = time.time()
			del_time = t2-t1

			log_file.write('\trunning time: %.1f (s)\n' % del_time)
		except  Exception,data:
			print 'Failed to compare:',data
			log_file.write('\tFailed to compare: %s\n' % data)
			log_file.close()
			return
	# report
	try:
		sno = sno+1
		report_dir_path = os.path.join(output_dir_path,'%d.report' % sno)
		log_file.write('%s\n' % ('-'*30))
		log_file.write('\nstep %d: report\n' % sno)
		log_file.write('\tdirectory: %s\n' % report_dir_path)

		t1 = time.time()

		Report.Run(initparam,Sampleinfo_list,index_dir_path, compare_path,report_dir_path) # modified at 2013.12.23

		t2 = time.time()
		del_time = t2-t1

		log_file.write('\trunning time: %.1f (s)\n' % del_time)
	except  Exception,data:
		print 'Failed to report results:',data
		log_file.write('\tFailed to report results: %s\n' % data)
		log_file.close()
		return

	# end log
	end_log(log_file, log_path, start_time)
	return

def end_log(log_file, log_path, start_time):
	# end log
	end_time = time.time()
	log_file.write('\nend time: %s\n' % time.ctime(end_time))
	log_file.write('total running time: %.1f (s)\n' % (end_time-start_time))
	log_file.close()
	
	if platform.system() == 'Linux':
		cmdline = 'cat "' + log_path + '"'
		os.system(cmdline)
	if platform.system() == 'Windows':
		cmdline = 'type "' + log_path + '"'
		os.system(cmdline)

if __name__ == '__main__':
	main()
