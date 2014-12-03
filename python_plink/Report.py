#!/usr/bin/env python
#-*-coding:utf8 -*-
import os
import sys
import string
import time
import fnmatch
import shutil
import thread

import Img
import Template
import Path
import Sampleinfo
import PSMDraw

def _CopypIndexInfo(index_path, report_path):
	Path.Copy(index_path,'*.csv',report_path)
	

def _CopypNormalInfo(normal_path, SampleIDlist, report_path):
	raw_idx=0
	for SampleID in SampleIDlist:
		# copy FDR files
		srcpath = os.path.join(normal_path,SampleID,SampleID+'_Java','Charge')
		dir_list = os.listdir(srcpath)
		for f1 in dir_list:
			cur_path = os.path.join(srcpath,f1)
			if os.path.isfile(cur_path) and (len(f1)>14 and f1[0:14]=='Sample1_Charge'):
                                shutil.copy(cur_path,os.path.join(report_path, ('%d.' % raw_idx) + f1))
				#cmdline = 'copy '+cur_path+' '+ os.path.join(report_path, ('%d.' % raw_idx) + f1)
				#os.system(cmdline)

		# copy the errorfigure file
		#cmdline = 'copy '+os.path.join(normal_path,SampleID,SampleID+'_Java','ErrorFigure','Sample1_Da.errorfigure')+' '+ os.path.join(report_path, '%d_Da.errorfigure' % raw_idx)
		#os.system(cmdline)
                shutil.copy(os.path.join(normal_path,SampleID,SampleID+'_Java','ErrorFigure','Sample1_Da.errorfigure'),os.path.join(report_path, '%d_Da.errorfigure' % raw_idx))
		#cmdline = 'copy '+os.path.join(normal_path,SampleID,SampleID+'_Java','ErrorFigure','Sample1_ppm.errorfigure')+' '+ os.path.join(report_path, '%d_ppm.errorfigure' % raw_idx)
		#os.system(cmdline)
                shutil.copy(os.path.join(normal_path,SampleID,SampleID+'_Java','ErrorFigure','Sample1_ppm.errorfigure'),os.path.join(report_path, '%d_ppm.errorfigure' % raw_idx))
                
		# copy the pbuild result files
		#cmdline = 'copy '+os.path.join(normal_path,SampleID,'*.xls')+' '+report_path
		#os.system(cmdline)
		Path.Copy(os.path.join(normal_path,SampleID),'*.xls',report_path)
		
		# copy statistics files
		srcpath = os.path.join(normal_path,SampleID)
		dir_list = os.listdir(srcpath)
		for f1 in dir_list:
			cur_path = os.path.join(srcpath,f1)
			if os.path.isfile(cur_path) and fnmatch.fnmatch(f1,'*.statistics.txt'):
				#cmdline = 'copy '+cur_path+' '+ os.path.join(report_path, ('spec_u1_'+ f1))
				#os.system(cmdline)
				shutil.copy(cur_path,os.path.join(report_path, ('spec_u1_'+ f1)))

		raw_idx = raw_idx+1

def _CopypCompareResult(compare_path, report_path):
	venn_path = os.path.join(compare_path,'compare.csv')
	#cmdline = 'copy '+venn_path+' '+report_path
	#os.system(cmdline)
	shutil.copy(venn_path,report_path)

def _GetTime(file_path):
	mytime = 0.0
	if False==os.path.isfile( file_path ):
		raise Exception,file_path
	else:
		file = open(file_path,'r')
		str = file.readline()
		mytime = string.atof(str)
		file.close()

	return mytime

def _GetRunningTime(index_path,normal_path):
	running_time_list = []
	time_file = 'running_time.txt'

	file_path = os.path.join(index_path,time_file)
	mytime = _GetTime(file_path)
	running_time_list.append(mytime)

	file_path = os.path.join(normal_path,time_file)
	mytime = _GetTime(file_path)
	running_time_list.append(mytime)

	return running_time_list


def Run(initparam,Sampleinfo_list,index_path,compare_path,report_path):
	output_path = initparam['output.path']
	bin_path = initparam['bin.path']
	show_precursor_deviation_type = initparam['show_peptide_tol_type']
	
	Path.Check(report_path)
	GeneralImgPath = os.path.join(report_path, 'img')
	Path.Check(GeneralImgPath)

	# get general information
	
	print 'Step : Copy Index Info';
	_CopypIndexInfo(index_path, report_path)
	if len(Sampleinfo_list) > 1:
		print 'Copy Compare Result';
		_CopypCompareResult(compare_path, report_path)


	for i in range(0,len(Sampleinfo_list)):
		IDlist = Sampleinfo_list[i].IDlist
		mgf_path = Sampleinfo_list[i].mgf_path
		normal_path = Sampleinfo_list[i].normal_path

		sample_path = os.path.join(report_path, 'sample%d' % (i+1))
		Path.Check(sample_path)

		img_path = os.path.join(sample_path, 'img')
		Path.Check(img_path)
		
		# copy the psm_analysis result files
		print 'Step: copy *.psm.output for Sample '+ Sampleinfo_list[i].spectra_title;
		#cmdline = 'copy '+os.path.join(normal_path,'*.psm.output') + ' ' + sample_path
		#os.system(cmdline)
		Path.Copy(normal_path,'*.psm.output',sample_path)
		
		# copy the *.plabel files
		print 'Step: copy *.plabel for Sample '+ Sampleinfo_list[i].spectra_title;
		#cmdline = 'copy '+os.path.join(normal_path,'*.plabel') + ' ' + sample_path
		#os.system(cmdline)
		Path.Copy(normal_path,'*.plabel',sample_path)

		print 'Step: Copy Normal Info'
		_CopypNormalInfo(normal_path, IDlist, sample_path)

	GeneralReportPath = os.path.join(report_path, 'general.html')
	print 'Step: Write General '+GeneralReportPath;
	Template.WriteGeneral(GeneralReportPath, initparam, Sampleinfo_list, report_path)
	report_time = 0
	
	# get information for each sample
	for i in range(0,len(Sampleinfo_list)):
		spectra_title = Sampleinfo_list[i].spectra_title
		spectra_path = Sampleinfo_list[i].spectra_path
		IDlist = Sampleinfo_list[i].IDlist
		spectra_list = Sampleinfo_list[i].mgf_path
		normal_path = Sampleinfo_list[i].normal_path

		sample_path = os.path.join(report_path, 'sample%d' % (i+1))
		img_path = os.path.join(sample_path, 'img')

		mgflist = spectra_list

		# AllPrecur
		AllPrecurImg = 'allprecur0.png'
		
		if show_precursor_deviation_type == 'ppm':
			error_file = '0_ppm.errorfigure'
		else:
			error_file = '0_Da.errorfigure'
			
		print 'Step: Draw precursor deviation for Sample ' + spectra_title;
		Img.WriteAllPrecurImg(img_path, AllPrecurImg,show_precursor_deviation_type, error_file, sample_path)
		for j in range(1,len(mgflist)+1):
			raw_img = 'allprecur%d.png' % j
			if show_precursor_deviation_type == 'ppm':
				error_file = '%d_ppm.errorfigure' % j
			else:
				error_file = '%d_Da.errorfigure' % j
			print 'Step: Draw precursor deviation for raw %d of Sample ' % j + spectra_title;
			Img.WriteAllPrecurImg(img_path, raw_img,show_precursor_deviation_type, error_file, sample_path)
			
		# FDR
		FDRImg = 'fdr0.png';
		
		Img.WriteFDRImg(img_path, FDRImg, 0, sample_path)
		print 'Step: Draw FDR Curve for Sample ' + spectra_title;
		for j in range(1,len(mgflist)+1):
			FDRImg = 'fdr%d.png' % j
			print 'Step: Draw FDR Curve for raw %d of Sample '% j + spectra_title;
			Img.WriteFDRImg(img_path, FDRImg, j, sample_path)
			
		# PSM Draw
		for j in range(1,len(mgflist)+1):
			psm_path = os.path.join(sample_path, 'psm' )
			Path.Check(psm_path);

			#added at 2014.2.11, psm draw control, fan
			if not 'drawpsm' in initparam:
				initparam['drawpsm'] = 'true'
			DrawPSM = initparam['drawpsm']
			if DrawPSM != 'true':
				print 'Step: PSM Draw has been closed.'
				break
				
			print 'Step: Draw PSM for raw %d of Sample ' % j + spectra_title;
			psm_analysis_file = os.path.join(sample_path,'%s%d.inter.psm.output' % (spectra_title,j))

                        spectra_name_list = Img.DrawPSM(psm_path,psm_analysis_file)

                        
                        #modified at 2013.5.7 by fan, multithread
                        #annotated at 2013.5.11 by fan, multithread not supported in matplot
##                        cpu_cores = 1 # wait for modification
##
##                	psm_info_list_total = PSMDraw.LoadPSMData(psm_analysis_file); #modified at 2013.5.7 by fan, multithread
##                        spectra_name_list = []
##
##                        if cpu_cores == 1:
##                                spectra_name_list = thread.start_new_thread(Img.DrawPSM, (psm_path,psm_analysis_file, 0, cpu_cores, psm_info_list_total, ))
##                        else:
##                                for index in range(0, cpu_cores-1):
##                                        spectra_name = thread.start_new_thread(Img.DrawPSM, (psm_path,psm_analysis_file, index, cpu_cores, psm_info_list_total, ))
##                                        spectra_name_list.extend(spectra_name)
						
		# RunningTime
		RunningTimeImg = 'running.time.png'
		running_time_list = _GetRunningTime(index_path,normal_path)
		
		Img.WriteRunningTimeImg(img_path, RunningTimeImg, running_time_list)

		ReportFilePath = os.path.join(sample_path, 'report.html')
		print 'Step: Write Sample ' + ReportFilePath ;
		Template.WriteSample(ReportFilePath, i+1, sample_path, spectra_title, mgflist, RunningTimeImg, initparam)
		
		report_time = os.stat(ReportFilePath).st_mtime

	return report_time
