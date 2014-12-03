#!/usr/bin/env python
# -*-coding:utf8 -*-
import os
import sys
import math
import string
import thread
import copy

import Path
import PSMDraw

import numpy
import matplotlib.pyplot as plt

### multithread drawing, annotated at 2013.5.11 since not supported in matplot
##def DrawPartPSM(psm_info_list, linecolor):
##        spectra_name_list = [];
##        for cur_psm in psm_info_list:
##		figfile = os.path.join(img_path, cur_psm.spectrum + '.png');
##		spectra_name_list.append(cur_psm.spectrum);
##		if os.path.isfile(figfile) == 0:
##			print 'Step : drawing %s' % cur_psm.spectrum
##			plt.clf();
##			# plt.getp(fig.patch)
##			mz_min = mz_max = 0;
##			if len(cur_psm.mass_list) > 0:
##				mz_min = cur_psm.mass_list[0].mz;
##				mz_max = mz_min;
##				for cur_peak in cur_psm.mass_list:
##					plt.vlines(cur_peak.mz, 0, cur_peak.intensity, color=linecolor[cur_peak.iontype - 1]);
##					if cur_peak.iontype != 1:
##						plt.text(cur_peak.mz, cur_peak.intensity + 0.025 * len(cur_peak.iontag), cur_peak.iontag , rotation='vertical', color=linecolor[cur_peak.iontype - 1]);
##					if mz_min > cur_peak.mz:
##						mz_min = cur_peak.mz;
##					if mz_max < cur_peak.mz:
##						mz_max = cur_peak.mz;
##					for cur_peak in cur_psm.mass_list:
##						if cur_peak.iontype != 1:
##							plt.vlines(cur_peak.mz, 0, cur_peak.intensity, color=linecolor[cur_peak.iontype - 1]);
##				titlestr = '\n'.join(['spectrum=%s' % cur_psm.spectrum, 'MH=%s' % cur_psm.MH, 'sequence=%s' % cur_psm.sequence, 'site=%s' % cur_psm.site, 'matchinfo=%s' % cur_psm.matchinfo]);
##				plt.title(titlestr)
##				plt.ylim(0, 1.2);
##				plt.xlim(mz_min - 10, mz_max + 10);
##				plt.xlabel('m/z')
##				plt.ylabel('intensity')
##				# plt.show()
##				plt.savefig(figfile)
##	return spectra_name_list


## def DrawPSM(img_path, psm_analysis_file, index, cpu_cores, psm_info_list_total):
def DrawPSM(img_path, psm_analysis_file):
	print 'Step: load psm info from ' + psm_analysis_file;
##	psm_info_list_total = PSMDraw.LoadPSMData(psm_analysis_file); #modified at 2013.5.7 by fan, multithread
	psm_info_list = PSMDraw.LoadPSMData(psm_analysis_file);

##        if cpu_cores == 1:
##                psm_info_list = copy.deepcopy(psm_info_list_total)
##	else:
##                if index == 0:
##                        psm_info_list = copy.deepcopy( psm_info_list_total[( len(psm_info_list_total)/(cpu_cores - 1)*(index) ): len(psm_info_list_total)/(cpu_cores - 1)*(index+1)] )#modified at 2013.5.7 by fan, multithread
##                else:
##                        psm_info_list = copy.deepcopy( psm_info_list_total[( len(psm_info_list_total)/(cpu_cores - 1)*(index) ) + 1: len(psm_info_list_total)/(cpu_cores - 1)*(index+1)] )#modified at 2013.5.7 by fan, multithread


	linecolor = ['k', 'g', 'b', 'r', 'm', 'y', 'c', 'c']
	spectra_name_list = [];
	fig = plt.figure(figsize=(12, 6));
	plt.subplots_adjust(left=0.05, bottom=0.12, right=0.98, top=0.8)
	print 'Step : begin draw psm ';

	for cur_psm in psm_info_list:
		figfile = os.path.join(img_path, cur_psm.spectrum + '.png');
		spectra_name_list.append(cur_psm.spectrum);
		if os.path.isfile(figfile) == 0:
			print 'Step : drawing %s' % cur_psm.spectrum
			plt.clf();
			# plt.getp(fig.patch)
			mz_min = mz_max = 0;
			if len(cur_psm.mass_list) > 0:
				mz_min = cur_psm.mass_list[0].mz;
				mz_max = mz_min;

				for cur_peak in cur_psm.mass_list:
					plt.vlines(cur_peak.mz, 0, cur_peak.intensity, color=linecolor[cur_peak.iontype - 1]);
					if cur_peak.iontype != 1:
						plt.text(cur_peak.mz, cur_peak.intensity + 0.025 * len(cur_peak.iontag), cur_peak.iontag , rotation='vertical', color=linecolor[cur_peak.iontype - 1]);
					if mz_min > cur_peak.mz:
						mz_min = cur_peak.mz;
					if mz_max < cur_peak.mz:
						mz_max = cur_peak.mz;
					for cur_peak in cur_psm.mass_list:
						if cur_peak.iontype != 1:
							plt.vlines(cur_peak.mz, 0, cur_peak.intensity, color=linecolor[cur_peak.iontype - 1]);
							
				titlestr = '\n'.join(['spectrum=%s' % cur_psm.spectrum, 'MH=%s' % cur_psm.MH, 'sequence=%s' % cur_psm.sequence, 'site=%s' % cur_psm.site, 'matchinfo=%s' % cur_psm.matchinfo]);
				plt.title(titlestr)
				plt.ylim(0, 1.2);
				plt.xlim(mz_min - 10, mz_max + 10);
				plt.xlabel('m/z')
				plt.ylabel('intensity')
				# plt.show()
				plt.savefig(figfile)

	return spectra_name_list;
	
def WriteAllPrecurImg(img_path, AllPrecurImg, error_type, error_file, report_path):
	plt.clf()
	fig = plt.figure()

	filename = os.path.join(report_path, error_file)
	if False == os.path.isfile(filename):
		raise Exception, filename

	file = open(filename, 'r')
	L = file.readlines()
	file.close()

	# F
	score1_list = []
	err1_list = []
	# T
	score0_list = []
	err0_list = []
	# U
	score2_list = []	
	err2_list = []

	for R in L:
		[str1, str2, str3, str4] = string.split(R, ' ')
		score = string.atof(str2)
		err = string.atof(str3)
		flag = string.atoi(str4)
		if 1 == flag:
			score1_list.append(score)
			err1_list.append(err)
		elif 0 == flag:
			score0_list.append(score)
			err0_list.append(err)
		else:
			score2_list.append(score)
			err2_list.append(err)

	logo = []
	if 0 == len(err0_list):
		plt.plot(err1_list, score1_list, linestyle='', marker='.', color='b')
		logo.append('target')
	else:
		plt.plot(err1_list, score1_list, linestyle='', marker='.', color='b')
		plt.plot(err2_list, score2_list, linestyle='', marker='+', color='g')
		plt.plot(err0_list, score0_list, linestyle='', marker='+', color='r')
		logo.append('T')
		logo.append('F')
		logo.append('U')

	plt.legend(logo, loc=1)
	plt.xlabel('precursor tolerance (%s)' % error_type)
	plt.ylabel('pFind score')

	AllPrecurImgPath = os.path.join(img_path, AllPrecurImg)
	plt.savefig(AllPrecurImgPath)
	return

def WriteFDRImg(img_path, FDRImg, raw_no, report_path):
	plt.clf()
	fig = plt.figure()
	q_lim = 1
	num_max = 0
	logo = []
	
	n2 = 81;
	for i in range(1, 7):
		n1 = 2 ** i;
		tmpn = n1 * n2;
		filename = os.path.join(report_path, '%d.Sample1_Charge%d.txt' % (raw_no, tmpn))
		if False == os.path.isfile(filename):
			continue

		logo.append('chg%d' % i)

		file = open(filename, 'r')
		L = file.readlines()
		file.close()

		no_list = []
		q_list = []
		for R in L:
			[str1, str2] = string.split(R, ' ')
			no = string.atoi(str1)
			qvalue = string.atof(str2)
			if qvalue < q_lim:
				q_list.append(qvalue * 100.0)
				no_list.append(no)
				if num_max < no:
					num_max = no
			else:
				break

		plt.plot(q_list, no_list, linewidth=2)

	plt.legend(logo, loc=4)
	plt.title('FDR curve')
	plt.xlabel('FDR (%)')
	plt.ylabel('Number of spectra')

	plt.plot([5, 5], [0, num_max], linestyle='-.', linewidth=2, color='k')

	FDRImgPath = os.path.join(img_path, FDRImg)
	plt.savefig(FDRImgPath)

	return

def WriteRunningTimeImg(img_path, RunningTimeFile, running_time_list):
	totaltime = 0.000001
	for time in running_time_list:
		totaltime = totaltime + time

	fracs = []
	for i in range(0, len(running_time_list)):
		fracs.append(running_time_list[i] * 100 / totaltime)

	plt.clf()
	plt.figure(figsize=(8, 8))
	explode = (0.03, 0.03)
	labels = '0.make index', '1.database search',

	plt.pie(fracs, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True)
	plt.title('Running time')
	RunningTimeImgPath = os.path.join(img_path, RunningTimeFile)
	plt.savefig(RunningTimeImgPath)

	return 
