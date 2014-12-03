#!/usr/bin/env python
#-*-coding:utf8 -*-
import os
import sys
import time
import string

class mass_list_info:
	intensity = 0.0
	mz = 0.0
	iontype = 0
	iontag = '';

class psm_info:
	spectrum=''
	MH=0.0
	sequence=''
	site=''
	matchinfo=''
	mass_list = [];

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


def LoadPSMData(psm_info_file):
	psm_info_list = [];
	file=open(psm_info_file, 'r');
	psm_info_list = [];
	cur_psm = psm_info()
	while 1:
		str = file.readline();	
		if not str:
			break;
		
		str = str.strip();
		cur_line = str.split('=');
		if cur_line[0] == 'spectrum':
			if cur_psm.spectrum != '':
				psm_info_list.append(cur_psm);
			cur_psm = psm_info();
			cur_psm.mass_list = [];
			#cur_psm.spectrum = cur_line[1]; # modified at 2013.5.3
			cur_psm.spectrum = RmIllegalCharsForFile(cur_line[1]);
		elif cur_line[0] == 'MH':
			cur_psm.MH = cur_line[1];
		elif cur_line[0] == 'sequence':
			cur_psm.sequence = cur_line[1];
		elif cur_line[0] == 'site':
			cur_psm.site = cur_line[1];
		elif cur_line[0] == 'matchinfo':
			cur_psm.matchinfo = cur_line[1];
		else:
			cur_line = str.split('	');
			cur_mass_list = mass_list_info();
			if len(cur_line) >= 2:
				cur_mass_list.intensity = string.atof(cur_line[0])
				cur_mass_list.mz = string.atof(cur_line[1])
				
				if len(cur_line) == 2:
					cur_mass_list.iontag = '';
					cur_mass_list.iontype = 1;
				else:
					cur_mass_list.iontag = '|'.join(cur_line[2:len(cur_line)]);
					if cur_mass_list.iontag[0] == 'N':
						cur_mass_list.iontype = 2;
					elif cur_mass_list.iontag[0] == 'C':
						cur_mass_list.iontype = 4;
					elif cur_mass_list.iontag[0] == 'n':
						cur_mass_list.iontype = 3;
					elif cur_mass_list.iontag[0] == 'c':
						cur_mass_list.iontype = 5;
					elif cur_mass_list.iontag[0] == 'M':
						cur_mass_list.iontype = 6;
					elif cur_mass_list.iontag[0] == 'P':
						cur_mass_list.iontype = 7;
					elif cur_mass_list.iontag[0] == 'Q':
						cur_mass_list.iontype = 8;
					else:
						cur_mass_list.iontype = 1;
						cur_mass_list.iontag = '';
				if len(cur_mass_list.iontag) > 2:
					cur_mass_list.iontag = cur_mass_list.iontag[2:len(cur_mass_list.iontag)];
				cur_psm.mass_list.append(cur_mass_list);
				
	if cur_psm.spectrum != '':
		psm_info_list.append(cur_psm);

	file.close();
	return psm_info_list;
	
