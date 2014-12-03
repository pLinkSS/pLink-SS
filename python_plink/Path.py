#!/usr/bin/env python
#-*-coding:utf8 -*-
import os
import fnmatch
import shutil

def Copy(src_path,extend,des_path):
	if False==os.path.isdir(src_path):
		raise Exception,src_path
	if False == os.path.isdir(des_path):
		raise Exception,des_path
	dir_list = os.listdir(src_path)
	for d in dir_list:
		cur_path = os.path.join(src_path,d)
		if os.path.isfile(cur_path):
			if fnmatch.fnmatch(d.lower(),extend.lower()):
				shutil.copy(cur_path,des_path)
	return

def Check(output_path):
	if False==os.path.isdir(output_path):
		os.mkdir(output_path)
		if False==os.path.isdir(output_path):
			raise Exception,output_path
	return

def FileCheck(output_path):
	if False==os.path.isfile(output_path):
		raise Exception, "Error: failed to find %s" %output_path
	return

def CheckTime(output_path,input_time):
	run = False
	if os.path.isfile(output_path):
		t = os.stat(output_path).st_mtime
		if t < input_time:
			run = True
	else:
		run = True
	return run

def Walk(SourcePath,extend,filelist):
	if False==os.path.isdir(SourcePath):
		raise Exception,SourcePath

	dir_list = os.listdir( SourcePath )
	for d in dir_list:
		cur_path = os.path.join(SourcePath,d)
		if os.path.isdir(cur_path):
			Walk(cur_path,extend,filelist)
		else:
			if fnmatch.fnmatch(d,extend):
				filelist.append((SourcePath,d))
	return

def WriteTime(time_file,input_time,del_time):
	if CheckTime(time_file,input_time):
		tfile = open(time_file,'w')
		tfile.write('%f' % del_time)
		tfile.close()
	return 
