#!/usr/bin/env python
#-*-coding:utf8 -*-
import os
import sys
import time
import string

class Sampleinfo:
	spectra_title=''
	spectra_path=''
	IDlist=[]
	mgf_path=''
	normal_path=''

	def __init__(self,spectra_title,spectra_path,IDlist,mgf_path,normal_path):
		self.spectra_title=spectra_title
		self.spectra_path=spectra_path
		self.IDlist=IDlist
		self.mgf_path=mgf_path
		self.normal_path=normal_path
		return 