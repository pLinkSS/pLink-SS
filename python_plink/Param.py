#!/usr/bin/env python
#-*-coding:utf8 -*-
import os
import string 
import ConfigParser

def ReadAA(initparam):
	bin_path = initparam['bin.path']
	aa_ini = os.path.join(bin_path,'aa.ini')
	if False==os.path.isfile(aa_ini):
		raise Exception,aa_ini

	cfg = ConfigParser.ConfigParser()
	cfg.read(aa_ini)
	v = cfg.items('aa')

	dic={}
	for x in v:
		dic[x[0]] = x[1]
	return dic

def ReadModify(initparam):
	bin_path = initparam['bin.path']
	modify_ini = os.path.join(bin_path,'modify.ini')
	if False==os.path.isfile(modify_ini):
		raise Exception,modify_ini

	cfg = ConfigParser.ConfigParser()
	cfg.read(modify_ini)
	v = cfg.items('modify')

	tmpdic = {}
	for x in v:
		tmpdic[x[0]] = x[1]

	modify_dic={}
	total = string.atoi(tmpdic['total'])
	for i in range (1, total+1):
		key = 'name%d' % i
		modifyname = tmpdic[key]
		modifyline = tmpdic[string.lower(modifyname)]
		modify_dic[modifyname] = modifyline
	return modify_dic

def ReadLinker(initparam):
	bin_path = initparam['bin.path']
	linker_ini = os.path.join(bin_path,'xlink.ini')
	if False==os.path.isfile(linker_ini):
		raise Exception,linker_ini

	cfg = ConfigParser.ConfigParser()
	cfg.read(linker_ini)
	v = cfg.items('xlink')

	tmpdic = {}
	for x in v:
		tmpdic[x[0]] = x[1]

	linker_dic={}
	total = string.atoi(tmpdic['total'])
	for i in range (1, total+1):
		key = 'name%d' % i
		linkername = tmpdic[key]
		linkerline = tmpdic[string.lower(linkername)]
		linker_dic[linkername] = linkerline
	return linker_dic
	
def ReadInstrument(initparam):
	spectra_instrument = initparam['spectra.instrument']
	bin_path = initparam['bin.path']
	instrument_ini = os.path.join(bin_path,'instrument.ini')

	if False==os.path.isfile(instrument_ini):
		raise Exception,instrument_ini

	cfg = ConfigParser.ConfigParser()
	cfg.read(instrument_ini)
	v = cfg.items(spectra_instrument)

	return v

def Read(parafile):
	if False==os.path.isfile(parafile):
		raise Exception,parafile

	cfg = ConfigParser.ConfigParser()
	cfg.read(parafile)
	v = cfg.items('pLink')

	dic = {'para.path':parafile}
	for x in v:
		dic[x[0]] = x[1]

	return dic

if __name__ == '__main__':
	Read('test.ini')
