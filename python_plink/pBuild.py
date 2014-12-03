#!/usr/bin/env python
#-*-coding:utf8 -*-

import os
import sys
import string

def Write(pbuild_file,pBuild_Param,Sample_List):
        #2012.2.7 add the default param
        if not 'UseFDR' in pBuild_Param:
                pBuild_Param['UseFDR']=1
        if not 'ScoreMax' in pBuild_Param:
                pBuild_Param['ScoreMax']=10000
        if not 'ScoreMin' in pBuild_Param:
                pBuild_Param['ScoreMin']=0
        if not 'CrossLink' in pBuild_Param:
                pBuild_Param['CrossLink']=1
        
	file = open(pbuild_file,'w')
	file.write('SampleID=%s\n' % pBuild_Param['SampleID'])
	file.write('Total_Samples=%d\n' % len(Sample_List))
	sample_no = 0
	for SubItems in Sample_List:
		sample_no = sample_no + 1
		file.write('[Sample%d]\n' % sample_no)
		file.write('EngineType=pFind\n')
		file.write('SubItems=%d\n' % len(SubItems))
		for i in range(0,len(SubItems)):
			file.write('SubItem%d=%s,1\n' % (i+1,SubItems[i]))
	file.write('[END]\n')
	file.write('[DATABASE INFORMATION]\n')
	file.write('Fasta_File_Path=1,%s\n' % pBuild_Param['Fasta_File_Path'])
	file.write('Decoy_Tags=1,REVERSE_\n')
	file.write('[END]\n')
	file.write('[filtration]\n')
	file.write('UseFilter=1\n')
	file.write('PepMassLower=0\n')
	file.write('PepMassUpper=600000\n')
	file.write('LengthLower=0\n')
	file.write('LengthUpper=6000\n')
	file.write('Separate=0\n')
	file.write('PepTolBase=%s\n' % pBuild_Param['PepTolBase'])
	file.write('PepTolLower=%s\n' % pBuild_Param['PepTolLower'])
	file.write('PepTolUpper=%s\n' % pBuild_Param['PepTolUpper'])
	file.write('PepTolType=%s\n' % pBuild_Param['PepTolType'])
	file.write('ChargeState=1,2,3,4,5,6\n')
	file.write('RankLimit=1\n')
	file.write('Redundant=%s\n' % pBuild_Param['Redundant'])
	file.write('FixedDeltCn=1\n')
	file.write('DeltCnLower=0.1\n')
	file.write('UseFDR=%s\n' % pBuild_Param['UseFDR'])#2012.2.7 Modify
	file.write('FDR=0.05\n')
	file.write('FDRFormula=2\n')
	file.write('ScoreMin=%s\n' % pBuild_Param['ScoreMin'])#2012.2.7 Modify
	file.write('ScoreMax=%s\n' % pBuild_Param['ScoreMax'])#2012.2.7 Modify
	file.write('DistinctPep=%s\n' % pBuild_Param['DistinctPep'])
	file.write('DistinctSpec=1\n')
	file.write('CrossLink=%s\n' % pBuild_Param['CrossLink']) #2012.2.8 Modify
	file.write('LinkerID=0,1\n')
	file.write('LinkerType=4\n')
	file.write('ModSites=\n')
	file.write('CTerminal=\n')
	file.write('NTerminal=\n')
	file.write('[END]\n')
	file.write('OutPutPath=%s\n' % pBuild_Param['OutPutPath'])
	file.write('OutPutName=%s\n' % pBuild_Param['OutPutName'])
	file.write('ExportPath=%s\n' % pBuild_Param['ExportPath'])
	file.write('ExportFormat=%s\n' % pBuild_Param['ExportFormat'])
	file.write('ExportFile=0,1,2\n')
	file.write('[END]\n')
	file.close()
	return 
