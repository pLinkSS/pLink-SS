#!/usr/bin/env python
# -*-coding:utf8 -*-
'''
Created on 23 Jan 2013

@author: Jeremy
'''
import os
import Path
import platform

def _ConstructFilterParam(pfind_file_path, input_file_path, output_dir_path, extra_param):
    filterParam = {}
    # IO
    _io = []
    _io.append(('pFindFile', pfind_file_path))
    _io.append(('InputReportFile', input_file_path))
    _io.append(('OutputPath', output_dir_path))
    _io.append(('title', extra_param['title']))
    filterParam['IO'] = _io
    
    # filter
    _filter = []
    _filter.append(('FilterMethod', '0'))
    _filter.append(('XLinkType', extra_param['XLinkType']))
    _filter.append(('LinkerId', '-1'))
    _filter.append(('PepTolType', 'Da'))
    _filter.append(('PepTol', '100'))
    _filter.append(('PepTolBaseType', 'Da'))
    _filter.append(('PepTolBaseTotal', '1'))
    _filter.append(('PepTolBase0', '0'))
    _filter.append(('MaxEvalue', '1'))
    _filter.append(('MinScore', '0'))
    _filter.append(('MinPepLength', '0'))
    _filter.append(('MinContiAAnum', '0'))
    _filter.append(('AAconfLevel', '1'))
    _filter.append(('MinMatchOverUnMatch', '0.000000'))
    _filter.append(('SaveCTerm', '1'))
    _filter.append(('FDR', '-1'))
    _filter.append(('ReverseTag', 'REVERSE_'))
    filterParam['filter'] = _filter
    return filterParam

def _WriteFilter(filterParam, filter_file_path):
    pfile = open(filter_file_path,'w')
    keys = ['IO', 'filter']
    for key in keys:
        pfile.write('['+key+']\n')
        for elem in filterParam[key]:
            pfile.write(elem[0]+'='+elem[1]+'\n')
    pfile.close()
    return

def _ConstructConvertParam(pfind_file_path, report_file_path, output_dir_path, _title, _format):
    convertParam = {}
    # IO
    _io = []
    _io.append(('pFindFile', pfind_file_path))
    _io.append(('ReportFile', report_file_path))
    _io.append(('OutputPath', output_dir_path))
    _io.append(('title', _title))
    _io.append(('format', _format))
    convertParam['IO'] = _io
    return convertParam

def _WriteConvert(convertParam, convert_file_path):
    pfile = open(convert_file_path,'w')
    keys = ['IO']
    for key in keys:
        pfile.write('['+key+']\n')
        for elem in convertParam[key]:
            pfile.write(elem[0]+'='+elem[1]+'\n')
    pfile.close()
    return

def _ConstructpBuildParam(initparam, sample_list, output_dir_path, extra_param):
    pBuildParam = {}

    #2012.2.7 add the default param
    if not 'UseFDR' in extra_param:
            extra_param['UseFDR']='1'
    if not 'ScoreMax' in extra_param:
            extra_param['ScoreMax']='10000'
    if not 'ScoreMin' in extra_param:
            extra_param['ScoreMin']='0'
    if not 'CrossLink' in extra_param:
            extra_param['CrossLink']='1'

    # samplepBuild_Param
    samples = []
    samples.append(('SampleID', extra_param['SampleID']))
    samples.append(('Total_Samples', '%d' % len(sample_list)))
    for subItems in sample_list:
        samples.append(('EngineType', 'pFind'))
        samples.append(('SubItems', '%d' % len(subItems)))
        for i in range(0, len(subItems)):
            samples.append(('SubItem%d' % (i+1), '%s,1' % subItems[i]))
    pBuildParam['samples'] = samples
    
    # database information
    databaseInformation = []
    databaseInformation.append(('Fasta_File_Path', '1,%s' % initparam['database.path']))
    databaseInformation.append(('Decoy_Tags', '1,REVERSE_'))
    pBuildParam['DATABASE INFORMATION'] = databaseInformation
    
    # filtration
    filtration = []
    filtration.append(('UseFilter', '1'))
    filtration.append(('PepMassLower', '0'))
    filtration.append(('PepMassUpper', '600000'))
    filtration.append(('LengthLower', '0'))
    filtration.append(('LengthUpper', '6000'))
    filtration.append(('Separate', '0'))
    filtration.append(('PepTolBase', extra_param['PepTolBase']))
    filtration.append(('PepTolLower', extra_param['PepTolLower']))
    filtration.append(('PepTolUpper', extra_param['PepTolUpper']))
    filtration.append(('PepTolType', extra_param['PepTolType']))
    filtration.append(('ChargeState', '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15'))
    filtration.append(('RankLimit', '1'))
    filtration.append(('Redundant', '1'))
    filtration.append(('FixedDeltCn', '1'))
    filtration.append(('DeltCnLower', '0.1'))
    if not 'usefdr' in extra_param:
        filtration.append(('UseFDR', '1'))
    else:
        filtration.append(('UseFDR', '0'))
    if not 'fdr' in initparam: # added at 2013.10.12
        filtration.append(('FDR', '0.05'))
    else:
        filtration.append(('FDR', initparam['fdr']))
#     filtration.append(('FDR', '0.05'))
    filtration.append(('FDRFormula', '2'))
    filtration.append(('ScoreMin', '0'))
    if not 'evalue_max' in initparam:
        filtration.append(('ScoreMax', '10000'))
    else:
        filtration.append(('ScoreMax', initparam['evalue_max']))
    filtration.append(('DistinctPep', '1'))
    filtration.append(('DistinctSpec', '1'))
    filtration.append(('CrossLink', extra_param['CrossLink']))
    filtration.append(('LinkerID', '0,1'))
    filtration.append(('LinkerType', '4'))
    filtration.append(('ModSites', ''))
    filtration.append(('CTerminal', ''))
    filtration.append(('NTerminal', ''))
    pBuildParam['filtration'] = filtration
    
    # output
    output = []
    output.append(('OutPutPath', output_dir_path))
    output.append(('OutPutName', extra_param['OutPutName']))
    output.append(('ExportPath', output_dir_path))
    output.append(('ExportFormat', '0,2,3'))
    output.append(('ExportFile', '0,1,2'))
    pBuildParam['output'] = output
    return pBuildParam;

def _WritepBuild(pBuildParam, sample_list, pbuild_file_path):
    pfile = open(pbuild_file_path,'w')
    key = 'samples'
    vec = pBuildParam[key]
    for i in range(0, 2):
        elem = vec[i]
        pfile.write(elem[0]+'='+elem[1]+'\n')
    for i in range(0, len(sample_list)):
        pfile.write('[Sample%d]\n' % (i+1))
        for j in range(2, 4):
            elem = vec[j]
            pfile.write(elem[0]+'='+elem[1]+'\n')
        subItems = sample_list[i]
        for j in range(0, len(subItems)):
            elem = vec[j+4]
            pfile.write(elem[0]+'='+elem[1]+'\n')
    pfile.write('[END]\n')
    
    keys = ['DATABASE INFORMATION', 'filtration']
    for key in keys:
        pfile.write('['+key+']\n')
        for elem in pBuildParam[key]:
            pfile.write(elem[0]+'='+elem[1]+'\n')
        pfile.write('[END]\n')
    
    key = 'output'
    for elem in pBuildParam[key]:
        pfile.write(elem[0]+'='+elem[1]+'\n')
    pfile.write('[END]\n')
    pfile.close()
    return

def _ConstructPSMParam(initparam, pfind_file_path, report_file_path, output_dir_path):
    psmParam = {}
    # pfd
    pfd = []
    pfd.append(('pFindFile', pfind_file_path))
    pfd.append(('ReportFile', report_file_path))
    pfd.append(('InputType', 'pbuild'))
    pfd.append(('OutputPath', output_dir_path))
    psmParam['pfd'] = pfd
    
    # ion
    ion = []
    if initparam['spectra.instrument'] == 'HCD':
        ion.append(('tolerance', '20'))
        ion.append(('tolerance_type', 'ppm'))
    else:
        ion.append(('tolerance', '0.5'))
        ion.append(('tolerance_type', 'Da'))
    if initparam['spectra.instrument']=='ETD':
        ion.append(('ion_type_total', '8'))
        ion.append(('ion_type1', 'z 1 0 0 0'))
        ion.append(('ion_type2', 'c 1 0 0 0'))
        ion.append(('ion_type3', 'z 2 0 0 0'))
        ion.append(('ion_type4', 'c 2 0 0 0'))
        ion.append(('ion_type5', 'z 1 0 0 -1.007825'))
        ion.append(('ion_type6', 'c 1 0 0 1.007825'))
        ion.append(('ion_type7', 'z 2 0 0 -1.007825'))
        ion.append(('ion_type8', 'c 2 0 0 1.007825'))
    elif initparam['spectra.instrument']=='HCD':
        ion.append(('ion_type_total', '14'))
        ion.append(('ion_type1', 'b 1 0 0 0'))
        ion.append(('ion_type2', 'y 1 0 0 0'))
        ion.append(('ion_type3', 'b 2 0 0 0'))
        ion.append(('ion_type4', 'y 2 0 0 0'))
        ion.append(('ion_type5', 'b 3 0 0 0'))
        ion.append(('ion_type6', 'y 3 0 0 0'))
        ion.append(('ion_type7', 'a 1 0 0 0'))
        ion.append(('ion_type8', 'a 2 0 0 0'))
        ion.append(('ion_type9', 'a 3 0 0 0'))
        ion.append(('ion_type10', 'q 1 0 0 0'))                  
        ion.append(('ion_type11', 'q 2 0 0 0'))
        ion.append(('ion_type12', 'q 3 0 0 0'))
        ion.append(('ion_type13', 'q 4 0 0 0'))
        ion.append(('ion_type14', 'q 5 0 0 0'))
    else:
        ion.append(('ion_type_total', '6'))
        ion.append(('ion_type1', 'b 1 0 0 0'))
        ion.append(('ion_type2', 'y 1 0 0 0'))
        ion.append(('ion_type3', 'b 2 0 0 0'))
        ion.append(('ion_type4', 'y 2 0 0 0'))
        ion.append(('ion_type5', 'b 3 0 0 0'))
        ion.append(('ion_type6', 'y 3 0 0 0'))
    psmParam['ion'] = ion
    
    # gap
    gap = []
    gap.append(('mass_scope', '100'))
    psmParam['gap'] = gap
    
    # filter
    _filter = []
    _filter.append(('RemoveNoise', '0'))
    psmParam['filter'] = _filter
    return psmParam

def _WritePSM(psmParam, psm_file_path):
    pfile = open(psm_file_path,'w')
    keys = ['pfd', 'ion', 'gap', 'filter']
    for key in keys:
        pfile.write('['+key+']\n')
        for elem in psmParam[key]:
            pfile.write(elem[0]+'='+elem[1]+'\n')
    pfile.close()
    return

def RunForType(bfile, xlinktype, initparam, spectra_list, filter_dir_path, bEvalueMax): 
    if xlinktype == 'inter':
        typeno = '3'
    elif xlinktype == 'loop':
        typeno = '2'
    elif xlinktype == 'mono':
        typeno = '1'
    elif xlinktype == 'common':
        typeno = '0'
    else:
        print 'Error: not the right xlink type.\n'
        return
    
    Path.Check(filter_dir_path)
    bin_path = initparam['bin.path']
    spectra_title = initparam['spectra.title']
    
    # extract specific linked spectra by XLinkResultFilter, .filter
    print 'Step : Extract %s-linked spectra' % xlinktype
    program_path = os.path.join(bin_path, 'XLinkResultFilter')
    for i in range(0, len(spectra_list)):
        pfind_file_path = os.path.join(filter_dir_path, '%s%d.pfind' % (spectra_title, (i+1)))
        filter_file_path = os.path.join(filter_dir_path, '%s%d.%s.filter' % (spectra_title, (i+1), xlinktype))
        filter_input_file_path = os.path.join(filter_dir_path, '%s%d_qry.proteins.txt' % (spectra_title, (i+1)))
        filter_result_file_path = os.path.join(filter_dir_path, '%s%d_%s_qry.proteins.txt' % (spectra_title, (i+1), xlinktype))
        if os.path.isfile(filter_result_file_path):
            print filter_result_file_path + ' is existed, skip the step'
        else:
            print 'Filter of ' + filter_input_file_path
            extraParam = {}
            extraParam['title'] = xlinktype
            extraParam['XLinkType'] = typeno
            filterParam = _ConstructFilterParam(pfind_file_path, filter_input_file_path, filter_dir_path, extraParam)
            _WriteFilter(filterParam, filter_file_path)
            bfile.write('"%s" "%s"\n' % (program_path, filter_file_path))
            
    # convert inter-linked spectra by XLinkResultConverter, .convert, .pbuild.txt
    print 'Step : Convert to *.pbuild by XLinkResultConverter'
    program_path = os.path.join(bin_path, 'XLinkResultConverter')
    for i in range(0, len(spectra_list)):
        pfind_file_path = os.path.join(filter_dir_path, '%s%d.pfind' % (spectra_title, (i+1)))
        convert_file_path = os.path.join(filter_dir_path, '%s%d.%s.convert' % (spectra_title, (i+1), xlinktype))
        convert_input_file_path = os.path.join(filter_dir_path, '%s%d_%s_qry.proteins.txt' % (spectra_title, (i+1), xlinktype))
        convert_result_file_path = os.path.join(filter_dir_path, '%s%d.%s.pbuild.txt' % (spectra_title, (i+1), xlinktype))
        if os.path.isfile(filter_result_file_path):
            print convert_result_file_path + ' is existed, skip the step'
        else:
            print 'Filter of ' + convert_input_file_path
            filterParam = _ConstructConvertParam(pfind_file_path, convert_input_file_path, filter_dir_path, xlinktype, 'pbuild')
            _WriteConvert(filterParam, convert_file_path)
            bfile.write('"%s" "%s"\n' % (program_path, convert_file_path))
    
    # convert inter-linked spectra by XLinkResultConverter, .convert, .pxbuild.txt
    print 'Step : Convert to *.pXBuild by XLinkResultConverter'
    for i in range(0, len(spectra_list)):
        pfind_file_path = os.path.join(filter_dir_path, '%s%d.pfind' % (spectra_title, (i+1)))
        convert_file_path = os.path.join(filter_dir_path, '%s%dx.%s.convert' % (spectra_title, (i+1), xlinktype))
        convert_input_file_path = os.path.join(filter_dir_path, '%s%d_%s_qry.proteins.txt' % (spectra_title, (i+1), xlinktype))
        convert_result_file_path = os.path.join(filter_dir_path, '%s%d.%s.pXbuild' % (spectra_title, (i+1), xlinktype))
        if os.path.isfile(filter_result_file_path):
            print convert_result_file_path + ' is existed, skip the step'
        else:
            print 'Convert of ' + convert_input_file_path
            filterParam = _ConstructConvertParam(pfind_file_path, convert_input_file_path, filter_dir_path, xlinktype, 'pxbuild')
            _WriteConvert(filterParam, convert_file_path)
            bfile.write('"%s" "%s"\n' % (program_path, convert_file_path))
            
    # filter inter-linked spectra by builder, .pbuild, one input one output
    print 'Step : Filter by pBuild '
    program_path = os.path.join(bin_path, 'builder')
    samples = []
    for i in range(0, len(spectra_list)):
        builder_input_file_path = os.path.join(filter_dir_path, '%s%d_%s.pbuild.txt' % (spectra_title, (i+1), xlinktype))
        samples.append(builder_input_file_path)
        builder_result_file_path = os.path.join(filter_dir_path, '%s_%d_%s' % (spectra_title, (i+1), xlinktype), '%s_%d_%s_DataSet1_pfind.txt' % (spectra_title,(i+1), xlinktype))
        pbuild_file_path = os.path.join(filter_dir_path, '%d_%s.filter.pbuild' % ((i+1), xlinktype))
        if os.path.isfile(builder_result_file_path):
            print builder_result_file_path + ' is existed, skip the step'
        else:
            print 'Filter of ' + builder_input_file_path
            extraParam = {}
            extraParam['SampleID'] = '%s_%d_%s' % (spectra_title, (i+1), xlinktype)
            extraParam['OutPutName'] = '%s_%d_%s' % (spectra_title, (i+1), xlinktype)
            extraParam['PepTolBase'] = initparam['filter_peptide_tol_base']
            extraParam['PepTolLower'] = initparam['filter_peptide_tol_lb']
            extraParam['PepTolUpper'] = initparam['filter_peptide_tol_ub']
            extraParam['PepTolType'] = initparam['filter_peptide_tol_type']
            if xlinktype == 'inter':
                if not initparam['index_content'] == "PEPTIDE_TRI_ALL":
                    extraParam['CrossLink'] = '1'
                else:
                    extraParam['CrossLink'] = '2'
            else:
                extraParam['CrossLink'] = '0'
            pBuildParam = _ConstructpBuildParam(initparam, [[builder_input_file_path]], filter_dir_path, extraParam)
            _WritepBuild(pBuildParam, [[builder_input_file_path]], pbuild_file_path)
            bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path))
            
            #Evalue cut, refilter, 2013.10.22
            if (bEvalueMax):
                if platform.system() == 'Linux':
                        bfile.write('cp -b "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
                        bfile.write('rm "%s"\n' % builder_result_file_path)
                elif platform.system() == 'Windows':
                        bfile.write('move "%s" "%s"\n' % (builder_input_file_path, builder_input_file_path+"~1"))
                        bfile.write('move "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
                print 'Evalue filter of ' + builder_input_file_path

                extraParam['usefdr'] = 0
                pbuild_file_path_two = pbuild_file_path[0:-7]+"2.pbuild"
                pBuildParam = _ConstructpBuildParam(initparam, [[builder_input_file_path]], filter_dir_path, extraParam)
                _WritepBuild(pBuildParam, [[builder_input_file_path]], pbuild_file_path_two)
                bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path_two))

    
    # filter inter-linked spectra by builder, .pbuild, all inputs one output
    print 'Step: Filter all by pBuild '
    pbuild_file_path = os.path.join(filter_dir_path, 'total_%s.filter.pbuild' % xlinktype)
    builder_result_file_path = os.path.join(filter_dir_path, spectra_title+'_'+xlinktype, '%s_%s_DataSet1_pfind.txt' % (spectra_title, xlinktype))
    if os.path.isfile(builder_result_file_path):
        print builder_result_file_path + ' is existed, skip the step'
    else:
        print 'Filter of ' + pbuild_file_path
        extraParam = {}
        extraParam['SampleID'] = '%s_%s' % (spectra_title, xlinktype)
        extraParam['OutPutName'] = '%s_%s' % (spectra_title, xlinktype)
        extraParam['PepTolBase'] = initparam['filter_peptide_tol_base']
        extraParam['PepTolLower'] = initparam['filter_peptide_tol_lb']
        extraParam['PepTolUpper'] = initparam['filter_peptide_tol_ub']
        extraParam['PepTolType'] = initparam['filter_peptide_tol_type']
        if xlinktype == 'inter':
            extraParam['CrossLink'] = '1'
        else:
            extraParam['CrossLink'] = '0'
        pBuildParam = _ConstructpBuildParam(initparam, [samples], filter_dir_path, extraParam)
        _WritepBuild(pBuildParam, [samples], pbuild_file_path)
        bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path))
        #Evalue cut, refilter, 2013.10.22
        if (bEvalueMax):
            if platform.system() == 'Linux':
                    bfile.write('cp -b "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
                    bfile.write('rm "%s"\n' % builder_result_file_path)
            elif platform.system() == 'Windows':
                    bfile.write('move "%s" "%s"\n' % (builder_input_file_path, builder_input_file_path+"~2"))
                    bfile.write('move "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
            print 'Evalue filter of ' + builder_input_file_path

            extraParam['usefdr'] = 0
            pbuild_file_path_two = pbuild_file_path[0:-7]+"2.pbuild"
            pBuildParam = _ConstructpBuildParam(initparam, [[builder_input_file_path]], filter_dir_path, extraParam)
            _WritepBuild(pBuildParam, [[builder_input_file_path]], pbuild_file_path_two)
            bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path_two))
    
    # filter inter-linked spectra by builder, .pbuild, all inputs one output
    print 'Step: Combine by pBuild '
    samples = []
    for i in range(0, len(spectra_list)):
        middle_file_path = os.path.join(filter_dir_path, '%s_%d_%s' % (spectra_title, (i+1), xlinktype), '%s_%d_%s_DataSet1_pfind.txt' % (spectra_title,(i+1), xlinktype))
        samples.append(middle_file_path)
    pbuild_file_path = os.path.join(filter_dir_path, 'total_%s.combine.pbuild' % xlinktype)
    pbuild_result_file_path = os.path.join(filter_dir_path, '%s_%s_combine' % (spectra_title, xlinktype), '%s_%s_combine_DataSet1_pfind.txt' % (spectra_title, xlinktype))
    if os.path.isfile(pbuild_result_file_path):
        print builder_result_file_path + ' is existed, skip the step'
    else:
        print 'Filter of ' + pbuild_file_path
        extraParam = {}
        extraParam['SampleID'] = '%s_%s_combine' % (spectra_title, xlinktype)
        extraParam['OutPutName'] = '%s_%s_combine' % (spectra_title, xlinktype)
        extraParam['PepTolBase'] = '0' 
        extraParam['PepTolLower'] = '-100'
        extraParam['PepTolUpper'] = '100'
        extraParam['PepTolType'] = 'Da'
        if xlinktype == 'inter':
            extraParam['CrossLink'] = '1'
        else:
            extraParam['CrossLink'] = '0'
        pBuildParam = _ConstructpBuildParam(initparam, [samples], filter_dir_path, extraParam)
        _WritepBuild(pBuildParam, [samples], pbuild_file_path)
        bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path))
        #Evalue cut, refilter, 2013.10.22
        if (bEvalueMax):
            if platform.system() == 'Linux':
                    bfile.write('cp -b "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
                    bfile.write('rm "%s"\n' % builder_result_file_path)
            elif platform.system() == 'Windows':
                    bfile.write('move "%s" "%s"\n' % (builder_input_file_path, builder_input_file_path+"~3"))
                    bfile.write('move "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
            print 'Evalue filter of ' + builder_input_file_path

            extraParam['usefdr'] = 0
            pbuild_file_path_two = pbuild_file_path[0:-7]+"2.pbuild"
            pBuildParam = _ConstructpBuildParam(initparam, [[builder_input_file_path]], filter_dir_path, extraParam)
            _WritepBuild(pBuildParam, [[builder_input_file_path]], pbuild_file_path_two)
            bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path_two))
        
    
    
    # PSMAnalysis
    print 'Step : Generate PSM file by PSMAnalysis '
    program_path = os.path.join(bin_path, 'PSMAnalysis')
    for i in range(0, len(spectra_list)):
        pfind_file_path = os.path.join(filter_dir_path, '%s%d.pfind' % (spectra_title, (i+1)))
        psm_file_path = os.path.join(filter_dir_path, '%s%d.%s.psm' % (spectra_title, (i+1), xlinktype))
        psm_input_file_path = os.path.join(filter_dir_path, '%s_%d_%s' % (spectra_title, (i+1), xlinktype), '%s_%d_%s_DataSet1_pfind.txt' % (spectra_title, (i+1), xlinktype))
        psm_result_file_path = os.path.join(filter_dir_path, '%s%d.%s.psm.output' % (spectra_title, (i+1), xlinktype))
        if os.path.isfile(psm_result_file_path):
            print psm_result_file_path + ' is existed, skip the step'
        else:
            print 'Gen PSM file of ' + psm_file_path
            psmParam = _ConstructPSMParam(initparam, pfind_file_path, psm_input_file_path, filter_dir_path)
            _WritePSM(psmParam, psm_file_path)
            bfile.write('"%s" "%s"\n' % (program_path, psm_file_path))
    
    # plabel
    print 'Step : Convert to *.plabel by XLinkResultConverter'
    program_path = os.path.join(bin_path, 'XLinkResultConverter')
    for i in range(0, len(spectra_list)):
        pfind_file_path = os.path.join(filter_dir_path, '%s%d.pfind' % (spectra_title, (i+1)))
        convert_file_path = os.path.join(filter_dir_path, '%s%d.%s.2plabel.convert' % (spectra_title, (i+1), xlinktype))
        convert_input_file_path = os.path.join(filter_dir_path, '%s_%d_%s' % (spectra_title, (i+1), xlinktype), '%s_%d_%s_DataSet1_pfind.txt' % (spectra_title, (i+1), xlinktype))
        convert_result_file_path = os.path.join(filter_dir_path, '%s%d.%s.plabel' % (spectra_title, (i+1), xlinktype))
        if os.path.isfile(convert_result_file_path):
            print convert_result_file_path + ' is existed, skip the step'
        else:
            print 'Convert of ' + convert_file_path
            convertParam = _ConstructConvertParam(pfind_file_path, convert_input_file_path, filter_dir_path, xlinktype, 'plabel')
            convertParam['IO'].append(('InputType', 'pbuild'))
            _WriteConvert(convertParam, convert_file_path)
            bfile.write('"%s" "%s"\n' % (program_path, convert_file_path))
    
    return

def RunForTri(bfile, initparam, spectra_list, filter_dir_path, bEvalueMax): 
    xlinktype = 'tri'
    
    Path.Check(filter_dir_path)
    bin_path = initparam['bin.path']
    spectra_title = initparam['spectra.title']
    
    # filter tri-linked spectra by builder, .pbuild, one input one output
    print 'Step : Filter by pBuild '
    program_path = os.path.join(bin_path, 'builder')
    samples = []
    for i in range(0, len(spectra_list)):
        builder_input_file_path = os.path.join(filter_dir_path, '%s%d_qry.proteins.txt' % (spectra_title, (i+1)))
        samples.append(builder_input_file_path)
        builder_result_file_path = os.path.join(filter_dir_path, '%s_%d_%s' % (spectra_title, (i+1), xlinktype), '%s_%d_%s_DataSet1_pfind.txt' % (spectra_title,(i+1), xlinktype))
        pbuild_file_path = os.path.join(filter_dir_path, '%d_%s.filter.pbuild' % ((i+1), xlinktype))
        if os.path.isfile(builder_result_file_path):
            print builder_result_file_path + ' is existed, skip the step'
        else:
            print 'Filter of ' + builder_input_file_path
            extraParam = {}
            extraParam['SampleID'] = '%s_%d_%s' % (spectra_title, (i+1), xlinktype)
            extraParam['OutPutName'] = '%s_%d_%s' % (spectra_title, (i+1), xlinktype)
            extraParam['PepTolBase'] = initparam['filter_peptide_tol_base']
            extraParam['PepTolLower'] = initparam['filter_peptide_tol_lb']
            extraParam['PepTolUpper'] = initparam['filter_peptide_tol_ub']
            extraParam['PepTolType'] = initparam['filter_peptide_tol_type']
            extraParam['CrossLink'] = '2'
            
            pBuildParam = _ConstructpBuildParam(initparam, [[builder_input_file_path]], filter_dir_path, extraParam)
            _WritepBuild(pBuildParam, [[builder_input_file_path]], pbuild_file_path)
            bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path))
            #Evalue cut, refilter, 2013.10.22
            if (bEvalueMax):
                if platform.system() == 'Linux':
                        bfile.write('cp -b "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
                        bfile.write('rm "%s"\n' % builder_result_file_path)
                elif platform.system() == 'Windows':
                        bfile.write('move "%s" "%s"\n' % (builder_input_file_path, builder_input_file_path+"~1"))
                        bfile.write('move "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
                print 'Evalue filter of ' + builder_input_file_path
    
                extraParam['usefdr'] = 0
                pbuild_file_path_two = pbuild_file_path[0:-7]+"2.pbuild"
                pBuildParam = _ConstructpBuildParam(initparam, [[builder_input_file_path]], filter_dir_path, extraParam)
                _WritepBuild(pBuildParam, [[builder_input_file_path]], pbuild_file_path_two)
                bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path_two))
    
    # filter tri-linked spectra by builder, .pbuild, all inputs one output
    print 'Step: Filter all by pBuild '
    pbuild_file_path = os.path.join(filter_dir_path, 'total_%s.filter.pbuild' % xlinktype)
    builder_result_file_path = os.path.join(filter_dir_path, spectra_title+'_'+xlinktype, '%s_%s_DataSet1_pfind.txt' % (spectra_title, xlinktype))
    if os.path.isfile(builder_result_file_path):
        print builder_result_file_path + ' is existed, skip the step'
    else:
        print 'Filter of ' + pbuild_file_path
        extraParam = {}
        extraParam['SampleID'] = '%s_%s' % (spectra_title, xlinktype)
        extraParam['OutPutName'] = '%s_%s' % (spectra_title, xlinktype)
        extraParam['PepTolBase'] = initparam['filter_peptide_tol_base']
        extraParam['PepTolLower'] = initparam['filter_peptide_tol_lb']
        extraParam['PepTolUpper'] = initparam['filter_peptide_tol_ub']
        extraParam['PepTolType'] = initparam['filter_peptide_tol_type']
        extraParam['CrossLink'] = '2'
        pBuildParam = _ConstructpBuildParam(initparam, [samples], filter_dir_path, extraParam)
        _WritepBuild(pBuildParam, [samples], pbuild_file_path)
        bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path))
        #Evalue cut, refilter, 2013.10.22
        if (bEvalueMax):
            if platform.system() == 'Linux':
                    bfile.write('cp -b "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
                    bfile.write('rm "%s"\n' % builder_result_file_path)
            elif platform.system() == 'Windows':
                    bfile.write('move "%s" "%s"\n' % (builder_input_file_path, builder_input_file_path+"~2"))
                    bfile.write('move "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
            print 'Evalue filter of ' + builder_input_file_path

            extraParam['usefdr'] = 0
            pbuild_file_path_two = pbuild_file_path[0:-7]+"2.pbuild"
            pBuildParam = _ConstructpBuildParam(initparam, [[builder_input_file_path]], filter_dir_path, extraParam)
            _WritepBuild(pBuildParam, [[builder_input_file_path]], pbuild_file_path_two)
            bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path_two))
    
    # filter tri-linked spectra by builder, .pbuild, all inputs one output
    print 'Step: Combine by pBuild '
    samples = []
    for i in range(0, len(spectra_list)):
        middle_file_path = os.path.join(filter_dir_path, '%s_%d_%s' % (spectra_title, (i+1), xlinktype), '%s_%d_%s_DataSet1_pfind.txt' % (spectra_title,(i+1), xlinktype))
        samples.append(middle_file_path)
    pbuild_file_path = os.path.join(filter_dir_path, 'total_%s.combine.pbuild' % xlinktype)
    pbuild_result_file_path = os.path.join(filter_dir_path, '%s_%s_combine' % (spectra_title, xlinktype), '%s_%s_combine_DataSet1_pfind.txt' % (spectra_title, xlinktype))
    if os.path.isfile(pbuild_result_file_path):
        print builder_result_file_path + ' is existed, skip the step'
    else:
        print 'Filter of ' + pbuild_file_path
        extraParam = {}
        extraParam['SampleID'] = '%s_%s_combine' % (spectra_title, xlinktype)
        extraParam['OutPutName'] = '%s_%s_combine' % (spectra_title, xlinktype)
        extraParam['PepTolBase'] = '0' 
        extraParam['PepTolLower'] = '-100'
        extraParam['PepTolUpper'] = '100'
        extraParam['PepTolType'] = 'Da'
        extraParam['CrossLink'] = '2'
        
        pBuildParam = _ConstructpBuildParam(initparam, [samples], filter_dir_path, extraParam)
        _WritepBuild(pBuildParam, [samples], pbuild_file_path)
        bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path))
        #Evalue cut, refilter, 2013.10.22
        if (bEvalueMax):
            if platform.system() == 'Linux':
                    bfile.write('cp -b "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
                    bfile.write('rm "%s"\n' % builder_result_file_path)
            elif platform.system() == 'Windows':
                    bfile.write('move "%s" "%s"\n' % (builder_input_file_path, builder_input_file_path+"~3"))
                    bfile.write('move "%s" "%s"\n' % (builder_result_file_path, builder_input_file_path))
            print 'Evalue filter of ' + builder_input_file_path

            extraParam['usefdr'] = 0
            pbuild_file_path_two = pbuild_file_path[0:-7]+"2.pbuild"
            pBuildParam = _ConstructpBuildParam(initparam, [[builder_input_file_path]], filter_dir_path, extraParam)
            _WritepBuild(pBuildParam, [[builder_input_file_path]], pbuild_file_path_two)
            bfile.write('"%s" "%s"\n' % (program_path, pbuild_file_path_two))
    
#     # PSMAnalysis
#     print 'Step : Generate PSM file by PSMAnalysis '
#     program_path = os.path.join(bin_path, 'PSMAnalysis')
#     for i in range(0, len(spectra_list)):
#         pfind_file_path = os.path.join(filter_dir_path, '%s%d.pfind' % (spectra_title, (i+1)))
#         psm_file_path = os.path.join(filter_dir_path, '%s%d.%s.psm' % (spectra_title, (i+1), xlinktype))
#         psm_input_file_path = os.path.join(filter_dir_path, '%s_%d_%s' % (spectra_title, (i+1), xlinktype), '%s_%d_%s_DataSet1_pfind.txt' % (spectra_title, (i+1), xlinktype))
#         psm_result_file_path = os.path.join(filter_dir_path, '%s%d.%s.psm.output' % (spectra_title, (i+1), xlinktype))
#         if os.path.isfile(psm_result_file_path):
#             print psm_result_file_path + ' is existed, skip the step'
#         else:
#             print 'Gen PSM file of ' + psm_file_path
#             psmParam = _ConstructPSMParam(initparam, pfind_file_path, psm_input_file_path, filter_dir_path)
#             _WritePSM(psmParam, psm_file_path)
#             bfile.write('"%s" "%s"\n' % (program_path, psm_file_path))
#     
#     # plabel
#     print 'Step : Convert to *.plabel by XLinkResultConverter'
#     program_path = os.path.join(bin_path, 'XLinkResultConverter')
#     for i in range(0, len(spectra_list)):
#         pfind_file_path = os.path.join(filter_dir_path, '%s%d.pfind' % (spectra_title, (i+1)))
#         convert_file_path = os.path.join(filter_dir_path, '%s%d.%s.2plabel.convert' % (spectra_title, (i+1), xlinktype))
#         convert_input_file_path = os.path.join(filter_dir_path, '%s_%d_%s' % (spectra_title, (i+1), xlinktype), '%s_%d_%s_DataSet1_pfind.txt' % (spectra_title, (i+1), xlinktype))
#         convert_result_file_path = os.path.join(filter_dir_path, '%s%d.%s.plabel' % (spectra_title, (i+1), xlinktype))
#         if os.path.isfile(convert_result_file_path):
#             print convert_result_file_path + ' is existed, skip the step'
#         else:
#             print 'Convert of ' + convert_file_path
#             convertParam = _ConstructConvertParam(pfind_file_path, convert_input_file_path, filter_dir_path, xlinktype, 'plabel')
#             convertParam['IO'].append(('InputType', 'pbuild'))
#             _WriteConvert(convertParam, convert_file_path)
#             bfile.write('"%s" "%s"\n' % (program_path, convert_file_path))
    
    return

    
