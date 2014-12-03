#!/usr/bin/env python
# -*-coding:utf8 -*-
'''
Created on 22 Jan 2013

@author: Jeremy
'''
import os
import platform
import string
import Param
import ClusterSetting

def _ConstructpFindParam(initparam, search_mode, normal_path):
    pFind_Param = {}
    # Version
    version = []
    version.append(('pLink_Version', 'EVA.0.0.1.2011.5.20'))
    pFind_Param['version'] = version;
    
    # database
    database = []
    database.append(('xlink_List', os.path.join(initparam['bin.path'],'xlink.ini')))

    index_content = 'PEPTIDE_PAIR'
    if not 'index_content' in initparam:
        if search_mode == 0:
            database.append(('index_content', 'PEPTIDE_PAIR'))
        elif search_mode == 1:
            database.append(('index_content', 'PEPTIDE_PAIR_OPEN'))
        elif search_mode == 2:
            database.append(('index_content', 'PEPTIDE_SUMO'))
        elif search_mode == 3:
            database.append(('index_content', 'PEPTIDE_TRI_ALL'))
        elif search_mode == 4:
            database.append(('index_content', 'PEPTIDE_TRI_ION'))
    else:
        index_content = initparam['index_content']
        database.append(('index_content', index_content))
            
    database.append(('db_type', 'DB_INDEXED_FILE_ALL'))
    database.append(('index_type', 'MEMORY_MAPPED'))
    database.append(('db_total', '1'))
    database.append(('db_name1', initparam['database.name']))
    pFind_Param['database'] = database
    
    # enzyme
    enzyme = []
    enzyme.append(('enzyme_name', initparam['enzyme.name']))
    pFind_Param['enzyme'] = enzyme
    
    #modify
    modify = []
    mod_fixed_list = []
    mod_fixed_num = string.atoi(initparam['mod.fixed.total'])
    for i in range(1,mod_fixed_num+1):
        mod_fixed_list.append(initparam['mod.fixed.%d' % i])
    modify.append(('fix_total', '%d' %mod_fixed_num))
    for i in range(0, mod_fixed_num):
        modify.append(('fix_mod%d' %(i+1), mod_fixed_list[i]))
    mod_variable_list = []
    mod_var_num = string.atoi(initparam['mod.variable.total'])
    for i in range(1,mod_var_num+1):
        mod_variable_list.append(initparam['mod.variable.%d' % i])
    modify.append(('var_total', '%d' %mod_var_num))
    for i in range(0, mod_var_num):
        modify.append(('var_mod%d' %(i+1), mod_variable_list[i]))
    modify.append(('max_number', '3'))
    pFind_Param['modify'] = modify
    
    # xlink
    xlink = []
    linkerNum = string.atoi(initparam['linker.total']);
    if linkerNum <= 0:
        raise Exception,'no linker';
    xlink.append(('linker_total', '%d' %linkerNum))
    xlink.append(('linker', initparam['linker.name1']))
    for i in range(1,linkerNum):
        xlink.append(('linker%d' %i, initparam['linker.name%d' %(i+1)] ))
    pFind_Param['xlink'] = xlink
    
    #ions
    ions = []
    InstrInfo = Param.ReadInstrument(initparam)
    for elem in InstrInfo:
        ions.append(elem)
    if search_mode == 0:
        initparam['peptide_tol_total'] = '1'
    tolNum = string.atoi(initparam['peptide_tol_total'])
    if tolNum <= 0:
        raise Exception,'no tolerance window'
    ions.append(('peptide_tol_total', '%d' %tolNum))
    ions.append(('peptide_tol', initparam['peptide_tol1']))
    ions.append(('peptide_tol_type', initparam['peptide_tol_type1']))
    ions.append(('peptide_tol_base', initparam['peptide_tol_base1']))
    ions.append(('peptide_tol_base_type', initparam['peptide_tol_base_type1']))
    for i in range(1, tolNum):
        InstrInfo.append(('peptide_tol%d' %i, initparam['peptide_tol%d' %(i+1)]))
        InstrInfo.append(('peptide_tol_type%d' %i, initparam['peptide_tol_type%d' %(i+1)]))
        InstrInfo.append(('peptide_tol_base%d' %i, initparam['peptide_tol_base%d' %(i+1)]))
        InstrInfo.append(('peptide_tol_base_type%d' %i, initparam['peptide_tol_base_type%d' %(i+1)]))
    ions.append(('min_peptide_tol', '1'))
    ions.append(('min_peptide_tol_type', 'ppm'))

    # max_miss_site added at 2013.11.17
    if not 'max_miss_site' in initparam:
            max_miss_site = '2'
    else:
            max_miss_site = initparam['max_miss_site']
    ions.append(('max_miss', max_miss_site))
    
    pFind_Param['ions'] = ions
    
    #### used for simplescore and flow
    instru = initparam['spectra.instrument']
    preproc = 'PRE_PROC_XLINK_HCD'
    if not 'preproc' in initparam:
        if instru == 'HCD':
            preproc = 'PRE_PROC_XLINK_HCD';
        elif instru == 'ETD':
            preproc = 'PRE_PROC_ETD_SIMPLE';
        else:
            preproc = 'PRE_PROC_XLINK_HCD';
##            preproc = 'PRE_PROC_DEFAULT'; #Default changed at 2014.7.17
    else:
        preproc = initparam['preproc']

    evaluat = 'EV_DEFAULT'
    if not 'evaluat' in initparam:
        evaluat = 'EV_DEFAULT'
    else:
        evaluat = initparam['evaluat']

    score = 'SCORE_KSDP'
    if not 'score' in initparam:
        score = 'SCORE_KSDP'
    else:
        score = initparam['score']

    if not 'pepgen' in initparam:
        pepgen = 'PEPGENE_DEFAULT'
    else:
        pepgen = initparam['pepgen']
    
    #simplescore
    if search_mode == 1 or search_mode == 4: #add tri ion mode loading at 2014.4.18
        simplescore = []
        simplescore.append(('preproc', preproc))
        simplescore.append(('simple_score_cutoff', '3'))
        simplescore.append(('report_peptide_number', '500'))
        if instru == 'ETD':
            simplescore.append(('ion_type_total', '4'))
            simplescore.append(('ion_type1', 'c 1 0 0 0'))        
            simplescore.append(('ion_type2', 'c 2 0 0 0'))
            simplescore.append(('ion_type3', 'z 1 0 0 0'))
            simplescore.append(('ion_type4', 'z 2 0 0 0'))
        else:
            simplescore.append(('ion_type_total', '4'))
            simplescore.append(('ion_type1', 'b 1 0 0 0'))            
            simplescore.append(('ion_type2', 'b 2 0 0 0'))
            simplescore.append(('ion_type3', 'y 1 0 0 0'))
            simplescore.append(('ion_type4', 'y 2 0 0 0'))
        pFind_Param['simplescore'] = simplescore
    elif search_mode == 2:
        sumo = []
        print 'here in SUMO'
        os.system("pause")
        if not 'sumoseq' in initparam:
            sumoseq = "MLGG"
        else:
            sumoseq = initparam['sumoseq']
        if not 'targetsite' in initparam:
            targetsite = "K"
        else:
            targetsite = initparam['targetsite']
        if not 'fixsite' in initparam:
            fixsite = "G"
        else:
            fixsite = initparam['fixsite']
        sumo.append(('sumoseq', sumoseq))
        sumo.append(('fixsite', fixsite))
        sumo.append(('targetsite', targetsite))
        pFind_Param['sumo'] = sumo
    
    procnum = '1'
    if not 'processor_num' in initparam:
        procnum = '1'
    else:
        procnum = initparam['processor_num']
    
    max_ev = '1.0'
    if not 'max_ev' in initparam:
        max_ev = '1.0'
    else:
        max_ev = initparam['max_ev']
    
    #flow
    flow = []
    flow.append(('processor_num', procnum))
    flow.append(('max_score_num', '0'))
    flow.append(('min_score_num', '5000'))
    if procnum == '1':
        if initparam['search_mode'] == '4':
            flow.append(('salvo_batch_size', '8000'))
        else:
            flow.append(('salvo_batch_size', '4000'))
    else:
        flow.append(('salvo_batch_size', str(4000/string.atoi(procnum))))
    flow.append(('log_rank', 'LOG_RANK_INFO'))
    flow.append(('flow', 'FLOW_XLINK'))
    flow.append(('instrument', instru))
    flow.append(('preproc', preproc))
    flow.append(('score', score))
    flow.append(('evaluat', evaluat))
    flow.append(('pepgen', pepgen))
    flow.append(('max_ev', max_ev))
    if initparam['search_mode'] == '4':
        flow.append(('report_peptide_number', '10'))
    else:
        flow.append(('report_peptide_number', '3'))
    flow.append(('report_peptide_type', 'PEPTIDE_REPORT_XML'))
    flow.append(('report_protein_number', '2000'))
    flow.append(('report_protein_type', 'PROTEIN_REPORT_TXT'))
    flow.append(('false_positive_rate', '1'))
    flow.append(('false_positive_sign', 'REVERSE'))
    flow.append(('output_path', normal_path))
    pFind_Param['flow'] = flow
    
    # spectrum
    spectrum = []
    pFind_Param['spectrum'] = spectrum
    
    # cluster
    if platform.system() == 'Linux':
        cluster = []
        cluster.append(('block_path', ClusterSetting.TempFilePath))
        cluster.append(('block_format_type', 'SPECTRA_INDEX_SIMPLE'))
        cluster.append(('load_balance_type', 'LOAD_BALANCE_MASS_DYNAMIC_MUL4'))
        cluster.append(('processor_num', '60'))
        cluster.append(('block_max_num', '100'))
        pFind_Param['cluster'] = cluster
        
    #triple peptide
    if not 'sidechain_singlesite' in initparam:
        sidechain_singlesite = 'true'
    else:
        sidechain_singlesite = initparam['sidechain_singlesite']
    
    triplepeptide = []
    triplepeptide.append(('SingleCAsSitePep', sidechain_singlesite))

    if not 'reversetag' in initparam:
        reversetag = 'REVERSE_'
    else:
        reversetag = initparam['reversetag']
        
    triplepeptide.append(("ReverseTag", reversetag))
    
    pFind_Param['triplepeptide'] = triplepeptide
    return pFind_Param

def _WritepFind(pFind_Param, pfind_file, search_mode):
    pfile = open(pfind_file,'w')
    keys = ['version', 'database', 'enzyme', 'modify', 'xlink', 'ions', 'simplescore', 'flow', 'sumo', 'spectrum', 'cluster', 'triplepeptide']
    for key in keys:
        if key == 'simplescore' and search_mode != 1 and search_mode != 4: #add tri ion mode loading at 2014.4.18
            continue
        if key == 'cluster' and platform.system() == 'Windows':
            continue
        if key == 'sumo' and search_mode != 2:
            continue
        pfile.write('['+key+']\n')
        for elem in pFind_Param[key]:
            pfile.write(elem[0]+'='+elem[1]+'\n')
    pfile.close()
    return 
