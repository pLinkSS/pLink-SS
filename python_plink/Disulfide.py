#!/usr/bin/env python
#-*-coding:utf8 -*-
import os
import sys
import platform
import time
import string
import Path
import ClusterSetting
import Search
import ConfigParser

class Spectrum:
    """ a spectrum """
    
    def __init__(self):
        self.title = ''
        self.charge = 1
        self.mh = 0.0
        self.mz = 0.0
    
    def __str__(self):
        return 'title=%s\ncharge=%d\nMH=%.6f\nMZ=%.6f' % (self.title, self.charge, self.mh, self.mz)

class Protein:
    """ a protein item """
    
    def __init__(self, name='', id=-1, site=-1):
        self.name = name
        self.id = id
        self.site = site
        
    def __str__(self):
        return '(%s, %d, %d)' % (self.name, self.id, self.site)

class ModificationSite:
    """ a modification item """
    
    def __init__(self, name='', site=-1):
        self.name = name
        self.site = site
        
    def __str__(self):
        return '(%s, %d)' % (self.name, self.site)
            
class Peptide:
    """ a peptide """
    
    def __init__(self):
        self.sq = ''
        self.mass = 0.0
        self.proteins = list()
        self.linker_site = ModificationSite()
        self.modification = list()
        
    def __str__(self):
        ss = 'sq=%s\nmass=%.6f\nlinker_site=%s\n' % (self.sq, self.mass, str(self.linker_site))
        ss += 'modifications='
        for mod in self.modification:
            ss += str(mod)
        ss += '\n'
        ss += 'proteins='
        for pro in self.proteins:
            ss += str(pro)
        return ss
        
class PeptidePair:
    """ a peptide pair """
    
    def __init__(self):
        self.mass = 0.0
        self.score = 0.0
        self.evalue = 0.0
        self.qvalue = 100.0
        self.peptide_type = 0 # 0 - common, 1 - mono, 2 - loop, 3 - inter
        self.protein_type = 0 # 0 - none, 1 - intra, 2 - inter, 3 - identical site
        self.td_type = 0 # 0 - FF, 1 - FT/TF, 2 - TT
        self.dalton_error = 0.0
        self.ppm_error = 0.0
        self.alpha = Peptide()
        self.beta = Peptide()
        
    def __str__(self):
        pep_type = ['common', 'mono', 'loop', 'inter']
        pro_type = ['none', 'intra', 'inter', 'identical_site']
        td_type = ['FF', 'FT/TF', 'TT']
        return 'pair_mass=%.6f\nscore=%.6f\nevalue=%.6e\npeptide_type=%s\nprotein_type=%s\ntd_type=%s\ndalton_error=%.6f\nppm_error=%.6f\n%s\n%s' % \
            (self.mass, self.score, self.evalue, \
             pep_type[self.peptide_type], pro_type[self.protein_type], td_type[self.td_type], \
             self.dalton_error, self.ppm_error, str(self.alpha), str(self.beta))
          
class PSM:
    """ only save top 1, easy to expand to top K """
    
    def __init__(self):
        self.spec = Spectrum()
        self.pair = PeptidePair()
        
    def __str__(self):
        return '%s\n%s' % (str(self.spec), str(self.pair))
        
def set_current_spectra_info(sno, init_param):
    # set current spectra basic information
    init_param['spectra.instrument'] = init_param['sample%d.spectra.instrument' % sno]
    init_param['spectra.format'] = init_param['sample%d.spectra.format' % sno]
    init_param['spectra.path'] = init_param['sample%d.spectra.path' % sno]
    init_param['spectra.title'] = init_param['sample%d.spectra.title' % sno]
    
    Path.Check(os.path.join(init_param['origin.output.path'], '%d.sample' % sno) )
    init_param['output.path'] = os.path.join(init_param['origin.output.path'], '%d.sample' % sno, 'search')
    Path.Check(init_param['output.path'])
    
    # set current file list, walk through the directories recursively
    init_param['file_list'] = list()
    if platform.system() == 'Linux':
        spec_format = ('*.%s' % init_param['spectra.format']).upper()
        Path.Walk(init_param['spectra.path'], spec_format, init_param['file_list']);
    spec_format = ('*.' + init_param['spectra.format']).lower();
    Path.Walk(init_param['spectra.path'], spec_format, init_param['file_list']);

    # join current file path
    init_param['spectra_list'] = list()
    for cur_file in init_param['file_list']:
        init_param['spectra_list'].append(os.path.join(cur_file[0], cur_file[1]))

def search(init_param):
    try:
        Path.Check(init_param['output.path'])
        
        if platform.system() == 'Linux':
            bat_file = os.path.join(init_param['output.path'], 'normal.bash')
            bat_fp = open(bat_file, 'w')
            bat_fp.write('export PATH=%s:$PATH\n' % ClusterSetting.MPIPath) #modified 2012.6.11
            bat_fp.write('export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n' % ClusterSetting.pLinkBinPath)
        elif platform.system() == 'Windows':
            bat_file = os.path.join(init_param['output.path'], 'normal.bat')
            bat_fp = open(bat_file, 'w')
            bat_fp.write('@echo off\n')
            bat_fp.write('%s\n' % init_param['bin.path'][0:2])
        else:
            raise Exception('search_and_filter', 'unknown platform, only support Windows and Linux')
        
        bat_fp.write('cd "%s"\n' % init_param['bin.path'])
        
        search_mode = string.atoi(init_param['search_mode'])
        pfind_param = Search._ConstructpFindParam(init_param, search_mode, init_param['output.path'])
        print 'Step : Search by Searcher'
        spectra_list = init_param['spectra_list']
        spectra_title = init_param['spectra.title']
        bin_path = init_param['bin.path']
        
        for i in range(0, len(init_param['spectra_list'])):
            pfind_file = os.path.join(init_param['output.path'], '%s%d.pfind' % (init_param['spectra.title'], i+1))
            pfind_result_file = os.path.join(init_param['output.path'], '%s%d_qry.proteins.txt' % (init_param['spectra.title'], i+1))
    
            if os.path.isfile(pfind_result_file):
                print os.path.split(pfind_result_file)[-1] + ' did exist, skip the step';
            else:
                print 'Searcher of '+ spectra_list[i];
                spectrum = []
                spectrum.append(('spec_title', spectra_title+'%d' %(i+1)))
                spectrum.append(('spec_type', '%s' % init_param['spectra.format'].upper()))
                spectrum.append(('spec_path', spectra_list[i]))
                pfind_param['spectrum'] = spectrum
                
                Search._WritepFind(pfind_param, pfind_file, search_mode)
                if platform.system() == 'Windows':
                    bat_fp.write('"%s" "%s"\n' % (os.path.join(bin_path,'Searcher'), pfind_file))
                else:
                    bat_fp.write('"%s" "%s"\n' % (os.path.join(bin_path,'Importer'), pfind_file))
                    if 'mpicores' in init_param:
                        mpicores = init_param['mpicores']
                    else:
                        mpicores = ClusterSetting.DefaultCores
                    if 'mpihosts' in init_param:
                        mpihosts = init_param['mpihosts']
                    else:
                        mpihosts = ClusterSetting.DefaultHosts
                    bat_fp.write('mpirun -np %s -host %s "%s" "%s"\n' %(mpicores, mpihosts, os.path.join(bin_path,'MPISearcher'), pfind_file))
        bat_fp.close()

    except Exception, e:
        print Exception + ": " + e

    if platform.system() == 'Linux':
        os.system('chmod 766 %s' % bat_file)
        
    os.system('"%s"' % bat_file)
    return os.stat(pfind_result_file).st_mtime

def set_peptide(conf, peptide, section, subpeptide):
    try:
        if not isinstance(peptide, Peptide):
            raise Exception('set_peptide', 'unknown type: peptide')
        if not isinstance(section, basestring):
            raise Exception('set_peptide', 'unknown type: section')
        if not isinstance(subpeptide, basestring):
            raise Exception('set_peptide', 'unknown type: subpeptide')
        
        peptide.sq = conf.get(section, 'NO1_%s_SQ' % subpeptide)
        ss_v = conf.get(section, 'NO1_%s_Proteins' % subpeptide)
        ll_v = ss_v.split(',') 
        protein_list = ll_v[1:]
        ss_v = conf.get(section, 'NO1_%s_ProteinIDs' % subpeptide)
        ll_v = ss_v.split(',')
        protein_id_list = [string.atoi(s) for s in ll_v[1:]]
        ss_v = conf.get(section, 'NO1_%s_ProteinSites' % subpeptide)
        ll_v = ss_v.split(',')
        protein_site_list = [string.atoi(s) for s in ll_v[1:]]
        for item in zip(protein_list, protein_id_list, protein_site_list):
            protein = Protein()
            protein.name = item[0]
            protein.id = item[1]
            protein.site = item[2]
            peptide.proteins.append(protein)

        ss_v = conf.get(section, 'NO1_%s_Modify_Pos' % subpeptide)
        ll_v = ss_v.split(',')
        mod_pos_list = [string.atoi(s) for s in ll_v[1:]]
        ss_v = conf.get(section, 'NO1_%s_Modify_Name' % subpeptide)
        ll_v = ss_v.split(',')
        mod_list = ll_v[1:]
        for item in zip(mod_list, mod_pos_list):
            modification = ModificationSite(name=item[0], site=item[1])
            peptide.modification.append(modification)
    except Exception, e:
        print e

def target_or_decoy(peptide):
    for protein in peptide.proteins:
        if protein.name[0:8] != 'REVERSE_':
            return 1
    return 0

def protein_type_of(psm_pair):
    if psm_pair.peptide_type < 3:
        return 0 # None
    
    alpha_id_list = list()
    beta_id_list = list()
    if psm_pair.td_type == 2:
        for protein in psm_pair.alpha.proteins:
            if protein.name[0:8] != 'REVERSE_':
                alpha_id_list.append(protein.id)
        
        for protein in psm_pair.beta.proteins:
            if protein.name[0:8] != 'REVERSE_':
                beta_id_list.append(protein.id)
    elif psm_pair.td_type == 0:
        for protein in psm_pair.alpha.proteins:
            alpha_id_list.append(protein.id)
        
        for protein in psm_pair.beta.proteins:
            beta_id_list.append(protein.id)
    else:
        return 2 # inter
    
    alpha_id_list.sort()
    beta_id_list.sort()
    i = j = 0
    while i < len(alpha_id_list) and j < len(beta_id_list):
        if alpha_id_list[i] == beta_id_list[j]:
            blen = len(psm_pair.beta.sq)
            if psm_pair.beta.sq == psm_pair.alpha.sq[:blen] or psm_pair.beta.sq == psm_pair.alpha.sq[-blen:]:
                return 3 # identical site
            return 1 # intra
        elif alpha_id_list[i] < beta_id_list[j]:
            i = i + 1
        else:
            j = j + 1
    return 2 # inter

def get_dalton_error(thr_m, exp_m):
    mass_base = [0, 1.003, 2.006, 3.009, 4.012]
    dalton_e = abs(exp_m - thr_m)
    base_i = 0
    min_e = sys.float_info.max
    for base_m in mass_base:
        cur_e = abs(dalton_e - base_m)
        if cur_e < min_e:
            min_e = cur_e
            base_i = base_m
            
    if exp_m > thr_m:
        exp_m = exp_m - base_i
    else:
        exp_m = exp_m + base_i
        
    return exp_m - thr_m

def get_ppm_error(thr_m, exp_m):
    mass_base = [0, 1.003, 2.006, 3.009, 4.012]
    dalton_e = abs(exp_m - thr_m)
    base_i = 0
    min_e = sys.float_info.max
    for base_m in mass_base:
        cur_e = abs(dalton_e - base_m)
        if cur_e < min_e:
            min_e = cur_e
            base_i = base_m
            
    if exp_m > thr_m:
        exp_m = exp_m - base_i
    else:
        exp_m = exp_m + base_i
        
    return (exp_m - thr_m) * 1.0e+6 / thr_m
         
def set_pro_td_type(psm_pair):
    try:
        if not isinstance(psm_pair, PeptidePair):
            raise Exception('set_pro_td_type', 'unknown type: psm_pair')
        
        if psm_pair.peptide_type == 3:
            alpha_tag = target_or_decoy(psm_pair.alpha)
            beta_tag = target_or_decoy(psm_pair.beta)
            psm_pair.td_type = alpha_tag + beta_tag
            psm_pair.protein_type = protein_type_of(psm_pair)
        else:
            psm_pair.td_type = target_or_decoy(psm_pair.alpha)
            if psm_pair.td_type == 1:
                psm_pair.td_type = 2
            psm_pair.protein_type = 0
    except Exception, e:
        print e

def seperate_psms(category_psms, psm, init_param):
    if psm.pair.peptide_type == 0:
        category_psms['common'].append(psm)
    elif psm.pair.peptide_type == 2:
        category_psms['loop'].append(psm)
    elif psm.pair.peptide_type == 3:
        if init_param['fdr_strategy'] == 1:
            category_psms['intra'].append(psm)
        else:
            if psm.pair.protein_type == 1:
                category_psms['intra'].append(psm)
            elif psm.pair.protein_type == 2:
                category_psms['inter'].append(psm)
            elif psm.pair.protein_type == 3:
                category_psms['identical'].append(psm)
            else:
                raise Exception('seperate_psms', 'unknown type')
    else:
        raise Exception('seperate_psms', 'unknown type')
        
def load_psms(category_psms, init_param):
    for i in range(0, len(init_param['spectra_list'])):
#        if i > 0:
#            break
        pfind_result_file = os.path.join(init_param['output.path'], '%s%d_qry.proteins.txt' % (init_param['spectra.title'], i+1))
        if not os.path.isfile(pfind_result_file):
            raise Exception('filter', 'failed to find result file %s' % pfind_result_file)
        
        conf = ConfigParser.ConfigParser()
        conf.read(pfind_result_file)
        
        spectra_num = conf.getint('Total', 'Spectra')
        for num in range(1, spectra_num+1):
            section = 'Spectrum%d' % num
            valid_num = conf.getint(section, 'ValidCandidate')
            if valid_num > 0 and conf.getfloat(section, 'NO1_EValue') < 1.0:
                psm = PSM()
                psm.spec.title = conf.get(section, 'Input')
                psm.spec.charge = conf.getint(section, 'Charge')
                psm.spec.mh = conf.getfloat(section, 'MH')
                psm.spec.mz = conf.getfloat(section, 'MZ')
                psm.pair.mass = conf.getfloat(section, 'NO1_Mass')
                psm.pair.score = conf.getfloat(section, 'NO1_Score')
                psm.pair.evalue = conf.getfloat(section, 'NO1_EValue')
                psm.pair.peptide_type = conf.getint(section, 'NO1_XLink_Type')
                if psm.pair.peptide_type >= 1:
                    psm.pair.alpha.linker_site.id = conf.getint(section, 'NO1_Linker_Id')
                    psm.pair.alpha.linker_site.site = conf.getint(section, 'NO1_XLink_Pos1')
                if psm.pair.peptide_type >= 2:
                    psm.pair.beta.linker_site.id = conf.getint(section, 'NO1_Linker_Id')
                    psm.pair.beta.linker_site.site = conf.getint(section, 'NO1_XLink_Pos2')
                set_peptide(conf, psm.pair.alpha, section, 'Alpha')
                if psm.pair.peptide_type >= 3:
                    set_peptide(conf, psm.pair.beta, section, 'Beta')
                psm.pair.dalton_error = psm.spec.mh - psm.pair.mass
                psm.pair.ppm_error = get_ppm_error(psm.pair.mass, psm.spec.mh)
                set_pro_td_type(psm.pair)
                bases = init_param['filter_peptide_tol_base']
                bases_v = bases.split(',')
                if psm.pair.peptide_type != 1 and abs(psm.pair.dalton_error) < (len(bases_v) - 0.5):
                    seperate_psms(category_psms, psm, init_param)
            print 'complete %.2f%%\r' % (100.0*num/(spectra_num+1)),
        print 'complete 100.00%'
        print 'loaded file: %s' % os.path.split(pfind_result_file)[-1]        

def compute_qvalue(category_psms, psm_types, init_param):
    thr_e = string.atoi(init_param['filter_peptide_tol_ub'].split(',')[0])
    for psm_type in psm_types:
        print 'Computing %s q-values...' % psm_type
        td_cnt = [0, 0, 0]
        category_psms[psm_type].sort(key=lambda x: x.pair.evalue)
        
        for psm in category_psms[psm_type]:
            # precursor error, length and mass filter
            if init_param['filter_peptide_tol_type'] == 'ppm':
                if abs(psm.pair.ppm_error) > thr_e:
                    continue
            else:
                if abs(get_dalton_error(psm.pair.mass, psm.spec.mh)) > thr_e:
                    continue
                
            td_cnt[psm.pair.td_type] += 1
            
            if td_cnt[2] == 0 and (td_cnt[0] > 0 or td_cnt[1] > 0):
                psm.pair.qvalue = 1.0
            else:
                target = td_cnt[2]
                if psm_type == 'intra':
                    random = td_cnt[0]
                elif td_cnt[1] > td_cnt[0]:
                    random = td_cnt[1] - td_cnt[0]
                else:
                    random = td_cnt[0]
                if random > target:
                    psm.pair.qvalue = 1.0
                else:
                    psm.pair.qvalue = 1.0 * random / target
        
        cur_qvalue = 100.0
        for psm in category_psms[psm_type][::-1]:
            if psm.pair.qvalue == 100.0:
                continue
            elif psm.pair.qvalue < cur_qvalue:
                cur_qvalue = psm.pair.qvalue
            else:
                psm.pair.qvalue = cur_qvalue
                
        category_psms[psm_type].sort(key=lambda x: (x.pair.qvalue, x.pair.evalue))
          
def filter(init_param):
    try:
        print 'Step : Filtering'
        psm_types = ['common', 'loop', 'intra', 'inter', 'identical']
        category_psms = dict()
        for psm_type in psm_types:
            category_psms[psm_type] = list()
        load_psms(category_psms, init_param)
        compute_qvalue(category_psms, psm_types, init_param)
        
        if 'fdr' not in init_param:
            init_param['fdr'] = 0.05
        if 'evalue_max' not in init_param:
            init_param['evalue_max'] = 1.0
        filtered_psms = dict()
        for psm_type in psm_types:
            filtered_psms[psm_type] = list()
        for psm_type in psm_types:
            for psm in category_psms[psm_type]:
                if psm.pair.qvalue > init_param['fdr'] or psm.pair.evalue > init_param['evalue_max']:
                    break
                if psm.pair.td_type == 2:
                    filtered_psms[psm_type].append(psm)
            print '%s: %d' % (psm_type, len(filtered_psms[psm_type]))
        return filtered_psms
    except Exception, e:
        print e
        return None

def count_c_num(sq):
    num = 0
    for c in sq:
        if c == 'C':
            num = num + 1
    return num

def get_detail_type(psm, psm_type):
    if psm_type == 'common':
        alpha_c_num = count_c_num(psm.pair.alpha.sq)
        if alpha_c_num == 0:
            detail_type = 0
        elif alpha_c_num % 2 == 0:
            detail_type = 1
        else:
            detail_type = -1
    elif psm_type == 'loop':
        alpha_c_num = count_c_num(psm.pair.alpha.sq)
        if alpha_c_num % 2 == 0:
            detail_type = 1
        else:
            detail_type = -1
    elif psm_type == 'intra':
        alpha_c_num = count_c_num(psm.pair.alpha.sq)
        beta_c_num = count_c_num(psm.pair.beta.sq)
        c_num = alpha_c_num + beta_c_num
        if c_num == 2:
            detail_type = 2
        elif c_num % 2 == 0:
            detail_type = 3
        else:
            detail_type = 4
    elif psm_type == 'inter':
        alpha_c_num = count_c_num(psm.pair.alpha.sq)
        beta_c_num = count_c_num(psm.pair.beta.sq)
        c_num = alpha_c_num + beta_c_num
        if c_num == 2:
            detail_type = 5
        elif c_num % 2 == 0:
            detail_type = 6
        else:
            detail_type = 7
    elif psm_type == 'identical':
        alpha_c_num = count_c_num(psm.pair.alpha.sq)
        beta_c_num = count_c_num(psm.pair.beta.sq)
        c_num = alpha_c_num + beta_c_num
        if c_num % 2 == 0:
            detail_type = 8
        else:
            detail_type = -1
    else:
        detail_type = -1
    return detail_type

def get_sq(psm, detail_type):
    if detail_type == 0:
        return psm.pair.alpha.sq
    elif detail_type == 1:
        alpha_site = list()
        for mod in psm.pair.alpha.modification:
            if mod.name == 'C-1':
                alpha_site.append(mod.site)
        alpha_site.sort()
        return '%s(%s)' % (psm.pair.alpha.sq, ','.join([str(c) for c in alpha_site]))
    elif detail_type >=2 and detail_type <= 8:
        alpha_site = list()
        for mod in psm.pair.alpha.modification:
            if mod.name == 'C-1':
                alpha_site.append(mod.site)
        alpha_site.sort()
        beta_site = list()
        for mod in psm.pair.beta.modification:
            if mod.name == 'C-1':
                beta_site.append(mod.site)
        beta_site.sort()
        return '%s(%s)-%s(%s)' % (psm.pair.alpha.sq, ','.join([str(c) for c in alpha_site]), psm.pair.beta.sq, ','.join([str(c) for c in beta_site]))
    else:
        return ''

def get_out_type(psm, detail_type):
    if detail_type == 1:
        return 0
    elif detail_type == 2 or detail_type == 5:
        return 1
    elif detail_type == 3 or detail_type == 6:
        return 2
    elif detail_type == 8:
        alpha_c_num = count_c_num(psm.pair.alpha.sq)
        beta_c_num = count_c_num(psm.pair.beta.sq)
        c_num = alpha_c_num + beta_c_num
        if c_num == 2:
            return 1
        else:
            return 2
    else:
        return -1

def get_mod_str(psm, detail_type):
    if detail_type == 0 or detail_type == 1:
        mods = list()
        for mod in psm.pair.alpha.modification:
            if mod.site == 0:
                mods.append('%d,%c(%s)' % (mod.site, psm.pair.alpha.sq[mod.site], mod.name)) # N terminal
            else:
                mods.append('%d,%c(%s)' % (mod.site, psm.pair.alpha.sq[mod.site-1], mod.name))
        if len(mods) == 0:
            return 'null'
        return ';'.join(mods)
    elif detail_type >=2 and detail_type <= 8:
        mods = list()
        for mod in psm.pair.alpha.modification:
            if mod.site == 0:
                mods.append('%d,%c(%s)' % (mod.site, psm.pair.alpha.sq[mod.site], mod.name))
            else:
                mods.append('%d,%c(%s)' % (mod.site, psm.pair.alpha.sq[mod.site-1], mod.name))
        alpha_len = len(psm.pair.alpha.sq)
        for mod in psm.pair.beta.modification:
            if mod.site == 0:
                mods.append('%d,%c(%s)' % (alpha_len+mod.site+3, psm.pair.beta.sq[mod.site], mod.name))
            else:
                mods.append('%d,%c(%s)' % (alpha_len+mod.site+3, psm.pair.beta.sq[mod.site-1], mod.name))
        if len(mods) == 0:
            return 'null'
        return ';'.join(mods)
    else:
        return 'null'

def get_pro_str(psm, detail_type):
    if detail_type == 0:
        pros = list()
        for pro in psm.pair.alpha.proteins:
            pros.append('%s(%d)' % (pro.name, pro.site))
        return '/'.join(pros)
    elif detail_type == 1:
        pros = list()
        for pro in psm.pair.alpha.proteins:
            sites = list()
            for mod in psm.pair.alpha.modification:
                if mod.name == 'C-1':
                    sites.append(mod.site+pro.site)
            sites.sort()
            pros.append('%s(%s)' % (pro.name, ','.join([str(c) for c in sites])))
        return '/'.join(pros)
    elif detail_type >=2 and detail_type <= 8:
        alpha_pros = list()
        for pro in psm.pair.alpha.proteins:
            sites = list()
            for mod in psm.pair.alpha.modification:
                if mod.name == 'C-1':
                    sites.append(mod.site+pro.site)
            sites.sort()
            alpha_pros.append('%s(%s)' % (pro.name, ','.join([str(c) for c in sites])))
        
        beta_pros = list()
        for pro in psm.pair.beta.proteins:
            sites = list()
            for mod in psm.pair.beta.modification:
                if mod.name == 'C-1':
                    sites.append(str(mod.site+pro.site))
            sites.sort()
            beta_pros.append('%s(%s)' % (pro.name, ','.join([str(c) for c in sites])))
            
        pros = list()
        for a in alpha_pros:
            for b in beta_pros:
                pros.append('%s-%s' % (a, b))
        return '/'.join(pros)
    else:
        return ''

def get_output_line(psm, psm_type):
    detail_type = get_detail_type(psm, psm_type)
    
    sq = get_sq(psm, detail_type)
    
    out_type = get_out_type(psm, detail_type)
    
    out_str = '%s\t%.2e\t%.6f\t%.6f\t%.6f\t%s\t%s\t%s\t%s\t%s' % \
        (psm.spec.title, psm.pair.evalue, psm.pair.mass, \
         psm.pair.dalton_error, psm.pair.ppm_error, \
         get_mod_str(psm, detail_type), '1.pFind', '1.pFind', '1', \
         get_pro_str(psm, detail_type))
    return (sq, out_type, detail_type, out_str)

def report(filtered_psms, init_param):
    print 'Step : Report'
    
    out_num = 3 # 0 - loop, 1 - inter, 2 - complex
    out_psms = list()
    for i in range(0, out_num):
        out_psms.append(dict())
    
    detail_num = 9 # 0 - common, 1 - loop, 2 - intraM.interX, 3 - intraM.CplX, 4 - intraM.oddC, 5 - interM.interX, 6 - interM.CplX, 7 - interM.OddC, 8 - same-site
    detail_psms = list()
    for i in range(0, detail_num):
        detail_psms.append(dict())
    
    for psm_type in filtered_psms:
        for psm in filtered_psms[psm_type]:
            (sq, out_type, detail_type, out_str) = get_output_line(psm, psm_type)
            if out_type >= 0 and out_type < out_num:
                if sq not in out_psms[out_type]:
                    out_psms[out_type][sq] = [out_str]
                else:
                    out_psms[out_type][sq].append(out_str)
            if detail_type >= 0 and detail_type < detail_num:
                if sq not in detail_psms[detail_type]:
                    detail_psms[detail_type][sq] = [out_str]
                else:
                    detail_psms[detail_type][sq].append(out_str)
    
    sample_num = string.atoi(init_param['sample.num'])
    out_dir_path = os.path.join(init_param['origin.output.path'], '%d.report' % (sample_num+1))
    Path.Check(out_dir_path)
    detail_dir_path = os.path.join(init_param['origin.output.path'], '%d.report' % (sample_num+1), 'details')
    Path.Check(detail_dir_path)
    
    suffix = ['loopLink.xls', 'interLink.xls', 'complexLink.xls']
    
    for i in range(0, len(suffix)):
        fp = open(os.path.join(out_dir_path, '%s.%s' % (init_param['spectra.title'], suffix[i])), 'w')
        fp.write('Order\tSequence\n')
        fp.write('\tOrder\tSpectrum\tScore\tCalc_M\tDelta\tppm\tModification\tSampleID\tEngine\tRank\tProteins\n')
        
        pep_num = 1
        for pep in out_psms[i]:
            fp.write('%d\t%s\n' % (pep_num, pep))
            for spec_num in range(0, len(out_psms[i][pep])):
                fp.write('\t%d,%d\t%s\n' % (spec_num+1, pep_num, out_psms[i][pep][spec_num]))
            pep_num = pep_num + 1
        fp.close()

    if init_param['fdr_strategy'] == 1:
        suffix = ['linear.pep.xls', 'loop.pep.xls', \
                  'interX.pep.xls', 'CplX.pep.xls', 'OddC.pep.xls']
        for i in range(0, len(suffix)):
            fp = open(os.path.join(detail_dir_path, '%s.%s' % (init_param['spectra.title'], suffix[i])), 'w')
            fp.write('Order\tSequence\n')
            fp.write('\tOrder\tSpectrum\tScore\tCalc_M\tDelta\tppm\tModification\tSampleID\tEngine\tRank\tProteins\n')
            pep_num = 1
            for pep in detail_psms[i]:
                fp.write('%d\t%s\n' % (pep_num, pep))
                for spec_num in range(0, len(detail_psms[i][pep])):
                    fp.write('\t%d,%d\t%s\n' % (spec_num+1, pep_num, detail_psms[i][pep][spec_num]))
                pep_num = pep_num + 1
            fp.close()
    else:      
        suffix = ['linear.pep.xls', 'loop.pep.xls', \
                  'intraM.interX.pep.xls', 'intraM.CplX.pep.xls', 'intraM.OddC.pep.xls',\
                  'interM.interX.pep.xls', 'interM.CplX.pep.xls', 'interM.OddC.pep.xls', \
                  'bt-same-site.xls']
        
        for i in range(0, len(suffix)):
            fp = open(os.path.join(detail_dir_path, '%s.%s' % (init_param['spectra.title'], suffix[i])), 'w')
            fp.write('Order\tSequence\n')
            fp.write('\tOrder\tSpectrum\tScore\tCalc_M\tDelta\tppm\tModification\tSampleID\tEngine\tRank\tProteins\n')
            pep_num = 1
            for pep in detail_psms[i]:
                fp.write('%d\t%s\n' % (pep_num, pep))
                for spec_num in range(0, len(detail_psms[i][pep])):
                    fp.write('\t%d,%d\t%s\n' % (spec_num+1, pep_num, detail_psms[i][pep][spec_num]))
                pep_num = pep_num + 1
            fp.close()
        
    return
    
def main(log_file, init_param):
    try:
        # samples
        sample_num = string.atoi(init_param['sample.num'])
        for sno in range(1, sample_num+1):
            log_file.write('%s\n' % ('-'*30))

            set_current_spectra_info(sno, init_param)
            
            log_info = 'spectra list in sample %d are:' % sno
            print log_info
            log_file.write(log_info)
            for cur_file in init_param['spectra_list']:
                print cur_file
                log_file.write('%s\n' % cur_file)
            
            # search
            log_file.write('\nstep %d: search\n' % sno)
            log_file.write('\tdirectory: %s\n' % init_param['output.path'])
            
            start_time = time.time()
            last_modified_time = search(init_param)
            filtered_psms = filter(init_param)
            report(filtered_psms, init_param)
            end_time = time.time()
            interval = end_time - start_time
            
            log_file.write('\tlast modified time: %s\n' % time.ctime(last_modified_time))
            log_file.write('\trunning time: %.1f (s)\n' % interval)
            time_file = os.path.join(init_param['output.path'], 'running_time.txt')
            Path.WriteTime(time_file, last_modified_time, interval)
            print 'Search Task Finished.'
    except Exception, data:
        print 'Failed to run normal search:', data
        log_file.write('\tFailed to run normal search: %s\n' % data)
        log_file.close()
