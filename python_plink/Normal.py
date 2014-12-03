'''
Created on 21 Mar 2013

@author: Jeremy
'''
import Path
import os
import platform
import Search
import Filter
import ClusterSetting

def Run(initparam, spectra_list, search_mode, normal_dir_path):
    Path.Check(normal_dir_path)
    bin_path = initparam['bin.path']
    spectra_title = initparam['spectra.title']
    
    if platform.system() == 'Linux':
        bat_file=os.path.join(normal_dir_path,'normal.bash')
        bfile=open(bat_file, 'w')
        bfile.write('export PATH=%s:$PATH\n' % ClusterSetting.MPIPath)#modified 2012.6.11
        bfile.write('export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH\n' % ClusterSetting.pLinkBinPath)

    if platform.system() == 'Windows':
        bat_file=os.path.join(normal_dir_path,'normal.bat')
        bfile=open(bat_file, 'w')
        bfile.write('%s\n' % bin_path[0:2])
    
    bfile.write('cd "%s"\n' % bin_path)
    
    pFind_Param = Search._ConstructpFindParam(initparam, search_mode, normal_dir_path)
    print 'Step : Search by Searcher'
    for i in range(0, len(spectra_list)):
        pfind_file = os.path.join(normal_dir_path, spectra_title+'%d.pfind' %(i+1))
        pfind_result_file = os.path.join(normal_dir_path, spectra_title+'%d_qry.proteins.txt' %(i+1));

        if os.path.isfile(pfind_result_file):
            print pfind_result_file + ' is existed , skip the step';
        else:
            print 'Searcher of '+ spectra_list[i];
            spectrum = []
            spectrum.append(('spec_title', spectra_title+'%d' %(i+1)))
            spectrum.append(('spec_type', '%s' %initparam['spectra.format'].upper()))
            spectrum.append(('spec_path', spectra_list[i]))
            pFind_Param['spectrum'] = spectrum
            
            Search._WritepFind(pFind_Param, pfind_file, search_mode)
            if platform.system() == 'Windows':
                bfile.write('"%s" "%s"\n' % (os.path.join(bin_path,'Searcher'), pfind_file))
            else:
                bfile.write('"%s" "%s"\n' % (os.path.join(bin_path,'Importer'), pfind_file))
                if 'mpicores' in initparam:
                    mpicores = initparam['mpicores']
                else:
                    mpicores = ClusterSetting.DefaultCores
                if 'mpihosts' in initparam:
                    mpihosts = initparam['mpihosts']
                else:
                    mpihosts = ClusterSetting.DefaultHosts
                bfile.write('mpirun -np %s -host %s "%s" "%s"\n' %(mpicores, mpihosts, os.path.join(bin_path,'MPISearcher'), pfind_file))
    
    if not "index_content" in initparam:
        initparam["index_content"] = "PEPTIDE_PAIR"
        
    if not 'noninterexport' in initparam:
        initparam['noninterexport'] = 'false'
        
    if 'evalue_max' in initparam and float(initparam['evalue_max']) < 1:#added at 2013.10.22, evalue cut, fan
        bEvalueMax = True
    else:
        bEvalueMax = False        
    
    if initparam['index_content'] != "PEPTIDE_TRI_ALL":
        Filter.RunForType(bfile, 'inter', initparam, spectra_list, normal_dir_path, bEvalueMax)
        
        if initparam['noninterexport'] == 'true':
            Filter.RunForType(bfile, 'loop', initparam, spectra_list, normal_dir_path, bEvalueMax)
            Filter.RunForType(bfile, 'mono', initparam, spectra_list, normal_dir_path, bEvalueMax)
            Filter.RunForType(bfile, 'common', initparam, spectra_list, normal_dir_path, bEvalueMax)
            
    else:
        Filter.RunForTri(bfile, initparam, spectra_list, normal_dir_path, bEvalueMax);

    bfile.close()
    if platform.system() == 'Linux':
        os.system('chmod 766 %s' % bat_file)
        
    os.system('"%s"' % bat_file)
    search_time = os.stat(pfind_result_file).st_mtime

    return search_time
