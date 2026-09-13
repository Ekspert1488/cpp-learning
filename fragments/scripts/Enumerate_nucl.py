#import mdna
import os
import subprocess
import csv
import time
import shutil
import multiprocessing as mp
import  glob

import re
def add_nucl(count_0 ,gene_0,a):

    count_0 += 1
    nucleotides = ['A','G','C','T']
    if (count_0 <= a):
        for i in range(4):
            print(gene_0 + nucleotides[i])
            add_nucl(count_0,gene_0 + nucleotides[i] , a)



def add_nucl_num (count_0,gene_0,a,list_0):
    count_0 += 1
    nucleotides = ['A', 'G', 'C', 'T']
    if (count_0 <= a):
        for i in range(4):
            gene_copy = gene_0 + nucleotides[i]
            add_nucl_num(count_0, gene_copy, a,list_0)
            if (len(gene_copy ) == a):
                list_0.append(gene_copy)

    return list_0

#len_nucl = 4
#add_nucl_num(0,'',len_nucl)

def list_print_like_str (list):
    str_0 = ''
    for i in range( len(list) ):
        str_0 += list[i]
    return str_0

def complem_chain (nucl):
    nucl = list(nucl)
    chain_len = len(nucl)
    for  i in range(chain_len):
        if (nucl[i] == 'A'):
            nucl[i] = 'T'
        elif (nucl[i] == 'T'):
            nucl[i] = 'A'
        elif (nucl[i] == 'G'):
            nucl[i] = 'C'
        elif (nucl[i] == 'C'):
            nucl[i] = 'G'
    reversed_chain = nucl[::-1]
    return  list_print_like_str(reversed_chain)

def result_chains (num):
    chains_complete_list = []
    list = add_nucl_num(0,'',num,[])
    for i in range(len(list)):
        chains_complete = list[i] + '&' + complem_chain(list[i])
        chains_complete_list.append(chains_complete)
    return chains_complete_list

full_library_for_mdna = add_nucl_num(0,'',3,[]) + add_nucl_num(0,'',4,[]) + add_nucl_num(0,'',5,[]) +add_nucl_num(0,'',6,[])
''' def mdna_chains (list):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_folder = os.path.join(script_dir, 'duplex_library')
    os.makedirs(output_folder, exist_ok=True)

    for i in range(len(list)):
        print('Загружаю цепь:', list[i])
        dna = mdna.make(sequence=list[i])
        print('Сохраняю в pdb файл')

        # Вот тут ключевая правка - собираем полный путь через output_folder
        file_path = os.path.join(output_folder, f'fragment_{list[i]}_duplex.pdb')
        dna.save_pdb(file_path)

        print("Готово! Последовательность:", dna.sequence)'''
#mdna_chains(full_library_for_mdna)

def unique_duplexes(seq_list):
    s = set()
    for i in range (len(seq_list)):
        comp = complem_chain(seq_list[i])
        if (seq_list[i] < comp):
            s.add(seq_list[i])

        else:
            s.add(comp)

    return sorted(s)
unique_chains = unique_duplexes(full_library_for_mdna)
os.path.expanduser('~/Scripts/duplex_library')
list_of_duplexes = sorted(os.listdir('/home/ekspert/Scripts/duplex_library'))

def find_duplexes(library_path):
    result = {}
    list_of_duplexes = sorted(os.listdir(library_path))
    for i in  range(len(list_of_duplexes)):
        file_name = list_of_duplexes[i]
        if file_name.startswith('fragment_'):
            position_end = file_name.find('_duplex')
            if (position_end != -1 ):
                file_name = file_name[9:position_end]
                result[file_name] = os.path.join(library_path, list_of_duplexes[i])

    return result

def build_worklist(fragments):
    result = []
    unique_chains = unique_duplexes(list(fragments.keys()))
    for i in range (len(unique_chains)):
        chain = unique_chains[i]
        path = fragments.get(chain)
        if path is None:
            print('ВНИМАНИЕ: нет файла для', chain, '- пропускаю')
            continue
        result.append( (chain, fragments[chain]))
    return result

def run_docking(seq, ligand_path, receptor_path, hdock_path,outdir, angle = None):
    if not os.path.isfile(receptor_path):
        print('Рецептор не найден:', receptor_path)
        return None
    if not os.path.isfile(ligand_path):
        print('Лиганд не найден:', ligand_path)
        return None
    workdir = os.path.join(outdir, seq)
    os.makedirs(workdir, exist_ok=True)
    out_file = os.path.join(workdir, 'job.out')
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return out_file
    shutil.copy(receptor_path, os.path.join(workdir, 'receptorA2A.pdb'))
    shutil.copy(ligand_path, os.path.join(workdir, f'ligand_{seq}.pdb'))
    cmd = [hdock_path, 'receptorA2A.pdb', f'ligand_{seq}.pdb']
    if angle is not None:
        cmd += ['-angle', str(angle)]

    cmd += ['-out', 'job.out']
    p = subprocess.run(cmd , cwd=workdir,
                       stdout=subprocess.DEVNULL)
    if (p.returncode != 0) :
        print('hdock упал на', seq)
        return None
    out_file = os.path.join(workdir, 'job.out')
    if   not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        return None

    return out_file

def parse_out(out_path):
    lines = open(out_path).readlines()
    if len(lines) < 6 :
        print(f'Файл ошибся в записи. Он записал ровно {len(lines)} строк')
        return None
    useful_lines = lines[5:]
    min_score = None
    for i in range (len(useful_lines)):
        parts = useful_lines[i].split()
        if len(parts) < 7 :
            continue
        try:
            score = float(parts[6])
        except ValueError:
            continue
        if min_score is None or score < min_score:
            min_score = score
    return min_score

def run_all(library_path, receptor_path, hdock_path, outdir, angle=None, limit=None):
    fragments = find_duplexes(library_path)
    worklists = build_worklist(fragments)
    seq_dict = {}
    if limit is not None :
        worklists = worklists[:limit]
    csv_path = os.path.join(outdir, 'scores.csv')
    if not os.path.exists(csv_path):
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['sequence', 'length', 'score', 'out_file'])  # шапка
    t_start = time.time()
    for i in range(len(worklists)):

        seq = worklists[i][0]
        ligand_path = worklists[i][1]
        complem_seq = complem_chain(seq)

        out_path = run_docking(seq, ligand_path,receptor_path,hdock_path,outdir,angle)
        if out_path is None:
            continue
        score = parse_out(out_path)
        if score is None:
            continue
        seq_dict[seq] = score
        if complem_seq != seq :
            with open(csv_path, 'a', newline='') as f:
                w = csv.writer(f)
                w.writerow([complem_seq, len(complem_seq), score, out_path])
        with open(csv_path, 'a', newline='') as f:
            w = csv.writer(f)
            w.writerow([seq, len(seq), score, out_path])
        print(f'[{i + 1}/{len(worklists)}] {seq} score={score}')
        elapsed = time.time() - t_start
        rate = (i + 1) / elapsed  # задач в секунду
        eta = (len(worklists) - i - 1) / rate  # секунд осталось
        if i % 10 == 0:

            print(f'... До окончания доккинга осталось {eta / 3600:.1f} ч')
    return seq_dict


def dock_one (data):
    seq, ligand_path, receptor_path, hdock_path, outdir, angle = data
    path = run_docking(seq, ligand_path, receptor_path, hdock_path, outdir, angle)
    if path is None : # Сначала название,потом score а затем путь
        return seq , None, None
    else:
        score = parse_out(path)
        return seq,score,path


def run_all_parallel(library_path, receptor_path, hdock_path, outdir, angle=None, limit=None, jobs=None):
    if jobs is None:
        seq_dict = run_all(library_path,receptor_path,hdock_path,outdir,angle,limit)
        return seq_dict
    seq_dict = {}
    fragments = find_duplexes(library_path)
    worklists = build_worklist(fragments)
    if limit is not None :
        worklists = worklists[:limit]
    list_of_data_for_proc = []
    os.makedirs(outdir, exist_ok=True)
    csv_path = os.path.join(outdir, 'scores.csv')
    if not os.path.exists(csv_path):
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['sequence', 'length', 'score', 'out_file'])  # шапка
    t_start = time.time()
    num_ready_dock = 0
    for i in range(len(worklists)):
        chain = worklists[i][0]
        ligand_path = worklists[i][1]
        list_of_data_for_proc.append((chain,ligand_path,receptor_path,hdock_path,outdir,angle))
    with mp.Pool(jobs) as pool:
        for answer in pool.imap_unordered(dock_one, list_of_data_for_proc):
            seq, score , path = answer
            seq_dict[seq] = score
            if score is None:
                continue
            complem_seq = complem_chain(seq)
            if complem_seq != seq:
                with open(csv_path, 'a', newline='') as f:
                    w = csv.writer(f)
                    w.writerow([complem_seq, len(complem_seq), score, path])
            with open(csv_path, 'a', newline='') as f:
                w = csv.writer(f)
                w.writerow([seq, len(seq), score, path])
            print(f'[{num_ready_dock + 1}/{len(worklists)}] {seq} score={score}')
            elapsed = time.time() - t_start
            rate = (num_ready_dock + 1) / elapsed  # задач в секунду
            eta = (len(worklists) - num_ready_dock - 1) / rate  # секунд осталось
            if num_ready_dock % 10 == 0:
                print(f'... До окончания доккинга осталось {eta / 3600:.1f} ч')
            num_ready_dock += 1

    return seq_dict


def min_dist_atoms (atom_xyz , atoms_recep_xyz_list):
    x_atom , y_atom ,z_atom = atom_xyz
    dist_list = []
    for i in range (len(atoms_recep_xyz_list)):
        x_atom_2, y_atom_2 , z_atom_2 = atoms_recep_xyz_list[i]
        dist = ( ( x_atom - x_atom_2   )**2 + (y_atom - y_atom_2 )**2 + (z_atom - z_atom_2)**2 )**0.5
        dist_list.append(dist)
    min_dist = min(dist_list)
    return  min_dist





def rescore_local(out_file , createpl_path, cutoff = 7.0, nmax = 100):
    EXTRA = {74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
             151, 152, 153, 154, 155, 156, 157, 158, 159, 160,
             161, 162, 163, 164, 165, 166, 167, 168, 169, 170,
             171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
             368, 369}  # 43 остатка, 216 атомов

    EPITOPE = {80, 81, 152, 154, 155, 156, 157,
                   161, 162, 163, 165, 166, 167, 168, 169}
    n_near_extra = 0
    n_near_epitop = 0
    workdir = os.path.dirname(out_file)
    for f in glob.glob(os.path.join(workdir, 'model_*.pdb')):
        os.remove(f)
    config = [createpl_path, 'job.out', 'top.pdb', '-nmax', '100', '-complex', '-models']
    p = subprocess.run(config,cwd=workdir)
    min_dist_extra_list = []
    min_dist_epitop_list = []
    for i_models in range(1, nmax+1):
        min_extra_local = []
        min_epitop_local = []
        out_path = os.path.join(workdir,f'model_{i_models}.pdb')
        if not os.path.exists(out_path):
            break
        lines = open(out_path).readlines()
        part_score = lines[3].split()
        score = float(part_score[2])
        start_recep = None
        end_recep = None
        start_ligand_A = None
        end_ligand_A = None
        start_ligand_B = None
        end_ligand_B     = None
        start = 0
        for i_molecula in range(3):
            while  not lines[start].startswith('ATOM') :
                start+=1
            end = start
            while lines[end+1].startswith('ATOM') :
                end+=1
            useful_lines = lines[start:end+1]
            part = useful_lines[0].split()
            if part[4] == 'C':
                start_recep = start
                end_recep = end
            elif part[4] == 'A':
                start_ligand_A = start
                end_ligand_A = end
            elif part[4] == 'B':
                start_ligand_B = start
                end_ligand_B = end
            start = end + 1
        atoms_recep_list = lines[start_recep:end_recep+1]
        atoms_recep_coord_extra = []
        atoms_recep_coord_epitop = []
        atoms_ligand_list = lines[start_ligand_A:end_ligand_A+1] + lines[start_ligand_B:end_ligand_B+1]
        for i_dist in range (len(atoms_recep_list)):
            part = atoms_recep_list[i_dist].split()
            if int(part[5]) in EXTRA:
                x = float(part[6])
                y = float(part[7])
                z = float(part[8])
                atoms_coord = (x, y, z)
                atoms_recep_coord_extra.append(atoms_coord)
            if int(part[5]) in EPITOPE:
                x = float(part[6])
                y = float(part[7])
                z = float(part[8])
                atoms_coord = (x, y, z)
                atoms_recep_coord_epitop.append(atoms_coord)
        for i_ligand in range(len(atoms_ligand_list)):
            part = atoms_ligand_list[i_ligand].split()
            x = float(part[6])
            y = float(part[7])
            z = float(part[8])
            coord_ligand = (x,y,z)
            dist_extra_local = min_dist_atoms(coord_ligand,atoms_recep_coord_extra)
            dist_epitop_local = min_dist_atoms(coord_ligand,atoms_recep_coord_epitop)
            if dist_extra_local <= cutoff :
                min_extra_local.append(dist_extra_local)
            if dist_epitop_local <= cutoff:
                min_epitop_local.append(dist_epitop_local)
        dist_extra = min(min_extra_local, default= None)
        dist_epitop = min(min_epitop_local, default=None)
        if len(min_extra_local) > 0 :
            n_near_extra+=1
            min_dist_extra_list.append((score,dist_extra))
        if len(min_epitop_local) > 0 :
            n_near_epitop+=1
            min_dist_epitop_list.append((score,dist_epitop))


    score_extra, dist_extra = min(min_dist_extra_list, key=lambda t: t[0], default=(None, None))
    score_epitop, dist_epitop = min(min_dist_epitop_list, key=lambda t: t[0], default=(None, None))
    for f in glob.glob(os.path.join(workdir, 'model_*.pdb')):
        os.remove(f)
    return score_extra,dist_extra,n_near_extra, score_epitop,dist_epitop,n_near_epitop

def run_all_sort (create_pl_path, outdir, limit=None, cutoff = None , nmax = None):
    list_dir = os.listdir(outdir)
    results_dict = {}
    worklist = []
    if cutoff is None:
        cutoff = 7
    if nmax is None:
        nmax = 100
    csv_path = os.path.join(outdir, 'scores_sort.csv')
    if not os.path.exists(csv_path):
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['sequence', 'length', 'score_extra', 'dist_extra', 'n_near_extra', 'score_epi', 'dist_epi', 'n_near_epi'])  # шапка

    for i in range (len(list_dir)):
        if os.path.isdir(os.path.join(outdir, list_dir[i])):
            workdir = os.path.join(outdir, list_dir[i])
            out_file = os.path.join(workdir, 'job.out')
            if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
                worklist.append((list_dir[i],out_file))
    if not limit is None:
        worklist=worklist[:limit]
    for j in range(len(worklist)):
        seq = worklist[j][0]
        complem_seq = complem_chain(seq)
        out_file = worklist[j][1]
        length = len(seq)
        results = rescore_local(out_file, create_pl_path, cutoff, nmax)
        score_extra = results[0]
        dist_extra = results[1]
        n_near_extra = results[2]
        score_epi = results[3]
        dist_epi = results[4]
        n_near_epi = results[5]
        if complem_seq != seq :
            with open(csv_path, 'a', newline='') as f:
                w = csv.writer(f)
                w.writerow([complem_seq, length, score_extra, dist_extra, n_near_extra, score_epi, dist_epi, n_near_epi])
        results_dict[seq] = results
        with open(csv_path, 'a', newline='') as f:
            w = csv.writer(f)
            w.writerow([seq, length, score_extra, dist_extra, n_near_extra, score_epi, dist_epi, n_near_epi])
        print(f'[{j + 1}/{len(worklist)}] {seq} score_extra={score_extra},dist_extra = {dist_extra},n_near_extra = {n_near_extra},'
                  f'score_epi = {score_epi}, dist_epi = {dist_epi} , n_near_epi = {n_near_epi}')
    return results_dict


def sort_one (data):
    seq,out_file, createpl_path , cutoff , nmax = data
    if out_file is None:
        return None, None,None,None , None, None, None
    results =  rescore_local(out_file,createpl_path,cutoff,nmax)
    score_extra = results[0]
    dist_extra = results[1]
    n_near_extra = results[2]
    score_epi = results[3]
    dist_epi = results[4]
    n_near_epi = results[5]
    return  seq,score_extra,dist_extra,n_near_extra,score_epi,dist_epi,n_near_epi


def sort_all_parallel (createpl_path, outdir,outfile_csv, limit = None, jobs = 7 , cutoff = None, nmax = None) :
    list_dir = os.listdir(outdir)
    results_dict = {}
    worklist = []
    list_for_proc = []
    num_ready_dock = 0
    if cutoff is None:
        cutoff = 7
    if nmax is None:
        nmax = 100
    csv_path = os.path.join(outdir, outfile_csv)
    if not os.path.exists(csv_path):
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['sequence', 'length', 'score_extra', 'dist_extra', 'n_near_extra', 'score_epi', 'dist_epi',
                        'n_near_epi'])  # шапка
    for i in range (len(list_dir)):
        if os.path.isdir(os.path.join(outdir, list_dir[i])):
            workdir = os.path.join(outdir, list_dir[i])
            out_file = os.path.join(workdir, 'job.out')
            if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
                worklist.append((list_dir[i],out_file))
    if not  limit is None:
        worklist = worklist[:limit]
    t_start = time.time()
    for j in range(len(worklist)):
        seq = worklist[j][0]
        out_file = worklist[j][1]
        list_for_proc.append((seq,out_file,createpl_path,cutoff,nmax))
    with mp.Pool(jobs) as pool:
        for answer in pool.imap_unordered(sort_one, list_for_proc):
            seq , score_extra , dist_extra , n_near_extra , score_epi,dist_epi, n_near_epi = answer
            if (score_extra is None) or (score_epi is None):
                continue
            length = len(seq)
            complem_seq = complem_chain(seq)
            if complem_seq != seq:
                with open(csv_path, 'a', newline='') as f:
                    w = csv.writer(f)
                    w.writerow(
                        [complem_seq, length, score_extra, dist_extra, n_near_extra, score_epi, dist_epi, n_near_epi])

            results_dict[seq] = answer
            with open(csv_path, 'a', newline='') as f:
                w = csv.writer(f)
                w.writerow([seq, length, score_extra, dist_extra, n_near_extra, score_epi, dist_epi, n_near_epi])
            print(
                f'[{num_ready_dock + 1}/{len(worklist)}] {seq} score_extra={score_extra},dist_extra = {dist_extra},n_near_extra = {n_near_extra},'
                f'score_epi = {score_epi}, dist_epi = {dist_epi} , n_near_epi = {n_near_epi}')
            elapsed = time.time() - t_start
            rate = (num_ready_dock + 1) / elapsed  # задач в секунду
            eta = (len(worklist) - num_ready_dock - 1) / rate  # секунд осталось
            if num_ready_dock % 10 == 0:
                print(f'... До окончания сортировки осталось {eta / 60:.1f} min')
            num_ready_dock += 1

    return results_dict

def clean_ter(in_path, out_path):

    if os.path.exists(in_path) and os.path.getsize(in_path) > 0:
        lines = open(in_path).readlines()

    num_of_ter_nec = 0
    num_of_ter = 0
    with open(out_path, 'w', newline='') as f:
        for i in range(len(lines)):
            if lines[i].startswith('ATOM'):
                last_chain = lines[i][21]
            if lines[i].startswith('TER'):
                new_chain = None
                num_of_ter +=1
                j = 1
                while i+j < len(lines) and not lines[i+j].startswith('ATOM'):
                    j+=1
                if i + j < len(lines) and lines[i+j].startswith('ATOM')   :
                    new_chain = lines[i + j][21]
                if last_chain != new_chain:
                    num_of_ter_nec += 1
                    f.write(lines[i])
            else:
                f.write(lines[i])
    print('Количество нужных TER:',num_of_ter_nec)
    print('КОличество TER всего:', num_of_ter)
    print('Удалено TER:',num_of_ter-num_of_ter_nec)
    return num_of_ter - num_of_ter_nec

def make_task_list():
    task_list = []
    list_of_name_ss = []
    list_second_sstr_ss = []

    list_of_sst_dup = []
    list_of_seq_ss = sorted(add_nucl_num(0,'',3,[]) + add_nucl_num(0,'',4,[])+add_nucl_num(0,'',5,[])+add_nucl_num(0,'',6,[]))
    list_of_name_dup = []
    seq_dup = []
    list_of_seq_dup = unique_duplexes(list_of_seq_ss)
    for i in range(len(list_of_seq_ss)):
        length = len(list_of_seq_ss[i])
        if i <  len(list_of_seq_dup):
            length_dup = len( list_of_seq_dup[i])
            dup = list_of_seq_dup[i] + '&' + complem_chain(list_of_seq_dup[i])
            seq_dup.append(dup)
            list_of_name_dup.append(str('ds_' + list_of_seq_dup[i]))
            list_of_sst_dup.append( str( '(' * length_dup + '&' +   ')' * length_dup ) )
            task_list.append((list_of_name_dup[i], seq_dup[i], list_of_sst_dup[i], 'ds'))
        list_second_sstr_ss.append(str('.'* length))
        list_of_name_ss.append(str('ss_' + list_of_seq_ss[i]))
        task_list.append((list_of_name_ss[i],list_of_seq_ss[i],list_second_sstr_ss[i],'ss'))

    return task_list










def build_one(task, out_root):
    RUN_SH = '/home/ekspert/3D_DNA/run_3dRNA_DNA/test_run_3dRNA/3dRNA/run.sh'
    JOB_PAR = '''_jobid {name}
    _version 2.0.0
    _lib_version 3dRNA-Lib1
    name {name}
    seq {seq}
    ss {sst}
    num 1
    num_samplings 1
    _minimize yes
    _pred_ss no
    loop_building partial_bires
    _mol_type dna
    _rna_type linear
    _routine Default
    '''
    name, seq, sst, kind = task
    workdir = os.path.join(out_root, name)
    os.makedirs(workdir, exist_ok=True)          # и корень, и папку; безопасно при повторе

    with open(os.path.join(workdir, 'job.par'), 'w') as f:
        f.write(JOB_PAR.format(name=name, seq=seq, sst=sst))

    cmd = f'source ~/nsp_env.sh && bash {RUN_SH} job.par > {name}.run.log 2>&1'  # окружение — внутри bash
    try:
        p = subprocess.run(['bash', '-c', cmd], cwd=workdir, timeout=600)
    except subprocess.TimeoutExpired:
        return {'name': name, 'ok': False, 'reason': 'timeout'}

    min_pdb = os.path.join(workdir, f'{name}.pred1.min.pdb')
    if not os.path.exists(min_pdb):
        return {'name': name, 'ok': False, 'reason': 'no min.pdb',
                'ret': p.returncode, 'files': sorted(os.listdir(workdir))[:20]}

    charge = None
    leap_log = os.path.join(workdir, f'{name}.pred1', 'leap.log')
    if os.path.exists(leap_log):
        with open(leap_log) as f:
            for line in f:
                m = re.search(r'unperturbed charge of the unit \((-?\d+\.?\d*)\)', line)
                if m:
                    charge = float(m.group(1))

    L = len(seq) if kind == 'ss' else len(seq.split('&')[0])
    expected = -(L - 1) if kind == 'ss' else -2 * (L - 1)

    return {'name': name, 'ok': charge is not None and abs(charge - expected) < 1e-6,
            'kind': kind, 'L': L, 'charge': charge, 'expected': expected,
            'ret': p.returncode, 'min_pdb': min_pdb}
task_test = make_task_list()
task_0 = task_test[0]
build_one(task_0 , '/home/ekspert/test_3D_DNA')
NT = {'A': 'DA', 'C': 'DC', 'G': 'DG', 'T': 'DT'}

LEAP_IN = '''source leaprc.DNA.OL15
model = sequence {{ {residues} }}
check model
saveamberparm model prmtop inpcrd
savepdb model init.pdb
quit
'''

MIN_IN = '''tleap-fallback: minimisation, те же настройки что в конвейере
&cntrl
  imin   = 1,
  ntxo   = 1,
  maxcyc = 1000,
  ncyc   = 500,
  ntb    = 0,
  igb    = 5,
  cut    = 12,
 /
'''


def _sh(cmd, cwd):
    return subprocess.run(['bash', '-c', cmd], cwd=cwd, timeout=600).returncode


def build_one_fallback(task, out_root):
    name, seq, sst, kind = task
    L = len(seq)
    workdir = os.path.join(out_root, name)

    # чистим следы неудачного NSP-прогона, начинаем с чистого листа
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)

    # ACGT -> DA5 DC DG DT3  (5'- и 3'-концевые варианты остатков)
    names = [NT[c] for c in seq]
    names[0] += '5'
    names[-1] += '3'

    with open(os.path.join(workdir, 'leap.in'), 'w') as f:
        f.write(LEAP_IN.format(residues=' '.join(names)))
    with open(os.path.join(workdir, 'min.in'), 'w') as f:
        f.write(MIN_IN)

    # 1) tleap: строит структуру и топологию (код выхода у tleap ненадёжен —
    #    проверяем артефакты)

    _sh('source ~/nsp_env.sh && tleap -f leap.in >> fallback.log 2>&1', workdir)
    if not (os.path.exists(os.path.join(workdir, 'prmtop'))
            and os.path.exists(os.path.join(workdir, 'inpcrd'))):
        return {'name': name, 'ok': False, 'reason': 'tleap: нет prmtop/inpcrd',
                'src': 'tleap'}

    # 2) msander: минимизация (LBFGS-ошибка в min.out, как мы знаем, безобидна)
    _sh('source ~/nsp_env.sh && msander -O -i min.in -o min.out '
        '-p prmtop -c inpcrd -r min.rst >> fallback.log 2>&1', workdir)
    if not os.path.exists(os.path.join(workdir, 'min.rst')):
        return {'name': name, 'ok': False, 'reason': 'msander: нет min.rst',
                'src': 'tleap'}

    # 3) ambpdb -> PDB -> нормализация (единый формат с NSP-веткой)
    _sh('source ~/nsp_env.sh && ambpdb -p prmtop < min.rst > min_raw.pdb 2>> fallback.log', workdir)    

    if not os.path.exists(os.path.join(workdir, 'min_raw.pdb')):
        return {'name': name, 'ok': False, 'reason': 'ambpdb: нет min_raw.pdb',
                'src': 'tleap'}

    min_pdb = os.path.join(workdir, f'{name}.pred1.min.pdb')
    _sh(f'python3 ~/normalize_pdb.py min_raw.pdb {name}.pred1.min.pdb >> fallback.log 2>&1', workdir)
    if not os.path.exists(min_pdb):
        return {'name': name, 'ok': False, 'reason': 'normalize: нет итога',
                'src': 'tleap'}

    # 4) заряд из НАШЕГО leap.log: ожидание -(L-1)
    charge = None
    with open(os.path.join(workdir, 'leap.log')) as f:
        for line in f:
            m = re.search(r'unperturbed charge of the unit \((-?\d+\.?\d*)\)', line)
            if m:
                charge = float(m.group(1))

    expected = -(L - 1)
    return {'name': name,
            'ok': charge is not None and abs(charge - expected) < 1e-6,
            'kind': kind, 'L': L, 'charge': charge, 'expected': expected,
            'min_pdb': min_pdb, 'src': 'tleap'}

'''run_all_sort('/home/ekspert/HDOCKlite/createpl','/home/ekspert/docking_results_FINAL')

'''
'''
sort_all_parallel('/home/ekspert/HDOCKlite/createpl','/home/ekspert/docking_results_FINAL','score_test_nmax_500.csv',
                  nmax = 500)
sort_all_parallel('/home/ekspert/HDOCKlite/createpl','/home/ekspert/docking_results_FINAL','score_test_nmax_1000.csv',
                  nmax = 1000)
sort_all_parallel('/home/ekspert/HDOCKlite/createpl','/home/ekspert/docking_results_FINAL','score_test_nmax_2000.csv',
                  nmax = 2000)

'''


'''run_all_parallel(
    '/home/ekspert/Scripts/duplex_library',
    '/home/ekspert/Receptor/Receptor_only_A2A.pdb',
    '/home/ekspert/HDOCKlite/hdock',
    '/home/ekspert/docking_results_FINAL', jobs = 7)
'''
'''run_all(library_path='/home/ekspert/Scripts/duplex_library',
        receptor_path='/home/ekspert/Receptor/Receptor_only_A2A.pdb',
        hdock_path='/home/ekspert/HDOCKlite/hdock',
        outdir='/home/ekspert/docking_test2',
         limit=3)'''

def run_all_sort_new(create_pl_path, outdir, limit=None, cutoff = None , nmax = None):
    list_dir = os.listdir(outdir)
    results_dict = {}
    worklist = []
    if cutoff is None:
        cutoff = 7
    if nmax is None:
        nmax = 100
    csv_path = os.path.join(outdir, 'scores_sort.csv')
    if not os.path.exists(csv_path):
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['sequence', 'length', 'score_extra', 'dist_extra', 'n_near_extra', 'score_epi', 'dist_epi',
                        'n_near_epi'])  # шапка

    for i in range(len(list_dir)):
        if os.path.isdir(os.path.join(outdir, list_dir[i])):
            workdir = os.path.join(outdir, list_dir[i])
            out_file = os.path.join(workdir, 'job.out')
            if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
                worklist.append((list_dir[i], out_file))
    if not limit is None:
        worklist = worklist[:limit]
    for j in range(len(worklist)):
        seq_full = worklist[j][0]
        seq = seq_full[3:]
        complem_seq = complem_chain(seq)
        out_file = worklist[j][1]
        length = len(seq)
        results = rescore_local(out_file, create_pl_path, cutoff, nmax)
        score_extra = results[0]
        dist_extra = results[1]
        n_near_extra = results[2]
        score_epi = results[3]
        dist_epi = results[4]
        n_near_epi = results[5]
        if (seq_full[:2] == 'ds') and (complem_seq != seq):
            with open(csv_path, 'a', newline='') as f:
                w = csv.writer(f)
                w.writerow(
                    ['ds_' + complem_seq, length, score_extra, dist_extra, n_near_extra, score_epi, dist_epi, n_near_epi])
        results_dict[seq_full] = results
        with open(csv_path, 'a', newline='') as f:
            w = csv.writer(f)
            w.writerow([seq_full, length, score_extra, dist_extra, n_near_extra, score_epi, dist_epi, n_near_epi])
        print(
            f'[{j + 1}/{len(worklist)}] {seq_full} score_extra={score_extra},dist_extra = {dist_extra},n_near_extra = {n_near_extra},'
            f'score_epi = {score_epi}, dist_epi = {dist_epi} , n_near_epi = {n_near_epi}')
    return results_dict


'''sort_all_parallel('/home/ekspert/HDOCKlite/createpl','/home/ekspert/docking_H','score_final_with_H_nmax_100.csv',
                  nmax = 100)'''







