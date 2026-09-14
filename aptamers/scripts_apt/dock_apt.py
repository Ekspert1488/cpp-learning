#!/usr/bin/env python3
# Этап 15: HDOCK всех аптамеров из ~/3D_DNA/apt_build + скоринг с эпитопным фильтром.
# Это адаптация боевого dock_all_chunk.py под аптамерную кампанию:
#   - источник задач: apt3d_results.csv (ok=True) вместо frag_build/build_results.csv
#   - без чанкования (одна машина), JOBS = 12
#   - те же EPITOPE/EXTRA, тот же порог 7.0 A, тот же createpl -nmax 100
# Запуск:  python3 dock_apt.py dock [limit]     -> ~/docking_apt/<id>/job.out
#          python3 dock_apt.py rescore [limit]  -> ~/docking_apt/scores_apt.csv
import csv
import glob
import os
import subprocess
import sys
import time
import multiprocessing as mp

HDOCK    = os.path.expanduser('~/HDOCKlite/hdock')
CREATEPL = os.path.expanduser('~/HDOCKlite/createpl')
RECEPTOR = os.path.expanduser('~/Receptor/Receptor_only_A2A.pdb')
BUILD    = os.path.expanduser('~/3D_DNA/apt_build')
RESULTS  = os.path.join(BUILD, 'apt3d_results.csv')
TASKS    = os.path.join(BUILD, 'apt_tasks.csv')
OUTDIR   = os.path.expanduser('~/docking_apt')
CSV_OUT  = os.path.join(OUTDIR, 'scores_apt.csv')
JOBS     = int(os.environ.get('DOCK_JOBS', '12'))

# те же списки остатков рецептора, что в фрагментной кампании (9XQB chain C)
EXTRA = {74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
         151, 152, 153, 154, 155, 156, 157, 158, 159, 160,
         161, 162, 163, 164, 165, 166, 167, 168, 169, 170,
         171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
         368, 369}
EPITOPE = {80, 81, 152, 154, 155, 156, 157,
           161, 162, 163, 165, 166, 167, 168, 169}


def load_tasks():
    """ok-аптамеры с реально существующим min.pdb."""
    seqs = {}
    if os.path.exists(TASKS):
        for row in csv.DictReader(open(TASKS)):
            seqs[row['id']] = row['seq']
    ok = {}
    with open(RESULTS) as f:
        for row in csv.DictReader(f):
            if str(row.get('ok')) == 'True':
                ok[row['name']] = row
    tasks, missing = [], 0
    for name in sorted(ok):
        pdb = os.path.join(BUILD, name, f'{name}.pred1.min.pdb')
        if os.path.exists(pdb):
            L = int(ok[name].get('L') or len(seqs.get(name, '')))
            tasks.append((name, 'ss', L, 'mcts', pdb))
        else:
            missing += 1
    if missing:
        print(f'ВНИМАНИЕ: у {missing} ok-задач нет min.pdb (пропущены)')
    return tasks


# ------------------ фаза 1: докинг ------------------

def dock_one(t):
    name, kind, L, src, pdb = t
    workdir = os.path.join(OUTDIR, name)
    out_file = os.path.join(workdir, 'job.out')
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return name, True, 'уже сделан'
    os.makedirs(workdir, exist_ok=True)
    rec = os.path.join(workdir, 'receptor.pdb')
    lig = os.path.join(workdir, 'ligand.pdb')
    if not os.path.exists(rec):
        os.symlink(RECEPTOR, rec)
    if not os.path.exists(lig):
        os.symlink(pdb, lig)
    try:
        p = subprocess.run([HDOCK, 'receptor.pdb', 'ligand.pdb', '-out', 'job.out'],
                           cwd=workdir, stdout=subprocess.DEVNULL, timeout=3600,
                           start_new_session=True)
    except subprocess.TimeoutExpired:
        return name, False, 'timeout'
    if p.returncode != 0:
        return name, False, f'hdock exit {p.returncode}'
    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        return name, False, 'пустой job.out'
    return name, True, ''


def run_docking(limit=None):
    tasks = load_tasks()
    if limit:
        tasks = tasks[:limit]
    print('задач на докинг:', len(tasks), f'(воркеров {JOBS})', flush=True)
    t0 = time.time()
    n_ok = n_bad = 0
    with mp.Pool(JOBS) as pool:
        for i, (name, ok, why) in enumerate(pool.imap_unordered(dock_one, tasks), 1):
            if ok:
                n_ok += 1
            else:
                n_bad += 1
                print('ПРОБЛЕМА:', name, '|', why, flush=True)
            if i % 25 == 0 or i == len(tasks):
                sec = time.time() - t0
                eta = sec / i * (len(tasks) - i)
                print(f'{i}/{len(tasks)} ok={n_ok} bad={n_bad} '
                      f'прошло {sec/60:.0f} мин, ETA {eta/3600:.1f} ч', flush=True)
    print('ДОКИНГ ГОТОВ:', n_ok, 'ok,', n_bad, 'проблемных')


# ------------------ фаза 2: скоринг + эпитоп ------------------

def _blocks(lines):
    blocks, cur = [], []
    for l in lines:
        if l.startswith('ATOM'):
            cur.append(l)
        elif cur:
            blocks.append(cur)
            cur = []
    if cur:
        blocks.append(cur)
    return blocks


def _mindist(coord, xyz_list):
    x, y, z = coord
    return min(((x - a) ** 2 + (y - b) ** 2 + (z - c) ** 2) ** 0.5
               for a, b, c in xyz_list)


def rescore_one(t):
    name, kind, L, src, pdb = t
    workdir = os.path.join(OUTDIR, name)
    out_file = os.path.join(workdir, 'job.out')
    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        return None
    for f in glob.glob(os.path.join(workdir, 'model_*.pdb')):
        os.remove(f)
    subprocess.run([CREATEPL, 'job.out', 'top.pdb', '-nmax', '100',
                    '-complex', '-models'], cwd=workdir,
                   stdout=subprocess.DEVNULL)
    list_extra, list_epi = [], []
    n_near_extra = n_near_epi = 0
    for i in range(1, 101):
        path = os.path.join(workdir, f'model_{i}.pdb')
        if not os.path.exists(path):
            break
        lines = open(path).readlines()
        try:
            score = float(lines[3].split()[2])
        except (IndexError, ValueError):
            continue
        recep_res, lig_atoms = [], []
        for b in _blocks(lines):
            chain = b[0].split()[4]
            for l in b:
                p = l.split()
                try:
                    if chain == 'C':
                        r = int(p[5])
                        if r in EXTRA or r in EPITOPE:
                            recep_res.append((r, float(p[6]), float(p[7]), float(p[8])))
                    else:
                        lig_atoms.append((float(p[6]), float(p[7]), float(p[8])))
                except (IndexError, ValueError):
                    continue
        ex_xyz = [(x, y, z) for r, x, y, z in recep_res if r in EXTRA]
        ep_xyz = [(x, y, z) for r, x, y, z in recep_res if r in EPITOPE]
        if not ex_xyz or not lig_atoms:
            continue
        d_ex = [d for d in (_mindist(c, ex_xyz) for c in lig_atoms) if d <= 7.0]
        d_ep = [d for d in (_mindist(c, ep_xyz) for c in lig_atoms) if d <= 7.0]
        if d_ex:
            n_near_extra += 1
            list_extra.append((score, min(d_ex)))
        if d_ep:
            n_near_epi += 1
            list_epi.append((score, min(d_ep)))
    for f in glob.glob(os.path.join(workdir, 'model_*.pdb')):
        os.remove(f)
    s_ex, d_ex = min(list_extra, key=lambda v: v[0], default=(None, None))
    s_ep, d_ep = min(list_epi,   key=lambda v: v[0], default=(None, None))
    return (name, kind, L, src, s_ex, d_ex, n_near_extra, s_ep, d_ep, n_near_epi)


def run_rescore(limit=None):
    tasks = load_tasks()
    if limit:
        tasks = tasks[:limit]
    done = set()
    if os.path.exists(CSV_OUT):
        with open(CSV_OUT) as f:
            done = {row['name'] for row in csv.DictReader(f)}
    todo = [t for t in tasks if t[0] not in done and
            os.path.exists(os.path.join(OUTDIR, t[0], 'job.out'))]
    print('задач на скоринг:', len(todo), flush=True)
    new_file = not os.path.exists(CSV_OUT) or os.path.getsize(CSV_OUT) == 0
    out = open(CSV_OUT, 'a', newline='')
    w = csv.writer(out)
    if new_file:
        w.writerow(['name', 'kind', 'L', 'src',
                    'score_extra', 'dist_extra', 'n_near_extra',
                    'score_epi', 'dist_epi', 'n_near_epi'])
        out.flush()
    t0 = time.time()
    with mp.Pool(JOBS) as pool:
        for i, res in enumerate(pool.imap_unordered(rescore_one, todo), 1):
            if res is not None:
                w.writerow(res)
                out.flush()
            if i % 50 == 0 or i == len(todo):
                sec = time.time() - t0
                eta = sec / i * (len(todo) - i)
                print(f'{i}/{len(todo)} прошло {sec/60:.0f} мин, '
                      f'ETA {eta/60:.0f} мин', flush=True)
    out.close()
    print('СКОРИНГ ГОТОВ ->', CSV_OUT)


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else ''
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else None
    os.makedirs(OUTDIR, exist_ok=True)
    if mode == 'dock':
        run_docking(limit)
    elif mode == 'rescore':
        run_rescore(limit)
    else:
        print('использование: python3 dock_apt.py dock [limit] | rescore [limit]')
