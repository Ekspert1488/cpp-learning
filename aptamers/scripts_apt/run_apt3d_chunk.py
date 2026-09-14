#!/usr/bin/env python3
# run_apt3d_chunk.py — массовое построение 3D полноразмерных аптамеров через NSP.
# Вход: ~/3D_DNA/apt_build/apt_tasks.csv (сделан apt3d_prep.py)
# Выход: ~/3D_DNA/apt_build/<id>/… + apt3d_results.csv
# Чанки на несколько компов:  PIPE_CHUNK=0 PIPE_NCHUNK=3 python3 run_apt3d_chunk.py
# Докачка: повторный запуск пропускает задачи с ok=True в apt3d_results.csv.
import csv
import os
import re
import subprocess
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

OUT_ROOT  = os.path.expanduser('~/3D_DNA/apt_build')
TASKS_CSV = os.path.join(OUT_ROOT, 'apt_tasks.csv')
RESULTS   = os.path.join(OUT_ROOT, 'apt3d_results.csv')
RUN_SH    = '/home/ekspert/3D_DNA/run_3dRNA_DNA/test_run_3dRNA/3dRNA/run_ss.sh'
MAX_WORKERS = 12                    # 8-ядерный комп: 7 в счёт, 1 системе
TIMEOUT   = 3600                   # сек на одну структуру (уточним по калибровке!)

CHUNK  = int(os.environ.get('PIPE_CHUNK', '0'))
NCHUNK = int(os.environ.get('PIPE_NCHUNK', '1'))

JOB_PAR = '''_jobid {name}
_version 2.0.0
_lib_version 3dRNA-Lib1
name {name}
seq {seq}
ss {ss}
num 1
_minimize yes
_pred_ss no
loop_building partial_bires
_mol_type dna
_rna_type linear
_routine Default
'''
# ВАЖНО: в файле после последней строки остаётся \n — пустая строка обязательна
# (боевой опыт: без неё run.sh может потерять последний параметр).


def my_task(name):
    """Стабильное чанкование: crc32 не зависит от перезапуска/машины."""
    return zlib.crc32(name.encode()) % NCHUNK == CHUNK


def build_apt(task):
    name, seq, ss, L = task
    workdir = os.path.join(OUT_ROOT, name)
    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, 'job.par'), 'w') as f:
        f.write(JOB_PAR.format(name=name, seq=seq, ss=ss))
    cmd = f'source ~/nsp_env.sh && bash {RUN_SH} job.par > {name}.run.log 2>&1'
    t0 = time.time()
    try:
        # start_new_session: при таймауте убиваем ВСЮ группу (run.sh + дети),
        # иначе брошенные perl/ambpdb копятся и съедают CPU
        p = subprocess.run(['bash', '-c', cmd], cwd=workdir, timeout=TIMEOUT,
                           start_new_session=True)
    except subprocess.TimeoutExpired:
        subprocess.run(['bash', '-c', f'pkill -9 -f "{name}.pred" 2>/dev/null; true'],
                       timeout=10)
        return {'name': name, 'ok': False, 'reason': f'timeout(>{TIMEOUT}c)',
                'seconds': TIMEOUT}
    sec = round(time.time() - t0, 1)
    min_pdb = os.path.join(workdir, f'{name}.pred1.min.pdb')
    if not os.path.exists(min_pdb):
        return {'name': name, 'ok': False, 'reason': 'no min.pdb',
                'ret': p.returncode, 'seconds': sec}
    with open(min_pdb) as f:
        atoms = sum(1 for l in f if l.startswith('ATOM'))
    if atoms < 50:
        return {'name': name, 'ok': False, 'reason': f'мало атомов: {atoms}',
                'ret': p.returncode, 'seconds': sec}
    charge = None
    leap_log = os.path.join(workdir, f'{name}.pred1', 'leap.log')
    if os.path.exists(leap_log):
        with open(leap_log, errors='replace') as f:
            for line in f:
                # терпим ОБА формата tleap: "unperturbed charge of the unit (-49.0)"
                # и "total unperturbed charge: -49.0"
                m = (re.search(r'unperturbed charge of the unit \((-?\d+\.?\d*)\)', line)
                     or re.search(r'[Uu]nperturbed charge[:\s]+(-?\d+\.?\d*)', line))
                if m:
                    charge = float(m.group(1))
    expected = -(L - 1)                        # одиночная цепь: -(L-1)
    ok = charge is not None and abs(charge - expected) < 1e-6
    if ok:
        # траектория сборки ~50 МБ на задачу — после успеха она мусор
        # (1451 × 50 МБ ≈ 75 ГБ; для докинга нужен только min.pdb)
        traj = os.path.join(workdir, f'{name}.traj.p1.pdb')
        try:
            if os.path.exists(traj):
                os.remove(traj)
        except OSError:
            pass
    return {'name': name, 'ok': ok, 'L': L, 'atoms': atoms, 'charge': charge,
            'expected': expected, 'ret': p.returncode, 'seconds': sec,
            'reason': '' if ok else f'заряд {charge} != {expected}'}


def worker(task):
    try:
        return build_apt(task)
    except Exception as e:
        return {'name': task[0], 'ok': False, 'reason': f'exception: {e!r}'}


def main():
    if not os.path.exists(TASKS_CSV):
        print('НЕТ задач:', TASKS_CSV, '— сначала запусти apt3d_prep.py')
        sys.exit(1)
    # опциональный лимит: `python3 run_apt3d_chunk.py 3` — дым/калибровка
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    tasks = []
    with open(TASKS_CSV) as f:
        for row in csv.DictReader(f):
            if my_task(row['id']):
                tasks.append((row['id'], row['seq'], row['ss'], int(row['L'])))
    if limit:
        # для калибровки берём короткую, среднюю и длинную
        by_len = sorted(tasks, key=lambda t: t[3])
        picks = [by_len[0], by_len[len(by_len) // 2], by_len[-1]]
        tasks = picks[:min(limit, len(picks))]
    print(f'ЧАНК {CHUNK}/{NCHUNK}: моих задач {len(tasks)}', flush=True)

    done = set()
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            for row in csv.DictReader(f):
                if row['ok'] == 'True':
                    done.add(row['name'])
    todo = [t for t in tasks if t[0] not in done]
    print(f'готово {len(done)} | осталось {len(todo)}', flush=True)
    if not todo:
        print('ВСЁ УЖЕ СДЕЛАНО')
        return

    new_file = not os.path.exists(RESULTS) or os.path.getsize(RESULTS) == 0
    with open(RESULTS, 'a', newline='') as out:
        w = csv.writer(out)
        if new_file:
            w.writerow(['name', 'ok', 'L', 'atoms', 'charge', 'expected',
                        'ret', 'reason', 'seconds'])
            out.flush()
        t0 = time.time()
        n_ok = n_bad = 0
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(worker, t): t for t in todo}
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                if r.get('ok'):
                    n_ok += 1
                else:
                    n_bad += 1
                    print('ПРОБЛЕМА:', r.get('name'), '|', r.get('reason'), flush=True)
                w.writerow([r.get('name'), r.get('ok'), r.get('L', ''),
                            r.get('atoms', ''), r.get('charge', ''),
                            r.get('expected', ''), r.get('ret', ''),
                            r.get('reason', ''), r.get('seconds', '')])
                out.flush()
                if i % 10 == 0 or i == len(todo):
                    sec = time.time() - t0
                    eta = sec / i * (len(todo) - i)
                    print(f'{i}/{len(todo)} ok={n_ok} bad={n_bad} '
                          f'прошло {sec/60:.0f} мин, ETA {eta/3600:.1f} ч',
                          flush=True)
    print('ГОТОВО:', n_ok, 'ok,', n_bad, 'проблемных', flush=True)


if __name__ == '__main__':
    main()
