#!/usr/bin/env python3
# repair_ss.py — ремонт одиночных цепей, разорванных сборкой NSP на куски.
# Аналог проверенного repair_ds (Задание 20-E), но для ss-аптамеров:
#   fit/pred pdb -> переписать ВСЕ остатки как ОДНУ цепь A 1..L (без TER внутри)
#   -> tleap сам достраивает фосфодиэфирные мостики через растянутые стыки
#   -> усиленная минимизация (maxcyc 2000, ncyc 1000 — как в repair_ds)
#   -> приёмка: заряд == -(L-1) И каждый стык O3'..O5' <= 3.05 A.
#
# Режим 1 (тест одного):   python3 repair_ss.py --one <pdb> <seq>
# Режим 2 (массовый ремонт провалов из apt3d_results.csv, с чанками):
#   PIPE_CHUNK=0 PIPE_NCHUNK=3 python3 repair_ss.py
import csv
import os
import re
import shutil
import subprocess
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

APT_ROOT = os.path.expanduser('~/3D_DNA/apt_build')
TASKS_CSV = os.path.join(APT_ROOT, 'apt_tasks.csv')
RESULTS   = os.path.join(APT_ROOT, 'apt3d_results.csv')

MAX_WORKERS = 12                    # как в основном раннере: поднять под nproc
NSP_BIN = ('/home/ekspert/3D_DNA/home/yiduo/3dRNA/NSP/cirRNA_and_DNA/nsp')

CHUNK  = int(os.environ.get('PIPE_CHUNK', '0'))
NCHUNK = int(os.environ.get('PIPE_NCHUNK', '1'))

NT = {'A': 'DA', 'C': 'DC', 'G': 'DG', 'T': 'DT'}

LEAP = '''source leaprc.DNA.OL15
model = loadpdb repaired.pdb
check model
saveamberparm model prmtop inpcrd
quit
'''

MIN = '''repair ss: minimisation (усиленная — стягиваем растянутые стыки)
&cntrl
  imin   = 1,
  ntxo   = 1,
  maxcyc = 2000,
  ncyc   = 1000,
  ntb    = 0,
  igb    = 5,
  cut    = 12,
 /
'''


def sh(cmd, cwd, timeout=900):
    subprocess.run(['bash', '-c', cmd], cwd=cwd, timeout=timeout)


def dist(a, b):
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


def residues_of(pdb):
    """Остатки в порядке следования: [(raw_key, resname, [строки])]."""
    res = []
    with open(pdb) as f:
        for l in f:
            if not l.startswith('ATOM'):
                continue
            key = l[17:27]
            if not res or res[-1][0] != key:
                res.append([key, l[17:20].strip(), []])
            res[-1][2].append(l)
    return res


def junctions_ok(pdb, cutoff=3.05):
    """Все стыки цепи не длиннее cutoff; возвращает (ok, описание, худший)."""
    res = {}
    with open(pdb) as f:
        for l in f:
            if not l.startswith('ATOM'):
                continue
            res.setdefault(int(l[22:26]), {})[l[12:16].strip()] = \
                (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    nums = sorted(res)
    worst, worst_at = 0.0, ''
    for n1, n2 in zip(nums, nums[1:]):
        if n2 == n1 + 1 and "O3'" in res[n1] and "O5'" in res[n2]:
            d = dist(res[n1]["O3'"], res[n2]["O5'"])
            if d > worst:
                worst, worst_at = d, f'{n1}->{n2}'
    ok = worst <= cutoff
    return ok, f'худший стык {worst_at}: {worst:.2f} A', worst


def repair_one(name, seq, src_pdb, workdir, out_name=None):
    """Ремонт одной цепи. Возвращает dict с ok/reason/charge/worst/seconds."""
    t0 = time.time()
    L = len(seq)
    out_name = out_name or f'{name}.pred1.min.pdb'
    rdir = os.path.join(workdir, 'repair')
    # КРИТИЧНО: чистим папку целиком — иначе повторный запуск подхватывает
    # вчерашние prmtop/min.rst и «проходит» за 0.3 с на мусоре (боевой баг)
    shutil.rmtree(rdir, ignore_errors=True)
    os.makedirs(rdir, exist_ok=True)

    # Сырой pred1.pdb из assemble несёт 5'-фосфат на первом остатке — tleap
    # на нём падает ("DG5 ... P does not have a type"). Боевой шаг `nsp fit`
    # этот фосфат снимает. Если вход уже fit.pdb — повторно не прогоняем.
    if os.path.basename(src_pdb) != 'fit.pdb':
        fitted = os.path.join(rdir, 'fit.pdb')
        sh(f'source ~/nsp_env.sh && {NSP_BIN} fit "{src_pdb}" -o "{fitted}" '
           f'>> {os.path.join(rdir, "tleap_repair.log")} 2>&1', rdir)
        if os.path.exists(fitted):
            src_pdb = fitted

    res = residues_of(src_pdb)
    if len(res) != L:
        return {'name': name, 'ok': False,
                'reason': f'остатков {len(res)} != L={L}'}
    got = [r[1] for r in res]
    exp = [NT[c] for c in seq]
    if got != exp:
        return {'name': name, 'ok': False,
                'reason': 'порядок остатков не совпал с seq (куски перепутаны?)'}

    # переписываем «по замыслу»: ОДНА цепь A, сплошная нумерация 1..L
    rep = os.path.join(rdir, 'repaired.pdb')
    with open(rep, 'w') as f:
        for i, (_, _, lines) in enumerate(res, 1):
            for l in lines:
                f.write(l[:21] + 'A' + f'{i:4d}' + l[26:])
        f.write('TER\nEND\n')

    with open(os.path.join(rdir, 'leap.in'), 'w') as f:
        f.write(LEAP)
    with open(os.path.join(rdir, 'min.in'), 'w') as f:
        f.write(MIN)

    sh('source ~/nsp_env.sh && tleap -f leap.in > tleap_repair.log 2>&1', rdir)
    if not (os.path.exists(os.path.join(rdir, 'prmtop'))
            and os.path.exists(os.path.join(rdir, 'inpcrd'))):
        return {'name': name, 'ok': False, 'reason': 'tleap: нет prmtop/inpcrd'}

    sh('source ~/nsp_env.sh && msander -O -i min.in -o min.out '
       '-p prmtop -c inpcrd -r min.rst >> tleap_repair.log 2>&1', rdir)
    if not os.path.exists(os.path.join(rdir, 'min.rst')):
        return {'name': name, 'ok': False, 'reason': 'msander: нет min.rst'}

    sh('source ~/nsp_env.sh && ambpdb -p prmtop < min.rst > min_raw.pdb '
       '2>> tleap_repair.log', rdir)
    min_pdb = os.path.join(workdir, out_name)
    sh(f'python3 ~/normalize_pdb.py min_raw.pdb "{min_pdb}" '
       '>> tleap_repair.log 2>&1', rdir)
    if not os.path.exists(min_pdb):
        return {'name': name, 'ok': False, 'reason': 'normalize: нет итога'}

    atoms = 0
    with open(min_pdb) as f:
        atoms = sum(1 for l in f if l.startswith('ATOM'))

    charge = None
    with open(os.path.join(rdir, 'tleap_repair.log'), errors='replace') as f:
        for line in f:
            m = (re.search(r'unperturbed charge of the unit \((-?\d+\.?\d*)\)', line)
                 or re.search(r'[Uu]nperturbed charge[:\s]+(-?\d+\.?\d*)', line))
            if m:
                charge = float(m.group(1))

    expected = -(L - 1)
    j_ok, j_why, worst = junctions_ok(min_pdb)
    ok = (charge is not None and abs(charge - expected) < 1e-6
          and j_ok and atoms >= 50)
    reason = '' if ok else (j_why if not j_ok
                            else (f'заряд {charge} != {expected}'
                                  if charge is None or abs(charge - expected) >= 1e-6
                                  else f'мало атомов: {atoms}'))
    return {'name': name, 'ok': ok, 'L': L, 'atoms': atoms, 'charge': charge,
            'expected': expected, 'worst': round(worst, 2),
            'seconds': round(time.time() - t0, 1), 'reason': reason,
            'min_pdb': min_pdb}


def task_pdb(workdir, name):
    """Откуда брать исходник: сперва fit.pdb (как в боевом run_min), потом pred1.pdb."""
    for cand in (os.path.join(workdir, f'{name}.pred1', 'fit.pdb'),
                 os.path.join(workdir, f'{name}.pred1.pdb')):
        if os.path.exists(cand):
            return cand
    return None


def my_task(name):
    return zlib.crc32(name.encode()) % NCHUNK == CHUNK


def main_one():
    pdb, seq = sys.argv[2], sys.argv[3].upper()
    workdir = os.path.dirname(os.path.abspath(pdb)) \
        if os.path.dirname(os.path.abspath(pdb)) else os.getcwd()
    r = repair_one('repair_test', seq, pdb, workdir,
                   out_name='repair_test.min.pdb')
    print('РЕЗУЛЬТАТ РЕМОНТА:', r)
    if r['ok']:
        print('итог:', r['min_pdb'])
        print('проверь стыки: python3 ~/Scripts/junction.py', r['min_pdb'])
    else:
        print('Не вышло:', r['reason'], '— пришли мне этот вывод.')


def _repair_worker(item):
    name, seq, src, wd = item
    try:
        return repair_one(name, seq, src, wd)
    except Exception as e:
        return {'name': name, 'ok': False, 'reason': f'exception: {e!r}'}


def main_mass():
    results = []
    if os.path.exists(RESULTS):
        results = list(csv.DictReader(open(RESULTS)))
    ok_names = {r['name'] for r in results if r['ok'] == 'True'}
    tasks = list(csv.DictReader(open(TASKS_CSV)))
    todo = []
    for t in tasks:
        if not my_task(t['id']) or t['id'] in ok_names:
            continue
        wd = os.path.join(APT_ROOT, t['id'])
        src = task_pdb(wd, t['id'])
        if src:
            todo.append((t['id'], t['seq'], src, wd))
    print(f'ЧАНК {CHUNK}/{NCHUNK}: ремонтировать {len(todo)} воркеров={MAX_WORKERS}',
          flush=True)
    if not todo:
        print('НЕЧЕГО РЕМОНТИРОВАТЬ (всё ок или нет исходников)')
        return
    new_file = not os.path.exists(RESULTS) or os.path.getsize(RESULTS) == 0
    out = open(RESULTS, 'a', newline='')
    w = csv.writer(out)
    if new_file:
        w.writerow(['name', 'ok', 'L', 'atoms', 'charge', 'expected',
                    'ret', 'reason', 'seconds'])
        out.flush()
    t0 = time.time()
    n_ok = n_bad = 0
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = [ex.submit(_repair_worker, item) for item in todo]
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            if r.get('ok'):
                n_ok += 1
            else:
                n_bad += 1
                print('НЕ ПОЧИНИЛ:', r.get('name'), '|', r.get('reason'), flush=True)
            print(f"{i}/{len(todo)} {r.get('name')} ok={r.get('ok')} "
                  f"заряд={r.get('charge')} худший={r.get('worst')} "
                  f"за {r.get('seconds')}с", flush=True)
            w.writerow([r.get('name'), r.get('ok'), r.get('L', ''),
                        r.get('atoms', ''), r.get('charge', ''),
                        r.get('expected', ''), 'repair-ss',
                        r.get('reason', ''), r.get('seconds', '')])
            out.flush()
            if i % 20 == 0 or i == len(todo):
                sec = time.time() - t0
                eta = sec / i * (len(todo) - i)
                print(f'== {i}/{len(todo)} ok={n_ok} bad={n_bad} '
                      f'прошло {sec/60:.0f} мин, ETA {eta/3600:.1f} ч', flush=True)
    out.close()
    print(f'РЕМОНТ ГОТОВ: {n_ok} починено, {n_bad} не поддались', flush=True)


if __name__ == '__main__':
    if len(sys.argv) >= 2 and sys.argv[1] == '--one':
        main_one()
    else:
        main_mass()
