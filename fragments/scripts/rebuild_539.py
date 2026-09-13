#!/usr/bin/env python3
# Задание 22: пересборка пустых pred1.min.pdb (запускать после патча ambpdb!)
import os, sys, csv, time, traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Enumerate_nucl import make_task_list, build_one, build_one_fallback
try:
    from repair_ds import repair_one
except ImportError:
    repair_one = None

FRAG    = os.path.expanduser('~/3D_DNA/frag_build')
OUT_CSV = os.path.join(FRAG, 'rebuild_results.csv')

def n_atoms(pdb):
    try:
        with open(pdb) as f:
            return sum(1 for L in f if L.startswith('ATOM'))
    except OSError:
        return 0

def pred_path(name):
    return os.path.join(FRAG, name, name + '.pred1.min.pdb')

tasks = {t[0]: t for t in make_task_list()}
names = sorted(n for n in os.listdir(FRAG)
               if os.path.isdir(os.path.join(FRAG, n)))
bad = [n for n in names if n_atoms(pred_path(n)) < 50]
print('задач в make_task_list:', len(tasks))
print('пустых/битых pred1.min.pdb:', len(bad))
miss = [n for n in bad if n not in tasks]
if miss:
    print('ВНИМАНИЕ, нет в make_task_list:', miss[:10])

new = open(OUT_CSV, 'a', newline='')
w = csv.writer(new)
if new.tell() == 0:
    w.writerow(['name', 'ok', 'src', 'atoms'])

t0 = time.time(); ok_cnt = 0
for i, name in enumerate(bad, 1):
    if name not in tasks:
        continue
    t = tasks[name]
    L = len(t[1].split('&')[0])
    src = ''
    try:
        if t[3] == 'ds':
            src = 'rebuild-ds'
            ok = repair_one(t) if repair_one else False
            if not ok:
                src = 'rebuild-ds-nsp'
                ok = build_one(t, FRAG)
        elif L <= 4:
            src = 'rebuild-ss-fb'
            ok = build_one_fallback(t, FRAG)
        else:
            src = 'rebuild-ss-nsp'
            ok = build_one(t, FRAG)
    except Exception:
        print(traceback.format_exc(limit=1)); ok = False
    atoms = n_atoms(pred_path(name))
    ok = atoms >= 50                      # истина = файл, а не return
    ok_cnt += bool(ok)
    w.writerow([name, int(ok), src, atoms]); new.flush()
    dt = time.time() - t0
    eta = dt / i * (len(bad) - i) / 60
    print(f'{i}/{len(bad)} {name} ok={ok} src={src} ATOM={atoms} '
          f'| прошло {dt/60:.0f} мин, осталось ~{eta:.0f} мин', flush=True)

print(f'ИТОГ: ok {ok_cnt}/{len(bad)}; журнал: {OUT_CSV}')
print('Контроль: find ~/3D_DNA/frag_build -name pred1.min.pdb -size -50c | wc -l -> ждём 0')
PYEOF
python3 -c "import ast,os; ast.parse(open(os.path.expanduser('~/Scripts/rebuild_539.py')).read()); print('синтаксис ок')"