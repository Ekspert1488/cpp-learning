#!/usr/bin/env python3
# diag_shape.py — почему аптамеры "палки": смотрим 2D из задач и реальную форму 3D.
# Запуск: python3 diag_shape.py [сколько_моделей_мерить]   (по умолчанию 200)
# Ничего не меняет, только читает.
import csv, math, os, random, sys

BUILD = os.path.expanduser('~/3D_DNA/apt_build')
TASKS = os.path.join(BUILD, 'apt_tasks.csv')
RES   = os.path.join(BUILD, 'apt3d_results.csv')
N     = int(sys.argv[1]) if len(sys.argv) > 1 else 200

# ---------- часть 1: что вообще в наших 2D-структурах ----------
tasks = {}
for r in csv.DictReader(open(TASKS)):
    tasks[r['id']] = (r['seq'], r['ss'], int(r['L']))

n_flat = 0
pair_frac = []
for _id, (seq, ss, L) in tasks.items():
    npair = ss.count('(')
    pair_frac.append(2.0 * npair / L)
    if npair == 0:
        n_flat += 1
pair_frac.sort()
print('=== 2D (из apt_tasks.csv) ===')
print(f'задач всего            : {len(tasks)}')
print(f'БЕЗ единой пары "(" )" : {n_flat}   <-- такие NSP строит линейной ниткой')
if pair_frac:
    print(f'доля спаренных нукл.   : медиана {pair_frac[len(pair_frac)//2]:.2f} | '
          f'min {pair_frac[0]:.2f} | max {pair_frac[-1]:.2f}')
print()

# ---------- часть 2: реальная геометрия min.pdb ----------
def load_p_atoms(path):
    """координаты P (или C1') по остаткам — скелет цепи"""
    pts, last = [], None
    for l in open(path):
        if not l.startswith('ATOM'):
            continue
        name = l[12:16].strip()
        if name not in ("P", "C1'"):
            continue
        rid = l[22:27]
        if name == 'P' or (last != rid):
            pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54]), rid, name))
        last = rid
    # приоритет P; если P мало (первый остаток без P), добираем C1'
    p = [(x, y, z) for x, y, z, rid, nm in pts if nm == 'P']
    c = [(x, y, z) for x, y, z, rid, nm in pts if nm == "C1'"]
    return p if len(p) >= len(c) - 1 else c


def rg(pts):
    n = len(pts)
    cx = sum(p[0] for p in pts) / n
    cy = sum(p[1] for p in pts) / n
    cz = sum(p[2] for p in pts) / n
    return math.sqrt(sum((p[0]-cx)**2 + (p[1]-cy)**2 + (p[2]-cz)**2 for p in pts) / n)


ok = [r['name'] for r in csv.DictReader(open(RES)) if str(r.get('ok')) == 'True']
random.seed(0)
sample = random.sample(ok, min(N, len(ok)))
rows = []
for name in sample:
    pdb = os.path.join(BUILD, name, f'{name}.pred1.min.pdb')
    if not os.path.exists(pdb):
        continue
    pts = load_p_atoms(pdb)
    if len(pts) < 5:
        continue
    L = len(pts)
    e2e = math.dist(pts[0], pts[-1])
    r = rg(pts)
    contour = 4.4 * (L - 1)                      # теоретический максимум растяжения
    rows.append((name, L, e2e, r, e2e / contour if contour else 0))

if not rows:
    print('нет моделей для замера'); sys.exit(0)

rows.sort(key=lambda v: v[4])
ext = [v[4] for v in rows]
e2es = sorted(v[2] for v in rows)
rgs = sorted(v[3] for v in rows)
print(f'=== 3D (замерено моделей: {len(rows)}) ===')
print(f'расстояние 5\'-3\'   : медиана {e2es[len(e2es)//2]:6.1f} A | '
      f'min {e2es[0]:6.1f} | max {e2es[-1]:6.1f}')
print(f'радиус инерции Rg  : медиана {rgs[len(rgs)//2]:6.1f} A | '
      f'min {rgs[0]:6.1f} | max {rgs[-1]:6.1f}')
print(f'растянутость e2e/L : медиана {sorted(ext)[len(ext)//2]:.2f}  '
      f'(1.00 = прямая палка, <0.3 = свёрнуто)')
print()
n_stick = sum(1 for e in ext if e > 0.7)
print(f'ПАЛОК (e2e/contour > 0.7): {n_stick} из {len(rows)} '
      f'({100.0*n_stick/len(rows):.0f}%)')
print()
print('--- 5 самых свёрнутых ---')
for name, L, e2e, r, x in rows[:5]:
    print(f'  {name} L={L:3d} e2e={e2e:6.1f} Rg={r:5.1f} ext={x:.2f}')
print('--- 5 самых растянутых ---')
for name, L, e2e, r, x in rows[-5:]:
    print(f'  {name} L={L:3d} e2e={e2e:6.1f} Rg={r:5.1f} ext={x:.2f}')
print()
print('ОРИЕНТИРЫ: свёрнутый аптамер 50 нт -> e2e 40-70 A, Rg 18-25 A, ext 0.2-0.35')
print('           прямая нитка 50 нт      -> e2e ~215 A, Rg ~63 A, ext ~1.0')
