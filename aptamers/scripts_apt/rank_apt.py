#!/usr/bin/env python3
# Этап 15b: сортировка по скору + отсев поз ВНЕ эпитопа.
# Вход:  ~/docking_apt/scores_apt.csv (из dock_apt.py rescore)
# Выход: ~/docking_apt/ranked_epi.csv   — прошли эпитопный фильтр, по возрастанию score_epi
#        ~/docking_apt/rejected_epi.csv — отсеяны (нет позы у эпитопа)
#        ~/docking_apt/top500.csv       — буфер для этапа 16
# Запуск: python3 rank_apt.py [N_TOP]     (по умолчанию 500)
import csv, os, sys

S    = os.path.expanduser('~/docking_apt/scores_apt.csv')
OUT  = os.path.dirname(S)
NTOP = int(sys.argv[1]) if len(sys.argv) > 1 else 500
MIN_NEAR = int(os.environ.get('MIN_NEAR_EPI', '1'))   # сколько поз из 100 должны быть у эпитопа

rows = {}
for r in csv.DictReader(open(S)):
    rows[r['name']] = r          # дедуп на случай повторного rescore

passed, rejected = [], []
for name, r in rows.items():
    se = r.get('score_epi') or ''
    nn = int(r.get('n_near_epi') or 0)
    if se and nn >= MIN_NEAR:
        passed.append((float(se), name, r))
    else:
        sx = r.get('score_extra') or ''
        rejected.append((float(sx) if sx else 0.0, name, r))

passed.sort(key=lambda v: v[0])          # HDOCK: чем отрицательнее, тем лучше
rejected.sort(key=lambda v: v[0])

COLS = ['rank', 'name', 'kind', 'L', 'src',
        'score_extra', 'dist_extra', 'n_near_extra',
        'score_epi', 'dist_epi', 'n_near_epi']


def dump(path, items, with_rank=True):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(COLS)
        for i, (_, name, r) in enumerate(items, 1):
            w.writerow([i if with_rank else '', name, r.get('kind', 'ss'),
                        r.get('L', ''), r.get('src', ''),
                        r.get('score_extra', ''), r.get('dist_extra', ''),
                        r.get('n_near_extra', ''),
                        r.get('score_epi', ''), r.get('dist_epi', ''),
                        r.get('n_near_epi', '')])


dump(os.path.join(OUT, 'ranked_epi.csv'), passed)
dump(os.path.join(OUT, 'rejected_epi.csv'), rejected, with_rank=False)
dump(os.path.join(OUT, f'top{NTOP}.csv'), passed[:NTOP])

print(f'всего в scores_apt.csv : {len(rows)}')
print(f'прошли эпитопный фильтр: {len(passed)}  (n_near_epi >= {MIN_NEAR})')
print(f'отсеяны (вне эпитопа)  : {len(rejected)}')
print()
if passed:
    sc = [p[0] for p in passed]
    print(f'score_epi: лучший {sc[0]:.2f} | медиана {sc[len(sc)//2]:.2f} | худший {sc[-1]:.2f}')
    if len(passed) > NTOP:
        print(f'порог отсечки top{NTOP}: score_epi <= {passed[NTOP-1][0]:.2f}')
    print()
    print('=== ТОП-20 ===')
    print(f'{"#":>3} {"name":<11} {"L":>3} {"score_epi":>10} {"d_epi":>6} {"n_epi":>6}'
          f' {"score_extra":>12} {"d_ex":>6} {"n_ex":>5}')
    for i, (sc_i, name, r) in enumerate(passed[:20], 1):
        d = float(r['dist_epi']) if r.get('dist_epi') else 0.0
        sx = float(r['score_extra']) if r.get('score_extra') else 0.0
        dx = float(r['dist_extra']) if r.get('dist_extra') else 0.0
        print(f'{i:>3} {name:<11} {r.get("L",""):>3} {sc_i:>10.2f} {d:>6.2f} '
              f'{r.get("n_near_epi",""):>6} {sx:>12.2f} {dx:>6.2f} '
              f'{r.get("n_near_extra",""):>5}')
print()
print('файлы:', os.path.join(OUT, 'ranked_epi.csv'), '|',
      os.path.join(OUT, f'top{NTOP}.csv'), '|', os.path.join(OUT, 'rejected_epi.csv'))
