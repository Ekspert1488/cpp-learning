#!/usr/bin/env python3
# Задание 23: топ-5 на (тип x длина) по score_epi -> 44 фрагмента AiDTA
import os, sys, csv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Enumerate_nucl import make_task_list

S = os.path.expanduser('~/docking_H/scores_H.csv')
rows = {}
for r in csv.DictReader(open(S)):
    rows[r['name']] = r                     # дедуп: последняя строка побеждает
tasks = {t[0]: t for t in make_task_list()}  # name -> (name, seq, ss, kind)

def L_of(seq):
    return len(seq.split('&')[0])

# группы (kind, L) -> кандидаты с ненулевым score_epi
groups = {}
for name, r in rows.items():
    if name not in tasks:
        print('ВНИМАНИЕ: имени нет в make_task_list:', name); continue
    seq, kind = tasks[name][1], tasks[name][3]
    if not r['score_epi']:
        continue                             # не дотянулся до эпитопа — мимо
    groups.setdefault((kind, L_of(seq)), []).append(
        (float(r['score_epi']), name, seq, r))

print('=== группы: кандидантов с score_epi -> взято в топ')
top = []
for key in sorted(groups):
    cand = sorted(groups[key])[:5]           # по возрастанию скора (левее = лучше)
    print(f'{key[0]} L{key[1]}: {len(groups[key]):4d} -> {len(cand)}')
    for rank, (sc, name, seq, r) in enumerate(cand, 1):
        top.append((key, rank, name, seq, sc,
                    r.get('dist_epi', ''), r.get('n_near_epi', '')))

mono = [('A', '.'), ('C', '.'), ('G', '.'), ('T', '.')]
total = len(top) + len(mono)
print(f'итого: {len(top)} из докинга + {len(mono)} моно = {total} (ждём 44)')
if total != 44:
    print('!!! НЕ 44 — какая-то группа дала меньше 5, смотри таблицу выше')

# вторичная структура в формате авторов: ds -> (((&))), ss -> ..., моно -> .
list1, list2 = [], []
for (kind, L), rank, name, seq, sc, d, n in top:
    list1.append(seq)
    list2.append('('*L + '&' + ')'*L if kind == 'ds' else '.'*L)
for s, s2 in mono:
    list1.append(s); list2.append(s2)

with open(os.path.expanduser('~/top44.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['kind', 'L', 'rank', 'name', 'seq', 'ss2',
                'score_epi', 'dist_epi', 'n_near_epi'])
    for (kind, L), rank, name, seq, sc, d, n in top:
        s2 = '('*L + '&' + ')'*L if kind == 'ds' else '.'*L
        w.writerow([kind, L, rank, name, seq, s2, sc, d, n])
    for s, s2 in mono:
        w.writerow(['mono', 1, '-', s, s2, '', '', ''])

with open(os.path.expanduser('~/list1_list2.py'), 'w') as f:
    f.write('list1 = ' + repr(list1) + '\n')
    f.write('list2 = ' + repr(list2) + '\n')

print()
print('=== top-44 (лучшая пятёрка каждой группы, скор по возрастанию)')
for (kind, L), rank, name, seq, sc, d, n in top:
    print(f'{kind} L{L} #{rank}  {seq:15s} score_epi={sc:8.2f} '
          f'dist={d} n_near={n}  [{name}]')
print()
print('файлы: ~/top44.csv, ~/list1_list2.py')
