#!/usr/bin/env python3
# renum_ss.py — перенумератор «по замыслу» для одиночных ДНК-цепей (ss-аптамеры).
# Отличие от normalize_pdb: НИКАКИХ измерений расстояний. Для ss-задачи мы
# заранее знаем правильный ответ: все остатки = одна цепь A, нумерация 1..L
# в порядке следования в файле. Поэтому растянутые стыки NSP НЕ разрезаются —
# tleap получает одну цепь и сам достраивает мостиковые фосфаты, а усиленная
# минимизация (2000/1000) их стягивает до ковалентной геометрии.
#
# Контроль качества встроен: порядок имён остатков сверяется с ожидаемым seq;
# при несовпадении (N сборщик перепутал куски) выходим с ненулевым кодом —
# задача уйдёт в ремонтный проход repair_ss.py.
#
# Использование:  python3 renum_ss.py  in.pdb  out.pdb  <SEQ>
import sys

NT = {'A': 'DA', 'C': 'DC', 'G': 'DG', 'T': 'DT'}


def main():
    src, dst, seq = sys.argv[1], sys.argv[2], sys.argv[3].upper()
    seen = []          # [(idx, resname)]
    out = []
    prev_raw = None
    for line in open(src):
        if not line.startswith('ATOM'):
            continue
        raw = line[17:27]
        if raw != prev_raw:
            seen.append(line[17:20].strip())
            prev_raw = raw
    exp = [NT[c] for c in seq]
    if len(seen) != len(exp):
        print(f'renum_ss: остатков {len(seen)} != L={len(exp)}', file=sys.stderr)
        sys.exit(1)
    # допускаем терминальные суффиксы имён вида DA5/DT3 от старых нормализаторов
    seen_base = [r.rstrip('53') for r in seen]
    if seen_base != exp:
        print(f'renum_ss: порядок остатков не совпал с seq', file=sys.stderr)
        sys.exit(1)

    prev_raw = None
    n = 0
    for line in open(src):
        if not line.startswith('ATOM'):
            continue
        raw = line[17:27]
        if raw != prev_raw:
            n += 1
            prev_raw = raw
        out.append(line[:17] + line[17:20].strip().ljust(3)[:3]
                   + ' A' + f'{n:4d}' + line[26:])
    out.append('TER\nEND\n')
    with open(dst, 'w') as f:
        f.writelines(out)
    print(f'renum_ss: цепь A, остатков {n} (по замыслу)')


if __name__ == '__main__':
    main()
