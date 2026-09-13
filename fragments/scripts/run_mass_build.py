import csv
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

from Enumerate_nucl import make_task_list, build_one, build_one_fallback

OUT_ROOT = '/home/ekspert/3D_DNA/frag_build'
RESULTS  = os.path.join(OUT_ROOT, 'build_results.csv')
MAX_WORKERS = 8


def route(task):
    name, seq, sst, kind = task
    if kind == 'ss' and len(seq) <= 4:
        return build_one_fallback      # ss L3/L4 — NSP не умеет
    return build_one                    # ds и ss L5/L6 — NSP


def worker(task):
    try:
        return route(task)(task, OUT_ROOT)
    except Exception as e:
        return {'name': task[0], 'ok': False, 'reason': f'exception: {e!r}'}


def main():
    tasks = make_task_list()

    # докачка: успешные пропускаем (прогон можно перезапускать сколько угодно)
    done = set()
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            for row in csv.DictReader(f):
                if row['ok'] == 'True':
                    done.add(row['name'])
    todo = [t for t in tasks if t[0] not in done]
    print(f'всего {len(tasks)} | готово {len(done)} | осталось {len(todo)}')
    if not todo:
        return

    os.makedirs(OUT_ROOT, exist_ok=True)
    t0 = time.time()
    n_ok = n_bad = 0
    new_file = not os.path.exists(RESULTS) or os.path.getsize(RESULTS) == 0
    with open(RESULTS, 'a', newline='') as out:
        w = csv.writer(out)
        if new_file:
            w.writerow(['name', 'ok', 'kind', 'L', 'charge', 'expected',
                        'ret', 'reason', 'src', 'seconds'])
            out.flush()
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(worker, t): t for t in todo}
            for i, fut in enumerate(as_completed(futs), 1):
                r = fut.result()
                sec = time.time() - t0
                if r.get('ok'):
                    n_ok += 1
                else:
                    n_bad += 1
                w.writerow([r.get('name'), r.get('ok'), r.get('kind', ''),
                            r.get('L', ''), r.get('charge', ''),
                            r.get('expected', ''), r.get('ret', ''),
                            r.get('reason', ''), r.get('src', 'nsp'),
                            round(sec, 1)])
                out.flush()
                if i % 100 == 0 or i == len(todo):
                    eta = sec / i * (len(todo) - i)
                    print(f'{i}/{len(todo)} ok={n_ok} bad={n_bad} '
                          f'прошло {sec/60:.0f} мин, ETA {eta/60:.0f} мин',
                          flush=True)
    print('ГОТОВО:', n_ok, 'ok,', n_bad, 'проблемных')


if __name__ == '__main__':
    main()