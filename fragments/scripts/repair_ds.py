import csv
import os
import re
import subprocess

from Enumerate_nucl import make_task_list, NT

FRAG    = '/home/ekspert/3D_DNA/frag_build'
RESULTS = os.path.join(FRAG, 'build_results.csv')

LEAP = '''source leaprc.DNA.OL15
model = loadpdb repaired.pdb
check model
saveamberparm model prmtop inpcrd
quit
'''

MIN = '''repair ds: minimisation (усиленная — стягиваем стык)
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


def sh(cmd, cwd):
    subprocess.run(['bash', '-c', cmd], cwd=cwd, timeout=900)


def dist(a, b):
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


def residues_of(pdb):
    """Остатки в порядке следования в файле: [(ключ, resname, [строки])]."""
    res = []
    for l in open(pdb):
        if not l.startswith('ATOM'):
            continue
        key = l[17:27]                      # resname + chain + resnum (сырые)
        if not res or res[-1][0] != key:
            res.append([key, l[17:20].strip(), []])
        res[-1][2].append(l)
    return res


def junctions_ok(pdb, cutoff=3.05):
    """Все стыки внутри каждой цепи не длиннее cutoff."""
    chains = {}
    for l in open(pdb):
        if not l.startswith('ATOM'):
            continue
        num = int(l[22:26])
        chains.setdefault((l[21], num), {})[l[12:16].strip()] = \
            (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    by_chain = {}
    for (ch, num), atoms in chains.items():
        by_chain.setdefault(ch, []).append((num, atoms))
    for ch, lst in by_chain.items():
        lst.sort()
        for (n1, a1), (n2, a2) in zip(lst, lst[1:]):
            if n2 == n1 + 1 and "O3'" in a1 and "O5'" in a2:
                d = dist(a1["O3'"], a2["O5'"])
                if d > cutoff:
                    return False, f'цепь {ch} {n1}->{n2}: {d:.2f} A'
    return True, ''


def repair_one(task):
    name, seq, sst, kind = task
    s1, s2 = seq.split('&')
    L = len(s1)
    workdir = os.path.join(FRAG, name)
    fit = os.path.join(workdir, f'{name}.pred1', 'fit.pdb')
    if not os.path.exists(fit):
        return {'name': name, 'ok': False, 'reason': 'нет fit.pdb'}

    res = residues_of(fit)
    if len(res) != 2 * L:
        return {'name': name, 'ok': False,
                'reason': f'остатков {len(res)} != {2 * L}'}

    got  = [r[1] for r in res]
    exp1 = [NT[c] for c in s1]
    exp2 = [NT[c] for c in s2]
    if got[:L] == exp1 and got[L:] == exp2:
        first, second = res[:L], res[L:]
    elif got[:L] == exp2 and got[L:] == exp1:
        first, second = res[L:], res[:L]      # цепи в файле идут в обратном порядке
    else:
        return {'name': name, 'ok': False, 'reason': 'порядок остатков не распознан'}

    # переписываем цепи «по замыслу»: A = первая цепь, B = вторая
    with open(os.path.join(workdir, 'repaired.pdb'), 'w') as f:
        for letter, block in (('A', first), ('B', second)):
            for i, (_, _, lines) in enumerate(block, 1):
                for l in lines:
                    f.write(l[:21] + letter + f'{i:4d}' + l[26:])
            f.write('TER\n')
        f.write('END\n')

    with open(os.path.join(workdir, 'leap.in'), 'w') as f:
        f.write(LEAP)
    with open(os.path.join(workdir, 'min.in'), 'w') as f:
        f.write(MIN)

    sh('source ~/nsp_env.sh && tleap -f leap.in > tleap_repair.log 2>&1', workdir)
    if not (os.path.exists(os.path.join(workdir, 'prmtop'))
            and os.path.exists(os.path.join(workdir, 'inpcrd'))):
        return {'name': name, 'ok': False, 'reason': 'tleap: нет prmtop/inpcrd'}

    sh('source ~/nsp_env.sh && msander -O -i min.in -o min.out '
       '-p prmtop -c inpcrd -r min.rst >> tleap_repair.log 2>&1', workdir)
    if not os.path.exists(os.path.join(workdir, 'min.rst')):
        return {'name': name, 'ok': False, 'reason': 'msander: нет min.rst'}

    sh('source ~/nsp_env.sh && ambpdb -p prmtop < min.rst > min_raw.pdb '
       '2>> tleap_repair.log', workdir)
    min_pdb = os.path.join(workdir, f'{name}.pred1.min.pdb')
    sh(f'python3 ~/normalize_pdb.py min_raw.pdb {name}.pred1.min.pdb '
       '>> tleap_repair.log 2>&1', workdir)
    if not os.path.exists(min_pdb):
        return {'name': name, 'ok': False, 'reason': 'normalize: нет итога'}

    charge = None
    with open(os.path.join(workdir, 'leap.log')) as f:
        for line in f:
            m = re.search(r'unperturbed charge of the unit \((-?\d+\.?\d*)\)', line)
            if m:
                charge = float(m.group(1))

    expected = -2 * (L - 1)
    j_ok, j_why = junctions_ok(min_pdb)
    ok = charge is not None and abs(charge - expected) < 1e-6 and j_ok
    why = '' if ok else (f'стык: {j_why}' if not j_ok else 'заряд не совпал')
    return {'name': name, 'ok': ok, 'kind': 'ds', 'L': L,
            'charge': charge, 'expected': expected,
            'reason': why, 'min_pdb': min_pdb}


def main():
    rows = list(csv.DictReader(open(RESULTS)))
    ok = {r['name'] for r in rows if r['ok'] == 'True'}
    bad_ds = [t for t in make_task_list() if t[3] == 'ds' and t[0] not in ok]
    print('дуплексов на ремонт:', len(bad_ds))
    out = open(RESULTS, 'a', newline='')
    w = csv.writer(out)
    n_ok = 0
    for t in bad_ds:
        r = repair_one(t)
        if r['ok']:
            n_ok += 1
        print(r['name'], '| ok:', r['ok'], '| заряд:', r.get('charge'),
              '|', r.get('reason', ''), flush=True)
        w.writerow([r.get('name'), r.get('ok'), r.get('kind', ''), r.get('L', ''),
                    r.get('charge', ''), r.get('expected', ''), '',
                    r.get('reason', ''), 'repair-ds', 0])
        out.flush()
    out.close()
    print('РЕМОНТ ГОТОВ:', n_ok, 'из', len(bad_ds))


if __name__ == '__main__':
    main()