
import math
import sys


def _xyz(line):
    """Координаты атома из PDB-строки."""
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def _atom_name(line):
    """Имя атома без пробелов: P, O3', C5' ..."""
    return line[12:16].strip()


def _raw_resi(line):
    """
    Сырой идентификатор остатка: буква цепи + номер, как записано в файле.
    Нужен только чтобы заметить, что остаток сменился.
    """
    return line[21:26]


def normalize_pdb(in_path, out_path, cutoff=3.0):
    """
    Перенумеровывает остатки сквозняком внутри цепи, проставляет буквы
    цепей A, B, C... и расставляет TER на границах.

    Возвращает список длин цепей, например [6, 6].
    """
    with open(in_path) as fh:
        lines = fh.readlines()

    chain_idx = -1        # -1 значит "ещё ни одной цепи не начали"
    res_num = 0           # номер остатка внутри текущей цепи
    prev_raw = None       # сырой id предыдущего остатка
    prev_o3 = None        # координаты O3' предыдущего остатка
    chain_lens = []       # сколько остатков в каждой цепи

    out = []

    for i, line in enumerate(lines):
        if not line.startswith('ATOM'):
            continue      # TER, END, REMARK из входа выбрасываем

        raw = _raw_resi(line)

        # --- сменился ли остаток ---
        if raw != prev_raw:
            # ищем P или O5' среди строк ЭТОГО остатка
            head = None
            k = i
            while k < len(lines) and lines[k].startswith('ATOM') \
                    and _raw_resi(lines[k]) == raw:
                nm = _atom_name(lines[k])
                if nm == 'P':
                    head = _xyz(lines[k])
                    break                    # P приоритетнее
                if nm == "O5'" and head is None:
                    head = _xyz(lines[k])    # запасной вариант
                k += 1

            # --- новая цепь или продолжение? ---
            if prev_o3 is None or head is None:
                new_chain = True             # самый первый остаток
            else:
                new_chain = math.dist(head, prev_o3) > cutoff

            if new_chain:
                if chain_idx >= 0:
                    chain_lens.append(res_num)
                    out.append("TER\n")
                chain_idx += 1
                res_num = 1
            else:
                res_num += 1

            prev_raw = raw

            # запоминаем O3' этого остатка для следующего сравнения
            prev_o3 = None
            k = i
            while k < len(lines) and lines[k].startswith('ATOM') \
                    and _raw_resi(lines[k]) == raw:
                if _atom_name(lines[k]) == "O3'":
                    prev_o3 = _xyz(lines[k])
                    break
                k += 1

        # --- переписываем строку с новой разметкой ---
        chain_letter = chr(ord('A') + chain_idx)
        out.append(line[:21] + chain_letter + f"{res_num:>4d}" + line[26:])

    if chain_idx >= 0:
        chain_lens.append(res_num)
        out.append("TER\n")
    out.append("END\n")

    with open(out_path, 'w') as fh:
        fh.writelines(out)

    return chain_lens


if __name__ == '__main__':
    src = sys.argv[1] if len(sys.argv) > 1 else 'fit_dirty.pdb'
    dst = sys.argv[2] if len(sys.argv) > 2 else 'fit_norm.pdb'
    lens = normalize_pdb(src, dst)
    print(f"{src} -> {dst}")
    print(f"цепей: {len(lens)}, длины: {lens}")
