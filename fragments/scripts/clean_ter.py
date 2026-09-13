import os
import  math
def clean_ter(in_path, out_path,cutoff=3.0):

    if os.path.exists(in_path) and os.path.getsize(in_path) > 0:
        lines = open(in_path).readlines()

    num_of_ter_nec = 0
    num_of_ter = 0
    with open(out_path, 'w', newline='') as f:
        last_chain_coord = None
        for i in range(len(lines)):
            name = lines[i][12:16].strip()
            if lines[i].startswith('ATOM'):
                if name == "O3'":
                    x = float(lines[i][30:38])
                    y = float(lines[i][38:46])
                    z = float(lines[i][46:54])
                    last_chain_coord = (x,y,z)
            if lines[i].startswith('TER'):
                new_chain_coord = (None , None , None)
                name_after = None
                num_of_ter +=1
                j = 1
                while i+j < len(lines) and not lines[i+j].startswith('ATOM') :
                    j+=1
                if i + j < len(lines):
                    k = i + j
                    resi = lines[k][22:26]
                    while k < len(lines) and lines[k].startswith('ATOM' ) and lines[k][22:26] == resi:
                        name_after =lines[k][12:16].strip()
                        if name_after == 'P' or name_after == "O5'"   :
                            x = float(lines[k][30:38])
                            y = float(lines[k][38:46])
                            z = float(lines[k][46:54])
                            new_chain_coord = (x,y,z)
                            break
                        k += 1
                if  not (new_chain_coord == (None, None, None)) and not(last_chain_coord == None):
                    if math.dist(new_chain_coord , last_chain_coord) > cutoff:
                        num_of_ter_nec += 1
                        f.write(lines[i])
                else:
                    num_of_ter_nec += 1
                    f.write(lines[i])

            else:
                f.write(lines[i])
    print('Количество нужных TER:',num_of_ter_nec)
    print('КОличество TER всего:', num_of_ter)
    print('Удалено TER:',num_of_ter-num_of_ter_nec)
    return num_of_ter - num_of_ter_nec
n = clean_ter('fit_dirty.pdb', 'fit_clean.pdb')
print('удалено TER:', n)
