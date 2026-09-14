## Этапы
1. Нахождение структуры рецептора с лигандом в внеклеточной части. 
2. Подбор эпитопа связывания A2A.
3. Отбор 40 лучших фрагментов размера 3-6 нт.
4. Обучение AiDTA ; Получение 2804 аптамера длины 50-61 нт
5. Создание PDB файлов аптамеров.
6. Затем «черновой» докинг HDOCK. - **МЫ ЗДЕСЬ**
7. Отбор лучших 200 (500) и запуск «чистового» докинга HDOCK.
8. Отбор лучших 20 кандидатов и проверка через 3d_DNA.
9. Проверка с помощью HADDOCK + опционально ALphaFold3.
10. Разработка методики проверки кандидатов на связывание
11. Проверка кандидатов на афинность в лабораторных условиях.




## Используемые версии: 
    Python 3.9.25
    torch 2.8.0+cu128
    numpy 1.26.4
    pandas 2.3.1            
    scipy 1.13.1 
    RNAstructure Fold: Version 6.6 (April 2, 2026). Copyright Mathews Lab, University of Rochester.     
    Amber classic sander: Version 22.0: поля  - **leaprc.DNA.OL15**(использовали данное поле ) leaprc.DNA.OL21 leaprc.DNA.bsc1.
## HDOCK - Protein-protein and Protein-RNA/DNA Docking Huang Lab @ HUST, http://huanglab.phys.hust.edu.cn/ HDOCKlite v1.2 -- The ab initio docking engine of the HDOCK server (http://hdock.phys.hust.edu.cn/)
    Required Linux system: CentOS 6.0 or later
    This program performs the ab initio protein-protein and protein-RNA/DNA docking using
    the intrinsic scoring functions for protein-protein and protein-RNA interactions based
    on a hierarchical algorithm.
## Окружение 
    Linux Ekspert 6.18.33.2-microsoft-standard-WSL2 #1 SMP PREEMPT_DYNAMIC Thu Jun 18 21:54:43 UTC 2026 x86_64 x86_64 x86_64 GNU/Linux
    No LSB modules are available.
    Distributor ID: Ubuntu
    Description:    Ubuntu 24.04.3 LTS
    Release:        24.04
    Codename:       noble
    PRETTY_NAME="Ubuntu 24.04.3 LTS"
    NAME="Ubuntu"
    VERSION_ID="24.04"
## Мое вычислительные мощности:
        14
               total        used        free      shared  buff/cache   available
    Mem:           9.7Gi       639Mi       8.6Gi       3.5Mi       676Mi       9.1Gi
    CPU(s):                                  14
    Model name:                              AMD Ryzen 7 5800U with Radeon Graphics



