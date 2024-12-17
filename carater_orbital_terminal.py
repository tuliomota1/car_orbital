# -*- coding: utf-8 -*-
"""
@author: tulio
"""
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from tabulate import tabulate
from termcolor import colored

def parseArgs():
    parser = argparse.ArgumentParser(description='Script para facilitar visualização de caráter orbital e atomico das bandas projetadas a partir de outputs do VASP (PROCAR file).',
									epilog= "Codigo desenvolvido para facilitar a aplicação do método DFT-1/2 a partir do conhecimento do caráter das bandas.",
									prog="carater_orbital.py")
    
    parser.add_argument('-p', '--procar', type=str, default='PROCAR', help="Arquivo PROCAR padrao do calculo DFT --> pode ser alterado")
    
    parser.add_argument('-k', '--kpers', type=int, default=1, help="Ponto-k de escolha do usuario")
    
    parser.add_argument('-b', '--bands', type=int, default=2, help="Banda de escolha do usuario para a tabela de resultado")
    
    parser.add_argument('-c', '--graficc', type=int, default=0, help="Escolha a forma da matriz. Padrao é simples. Com -c 1 a matriz é colorida")
    
    return parser.parse_args()
    

arguments=parseArgs()

ispin=1
def get_outcar():
    arq= 'OUTCAR'
    outcars = [line for line in open(arq) if line.strip()]   
    for i, line in enumerate(outcars):
        if 'NKPTS =' in line:
            kpontos = int(line.split()[3])
            nbandas = int(line.split()[-1]) 
            break      
    return kpontos,nbandas

nkpts,nband=get_outcar()

POSFILE='POSCAR'
PROCFILE='PROCAR'

def get_poscar(poscars):
    poscar = [line for line in open(str(poscars)) if line.strip()]
    atom_name = poscar[5].split()
    atom_num_list = poscar[6].split()
    atom_num = []
    atom_total = 0
    for terms in atom_num_list:
        atom_num.append(int(terms))
        atom_total += int(terms)
    return atom_name,atom_num,atom_total   
atom_name,atom_num,atom_total=get_poscar(POSFILE)

def get_procar(procars):
    procar = [line for line in open(str(procars)) if line.strip()]

    orbit_list = procar[4].split()
    orbit_list.remove('ion')
    orbit_num = len(orbit_list)

    atom_name,atom_num,atom_total = get_poscar(POSFILE)

    atom_total += 1
    kpoint_line = []
    if ispin == 1:
        spin_list = ['s_tot','s_x','s_y','s_z']
        spin_num = 4
        band_project = np.zeros([nkpts,nband,spin_num,atom_total,orbit_num])
        for ii,line in enumerate(procar):
            if 'k-point' in line and 'k-points' not in line:
                kpoint_line.append(ii)
        for ii in range(len(kpoint_line)-1):
            k_index = ii
            band_perk = procar[kpoint_line[k_index]:kpoint_line[k_index+1]]
            for jj,band_line in enumerate(band_perk):
                if 'band' in band_line:
                    band_index = int(band_line.split()[1])-1
                    
                    for spin_index in range(spin_num):
                        for atom_index in range(atom_total):
                            for orbit_index in range(orbit_num):
                                to_add = band_perk[jj+2+(atom_total)*spin_index+atom_index].split()
                                to_add = to_add[1:]
                                to_add1 = []
                                for terms in to_add:
                                    to_add1.append(float(terms))
                                band_project[k_index,band_index,spin_index,atom_index,orbit_index]\
                                =to_add1[orbit_index]
        k_index = len(kpoint_line)-1
        band_perk = procar[kpoint_line[k_index]:]
        for jj,band_line in enumerate(band_perk):
            if 'band' in band_line:
                band_index = int(band_line.split()[1])-1
                for spin_index in range(spin_num):
                    for atom_index in range(atom_total):
                        for orbit_index in range(orbit_num):
                            to_add = band_perk[jj+2+(atom_total)*spin_index+atom_index].split()
                            to_add = to_add[1:]
                            to_add1 = []
                            for terms in to_add:
                                to_add1.append(float(terms))
                            band_project[k_index,band_index,spin_index,atom_index,orbit_index]\
                            =to_add1[orbit_index]
    
    return band_project,atom_name,atom_num,atom_total,orbit_list,orbit_num,spin_list,spin_num

band_project,atom_name,atom_num,atom_total,orbit_list,orbit_num,spin_list,spin_num=get_procar(PROCFILE)
############################################

## band_project[kponto][banda][carater-spin]
band=arguments.bands #entrada do usuario
kperso=arguments.kpers #entrada do usuario
a=band_project[kperso-1][band-1][0]
#a=[linha][coluna]
#print(a)
cototal=a[sum(atom_num)][9]
carat_matx=[]
contador=0
for i in range(len(atom_num)):
    s_previ=0
    py_previ=0
    pz_previ=0
    px_previ=0
    dxy_previ=0
    dyz_previ=0
    dz2_previ=0
    dxz_previ=0
    dx2_previ=0
    for k in range(atom_num[i]):        
        s_previ+=a[k+contador][0]
        py_previ+=a[k+contador][1]
        pz_previ+=a[k+contador][2]
        px_previ+=a[k+contador][3]
        dxy_previ+=a[k+contador][4]
        dyz_previ+=a[k+contador][5]
        dz2_previ+=a[k+contador][6]
        dxz_previ+=a[k+contador][7]
        dx2_previ+=a[k+contador][8]
        
    carat_matx+=[[str(atom_name[i]),((s_previ)/cototal)*100,((py_previ+pz_previ+px_previ)/cototal)*100,
                  ((dxy_previ+dyz_previ+dz2_previ+dxz_previ+dx2_previ)/cototal)*100]]
    contador+=atom_num[i]
    
###################################

def cabecalho():
    header=['Atom','s', 'p','d']
    return [colored(c, 'cyan', attrs=['bold']) for c in header]

def tabela():
    tabela=[]
    for i in range(len(atom_name)):
        e=str(atom_name[i])
        tabela+=[
        [e,carat_matx[i][1],carat_matx[i][2],carat_matx[i][3]],        
        ]
    return [
        [colored(d[0], 'yellow', attrs=['bold']), d[1],d[2],d[3]] for d in tabela
    ]
print('----------------------------------------------------')
print('Elementos no sistema:',atom_name)
print('Quantidade de atomos por elemento:',atom_num)
print('BANDA escolhida pelo usuario:', band)
print('PONTO-K:', kperso)
print('----------------------------------------------------')
print('----------------------------------------------------')
print('  _   _   _   _   _   _   _   _   _ ') 
print(' / \ / \ / \ / \ / \ / \ / \ / \ / \\')
print('( C | . | O | r | b | i | t | a | l )')
print(' \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/ ')
print('----------------------------------------------------')
print('----------------------------------------------------')
if arguments.graficc == 1:
    print(tabulate(tabela(),headers=cabecalho(), tablefmt='fancy_grid'))
cabecalhos=['Atom', 's','p','d']
tabelas=[]
for i in range(len(atom_name)):
    e=str(atom_name[i])
    tabelas+=[
    [e,carat_matx[i][1],carat_matx[i][2],carat_matx[i][3]],   
    ]
if arguments.graficc == 0:
    print(tabulate(tabelas,headers=cabecalhos, tablefmt='fancy_grid'))