#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#xyzoverlay

import sys                                                      #sys.exit
import os.path                                                  #filename split 
import argparse                                                 #argument parser
import pandas as pd                                             #pandas tables
import numpy as np                                              #calculations
from scipy.spatial.distance import pdist, squareform, cosine    #for the calculations of the distance matrix and angles (cosine)
import matplotlib.pyplot as plt                                 #for molecule display
from mpl_toolkits.mplot3d import Axes3D                         #for molecule display
from cycler import cycler                                       #generate color cycle
from itertools import cycle                                     #for color cycling
import io                                                       #IO for (easy) saving multi xyz
import re                                                       #regex (get atom x y z from multi xyz)

#var for check if xyz is multi xyz
is_trj = 0

###### start view settings
#dpi for figure
plt.rcParams['savefig.dpi'] = 150
#alpha atoms and bonds
alpha_atoms = 0.55
alpha_bonds = 0.55
#sphere radius for atom view, change exponent
atom_scaler = 1e60
#cylinder radius for bond view, change exponent
bond_scaler = 1e6
#default color cycle
color_cycle = cycler(c=['b', 'r', 'g', 'c', 'm', 'y'])
###### end view settings 

#covalent radii from Alvarez (2008)
#DOI: 10.1039/b801115j
covalent_radii = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28,
    'Be': 0.96, 'B': 0.84, 'C': 0.76, 
    'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
    'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
    'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
    'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
    'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
    'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
    'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
    'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
    'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
    'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
    'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
    'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
    'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
    'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
    'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
    'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
    'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
    'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
    'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
    'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
    'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
    'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
    'Am': 1.80, 'Cm': 1.69
}

#https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
    '''
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    
    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    # set to 0.38 --> bigger molecules 
    plot_radius = 0.35*max([x_range, y_range, z_range])
    
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


#Kabsch algorithm
#https://github.com/charnley/rmsd
def kabsch(P, Q):
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # create Rotation matrix U
    U = np.dot(V, W)
    return U

#get rotmatrix from Kabsch and transform coordinates
def align_xyz(vec1,vec2,coord):
    rotmatrix = kabsch(vec1,vec2)
    return np.dot(coord,rotmatrix)

#argument parser
parser = argparse.ArgumentParser(
        prog='xyzoverlay', 
        description = "Overlay and align two or more molecules from (multi) xyz files.")

#filename is required
parser.add_argument("filename", 
    nargs = '+',
    help = "filename(s), xyz; e.g. mol1.xyz or mol1.xyz mol2.xyz")

#atoms for superposition
parser.add_argument('-a', '--atoms',
    nargs= '+',
    type= int,
    action= 'append',
    help =  "select atoms for superposition, e.g. -a 1 2 -a 3 4")

#select the same atoms in all xyz for superposition
parser.add_argument('-sa', '--sameatoms',
    nargs= '+',
    type= int,
    help =  "select the same atoms in all xyz for superposition, e.g. -sa 1 2 3 ")

#select all atoms for superposition
parser.add_argument('-aa','--allatoms',
    default=0, 
    action='store_true',
    help='select all atoms for superposition, xyz files must have the same number of atoms')

#plot or view, color by atom
parser.add_argument('-vca','--viewcba',
    default=0, action='store_true',
    help='view xyz coordinates (atoms) and bonds, color by atom')

#plot or view, color by molecule
parser.add_argument('-vcm','--viewcbm',
    default=0, action='store_true',
    help='view xyz coordinates (atoms) and bonds, color by molecule')

#use Matplotlib colormap
parser.add_argument('-cm','--cmap',
    type= str,
    help='use Matplotlib colormap in view mode')

#add +x% to radius (for view mode)
parser.add_argument('-r','--radius',
    default=8,
    type=float,
    help='enlarge atomic radii by x %%, e.g. -r 15.2, default is 8 %%')

#exclude elements in view mode
parser.add_argument('-ee','--excludeEl',
    nargs="+",
    type=str,
    help='exclude specified elements in view mode; e.g. -ee C or -ee C N')

#exclude atoms in view mode
parser.add_argument('-ea','--excludeAt',
    nargs="+",
    type=int,
    help='exclude specified atoms in view mode; e.g. -ea 1 or -ea 1 2')

#save single xyz
parser.add_argument('-s','--save',
    default=0, action='store_true',
    help='save superimposed / aligned data as single xyz file(s)')

#save multi xyz
parser.add_argument('-st','--savetr',
    default=0, action='store_true',
    help='save superimposed / aligned data as (multi) xyz or xyz trajectory file')

#parse arguments
args = parser.parse_args()

#open xyz
try:
    #header lines list
    head_list=[]
    #get the two header lines from xyz file
    for file in args.filename:
        head = open(file).readlines()[0:2]
        #if the second line is just a line break
        if head[-1] == "\n":
            head[-1] = " \n"
        head= ''.join(head).rstrip('\n')
        head_list.append(head)
    #read element x y z from xyz file
    #create a list of data frames, molecule = entry in list
    #molecular coordinates are in the data frame
    xyz_df_list = list(pd.read_csv(file, 
            sep='\s+', 
            skiprows=2, 
            names=["element", "x", "y", "z"],on_bad_lines='skip') 
            for file in args.filename)
    #index starts at 1, first atom is atom 1
    for xyz_df in xyz_df_list:
        xyz_df.index +=1
    #file not found
except IOError:
    print(f"File(s) not found. Exit.")
    sys.exit(1)

#check if xyz is trj or multi xyz. trj xyz contains nan after pandas csv import
for xyz_df in xyz_df_list:
    index_list = xyz_df.loc[pd.isna(xyz_df[['x','y','z']]).any(axis=1), :].index
    
    if index_list.size > 0:
        #check if more than one multi xyz is opened or multi xyz and normal xyz together
        #exit if True
        if len(args.filename) > 1:
            print('Warning! Only one multi xyz file is supported. Exit.')
            sys.exit(1)
        #some list for processing multi xyz files
        xyz_trj_list = []
        xyz_df_list = []
        head_list = []
        datas = []
        #is multi xyz = True
        is_trj = 1
        #open the multi xyz file again
        with open(args.filename[0], 'r') as xyz_trj_file:
            for line in xyz_trj_file:
                #regex for 'element x y z' --> grouping
                data = re.search(r'([a-zA-Z]{1,2})\s{1,}'
                       r'([+-]?\d*\.?\d*[e]?[+-]?\d*)\s{1,}'
                       r'([+-]?\d*\.?\d*[e]?[+-]?\d*)\s{1,}'
                       r'([+-]?\d*\.?\d*[e]?[+-]?\d*)',line)
                if data:
                    #add if 'element x y z'
                    datas.append(data.groups())
                else:
                    #if this is not 'element x y z' than this must the two xyz header lines
                    #add header to list
                    head = line
                    head= ''.join(head).rstrip('\n')
                    head_list.append(head)
                    if datas:
                        #'element x y z' --> data frame
                        xyz_df = pd.DataFrame(datas,columns=['element','x','y','z'])
                        #index starts at 1
                        xyz_df.index +=1
                        #set data types for 'element x y z'
                        xyz_df = xyz_df.astype({'element': str, 'x': float, 'y': float, 'z': float})
                        #collect data frames (molecular coordinates) 
                        #for each molecule in list of data frames
                        xyz_df_list.append(xyz_df)
                        #empty regex groups
                        datas = []
            #xyz_df = pd.DataFrame(datas,columns=['element','x','y','z'])
            #xyz_df.index +=1
            #xyz_df = xyz_df.astype({'element': str, 'x': float, 'y': float, 'z': float})
            #xyz_df_list.append(xyz_df)

#atoms for overlay / superpositioning, e.g. -a 1 2
#must be repeated for all openend xyz files, e.g. mol1.xyz mol2.xyz -a 1 2 -a 3 4
if args.atoms:
    if args.sameatoms or args.allatoms:
        #do not use these options together
        print('Warning! Do not use -a, -aa and -sa options together.')
    try:
        #extract selected atoms from the first molecule in the list of molecules, 
        #every molecule is overlaid with the first molecule
        atoms_mol1 = xyz_df_list[0][['x','y','z']].loc[args.atoms[0]]
        #generate centroid from selected atoms 
        centroid0 = np.mean(xyz_df_list[0][['x','y','z']].loc[args.atoms[0]], axis=0)
        #center the molecule by substracting the centroid from the coordinates
        xyz_df_list[0][['x','y','z']]=xyz_df_list[0][['x','y','z']].apply(lambda x: x-centroid0, axis=1)
        
        #loop over all molecules in the list of data frames
        for count, atoms in enumerate(args.atoms):
            #check for equal number of selected atoms
            if len(args.atoms[0]) == len(atoms):
                #generate centroid from selected atoms 
                centroid = np.mean(xyz_df_list[count][['x','y','z']].loc[atoms], axis=0)
                #center the molecule by substracting the centroid from the coordinates
                xyz_df_list[count][['x','y','z']]=xyz_df_list[count][['x','y','z']].apply(lambda x: x-centroid, axis=1)
                #overlay this molecule with the first molecule in the list (atoms_mol1) using the selected atoms
                xyz_df_list[count][['x','y','z']]=align_xyz(xyz_df_list[count][['x','y','z']].loc[atoms], 
                    atoms_mol1, xyz_df_list[count][['x','y','z']])
            else:
                #exit if the number of selected atoms is not equal
                print('Warning! Number of atoms must be equal. Exit.')
                sys.exit(1)
                
    except IndexError:
        #more -a options than molecules
        print('Warning! More atom pairs defined than xyz files. Exit.')
        sys.exit(1)
        
    except KeyError:
        #exit if more atoms in input than atoms in xyz file
        print('Warning! Number of atoms exceeded. Exit.')
        sys.exit(1)
    
#atoms for overlay / superpositioning, e.g. -sa 1 2
#use the same atoms in all open xyz files
if args.sameatoms:
    if args.atoms or args.allatoms:
        #do not use these options together
        print('Warning! Do not use -a, -aa and -sa options together.')
    try:
        #extract selected atoms from the first molecule in the list of molecules, 
        #every molecule is overlaid with the first molecule
        atoms_mol1 = xyz_df_list[0][['x','y','z']].loc[args.sameatoms]
        #generate centroid from selected atoms 
        centroid0 = np.mean(xyz_df_list[0][['x','y','z']].loc[args.sameatoms], axis=0)
        #center the molecule by substracting the centroid from the coordinates
        xyz_df_list[0][['x','y','z']]=xyz_df_list[0][['x','y','z']].apply(lambda x: x-centroid0, axis=1)
        
        for xyz_df in xyz_df_list:
            #generate centroid from selected atoms 
            centroid = np.mean(xyz_df[['x','y','z']].loc[args.sameatoms], axis=0)
            #center the molecule by substracting the centroid from the coordinates
            xyz_df[['x','y','z']]=xyz_df[['x','y','z']].apply(lambda x: x-centroid, axis=1)
            #overlay this molecule with the first molecule in the list (atoms_mol1) using the selected atoms
            xyz_df[['x','y','z']]=align_xyz(xyz_df[['x','y','z']].loc[args.sameatoms], 
            atoms_mol1, xyz_df[['x','y','z']])

    except KeyError:
        #exit if more atoms in input than atoms in xyz file
        print('Warning! Number of atoms exceeded. Exit.')
        sys.exit(1)

#atoms for overlay / superpositioning
#use all atoms in all open xyz files
#atom number in all xyz files must be equal
if args.allatoms:
    try:
        #extract all atoms from the first molecule in the list of molecules, 
        #every molecule is overlaid with the first molecule
        atoms_mol1 = xyz_df_list[0][['x','y','z']]
        #generate centroid from all atoms 
        centroid0 = np.mean(xyz_df_list[0][['x','y','z']], axis=0)
        #center the molecule by substracting the centroid from the coordinates
        xyz_df_list[0][['x','y','z']]=xyz_df_list[0][['x','y','z']].apply(lambda x: x-centroid0, axis=1)
        
        for xyz_df in xyz_df_list:
            #generate centroid from all atoms 
            centroid = np.mean(xyz_df[['x','y','z']], axis=0)
            #center the molecule by substracting the centroid from the coordinates
            xyz_df[['x','y','z']]=xyz_df[['x','y','z']].apply(lambda x: x-centroid, axis=1)
            #overlay this molecule with the first molecule in the list (atoms_mol1) using the selected atoms
            xyz_df[['x','y','z']]=align_xyz(xyz_df[['x','y','z']], 
            atoms_mol1, xyz_df[['x','y','z']])
    except ValueError:
        #exit if xyz files have differnet number of atoms 
        print('Warning! Number of atoms must be equal in all xyz files. Exit.')
        sys.exit(1)

#view / show the overlaid molecules 
#color by molecule
if args.viewcbm:
    #generate some empty lists
    bond_mat_list = []
    a1_a2_bond_list = []
    num_of_xyz = len(xyz_df_list)

    #setup matplotlib figure
    fig = plt.figure(figsize=(10,8))
    ax = plt.axes(projection='3d')
    #otherwise molecule looks strange
    ax.set_box_aspect((1, 1, 1))
    
    #exclude elements from molecule plot, e.g. -ee H, exlude all hydrogen atoms
    #the same elements will be excluded in all molecules
    if args.excludeEl:
        xyz_df_list = [xyz_df[~xyz_df.element.isin(args.excludeEl)] for xyz_df in xyz_df_list]
        #start index at 1 in all data frames
        xyz_df_list = [xyz_df.set_index(np.arange(1,len(xyz_df.index)+1,dtype=int)) for xyz_df in xyz_df_list]

    #exclude atoms from molecule plot, e.g. -ea 4 5, exclude atoms 4 and 5 from plot
    #the same atoms will be excluded in all molecules
    if args.excludeAt:
        if args.excludeEl:
            print('Warning! The -ee option has been applied. Numbering of atoms may have changed.')
        try:
            xyz_df_list = [xyz_df.drop(args.excludeAt) for xyz_df in xyz_df_list]
            #start index at 1 in all data frames
            xyz_df_list = [xyz_df.set_index(np.arange(1,len(xyz_df.index)+1,dtype=int)) for xyz_df in xyz_df_list]
        except KeyError:
            #exit if atom number is not in xyz file
            print('Warning! Number of atoms exceeded. Exit.')
            sys.exit(1)
    
    #get the total number of atoms from all xyz files 
    #minus excluded atoms for the atom size in the plot
    num_atom_xyz = sum(len(xyz_df) for xyz_df in xyz_df_list)
    
    #use Matplotlib colormap
    if args.cmap:
        try:
            new_colors = [plt.get_cmap(args.cmap)(1. * i/num_of_xyz) for i in range(num_of_xyz)]
            color_cycle = cycler('color',new_colors)
        except ValueError:
            print('Warning! Not a valid Matplotlib colormap. Applying default colormap instead.')
    
    #add a column with covalent radii to data frame
    #add +x % to radius
    for xyz_df in xyz_df_list:
        xyz_df['cov_radius'] = xyz_df['element'].apply(lambda x: covalent_radii[x] \
                               + covalent_radii[x] * args.radius/100)
    
    #generate a list with distance matrices
    dist_mat_list = list(pd.DataFrame(squareform(pdist(xyz_df[['x','y','z']],'euclid'))) for xyz_df in xyz_df_list)
    
    #generate a list with matrices with the sum of the atomic radii
    radii_sum_list = list(pd.DataFrame([[x + y for x in xyz_df['cov_radius']] \
                     for y in xyz_df['cov_radius']]) for xyz_df in xyz_df_list)
    
    #zero the matrix diagonal
    for radii_sum_df in radii_sum_list:
        np.fill_diagonal(radii_sum_df.values, 0)
    
    #a distance is considered as bond if the distance is smaller than the sum of the radii 
    #generate a list of bond matrices
    for count, dist_mat in enumerate(dist_mat_list):
        bond_mat_list.append(dist_mat.combine(dist_mat < radii_sum_list[count],np.multiply))
    
    #remove the upper triangle of the matrix 
    #which contains the same information as the lower triangle
    #e.g. distance C1-C2 = distance C2-C1
    #generate list of data frames of bonding atoms 
    for bond_mat in bond_mat_list:
        bond_mat.values[np.triu_indices_from(bond_mat, k=1)] = np.nan
        pd.set_option("future.no_silent_downcasting", True)
        bond_mat = bond_mat.replace(0, np.nan)
        bond_mat = bond_mat.unstack().dropna()
        bond_mat = bond_mat.reset_index(level=1)
        bond_mat = bond_mat.reset_index(level=0)
        bond_mat.columns = ['atom_1', 'atom_2','bond_length']
        a1_a2_bond_list.append(bond_mat)
    
    #get the coordinates of the atoms and put them in an array
    #also apply the color cycle / colormap
    for count, (xyz_df, style) in enumerate(zip(xyz_df_list,cycle(color_cycle))):
        atoms1=a1_a2_bond_list[count]['atom_1']
        atoms2=a1_a2_bond_list[count]['atom_2']
        atom1_coord = xyz_df[['x','y','z']].loc[atoms1+1]
        atom2_coord = xyz_df[['x','y','z']].loc[atoms2+1]
        atom1_2_coord = np.array(list(zip(atom1_coord.to_numpy(),atom2_coord.to_numpy())))
        
        #scatter (atom) plot
        ax.scatter(*xyz_df[['x','y','z']].to_numpy().T,s=np.log10(atom_scaler/num_atom_xyz),alpha=alpha_atoms,**style)

        #bonds
        for bonds in atom1_2_coord:
            ax.plot(*bonds.T,linewidth=np.log10(bond_scaler/num_atom_xyz),alpha=alpha_bonds,**style)
    #no axes
    ax.set_axis_off()
    #tight layout 
    fig.tight_layout()
    #adjust 3d drawing behavior, otherwise molecules are not correctly displayes
    set_axes_equal(ax)
    #show the plot
    plt.show()
    
    #drop the covalent radius from the data frame
    #should not be part of the altered xyz file
    for xyz_df in xyz_df_list:
        xyz_df.drop(columns='cov_radius', inplace=True)

#view / show the overlaid molecules 
#color by atoms
if args.viewcba:
    #generate some empty lists
    bond_mat_list = []
    a1_a2_bond_list = []
    
    #assign elements to colors
    metals = ['Ac', 'Ag', 'Al', 'Am', 'Ar', 'As', 'At', 'Au', 'Ba', 'Be', \
                'Bi', 'Ca', 'Cd', 'Ce', 'Cf', 'Cm', 'Co', 'Cr', 'Cs', 'Cu', \
                'Db', 'Dy', 'Er', 'Es', 'Eu', 'Fe', 'Fm', 'Fr', 'Ga', 'Gd', \
                'Ge', 'Hf', 'Hg', 'Ho', 'Hs', 'In', 'Ir', 'K', 'La', 'Li', \
                'Lr', 'Lu', 'Md', 'Mg', 'Mn', 'Mo', 'Na', 'Nb', 'Nd', 'Ni', \
                'Np', 'Os', 'Pa', 'Pb', 'Pd', 'Pm', 'Po', 'Pr', 'Pt', 'Pu', \
                'Ra', 'Rb', 'Re', 'Rf', 'Rh', 'Rn', 'Ru', 'Sc', 'Sm', 'Sn', \
                'Sr', 'Ta', 'Tb', 'Tc', 'Te', 'Th', 'Ti', 'Tl', 'Tm', 'U', \
                'V', 'W', 'Y', 'Yb', 'Zn', 'Zr']
    green = ['F','Cl']
    brown = ['Br']
    purple = ['P','I']
    orange = ['Si']
    
    #setup matplotlib figure
    fig = plt.figure(figsize=(10,8))
    ax = plt.axes(projection='3d')
    #otherwise molecule looks strange
    ax.set_box_aspect((1, 1, 1))
    
    #send warning if colormaps are used in this view mode
    if args.cmap:
        print('Warning! Option -cm not applicable.')
    
    #exclude elements from molecule plot, e.g. -ee H, exlude all hydrogen atoms
    #the same elements will be excluded in all molecules
    if args.excludeEl:
        xyz_df_list = [xyz_df[~xyz_df.element.isin(args.excludeEl)] for xyz_df in xyz_df_list]
        #start index at 1 in all data frames
        xyz_df_list = [xyz_df.set_index(np.arange(1,len(xyz_df.index)+1,dtype=int)) for xyz_df in xyz_df_list]
        
    #exclude atoms from molecule plot, e.g. -ea 4 5, exclude atoms 4 and 5 from plot
    #the same atoms will be excluded in all molecules
    if args.excludeAt:
        if args.excludeEl:
            print('Warning! The -ee option has been applied. Numbering of atoms may have changed.')
        try:
            xyz_df_list = [xyz_df.drop(args.excludeAt) for xyz_df in xyz_df_list]
            #start index at 1 in all data frames
            xyz_df_list = [xyz_df.set_index(np.arange(1,len(xyz_df.index)+1,dtype=int)) for xyz_df in xyz_df_list]
        except KeyError:
            #exit if atom number is not in xyz file
            print('Warning! Number of atoms exceeded. Exit.')
            sys.exit(1)
            
    #start index from zero
    xyz_df_list = [xyz_df.set_index(np.arange(len(xyz_df.index),dtype=int)) for xyz_df in xyz_df_list]
    
    #get the total number of atoms from all xyz files 
    #minus excluded atoms for the atom size in the plot
    num_atom_xyz = sum(len(xyz_df) for xyz_df in xyz_df_list)
    
    #add a column with covalent radii to data frame
    #add +x % to radius
    for xyz_df in xyz_df_list:
        xyz_df['cov_radius'] = xyz_df['element'].apply(lambda x: covalent_radii[x] + covalent_radii[x] * args.radius/100)
        
    #generate a list with distance matrices    
    dist_mat_list = list(pd.DataFrame(squareform(pdist(xyz_df[['x','y','z']],'euclid'))) for xyz_df in xyz_df_list)
    
    #generate a list with matrices with the sum of the atomic radii
    radii_sum_list = list(pd.DataFrame([[x + y for x in xyz_df['cov_radius']] for y in xyz_df['cov_radius']]) for xyz_df in xyz_df_list)
    
    #zero the matrix diagonal
    for radii_sum_df in radii_sum_list:
        np.fill_diagonal(radii_sum_df.values, 0)
    
    #a distance is considered as bond if the distance is smaller than the sum of the radii 
    #generate a list of bond matrices
    for count, dist_mat in enumerate(dist_mat_list):
        bond_mat_list.append(dist_mat.combine(dist_mat < radii_sum_list[count],np.multiply))
    
    #remove the upper triangle of the matrix 
    #which contains the same information as the lower triangle
    #e.g. distance C1-C2 = distance C2-C1
    #generate list of data frames of bonding atoms 
    for bond_mat in bond_mat_list:
        bond_mat.values[np.triu_indices_from(bond_mat, k=1)] = np.nan
        pd.set_option('future.no_silent_downcasting', True)
        bond_mat = bond_mat.replace(0, np.nan)
        bond_mat = bond_mat.unstack().dropna()
        bond_mat = bond_mat.reset_index(level=1)
        bond_mat = bond_mat.reset_index(level=0)
        bond_mat.columns = ['atom_1', 'atom_2','bond_length']
        a1_a2_bond_list.append(bond_mat)
        
    #get the coordinates of the atoms and put them in an array
    #also apply the color cycle / colormap
    for count, (xyz_df, style) in enumerate(zip(xyz_df_list,cycle(color_cycle))):
        atoms1=a1_a2_bond_list[count]['atom_1']
        atoms2=a1_a2_bond_list[count]['atom_2']
        atom1_coord = xyz_df[['x','y','z']].loc[atoms1]
        atom2_coord = xyz_df[['x','y','z']].loc[atoms2]
        atom1_2_coord = np.array(list(zip(atom1_coord.to_numpy(),atom2_coord.to_numpy())))
        
        #clumsy but safe
        #for assigning colors to different elemments
        carr = xyz_df.index[xyz_df['element'].isin(['C'])].tolist()
        harr = xyz_df.index[xyz_df['element'].isin(['H'])].tolist()
        narr = xyz_df.index[xyz_df['element'].isin(['N'])].tolist()
        oarr = xyz_df.index[xyz_df['element'].isin(['O'])].tolist()
        sarr = xyz_df.index[xyz_df['element'].isin(['S'])].tolist()
        marr = xyz_df.index[xyz_df['element'].isin(metals)].tolist()
        grarr = xyz_df.index[xyz_df['element'].isin(green)].tolist()
        brarr = xyz_df.index[xyz_df['element'].isin(brown)].tolist()
        parr = xyz_df.index[xyz_df['element'].isin(purple)].tolist()
        orarr = xyz_df.index[xyz_df['element'].isin(orange)].tolist()
        restarr = [item for item in atoms1 if item not in carr]
        restarr = [item for item in restarr if item not in harr]
        restarr = [item for item in restarr if item not in narr]
        restarr = [item for item in restarr if item not in oarr]
        restarr = [item for item in restarr if item not in sarr]
        restarr = [item for item in restarr if item not in marr]
        restarr = [item for item in restarr if item not in grarr]
        restarr = [item for item in restarr if item not in brarr]
        restarr = [item for item in restarr if item not in parr]
        restarr = [item for item in restarr if item not in orarr]
        
        #scatter (atom) plot
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[carr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='black')
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[harr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='tan')
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[narr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='blue')
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[oarr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='green')
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[sarr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='yellow')
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[marr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='red')
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[grarr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='green')
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[brarr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='brown')
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[parr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='purple')
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[orarr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='orange')
        ax.scatter(*xyz_df[['x','y','z']].to_numpy()[restarr].T,s=np.log10(atom_scaler/num_atom_xyz),
            alpha=alpha_atoms,color='gray')
        
        #bonds
        for bonds in atom1_2_coord:
            ax.plot(*bonds.T,linewidth=np.log10(bond_scaler/num_atom_xyz),alpha=alpha_bonds,color='gray')
    
    #no axes
    ax.set_axis_off()
    #tight layout 
    fig.tight_layout()
    #adjust 3d drawing behavior, otherwise molecules are not correctly displayes
    set_axes_equal(ax)
    #show the plot
    plt.show()
    
    #drop the covalent radius from the data frame
    #should not be part of the altered xyz file
    for xyz_df in xyz_df_list:
        xyz_df.drop(columns='cov_radius', inplace=True)

#save the altered xyz files
if args.save:
    #in case a multi xyz file was processed
    if is_trj:
        head_list=['\n'.join(head) for head in zip(head_list[0::2], head_list[1::2])]
        try:
            for count, xyz_df in enumerate(xyz_df_list):
                my_numpy = xyz_df.to_numpy()
                np.savetxt(os.path.splitext(args.filename[0])[0] + '-' + str(count+1) + '-mod.xyz', 
                           my_numpy, fmt='%-2s  %12.8f  %12.8f  %12.8f', delimiter='', 
                           header=head_list[count], comments='')
        #write error -> exit here
        except IOError:
            print("Write error. Exit.")
            sys.exit(1)
    #in case several single xyz files have been processed
    else:
        try:
            for count, file in enumerate(args.filename):
                my_numpy = xyz_df_list[count].to_numpy()
                np.savetxt(os.path.splitext(file)[0] +'-mod.xyz', my_numpy,
                           fmt='%-2s  %12.8f  %12.8f  %12.8f', delimiter='', 
                           header=head_list[count], comments='')
        #write error -> exit here
        except IOError:
            print("Write error. Exit.")
            sys.exit(1)

#save the altered xyz files in a single multi xyz file
#uses io.StringIO, single xyz files will be combined and saved
if args.savetr:
    #in case a multi xyz file was processed
    if is_trj:
        all_xyz_io = io.StringIO()
        head_list=['\n'.join(head) for head in zip(head_list[0::2], head_list[1::2])]
        for count, xyz_df in enumerate(xyz_df_list):
            my_numpy = xyz_df.to_numpy()
            np.savetxt(all_xyz_io, my_numpy,
                       fmt='%-2s  %12.8f  %12.8f  %12.8f', delimiter='', 
                       header=head_list[count], comments='')
        try:
            with open(os.path.splitext(file)[0] +'-mod.xyz', 'w') as all_xyz_trj:
                all_xyz_trj.write(all_xyz_io.getvalue())
        #write error -> exit here
        except IOError:
            print("Write error. Exit.")
            sys.exit(1)
    #in case several single xyz files have been processed
    else:
        all_xyz_io = io.StringIO()
        for count, file in enumerate(args.filename):
            my_numpy = xyz_df_list[count].to_numpy()
            np.savetxt(all_xyz_io, my_numpy,
                       fmt='%-2s  %12.8f  %12.8f  %12.8f', delimiter='', 
                       header=head_list[count], comments='')
        try:
            with open('all_xyz_trj.xyz', 'w') as all_xyz_trj:
                all_xyz_trj.write(all_xyz_io.getvalue())
        #write error -> exit here
        except IOError:
            print("Write error. Exit.")
            sys.exit(1)
