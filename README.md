# xyzoverlay
<p align="center">
<img width="800" alt="fig10" src="/examples/figure_10.png">
</p>
A Python 3 script for overlaying or superimposing two or more molecules. The overlaid / superimposed molecules can be displayed and the modified coordinates can be saved as xyz files. Coloring by atom or molecule is possible. All subsequent molecules are superimposed on the first molecule.

## External modules
`pandas`, `numpy`, `scipy`, `matplotlib`

## Quick start
```console
python3 xyzoverlay.py mol-1.xyz mol-2.xyz -a 1 2 3 -a 3 4 5 -vcm -s
```
Open two xyz files (`mol-1.xyz` & `mol-2.xyz`) and overlay the atoms 1 2 3 (`-a 1 2 3`) of molecule 1 with the atoms 3 4 5 (`-a 3 4 5`) of molecule 2. View the superimposed molecules with different colors (`-vcm`) and save the modified xyz files (`-s`).

<p align="center">
<img width="400" alt="fig1" src="/examples/figure_1.png">
</p>

```console
python3 xyzoverlay.py mol-1.xyz mol-2.xyz mol-3.xyz -sa 1 2 3 4 -vca -ee H -st
```
Open three molecules and overlay the first four atoms in all opened xyz files (`-sa 1 2 3 4`). All subsequent molecules are superimposed on the first molecule (`mol-1.xyz`). View the superimposed molecules with differently colored atoms (`-vca`). Do not display hydrogen atoms (`-ee H`). Save the modified xyz files as a single (multi) xyz file (`-st`).

<p align="center">
<img width="400" alt="fig2" src="/examples/figure_2.png">
</p>

```console
python3 xyzoverlay.py multi-mol_trj.xyz -aa -vcm -cm coolwarm -ee H
```
Open a xyz file that contains several molecules. Use all atoms for superimposing (`-aa`). All subsequent molecules are superimposed on the first molecule in the xyz file. View the overlaid molecules with different colors (`-vcm`). Use the `Matplotlib` colormap `coolwarm` for coloring (`-cm coolwarm`). Do not display hydrogen atoms (`-ee H`).

<p align="center">
<img width="400" alt="fig3" src="/examples/figure_3.png">
</p>

---
**NOTE**

The script only performs something when either a view option (`-vcm` or `-vca`) or a save option (`-s` or `-st`) or both are specified.

The script only computes something when one of the following options is specified: `-a`, `-sa` or `-aa`.

---


## Command-line options
- `filename` , required: filename(s), e.g. `my_xyz.xyz` or `my_xyz-1.xyz my_xyz-2.xyz ...` or `multi_xyz.xyz`. The first two lines will be ignored, file format must be `element x y z`, cartesian coordinates, units in Å. In case of xyz files that contain several molecules (multi xyz files or xyz trajectory) the same format is expected with no extra separation between the molecules (empty lines for example). Only one of these multi xyz files can be opened at a time.    
- `-a` `atom(s)`, optional:  define atoms for superimposing molecules. The order corresponds to the order of the xyz files, e.g. `xyzoverlay.py mol1.xyz mol2.xyz -a 1 2 3 -a 4 5 6` (atoms `1 2 3` from `mol1.xyz` and atoms `4 5 6` from `mol2.xyz`).
- `-sa` `atom(s)`, optional: define the same atoms in all open xyz files for superimposing molecules, e.g. `xyzoverlay.py mol1.xyz mol2.xyz -sa 1 2 3` (atoms `1 2 3` from `mol1.xyz` and `mol2.xyz`).
-  `-aa`, optional: use all atoms in all open xyz files for superimposing molecules, e.g. `xyzoverlay.py mol1.xyz mol2.xyz -aa`. The molecules must have the same number of atoms.
-  `-vca`, optional: display the overlaid molecules. Color by atom type.
-  `-vcm`, optional: display the overlaid molecules. Molecules have different colors. With the `-cm` option a colormap can be defined, otherwise the standard color map will be applied. 
-  `-cm` `colormap`, optional: choose one of the many `Matplotlib` colormaps, e.g. `-cm plasma`. Only active with the `-vcm` option.
-   `-r` `radius`, optional: enlarge atomic radii by x%, e.g. `-r 10` enlarge atomic radii by 10%. Default is 8%. Atomic radii determine whether a bond between atoms is drawn or not. Only active with `-vca` or `-vcm`.
-  `-ee` `element(s)`, optional: exclude elements from the molecule plot, e.g. `-ee H N` exclude hydrogen and nitrogen atoms from the molecular representation. Only active with `-vca` or `-vcm`.
-  `-ea` `atom(s)`, optional: exclude atoms from the molecule plot, e.g. `-ea 4 7 8` exclude atoms 4,7 and 8 from the molecular representation. Only active with `-vca` or `-vcm`. The specified atoms are not displayed in all opened xyz files. Numbering of atoms can be different if `-ee` option was used. 
-  `-s`, save the superimposed / aligned data as single xyz file(s). The filename(s) of the modified xyz file(s) is `*-mod.xyz`. A multi xyz file is split into individual files.
- `-st`, save the superimposed / aligned data in one (multi) xyz file. The filename of the modified xyz file is `all_xyz_trj.xyz` if different xyz files were processed. In case of a multi xyz file the filename of the modified xyz file is `*-mod.xyz`.

## Remarks
- Only the standard XYZ file format is supported. Also xyz files with more than one molecule (multi xyz file or xyz trajectory) are supported. Only one multi xyz file can be processed at a time.
- If the script is opened with no further options the script does not compute or display anything. A view option (`-vcm` or `-vca`) or a save option (`-s` or `-st`) or both and/or `-a`, `-sa` or `-aa` options are necessary.
- The script can easily process several hundred small molecules in a reasonable time. The `Matplotlib`, however, is not intended to display so many data points or graphs. Zooming and rotating becomes very slow. `Matplotlib` also does not have full camera controls ("roll" is missing).
- There are no options to change the size or color of the displayed atoms or bonds. At the beginning of the script a few parameters are collected that can be adjusted quickly. There is also no labeling option.
- The output of the modified xyz-files depends on whether several single xyz files or one multi xyzfile have been processed. Please note the `-s` and `-st` options.

## Examples

### Example 1:
```console
python3 xyzoverlay.py cupor.xyz cutpp.xyz fepor.xyz -sa 1 2 3 -vcm
```
Open 3 xyz files, select atoms 1, 2 & 3 in all molecules for the superposition (`-sa 1 2 3`), display the result, color by molecule (`-vcm`).

<p align="center">
<img width="800" alt="example 1" src="/examples/show_use1.gif">
</p>
