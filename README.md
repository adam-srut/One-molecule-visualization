
# One-molecule-visualization
Script for visualization of small molecules.


2D plots of moderately size molecules can be easily created from .xyz files.
Script requires compiled Julia with packages: ArgParse, Interact, Luxor, Colors and Blink.
Run script from terminal specifying path to the .xyz file:
```
 $ one_mol_vis.jl /path/to/molecule.xyz 
 ```
 Interactive window will appear which allows you to rotate and resize your molecule.
 Molecule plot is automaticaly saved in the working directory every time you change the orientation of the molecule.
 Output format can be specified with option `-o [png, Default: svg]` (pdf format might not work properly).
 Try `$ one_mol_vis.jl --help` to list all the options.
 
*Check the dictionaries in the atom_types.jl, expand them accordingly to cover all atom types presented in your xyz file.*

Modification of the shebang line of *one_mol_vis.jl* might be needed. e.g. `#! /usr/bin/env julia`

#### Examples of different representation and visualization of vibrational modes:

![example](https://user-images.githubusercontent.com/43886886/208878571-ca1aee93-6704-40cd-81cd-aa646110f85d.png)



### Coming soon:

- [ ] More complete support of various atom types.
- [x] Visualization of normal modes from quantum chemistry calculations. 
    - [x] Support for Orca calculations.
	- [x] Support for ADF calculations (--adf flag reads directly the output file)
    - [x] Support for Turbomole calculations.
