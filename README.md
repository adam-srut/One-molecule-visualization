
# One-molecule-visualization
Script for visualization of small molecules.


2D plots of moderately size molecules can be easily created from .xyz files.
Script requires compiled Julia with packages: ArgParse, Interact, Luxor, Colors and Blink and PeriodicTable.
Run script from terminal specifying path to the .xyz file:
```
 one_mol_vis.jl /path/to/molecule.xyz 
 ```
 Interactive window will appear which allows you to rotate and resize your molecule.
 Molecule plot is automaticaly saved in the working directory every time you change the orientation of the molecule.
 Try `$ one_mol_vis.jl --help` to list all the options.
 

#### Examples of different representations and visualizations of vibrational modes:

![example](https://user-images.githubusercontent.com/43886886/208878571-ca1aee93-6704-40cd-81cd-aa646110f85d.png)



### Features:

- [x] Whole periodic table is supported through PeriodicTable package.
- [x] Visualization of normal modes from quantum chemistry calculations. 
    - [ ] Support for Orca calculations.
    - [x] Support for Turbomole calculations, now through Molden format (use `tm2molden`).
    - [x] Support for MOLDEN format.
    
