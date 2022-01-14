
# One-molecule-visualization
Script for visualization of small molecules.


Publication/presentation 2D plots of moderately size molecules can be easily created from .xyz files.
Script requires compiled Julia with packages: ArgParse, Interact, Luxor, Colors and Blink.
Run script from terminal specifying path to the .xyz file:
```
 $ ./one_mol_vis.jl -f /path/to/molecule.xyz 
 ```
 Interactive window will appear which allows you rotate and resize your molecule.
 Molecule plot is automaticaly saved in working directory every time you change the orientation.
 Output format can be specified with option `-o [png, svg::Default]` (pdf format might cause some issues).
 
*Check the dictionaries at the beginning of the script, expand them accordingly to cover all atom types presented in your xyz file.*

![one_mol_vis](https://user-images.githubusercontent.com/43886886/148700795-dbea7815-8d8a-49ed-a7b6-941573e8652b.png)

### Coming soon:

- [ ] More complete support of various atom types.
- [x] Visualization of normal modes from quantum chemistry calculations. 
    - [x] Support for Orca calculations.
    - [ ] Support for Turbomole calculations.
