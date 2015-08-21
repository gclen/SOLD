Level Diagrams
--------------

SOLD (Spatially Oriented Level Diagrams) produces diagrams that plot the
energies of single particle orbitals relative to one another.
Additionally, the lines show the spatial extent of the orbital with
darker colours indicating more weight in that region. The wavefunction
can be integrated along an arbitrary vector (given two points) which will yield different diagrams.

### Implementation

This code requires a Gaussian formatted checkpoint file (.fchk), G09 output file, number of occupied and virtual orbitals, and two points on the line of integration as input. SOLD is written as a plugin to ase-gui which is a GUI for the python module Atomic Simulation Environment. Fig. 1 shows a screenshot of the menu for the SOLD plugin. It can be used to set all the input required. The main advantage of using a GUI is one is able to select points by clicking on the desired atoms(Fig. 2). 

<figure>
<img src="Image_Files/SOLD_GUI.png" width="350">
<figcaption> <br>Figure 1: Screenshot of the SOLD tool menu in ase-gui. </figcaption>
</figure>
<br>
<figure>
<img src="Image_Files/ASE_GUI_VOLF.png" width="350">
<figcaption> <br>Figure 2: Screenshot of a molecule (VOLF) visualized in ase-gui. The two fluorine atoms are selected to use as input for SOLD. </figcaption>
</figure>
<br>

For the given number of orbitals, their corresponding cube files are generated using the cubegen utility. The points are used to find the line of integration using the method described in the previous section. For each orbital,  the square of every point in the cube file is projected onto the line and binned appropriately. The sums are then scaled to a value between 0-1 on a global basis. That is, they are scaled relative to the maximum and minimum sums over the set of all orbitals. After all orbitals have been projected and scaled, they are plotted at their respective energies using the scaled value as a colour map. If two orbitals are degenerate (within 0.05 eV) then they are shifted. This is indicated on the graph. The shift also prevents overlap of the two lines, which would distort the data. Additional space is left on the graph, on the left and right sides, to account for labels and bubble diagrams. The positions of selected atoms projected onto the vector are shown as dotted vertical lines.

### Test Cases

To test the integration and scaling, cube files containing spheres were
generated. Fig. 3 shows three spheres placed
equidistantly along the z-axis. The integration was done along the
z-axis. As expected this level diagram contains three identical lines
that are spaced appropriately. Similarly, Fig. 4 shows
three groups of spheres on the z-axis. The middle group contains five
spheres, while the outer groups contain only one. Integrating along the
z-axis yields more weight on the middle group which is the desired
result.

<figure>
<img src="Image_Files/Level_diagram_spheres_001.png" width="550">
<figcaption> <br>Figure 3: Three equidistant spheres aligned on the z-axis as a test case. The corresponding level diagram is shown on the right. </figcaption>
</figure>
<br>
<figure>
<img src="Image_Files/Level_diagram_spheres_151.png" width="550">
<figcaption> <br>Figure 4: Spheres at three points along the z-axis as a test case. The middle point contains more spheres which yields more weight. The corresponding level diagram is shown on the right. </figcaption>
</figure>
<br>

Lithium hydride (LiH) was also investigated as a benchmark. Fig. 5 shows the level diagram where the integration was done along the line connecting the two atoms. Atomic positions are denoted by vertical dotted lines. Comparing the levels with their corresponding bubble diagrams shows good agreement. Of particular note, there is an area of minimal weight on Li in LUMO+3 which can be seen on the scaled line. 

<figure>
<img src="Image_Files/Level_diagram_LiH_SOLD.png" width="550">
<figcaption> <br>Figure 5: Level diagram for LiH. Integration was done along the line connecting the two atoms. </figcaption>
</figure>
