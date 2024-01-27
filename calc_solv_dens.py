import MDAnalysis as mda
import sys
from MDAnalysis.analysis import density
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import matplotlib.pyplot as plt


top  = sys.argv[1]
traj = sys.argv[2]
out  = sys.argv[3]
####################################
u = mda.Universe(top, traj)

selection = "(name OH2)"

print("Selection = ", selection)

solv = u.select_atoms(selection, updating=True)
dens = density.DensityAnalysis(solv, delta=1.0)

print("\nRunning density analysis...")
dens.run()
dens.density.convert_density("A^{-3}") # inverse cubic Angstrom
#dens.density.convert_density('TIP3P')  # in relation to densities reported in
                   #W. Jorgensen, C. Jenson, J Comp Chem 19 (1998), 1179-1186

#grid = dens.results.density.grid
#print(grid)

dens.results.density.export(out, type="double")
print("Done!")


grid = dens.density.grid
print(grid.shape)

avg = grid.mean(axis=1)
print(avg.shape)
#exit()

fig, ax = plt.subplots()

B = np.einsum('ab->ba', avg)

im = ax.imshow(B, interpolation="bicubic")
cbar = plt.colorbar(im)
cbar.set_label('Mean density of water over TIP4P literature value')
plt.xlabel('X-axis ($\AA$)')
plt.ylabel('Z-axis ($\AA$)')
ax.invert_yaxis()

plt.show()
