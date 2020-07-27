import meshio
import numpy as np

larger_sphere = meshio.read("larger_sphere.msh")
smaller_sphere = meshio.read("smaller_sphere.msh")

large = larger_sphere.points
large = large[:1490]

small = smaller_sphere.points
small_x = small[:,0]
small_y = small[:,1]
small_z = small[:,2]

large_x = large[:,0]
large_y = large[:,1]
large_z = large[:,2]

displacements = np.zeros(np.shape(small))
for i in range(0,1490):
    dist = np.zeros(np.size(small_x))
    for j in range(0,1490):
        dist[j] = np.sqrt((small_x[j]-large_x[i])**2 + (small_y[j]-large_y[i])**2 + (small_z[j]-large_z[i])**2)

    index = np.where(dist == np.amin(dist))
    index = index[0][0]
    displacements[i,0] = (small_x[index]-large_x[i]) 
    displacements[i,1] = (small_y[index]-large_y[i]) 
    displacements[i,2] = (small_z[index]-large_z[i])

disps = displacements
print(len(disps))
displacements = np.hstack((large, displacements))
#print(displacements)

#norm_displacement = np.linalg.norm(displacements, axis=1)
np.savetxt('displacements.csv', displacements, delimiter=',')
np.savetxt('displacements_only.csv', disps, delimiter=',')
#np.savetxt('disps.csv', disps, delimiter=',')