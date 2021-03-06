import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, interp2d

tpc = [[0, 0], [256, 1036]]  # cm (x, z)

def get_line(point1, point2):
    point2 = [1056.8, 103.]
    m = (point2[1] - point1[1]) / (point2[0] - point1[0])
    b = point1[1] - m * point1[0]
    return m, b

def calc_residual(track):
    # these are uboone coordinates!
    x, y, z = track[1], track[2], track[3]
    middle_point = int(len(z)/2) - 150
    m, b = get_line([z[middle_point], x[middle_point]], [z[-1],x[-1]])

    #plot_line(m, b, z[0], z[-1])
    #plot_track(track)

    r = np.abs(m*z + b - x)
    #plot_residuals(z, x, r)
    #plt.show()

    return np.array([z, x, r])

def plot_track(track):
    x, y, z = track[1], track[2], track[3]
    plt.scatter(z, x)

def plot_line(m, b, start_z=0, stop_z=1036):
    x = np.linspace(0, 1056.8, 100)
    plt.plot(x, m*x + b)

def plot_residuals(z, x, residuals):
    plt.scatter(z, x, c=residuals*10)

filename = "laser_tracks.npy"
tracks = np.load(filename)
residuals = np.array([[0],[0],[0]])

for track in tracks:
    res = calc_residual(track)
    residuals = np.append(residuals, res, axis=1)

residuals = np.delete(residuals, 0, 1)

z = residuals[0]
x = residuals[1]
res = residuals[2]
grid_z, grid_x = np.mgrid[0:1057:8,0:257:8 ]


grid_z2 = griddata((z, x), res, (grid_z, grid_x), method='nearest')
#print grid_z2.T
print np.max(grid_z2)

v = np.linspace(0, 40, 9)
CS = plt.contourf(grid_z, grid_x, grid_z2, v, cmap=plt.viridis())


plt.colorbar(CS, label="Distrortion [cm]")

plt.scatter(z, x, alpha=.01, s=2)

plt.xlim([0,1036])
plt.ylim([0,250])
plt.xlabel("z [cm]")
plt.ylabel("x (drift) [cm]")
plt.show()