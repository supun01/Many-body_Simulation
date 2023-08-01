#Lab 2 - Many-body Simulation
#Hapuarachchi H.A.S.V - ET/2018/018 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig = plt.figure()
ax = plt.axes(projection='3d')

#mass of plantes (Sun, Earth, Venus)
m1 = 1.989e30  # mass of the Sun
m2 = 5.972e24  # mass of the Earth
m3 = 0.073e24 # mass of the Venus

r1 = np.array([0.0, 0.0,0.0])  # initial position of the Sun
r2 = np.array([149.0e9, 0.0,0.0])  # initial position of the Earth
r3 = np.array([149.984e9,0.0,0.0]) # initial position of the Venus

v1 = np.array([0.0,0.0,0.0]) # initial velocity of the Sun
v2 = np.array([0.0,29800.0,0.0])  # initial velocity of the Earth
v3 = np.array([0.0,30800.0,0.0]) # initial velocity of the Venus

G = 6.67430e-11 #Universal gravitational constant
dt = 50000 # time change
num_steps = 3600 #num. of time steps

# Initialize Zero arrays to store positions of the plants 
positions1 = np.zeros((num_steps+1,3))
positions2 = np.zeros((num_steps+1, 3))
positions3 = np.zeros((num_steps+1,3))

positions1[0]=r1
positions2[0]=r2
positions3[0]=r3

# initial acceleration calculations when time is 1
r = r2 - r1
r_mag = np.linalg.norm(r)
r_es=r3-r2
r_es_mag = np.linalg.norm(r_es)
r_ss=r3-r1
r_ss_mag = np.linalg.norm(r_ss)

F = (G * m1 * m2) / (r_mag ** 2)
direction = r / r_mag
F2 = -F * direction

F3= (G * m2 * m3) / (r_es_mag ** 2)
direction_es =r_es/r_es_mag
F_es= -F3*direction_es

F4= (G * m1 * m3) / (r_ss_mag ** 2)
direction_es =r_ss/r_ss_mag
F_ss= -F4*direction_es

# Calculate the acceleration at t=0
a1= -(F2/m1 + F_ss/m1)
a2 = F2 / m2 + F_es/m2
a3 = F_ss/m3 + (-F3/m3)

for step in range(num_steps):

    # New position of the Earth
    if step==0:
      # verlet algorithm for position calculations at t=0
      r1 = r1 + v1*dt + (0.5*a1*dt**2)
      r2 = r2 + v2*dt + (0.5*a2*dt**2) 
      r3 = r3 + v3*dt + (0.5*a3*dt**2)
    else:
      # verlet algorithm for position calculations at t>0
      r1 = 2*r1 - positions1[step-1] + a1*dt**2
      r2 = 2*r2 - positions2[step-1] + a2*dt**2 
      r3 = 2*r3 - positions3[step-1] + a3*dt**2

    # Calculate the vector r (distance between the Sun and the Earth)
    r = r2 - r1
    r_mag = np.linalg.norm(r) 

    r_es=r3-r2
    r_es_mag = np.linalg.norm(r_es)

    r_ss=r3-r1
    r_ss_mag = np.linalg.norm(r_ss)


    # Calculate the gravitational force
    F = (G * m1 * m2) / (r_mag ** 2)
    direction = r / r_mag
    F2 = -F * direction

    F3= (G * m2 * m3) / (r_es_mag ** 2)
    direction_es =r_es/r_es_mag
    F_es= -F3*direction_es

    F4= (G * m1 * m3) / (r_ss_mag ** 2)
    direction_es =r_ss/r_ss_mag
    F_ss= -F4*direction_es

    # Calculate the net accelerations from the gravitational forces
    a1 = -(F2 / m1 + F_ss/m1)
    a2 = F2 / m2 + F_es/m2
    a3 = F_ss/m3 + (-F3/m3)

    # New velocity of the Earth
    v1 = v1 + a1 * dt
    v1_mag=np.linalg.norm(v1)
    v2 = v2 + a2 * dt
    v2_mag=np.linalg.norm(v2)
    v3 = v3 + a3 * dt
    v3_mag=np.linalg.norm(v3)

    # Store positions for visualization
    positions1[step+1]=r1
    positions2[step+1]=r2 
    positions3[step+1]=r3

# Define the coordinates
object1_x =positions1[:,0]
object1_y =positions1[:,1]
object1_z =positions1[:,2]

object2_x= positions2[:,0]
object2_y= positions2[:,1]
object2_z= positions2[:,2]

object3_x =positions3[:,0]
object3_y =positions3[:,1]
object3_z =positions3[:,2]

# Plotting the positions of the Sun, Earth, Venus
plt.plot(r1[0],r1[1],"o", label="Sun")
plt.plot(positions2[:,0], positions2[:, 1],positions2[:,2],".",label="Earth")
plt.plot(positions3[:,0], positions3[:, 1],positions3[:,2],".",label="Venus")
plt.xlabel("X - Axis")
plt.ylabel("Y - Axis")
plt.title("Many-Body Simulation")
plt.legend()
plt.axis("equal")
plt.show()