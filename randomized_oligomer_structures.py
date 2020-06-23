import numpy as np
import os
import math
import random
from numpy import genfromtxt

geo = genfromtxt('mon_geo.xyz', delimiter=' ')

def yz_rotate(geom, yz_angle=np.pi/4): 
    """This will take the original xyz geometry, rotate by radians and displace by x,y,z"""   
    
    # builds new array with correct atom order
    new_geom = np.zeros(((int(len(geom[:,0]))), 4))
    k=1
    for i in geom[:,0]:
        new_geom[k-1:k,0] = i
        k += 1

    # rotates molecule by given yz-angle in radians



    Y = geom[:,2]
    Z = geom[:,3]
    theta_0 = np.zeros(int(len(geom[:,0])))
    hyp = np.zeros(int(len(geom[:,0])))
    k=1
    for atom_y, atom_z in np.nditer([Y, Z]):
        #theta_0[k-1:k] = math.atan(atom_y/atom_z)
        if atom_y > 0 and atom_z > 0:
            theta_0[k-1:k] = math.degrees(math.atan(atom_y/atom_z))
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_z**2)

        elif atom_y < 0 and atom_z > 0:
            theta_0[k-1:k] = math.degrees(math.atan(atom_y/atom_z))
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_z**2)
        
        elif atom_y > 0 and atom_z < 0:
            theta_0[k-1:k] = math.degrees(math.atan(atom_y/atom_z)) + 180
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_z**2)
        
        elif atom_y < 0 and atom_z < 0:
            theta_0[k-1:k] = math.degrees(math.atan(atom_y/atom_z)) + 180
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_z**2)
        
        elif atom_y == 0 and atom_z == 0:
            theta_0[k-1:k] = 0
            hyp[k-1:k] = 0
        
        elif atom_y > 0 and atom_z == 0:
            theta_0[k-1:k] = math.degrees(np.pi/2)
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_z**2)
        
        elif atom_y < 0 and atom_z == 0:
            theta_0[k-1:k] = math.degrees(3*np.pi/2)
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_z**2)

        elif atom_y == 0 and atom_z > 0:
            theta_0[k-1:k] = math.degrees(0)
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_z**2)

        elif atom_y == 0 and atom_z < 0:
            theta_0[k-1:k] = math.degrees(np.pi)
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_z**2)

        k += 1
    #print(theta_0)
    #print(hyp)


    theta_fin = np.zeros(int(len(geom[:,0])))
    k = 1
    
    for i in range(len(theta_0)):
        if math.isnan(theta_0[i]) == True:
            theta_0[i] = 0
    #print(theta_0[:])

    theta_fin = np.zeros(int(len(geom[:,0])))
    k = 1
    for theta in theta_0[:]:
        #theta_fin[k-1:k] = theta + angle
        theta_fin[k-1:k] = theta + math.degrees(yz_angle)
        k += 1
    #print(theta_fin)



    
    xs = geom[:,1]
    ys = np.zeros(int(len(geom[:,0])))
    zs = np.zeros(int(len(geom[:,0])))
    k = 1
    for theta, hypot in np.nditer([theta_fin, hyp]):
        ys[k-1:k] = round(hypot*math.sin(math.radians(theta)), 8)
        zs[k-1:k] = round(hypot*math.cos(math.radians(theta)), 8)
        k += 1
    #print(ys)
    #print(zs)
    

    k=1
    for x in xs[:]:
        new_geom[k-1:k,1] = x
        k += 1


    k=1
    for y in ys[:]:
        new_geom[k-1:k,2] = y
        k += 1

    k=1
    for z in zs[:]:
        new_geom[k-1:k,3] = z
        k += 1
    
  ####################################################  # rotates molecule by given yz-angle in radians




    return new_geom

def xy_rotate(geom, xy_angle=np.pi/4):
    
    new_geom = np.zeros(((int(len(geom[:,0]))), 4))
    k=1
    for i in geom[:,0]:
        new_geom[k-1:k,0] = i
        k += 1

    Y = geom[:,2]
    X = geom[:,1]
    theta_0 = np.zeros(int(len(geom[:,0])))
    hyp = np.zeros(int(len(geom[:,0])))
    k=1
    for atom_y, atom_x in np.nditer([Y, X]):
        #theta_0[k-1:k] = math.atan(atom_y/atom_z)
        if atom_y > 0 and atom_x > 0:
            theta_0[k-1:k] = math.degrees(math.atan(atom_y/atom_x))
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_x**2)

        elif atom_y < 0 and atom_x > 0:
            theta_0[k-1:k] = math.degrees(math.atan(atom_y/atom_x))
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_x**2)
        
        elif atom_y > 0 and atom_x < 0:
            theta_0[k-1:k] = math.degrees(math.atan(atom_y/atom_x)) + 180
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_x**2)
        
        elif atom_y < 0 and atom_x < 0:
            theta_0[k-1:k] = math.degrees(math.atan(atom_y/atom_x)) + 180
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_x**2)
        
        elif atom_y == 0 and atom_x == 0:
            theta_0[k-1:k] = 0
            hyp[k-1:k] = 0
        
        elif atom_y > 0 and atom_x == 0:
            theta_0[k-1:k] = math.degrees(np.pi/2)
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_x**2)
        
        elif atom_y < 0 and atom_x == 0:
            theta_0[k-1:k] = math.degrees(3*np.pi/2)
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_x**2)

        elif atom_y == 0 and atom_x > 0:
            theta_0[k-1:k] = math.degrees(0)
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_x**2)

        elif atom_y == 0 and atom_x < 0:
            theta_0[k-1:k] = math.degrees(np.pi)
            hyp[k-1:k] = math.sqrt(atom_y**2 + atom_x**2)

        k += 1
    #print(theta_0)
    #print(hyp)


    theta_fin = np.zeros(int(len(geom[:,0])))
    k = 1
    
    for i in range(len(theta_0)):
        if math.isnan(theta_0[i]) == True:
            theta_0[i] = 0
    #print(theta_0[:])

    theta_fin = np.zeros(int(len(geom[:,0])))
    k = 1
    for theta in theta_0[:]:
        #theta_fin[k-1:k] = theta + angle
        theta_fin[k-1:k] = theta + math.degrees(xy_angle)
        k += 1
    #print(theta_fin)



    
    zs = geom[:,3]
    ys = np.zeros(int(len(geom[:,0])))
    xs = np.zeros(int(len(geom[:,0])))
    k = 1
    for theta, hypot in np.nditer([theta_fin, hyp]):
        ys[k-1:k] = round(hypot*math.sin(math.radians(theta)), 8)
        xs[k-1:k] = round(hypot*math.cos(math.radians(theta)), 8)
        k += 1
    #print(ys)
    #print(zs)
    

    k=1
    for y in ys[:]:
        new_geom[k-1:k,2] = y
        k += 1

    k=1
    for x in xs[:]:
        new_geom[k-1:k,1] = x
        k += 1
    k=1
    for z in zs[:]:
        new_geom[k-1:k,3] = z
        k += 1

    return new_geom
 
def displacement(new_geom, x_dis=5, y_dis=0, z_dis=0):

    """ displacements for x, y, z coordinates for the whole monomer """

    cnt = 0
    for x in new_geom[:,1]:
        new_geom[cnt: cnt+1,1] = x + x_dis
        cnt += 1
    
    cnt = 0
    for y in new_geom[:,2]:
        new_geom[cnt: cnt+1,2] = y + y_dis
        cnt += 1
    
    cnt = 0
    for z in new_geom[:,3]:
        new_geom[cnt: cnt+1,3] = z + z_dis
        cnt += 1

    return new_geom

def ran_angle():
    return np.random.random_sample()*2*np.pi

def ran_dis():
    return np.random.random_sample()*15 - 7.5



final_geom = yz_rotate(geo, ran_angle())

final_geom = xy_rotate(final_geom, ran_angle())

last_geom = displacement(final_geom, ran_dis(), ran_dis(), ran_dis())
#print(last_geom)


def olig_random(geom, num):
    arching = geom[:,:]
    cnt = 0
    while (cnt < num - 1):
        yz = yz_rotate(geom, ran_angle())
        xy = xy_rotate(yz, ran_angle())
        dis = displacement(xy, ran_dis(), ran_dis(), ran_dis())
        arching = np.concatenate((arching, dis))
        cnt += 1
    return arching



ran_oligomer = olig_random(geo, 8)



out_file = "many.txt"

np.savetxt(out_file, ran_oligomer,
fmt="%s")




