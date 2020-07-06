import numpy as np
import os
import math
import random
from numpy import genfromtxt
import numpy as npimport 
from numpy import genfromtxt
import pandas as pd




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


    theta_fin = np.zeros(int(len(geom[:,0])))
    k = 1
    
    for i in range(len(theta_0)):
        if math.isnan(theta_0[i]) == True:
            theta_0[i] = 0


    theta_fin = np.zeros(int(len(geom[:,0])))
    k = 1
    for theta in theta_0[:]:

        theta_fin[k-1:k] = theta + math.degrees(yz_angle)
        k += 1




    
    xs = geom[:,1]
    ys = np.zeros(int(len(geom[:,0])))
    zs = np.zeros(int(len(geom[:,0])))
    k = 1
    for theta, hypot in np.nditer([theta_fin, hyp]):
        ys[k-1:k] = round(hypot*math.sin(math.radians(theta)), 8)
        zs[k-1:k] = round(hypot*math.cos(math.radians(theta)), 8)
        k += 1

    

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

    return new_geom

def xy_rotate(geom, xy_angle=np.pi/4):
    """ First, must initialize a np.array to store new geometry based on input geometry's size """
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
        """ Accounts for each quadrant that the inverse tangent function produces """


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

        

    theta_fin = np.zeros(int(len(geom[:,0])))
    k = 1
    
    for i in range(len(theta_0)):
        if math.isnan(theta_0[i]) == True:
            theta_0[i] = 0
    

    theta_fin = np.zeros(int(len(geom[:,0])))
    k = 1
    for theta in theta_0[:]:
        
        theta_fin[k-1:k] = theta + math.degrees(xy_angle)
        k += 1
    



    
    zs = geom[:,3]
    ys = np.zeros(int(len(geom[:,0])))
    xs = np.zeros(int(len(geom[:,0])))
    k = 1
    for theta, hypot in np.nditer([theta_fin, hyp]):
        ys[k-1:k] = round(hypot*math.sin(math.radians(theta)), 8)
        xs[k-1:k] = round(hypot*math.cos(math.radians(theta)), 8)
        k += 1
    
    """ Saves the rotated x and y values to new array"""

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
    """ Randomly selects a rotation from 0 to 2pi """
    return np.random.random_sample()*2*np.pi

def ran_dis():
    """ Change the coefficient and subtracted constant to build the rectangular box of your choosing """

    return np.random.random_sample()*15 - 7.5  ################change for rectangular box size


def distance(geom1, geom2):
    """ Three dimensional distance formula for evaluating the distance betweeen molecules """
    return math.sqrt( (geom1[1] - geom2[1]) **2 + (geom1[2] - geom2[2])**2 + ( geom1[3] - geom2[3])**2)




def random_arrangement(geom, num):
    arching = geom[:,:]
    cnt = 1
    molecule = [0]
    

    spacer = len(geom[:,0])

    while ( cnt < num ):
        yz = yz_rotate(geom, ran_angle())
        xy = xy_rotate(yz, ran_angle())
       
       



        dis = displacement(xy, ran_dis(), ran_dis(), ran_dis())
        check_tf = False

        for i in range(len(molecule)):
            dist_CC = distance(dis[0,:], arching[0 + i*spacer,:])

            print(dist_CC)
            print((dis[0,:]))
            print(arching[0 + i*spacer,:])
            
            if dist_CC < 4.0:  # change for the minimum distance between molecules
                check_tf = True

        if check_tf == False:
            molecule.append(cnt)
            arching = np.concatenate((arching, dis))
            cnt += 1

    arching = np.round_(arching, decimals=16)


    return arching, len(molecule), spacer





def clean_many_txt():

    """ This will replace the numerical forms of the elements as their letters numbered in order """

    f = open('many.txt','r')
    a = ['6.0 ', '8.0 ', '1.0 ']
    table = {
        '6.0 ': 'C', '8.0 ' : 'O', '1.0 ': 'H'
    }

    lst = []
    cnt2 = 0
    for line in f:
        cnt2 += 1
        for word in a: 
            if word in line:
                convert_wrd = table[word]
                line = line.replace(word, convert_wrd + str(cnt2) + " ")
                
        lst.append(line)
    f.close()
    f = open('many.txt','w')
    for line in lst:
        f.write(line)
    f.close()



def constraints(molecule_cnt, spacer):
    """ Must be called after random_arrangement() since it takes the value of len(molecule) from the function's output. Takes the bonds.txt and angles.txt of the monomer's geometry  """

    bsas = pd.read_csv('bonds.txt', sep=' ', header=None) 
    df_bds = bsas.replace(np.nan, ' ', regex=True)

    bsas = pd.read_csv('angles.txt', sep=' ', header=None)
    df_ang = bsas.replace(np.nan, ' ', regex=True)





    for row in df_bds.itertuples():
        if df_bds.columns[0]==int:
            df_bds.at[row.Index, 0] += 6
        if df_bds.columns[1]==int:
            df_bds.at[row.Index, 1] += 6

        df_bds.loc[row.Index, 0] += 6
        df_bds.loc[row.Index, 1] += 6



    for row in df_ang.itertuples():
        if df_ang.columns[0]==int:
            df_ang.at[row.Index, 0] += 6
        if df_ang.columns[1]==int:
            df_ang.at[row.Index, 1] += 6
        if df_ang.columns[2]==int:
            df_ang.at[row.Index, 2] += 6

        df_ang.loc[row.Index, 0] += 6
        df_ang.loc[row.Index, 1] += 6
        df_ang.loc[row.Index, 2] += 6


    num = molecule_cnt - 2

    fnames = [df_ang, df_bds]
    df = pd.concat(fnames, ignore_index=True)

    for i in range(num):

        for row in df_bds.itertuples():
            if df_bds.columns[0]==int:
                df_bds.at[row.Index, 0] += spacer
            if df_bds.columns[1]==int:
                df_bds.at[row.Index, 1] += spacer

            df_bds.loc[row.Index, 0] += spacer
            df_bds.loc[row.Index, 1] += spacer





        for row in df_ang.itertuples():
            if df_ang.columns[0]==int:
                df_ang.at[row.Index, 0] += spacer
            if df_ang.columns[1]==int:
                df_ang.at[row.Index, 1] += spacer
            if df_ang.columns[2]==int:
                df_ang.at[row.Index, 2] += spacer

            df_ang.loc[row.Index, 0] += spacer
            df_ang.loc[row.Index, 1] += spacer
            df_ang.loc[row.Index, 2] += spacer

        df = pd.concat([df, df_ang, df_bds], ignore_index=True)

    return df




def clean_dataframe(df):
    
    """ This cleans the output of the dataframe to remove blanks """
    df.to_csv('dataframe_test.csv', index=False, sep=" ")

    f = open('dataframe_test.csv','r')
    a = ['" "']
    lst = []
    for line in f:
        for word in a:
            if word in line:
                line = line.replace(word,'')
        lst.append(line)
    f.close()
    f = open('dataframe_test.csv','w')
    for line in lst:
        f.write(line)
    f.close()


def make_input_dir(dir_name_number):
    """ Combines the geometry output and the constrained output. Then makes the .com and .pbs files in a subdirectory """

    data = data2 = "" 
    
    with open('many.txt') as fp: 
        data = fp.read() 
    
    # Reading data from file2 
    with open('dataframe_test.csv') as fp: 
        data2 = fp.read() 
    

    data += "\n\n"
    data += data2 
    charges = "0 1"

    new_dir = "geom" + str(dir_name_number)
    os.mkdir(new_dir)
    with open (new_dir + '/mex.com', 'w') as fp: 
        fp.write("%mem=1600mb\n")
        fp.write("%nprocs=4\n")
        fp.write("#N wB97XD/6-31G(d) opt=ModRedundant FREQ\n")
        fp.write("\n")
        fp.write("Name ModRedundant\n")
        fp.write("\n")
        fp.write(charges + "\n")
        fp.write(data) 

    with open (new_dir + '/mex.pbs', 'w') as fp: 
        fp.write("#!/bin/sh\n")
        fp.write("#PBS -N mex\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m abe\n#PBS -l")
        fp.write("mem=15gb\n")
        fp.write("#PBS -l nodes=1:ppn=4\n#PBS -q gpu\n\nscrdir=/tmp/$USER.$PBS_JOBID\n\n")
        fp.write("mkdir -p $scrdir\nexport GAUSS_SCRDIR=$scrdir\nexport OMP_NUM_THREADS=1\n\n")
        fp.write("""echo "exec_host = $HOSTNAME"\n\nif [[ $HOSTNAME =~ cn([0-9]{3}) ]];\n""")
        fp.write("then\n")
        fp.write("  nodenum=${BASH_REMATCH[1]};\n  nodenum=$((10#$nodenum));\n  echo $nodenum\n\n")
        fp.write("""  if (( $nodenum <= 29 ))\n  then\n    echo "Using AVX version";\n""")
        fp.write("    export g16root=/usr/local/apps/gaussian/g16-b01-avx/\n  elif (( $nodenum > 29 ))\n")
        fp.write("""  then\n    echo "Using AVX2 version";\n    export g16root=/usr/local/apps/gaussian/g16-b01-avx2/\n  else\n""")
        fp.write("""    echo "Unexpected condition!"\n    exit 1;\n  fi\nelse\n""")
        fp.write("""  echo "Not on a compute node!"\n  exit 1;\nfi\n\n""")
        fp.write("cd $PBS_O_WORKDIR\n. $g16root/g16/bsd/g16.profile\ng16 mex.com mex.out\n\nrm -r $scrdir\n")


def make_input_files():
    """ Combines the geometry output and the constrained output. Then makes the .com and .pbs files in a subdirectory """

    data = data2 = "" 
    
    with open('many.txt') as fp: 
        data = fp.read() 
    
    # Reading data from file2 
    with open('dataframe_test.csv') as fp: 
        data2 = fp.read() 
    

    data += "\n\n"
    data += data2 
    charges = "0 1"

    
    with open ('mex.com', 'w') as fp: 
        fp.write("%mem=1600mb\n")
        fp.write("%nprocs=4\n")
        fp.write("#N wB97XD/6-31G(d) opt=ModRedundant FREQ\n")
        fp.write("\n")
        fp.write("CH2O3 ModRedundant - Minimalist working constrained optimisation\n")
        fp.write("\n")
        fp.write(charges + "\n")
        fp.write(data) 

    with open ('mex.pbs', 'w') as fp: 
        fp.write("#!/bin/sh\n")
        fp.write("#PBS -N mex\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -m abe\n#PBS -l")
        fp.write("mem=15gb\n")
        fp.write("#PBS -l nodes=1:ppn=4\n#PBS -q gpu\n\nscrdir=/tmp/$USER.$PBS_JOBID\n\n")
        fp.write("mkdir -p $scrdir\nexport GAUSS_SCRDIR=$scrdir\nexport OMP_NUM_THREADS=1\n\n")
        fp.write("""echo "exec_host = $HOSTNAME"\n\nif [[ $HOSTNAME =~ cn([0-9]{3}) ]];\n""")
        fp.write("then\n")
        fp.write("  nodenum=${BASH_REMATCH[1]};\n  nodenum=$((10#$nodenum));\n  echo $nodenum\n\n")
        fp.write("""  if (( $nodenum <= 29 ))\n  then\n    echo "Using AVX version";\n""")
        fp.write("    export g16root=/usr/local/apps/gaussian/g16-b01-avx/\n  elif (( $nodenum > 29 ))\n")
        fp.write("""  then\n    echo "Using AVX2 version";\n    export g16root=/usr/local/apps/gaussian/g16-b01-avx2/\n  else\n""")
        fp.write("""    echo "Unexpected condition!"\n    exit 1;\n  fi\nelse\n""")
        fp.write("""  echo "Not on a compute node!"\n  exit 1;\nfi\n\n""")
        fp.write("cd $PBS_O_WORKDIR\n. $g16root/g16/bsd/g16.profile\ng16 mex.com mex.out\n\nrm -r $scrdir\n")


for i in range(1, 6, 1):

    geo = genfromtxt('mon_geo.xyz', delimiter=' ')

    """ Takes array and saves it to file """

    final, mole, spacer = random_arrangement(geo, 10)

    out_file = "many.txt"

    np.savetxt(out_file, final,
    fmt="%s")
    """ end """

    clean_many_txt()

    df = constraints(mole, spacer)

    clean_dataframe(df)

    #make_input_files()
    make_input_dir(i)

