import matplotlib.pyplot as plt
import numpy as np
import os

file_name_1='4K.csv'
file_name_2='50K.csv'
file_name_3='100K.csv'

number_of_nodes=40
#path_length_range=4000e-9
os.chdir('Scattering analysis')

def scattering_analysis(file_name):
    with open(file_name,"r") as f:
        scattering_maps = np.genfromtxt(f, dtype='float', delimiter=",")
    
    plt.plot (scattering_maps[:,2], scattering_maps[:,3], 'o', color='#1adc00', markersize=0.5, alpha=0.05)
    plt.plot (scattering_maps[:,4], scattering_maps[:,5], 'o', color='#ff00eb', markersize=0.5, alpha=0.05)
    plt.plot (scattering_maps[:,0], scattering_maps[:,1], 'o', color='#0065fd', markersize=0.5, alpha=0.05)   #color='#0033e7'
    plt.xlabel('X (um)', fontsize=12)
    plt.ylabel('Y (um)', fontsize=12)
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig("100K.png",dpi=900, format = 'png', bbox_inches="tight")
    #plt.savefig('testplot.png',dpi=300, format = 'pdf', bbox_inches="tight")
    plt.show()
    return

def histogram_calculation(file_name):
    with open(file_name,"r") as f:
        scattering_maps = np.genfromtxt(f, dtype='float', delimiter=",")
    
    data_diffuse=[]
    data_specular=[]
    for i in range(len(scattering_maps[:,1])):
        if scattering_maps[i,1] > 5.98e-7 and scattering_maps[i,1] < 6.02e-7:   # Select only points at the oundary
           data_diffuse.append(scattering_maps[i,0]) 
    for i in range(len(scattering_maps[:,3])):
        if scattering_maps[i,3] > 5.98e-7 and scattering_maps[i,3] < 6.02e-7:
           data_specular.append(scattering_maps[i,2]) 

    x_min=min(scattering_maps[:,0])
    x_max=max(scattering_maps[:,0])
    dist=np.zeros((number_of_nodes,3))
    delta=(x_max-x_min)/number_of_nodes
    dist[:,0]=range(number_of_nodes)
    dist[:,0]=abs(x_max-x_min)*(dist[:,0]/number_of_nodes)-abs(x_min)
    dist[:,1]=[len(filter(lambda x: x >= j and x < j+delta and x!=0, data_diffuse)) for j in dist[:,0]]    # Diffusive scattering
    dist[:,2]=[len(filter(lambda x: x >= j and x < j+delta and x!=0, data_specular)) for j in dist[:,0]]      # Specular scattering
    plt.plot (dist[:,0], dist[:,1],'b')                                         # Diffusive scattering
    plt.plot (dist[:,0], dist[:,2],'g')                                         # Specular scattering
    plt.xlabel('X (um)', fontsize=12)
    plt.ylabel('Number of phonons', fontsize=12)
    plt.show()
    
    
    
    plt.hist(data_diffuse,bins=30, color='#1d68ff')
    plt.hist(data_specular,bins=30, color='#00ac0a')
    plt.xlabel('X (um)', fontsize=12)
    plt.ylabel('Number of phonons', fontsize=12)
    plt.savefig("100K_hist.png",dpi=900, format = 'png', bbox_inches="tight")
    plt.show()
    
    return dist

scattering_analysis(file_name_3)
#scattering_analysis(file_name_2)
#scattering_analysis(file_name_3)

#histogram_calculation(file_name_3)