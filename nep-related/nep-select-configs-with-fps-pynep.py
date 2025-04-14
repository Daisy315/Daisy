from pynep.calculate import NEP
from pynep.select import FarthestPointSample
from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from multiprocessing import Pool
import sys

def map_fun(frame):
    return np.mean(calc.get_property('descriptor', frame), axis=0)

if __name__=='__main__':
    # read old training data, new data and descriptor
    fxyz = "."
    data_ref=read(fxyz+"/train.xyz",index=':',format='extxyz')
    data_current=read(fxyz+"/added_sphere.xyz",index=':',format='extxyz')
    calc = NEP(fxyz+"/nep.txt")

    # select data
    des_ref,des_current=[],[]
    proc_n=int(sys.argv[1])
    with Pool(processes=proc_n) as pool:
        des_ref=np.array(pool.map(map_fun,data_ref))
        des_current=np.array(pool.map(map_fun,data_current))
    sampler = FarthestPointSample(min_distance=0.01)
    selected_i = sampler.select(des_current, des_ref, min_select=100,max_select=500)
    #selected_i = sampler.select(des_current, [])
    write('selected.xyz', [data_current[i] for  i in selected_i],format='extxyz')
    np.savetxt("selected_i",selected_i)

    # color for each phase
#    color={'bcc':'orange','fcc':'purple','Pbcm':'r','X':'b'}
    
    # fit PCA model
    reducer = PCA(n_components=2)
    reducer.fit(des_current)
    # reference data
    proj_ref=reducer.transform(des_ref)
    p0 = reducer.explained_variance_ratio_[0]
    p1 = reducer.explained_variance_ratio_[1]
    print(f'Explained variance for component 0: {p0:.2f}')
    print(f'Explained variance for component 1: {p1:.2f}')
    plt.scatter(proj_ref[:,0],proj_ref[:,1], label="init data",s=1,color="silver")
    # current data
    print(des_current[0])
    print(des_current.shape)
    print(des_ref[0])
    print(des_ref.shape)
    proj_current = reducer.transform(des_current)
    plt.scatter(proj_current[:,0], proj_current[:,1],label="init data",color="orange")
    # selected data
    proj_selected = reducer.transform(np.array([des_current[i] for i in selected_i]))
#    color_list=[color[data_current[i].info['phase']] for i in selected_i]
    #plt.scatter(proj_selected[:,0], proj_selected[:,1],s=1,label="selected data",color='blue')

    np.savetxt("proj_ref",proj_ref)
    np.savetxt("proj_current",proj_current)
    np.savetxt("proj_selected",proj_selected)

    plt.legend()
    plt.axis('off')
    plt.savefig('select.png')
