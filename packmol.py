#!/usr/bin/python
import math
import os
import numpy 
from subprocess import Popen,PIPE

#LixSi1
numSi=0
y = lambda x: 0.86196+1.5*math.pow(0.67145,x)
for i in numpy.linspace(0,20,80):
    density=y(i)
    numSi=int((70)/(i+1))
    numLi = max(round(i*numSi),1)
    volume=((28.0855*numSi+6.939*numLi)/(6.02214076*1E+23))/density/(1E-24)
    height=math.pow(volume,1/3.0) 
    #print(height)
    root = os.getcwd()
    name = '{0:3f}' .format(i)
    if not os.path.exists(name):
        os.mkdir(name)
    os.chdir(name)

    to_writen = ['tolerance 2.0']
    to_writen.append('filetype xyz')
    to_writen.append('output pack.xyz')
    if numLi != 0.0:
        to_writen.append('structure Li.xyz')
        to_writen.append(' number %d' %numLi)
        to_writen.append(' inside box 0. 0. 0. %f %f %f' % (height,height,height))
        to_writen.append('end structure')
    to_writen.append('')
    to_writen.append('structure Si.xyz')
    to_writen.append(' number %d' %numSi)
    to_writen.append(' inside box 0. 0. 0. %f %f %f' % (height,height,height))
    to_writen.append('end structure')
    with open('pack.inp','w') as writer:
        writer.writelines([j+'\n' for j in to_writen])
    with open('Si.xyz','w') as writer:
        writer.write('1\n\nSi  0 0 0\n')
    with open('Li.xyz','w') as writer:
        writer.write('1\n\nLi  0 0 0\n')
    p2 = Popen('/usr/bin/packmol < pack.inp',shell=True,stdout=PIPE)
    p2.communicate() #now wait
    with open('pack.xyz') as reader:
        number = int(reader.readline())
        reader.readline()
        coord = []
        for k in range(number):
            coord.append(' '.join(reader.readline().split()[1:]))

    with open('POSCAR','w') as writer:
        writer.write('nxu\n1.0\n')
        writer.write('%f 0.0 0.0\n' %(height+1))
        writer.write('0.0 %f 0.0\n' %(height+1))
        writer.write('0.0 0.0 %f\n' %(height+1))
        if numLi != 0.0:
            writer.write('Li Si\n')
            writer.write('%d %d\nCartesian\n' %(numLi,numSi))
        else:
            writer.write('Li Si\n')
            writer.write('0 %d\nCartesian\n' %(numSi))       
        for k in coord:
            writer.write(k+'\n')
   
    os.chdir(root)   
