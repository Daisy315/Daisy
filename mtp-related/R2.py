import os 
import numpy as np
import sys

if len(sys.argv) <=1:
    print("wrong")
    sys.exit(1)
data=np.loadtxt(sys.argv[1])
#data=np.loadtxt("E.txt")
#data[:,1]=data[:,1]-data[0,1]
#data[:,2]=data[:,2]-data[0,2]
a=int(data.shape[0])
data_1=[]
data_2=[]
for i in range(a):
    data_1.append((data[i,0]-data[0,0])/64*1000)
    data_2.append((data[i,1]-data[0,0])/64*1000)
x=data_1 #DFT，实际值
y=data_2 #MTP，预测值
print(len(y))
RSS=0
TSS=0
aver=data[:,0].mean()

for i in range(len(x)):
    RSS+=(y[i]-x[i])**2 #实际值-预测值
    TSS+=(aver-x[i])**2 #实际值-实际值平均
print("R2 is %f" %(1-RSS/TSS))
with open ("delta_E.txt","w") as writer:
    for i in range(a):
        writer.write( "%f %f \n"  %(float(data_1[i]),float(data_2[i])) )
    
