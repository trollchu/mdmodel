import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mdmodel import mdmodel as mm 

###### temperature

def Tplot(file):
    anp = mm()
    anp.read_lmp_dump(file)
    df = anp.data
    df.columns
    # df[df.type<3].z.hist()

    dz = 5
    input = df[(df.type<4)&(df.z>40)&(df.z<50)][['x','z','v_Tatom']].values
    vm = input[:,2].mean()*1.2
    xm = input[:,0].mean()
    ym = input[:,1].mean()
    box = [[xm-35,xm+35],[ym-35,ym+35]]
    counterplot(input,box=box,vmax=vm)
    plt.title('t={}ps'.format(i*2))


sys='raw'
v = '0.05'
angley = 0
number_of_dump = 50

inpath = r'D:\OneDrive - bit.edu.cn\Work\2020\0728-morphology\files\rawdata\{}-{}-{}-0'.format(sys,v,angley)
outpath = r'D:\test\counter\{}-{}-{}-0\T'.format(sys,v,angley)

for i in range(1,number_of_dump+1,2):
    file = os.path.join(inpath,'trj','shock.{}0000.dump'.format(i))
    Tplot(file)
    os.makedirs(outpath,exist_ok=True)
    plt.savefig(os.path.join(outpath, 'p{}.png'.format(i)))
    plt.close()    
