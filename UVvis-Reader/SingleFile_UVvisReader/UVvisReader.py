'''
Created on 14.09.2016

@author: jannik & matthias
'''
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


path = os.getcwd()
if not os.path.exists(path+'/IGOR/'):
    os.mkdir(path+'/IGOR/')


""" Smoothing"""
##################
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    return y_smooth


""" Extraction """
##################    
def convertFiles(name, start, stop, path, header=88, lineN=2137, Selector = 1):
    if i % Selector != 0:
        return
     
    myArray = []
    waveLength = []
    for i in range(int(start),int(stop)+1):
        partArray = []
        prozero = (i==0)*('0')+(i<10)*('0')+(i<100)*('0')+(i<1000)*('0')
        myfile = open(path+'/'+name+prozero+str(i)+".txt",'r')
        newfile = open(path+'/IGOR/'+name+"_"+prozero+str(i)+"_NEW.txt", 'w')

        line = "marc"
        j = 0
        while j < lineN:
            line = myfile.readline()
            if (j > header):
                newline = line.split(';')
                for i, element in enumerate(newline):
                    element2 = element.replace(',','.')
                    newline[i] = element2.replace('   ','0')
                                  
                writtenline =''
                for element in newline:  
                    writtenline += (element+' ')
                newfile.write(writtenline)
                if newline[1] == '0':
                  pass
                else:
                  if len(waveLength)==0 or (float(newline[1]) > waveLength[-1]):
                      waveLength.append(float(newline[1]))
                  partArray.append(float(newline[8]))
                
            j+=1
        
        myArray.append(partArray)
         
        newfile.close()
        myfile.close()
    return myArray, waveLength

""" Visualisation """
#####################
def visualisation(dataSet, name):
    fig, ax = plt.subplots()

    cax = ax.imshow(dataSet)
    ax.set_title(name)
    
    cbar = fig.colorbar(cax)#, ticks=[0, 100, 200])
    #cbar.ax.set_yticklabels(['0', '100', '200'])

    plt.show()





##################################
#########     CODE     ###########
##################################  
name = raw_input("File name without number postfix: ")
start = int(raw_input("First sequence number: "))
stop = int(raw_input("Last sequence number: "))
smoo = int(raw_input("Smoothing lenght: "))
selected = int(raw_input("Use only every ... file: "))


dataSet, wavelength = convertFiles(name, start, stop, path, Selected = selected)

print len(dataSet)
print len(dataSet[0])
print len(wavelength)
#print wavelength

originfile = open(path+'/Spectra.txt', 'w')
text = 'wavelength'

dataSetneo = []
wavelengthneo = []
wavelengthneo = smooth(wavelength, smoo)[:]

for sdata in dataSet:
    dataSetneo.append(smooth(sdata, smoo))

print len(dataSetneo[0])
print len(wavelengthneo)
                     
for k in range(start, stop):
    text += ' '+str(k)
text += '\n'
for j in range(len(wavelengthneo)):
    text += str(wavelengthneo[j])
    for i in range(stop-start):
      text +=' '+str(dataSetneo[i][j])
    text += '\n'
originfile.write(text)
originfile.close()

visualisation(dataSetneo, name)
# visualisation(np.clip(dataSetneo, 0, 200), name)



