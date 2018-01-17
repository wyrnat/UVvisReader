'''
Created on 14.09.2016

@author: jannik & matthias
'''
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

class dataArray(list):
     
    def __init__(self, filenumber):
        super(dataArray, self).__init__()
        self.filenumber = filenumber
     
    def boolZeroMeasure(self):
        mysum = 0
        for i in self:
            mysum +=  i
        meanval = (len(self)>0)*mysum/float(len(self))
        return meanval < 10.0
     
    def deleteEntries(self, minpos, maxpos):
        iterations = (len(self)-1) - maxpos
        for i in range(iterations):
            self.remove(-1)
        for i in range(minpos):
            self.remove(0)
            
    def getNumber(self):
        return self.filenumber
             


""" Smoothing"""
##################
def smooth(y, box_pts):
    """
    Smoothes data with numpy.convole package
    @param y: data set
    @param box_pts: number of measure points in both directions to smooth the value
    @return: y_smooth: smoothed data set
    """
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    return y_smooth

def smooth_files(dataSet, wavelength, smoo):
    """
    convoles a wavelength and absorbance
    @param dataSet: Absorbance dataset
    @param wavelength: measured wavelength data
    @param smoo: number of box points for smoothing
    @return: dataSetneo, wavelengthneo: smoothed datasets
    """
    dataSetneo = []
    wavelengthneo = []
    
    wavelengthneo = smooth(wavelength, smoo)[:]
    for sdata in dataSet:
        dataSetneo.append(smooth(sdata, smoo))
        
    return dataSetneo, wavelengthneo


""" Extraction """
##################    

def extractFiles(name):
    waveLength = dataArray()
    dataBib = {'TR': dataArray()}
    
    
    
    
     
    


def convertFiles(name, start, stop, path, header=88, lineN=2137, Selected = 1):
    """
    Loads the raw UVvis files and extracts the data.
    Models the data structure to import it with IGOR, parallel to modeling it saves
    the wavelength as well as the absorbance values in an array
    @param name: file name without indice or filetype
    @param start: first file indice to be extracted
    @param stop: last file indice to be extracted
    @param path
    """
     
    myArray = []
    waveLength = []
    filenumber = []
    selectCounter = 0
    startselection = False
    for i in range(int(start),int(stop)+1):
        if (selectCounter % Selected != 0) and (startselection == True):
            continue
        partArray = dataArray(i)
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
        
        if partArray.boolZeroMeasure() == False:
            filenumber.append(i)
            myArray.append(partArray)
            startselection = True
            selectCounter+=1
        else:
            startselection = False
         
        newfile.close()
        myfile.close()
    return myArray, waveLength, filenumber

def origin_safefile(dataSetneo, wavelengthneo, filenumber, path):
    """
    Creates a .txt File with the header wavelength and the extracted file indices.
    In the wavelength column there are all extracted wavelength.
    Under each extracted file indice there is the measured absorbance.
    @param dataSetneo: smoothed dataset of absorbance values for for every extracted file
    @param wavelengthneo: smoothed array of etracted wavelength
    @param filenumber the file indices of the extracted files
    """
    originfile = open(path+'/Spectra.txt', 'w')
    text = 'wavelength'

                
    for k in filenumber:
        text += ' '+str(k)
    text += '\n'
    for j in range(len(wavelengthneo)):
        text += str(wavelengthneo[j])
        for i in range(len(dataSetneo)):
            text +=' '+str(dataSetneo[i][j])
        text += '\n'
        
    originfile.write(text)
    originfile.close()


""" Visualisation """
#####################
def plot2d(dataSet, name):
    fig, ax = plt.subplots()

    cax = ax.imshow(dataSet)
    ax.set_title(name)
    ax.set_ylabel('wavelength (nm)')
    ax.set_xlabel('file number')
    
    cbar = fig.colorbar(cax)#, ticks=[0, 100, 200])
    #cbar.ax.set_yticklabels(['0', '100', '200'])

    plt.show()
    
def multipleLines2dPlot(dataSet, waveLength, pos, name, filenumber):
    for i in pos:
        print "dataset "+str(dataSet[i])
        plt.plot(waveLength, dataSet[i], color=np.random.rand(3,1), label=filenumber[i])
    plt.xlabel('wavelength (nm)')
    plt.ylabel('absorbance (a.u.)')
    plt.title(name)
    plt.legend()
    plt.show()

def plot3d(dataSet, waveLength, Name):
    X = np.arange(1,len(dataSet)+1)
    Y = np.array(waveLength)
    Y, X = np.meshgrid(Y, X)
    Z = dataSet
    fig = plt.figure()
    fig.canvas.set_window_title('UVvis Reader 3D Plot')
    ax = Axes3D(fig)
    ax.set_ylabel('wavelength (nm)')
    ax.set_xlabel('file number')
    ax.set_zlabel('absorbance (a.u.)')
    ax.set_title(Name)
    ax.plot_surface(X, Y, Z, rstride = 5, cstride = 5, cmap=plt.get_cmap('nipy_spectral'), linewidth=0, antialiased=True)
    plt.show()


""" Command Line """

def convertCommand(inputcommand):
    if type(inputcommand) != str:
        return []
    return inputcommand.split(' ')

def handleCommands(inputArray, dataSet, wavelength, filenumber, dataSetsmooth, ):
    if inputArray == 0:
        return False
    elif inputArray[0] == 'help':
        helpCommand()
    elif inputArray[0] == 'exit':
        return True
    elif inputArray[0] == 'plot':
        if inputArray[1] == '2d':
            plot2d(dataSetsmooth, name)
        elif inputArray[1] == '3d':
            plot3d(dataSetsmooth, wavelengthsmooth, name)
        elif inputArray[1] == 'n':
            multipleLines2dPlot(dataSetsmooth, wavelengthsmooth, range(len(dataSet)), name, filenumber)
        else:
            print "No valid plot object. [2d/3d/n]"
    elif inputArray[0] == 'cutoff':
        if inputArray[1] == 'wavelength':
            cutoffWavelengthCommand(inputArray[2])
        elif inputArray[1] == 'files':
            cutoffFilenumbersCommand(inputArray[2])
        else:
            print "No valid cutoff object. [files/wavelength]"
    elif inputArray[0] == 'selector':
        selectorCommand(inputArray[1])
    
    
    else:
        print "No valid command. Try help for list of commands" 
            
                


def helpCommand():
    print "UVvisReader commands"
    print "help -- this list"
    print "exit -- exit UVvisReader"
    print "plot [2d/3d/n] -- plot different figures. Example: UVvisReader>>plot 3d"
    print "cutoff [wavelength/files] [True/False/ [start,stop]]"
    
def cutoffWavelengthCommand(mydata):
    pass

def cutoffFilenumbersCommand(mydata):
    pass

def selectorCommand(selectvalue):
    pass
    
    
    
    
    
##################################
#########     CODE     ###########
##################################  
name = raw_input("File name without number postfix: ")
start = int(raw_input("First sequence number: "))
stop = int(raw_input("Last sequence number: "))
smoo = int(raw_input("Smoothing lenght: "))
selected = int(raw_input("Use only every ... file: "))

path = os.getcwd()
if not os.path.exists(path+'/IGOR/'):
    os.mkdir(path+'/IGOR/')


dataSet, wavelength, filenumber = convertFiles(name, start, stop, path, Selected = selected)

dataSetsmooth, wavelengthsmooth = smooth_files(dataSet, wavelength, smoo)

origin_safefile(dataSetsmooth, wavelengthsmooth, filenumber, path)

myexit = False
while myexit != True:
    commandLine = raw_input('UVvisReader >>')
    myexit = handleCommands(convertCommand(commandLine),
                            dataSet = dataSet, wavelength = wavelength, filenumber = filenumber,
                            dataSetsmooth = dataSetsmooth, wavelengthsmooth=wavelengthsmooth,
                            name=name, start=start, stop=stop, smoo=smoo, selected=selected
                            )
    
    
    
# plotBool = raw_input("2d or 3d Colorplot or normal plot? [2d/3d/n]: ")
# if plotBool == "2d":
#     plot2d(dataSetsmooth, name)
# elif plotBool == "3d":
#     plot3d(dataSetsmooth, wavelengthsmooth, name)
# elif plotBool == "n":
#     pos = range(len(dataSet))
#     multipleLines2dPlot(dataSetsmooth, wavelengthneo, pos, name, filenumber)
# else:
#     pass
# visualisation(np.clip(dataSetneo, 0, 200), name)



