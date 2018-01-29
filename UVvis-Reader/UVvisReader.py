'''
Created on 14.09.2016

@author: Jannik & Matthias Schwartzkopf
@version: 0.2 beta
'''
import os

try:
    import matplotlib.pyplot as plt
except:
    print 'matplotlib not installed.'
    raw_input('Close UVvisReader...')
try:
    import numpy as np
except:
    print 'numpy not installed.'
    raw_input('Close UVvisReader...')  

from mpl_toolkits.mplot3d import Axes3D

class dataArray(list):
    """
    specialised list object for data handling.
    contains the file number
    """
    
    def __init__(self, liste=[]):
        super(dataArray, self).__init__(liste)
        self.filenumber = None
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
    
    def setNumber(self, nr):
        self.filenumber = nr
             


""" Extraction """
#########################################################################

def extractFileName(fileNameString):
    """
    Splits the input string to extract
    - the preceding path
    - the clear file name
    - the file number
    - the data type
    """
    pathArray = fileNameString.rsplit('/', 1)
    nameAndType = pathArray[-1].rsplit('.', 1)
    nameList = list(nameAndType[0])
    revnameList = nameList[:]
    revnameList.reverse()
    numberlen = 0
    for chr in revnameList:
        try:
            int(chr)
            numberlen -= 1
        except:
            break
           
    while True:
        if nameList[numberlen] == '0':
            break
        else:
            numberlen +=1
            
    filenumber = ''.join(nameList[numberlen:])
    filename = ''.join(nameList[:numberlen])
    filetype = nameAndType[-1]
    if len(pathArray) == 1:
        additionalPath = ''
    else:
        additionalPath = '/'+pathArray[0]
    
    return additionalPath, filename, filenumber, filetype



def extractFiles(fileNameString, path, header=88, Cwl = 1, Ctr = 8):
    waveLength = []
    dataBib = {'TR': []}
    
    additionalPath, filename, filenumber, filetype = extractFileName(fileNameString)
    pzf = len(filenumber)   #prozero factor
    
    k = 1
    while True:
        partArray = dataArray()
        partArray.setNumber(k)
        prozero = ((k<10)*(pzf>1)*('0')+
                   (k<100)*(pzf>2)*('0')+
                   (k<1000)*(pzf>3)*('0')+
                   (k<10000)*(pzf>4)*('0')+
                   (k<100000)*(pzf>5)*('0')
                   )
        if not os.path.isfile(path+additionalPath+'/'+filename+prozero+str(k)+"."+filetype):
            print 'INFO: Finished Extraction at:'
            print path+additionalPath+'/'+filename+prozero+str(k)+"."+filetype
            break
        myfile = open(path+additionalPath+'/'+filename+prozero+str(k)+"."+filetype,'r')
        
        j=0
        while True:
            line = myfile.readline()
            if (j > header):
                if line == '':
                    break
                newline = line.split(';')
                for i, element in enumerate(newline):
                    element2 = element.replace(',','.')
                    newline[i] = element2.replace('   ','0')
                if newline[Cwl] != '0':
                    if k == 1:
                        waveLength.append(float(newline[Cwl]))
                    partArray.append(float(newline[Ctr]))
    
            j+=1
        k+=1
        
        if partArray.boolZeroMeasure() == False:
            dataBib['TR'].append(partArray)
        
        myfile.close()
    return dataBib['TR'], waveLength, filename



""" Handle Data"""
#######################################################################

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
    
    dataSetsmooth = []
    wavelengthsmooth = smooth(wavelength, smoo)[:]
    for sdata in dataSet:
        newdA = dataArray(smooth(sdata, smoo))
        newdA.setNumber(sdata.getNumber())
        dataSetsmooth.append(newdA)
        
    return dataSetsmooth, wavelengthsmooth

def prepareDataSet(dataSet, wavelength, dataSetProperties):
    """
    Cut off files, wavelength or skip files.
    This function provides the choosen dataSet for plotting and saving
    @param dataSet: list of dataArray objects of the measure
    @param wavelength: list of loaded wavelength
    @param dataSetProperties: dictionary with dataSet-changing values
    @return: prepared dataSet and wavelength with values from dataSetProperties
    """

    wlMin = dataSetProperties['wl'][0]
    wlMax = dataSetProperties['wl'][1]
    
    for i, wl in enumerate(wavelength):
        if wl >= wlMin:
            wlMinIndex = i
            break
    for i, wl in enumerate(wavelength[::-1]):
        if wl <= wlMax:
            wlMaxIndex = len(wavelength)-i-1
            break
        
    try:
        wlMinIndex
        wlMaxIndex
    except:
        wlMinIndex = 0
        wlMaxIndex = -1
        print 'ERROR: could not set wavelength cutoff. wavelength will not be cut off.'
    
    fileMin = dataSetProperties['file'][0]
    fileMax = dataSetProperties['file'][1]
    for i, dArray in enumerate(dataSet):
        if dArray.getNumber() >= fileMin:
            fileMinIndex = i
            break
    for i, dArray in enumerate(dataSet[::-1]):
        if dArray.getNumber() <= fileMax:
            fileMaxIndex = len(dataSet)-i-1
            break
        
    try:
        fileMinIndex
        fileMaxIndex
    except:
        fileMinIndex = 0
        fileMaxIndex = -1
        print 'ERROR: could not set file cutoff to dataSet. dataSet will not be cutted of.'
    
    slct = dataSetProperties['selector']
    cutFileArray = dataSet[fileMinIndex:fileMaxIndex:slct]
    preparedArray = []
    for dArray in cutFileArray:
        newArray = dataArray(dArray[wlMinIndex:wlMaxIndex])
        newArray.setNumber(dArray.getNumber())
        preparedArray.append(newArray)
        
    preparedwavelength = wavelength[wlMinIndex:wlMaxIndex]
    return preparedArray, preparedwavelength






""" Safe Files """
####################

def origin_safefile(dataSet, wavelength, dataSetProperties, path):
    """
    Creates a .txt File with the header wavelength and the extracted file indices.
    In the wavelength column there are all extracted wavelength.
    Under each extracted file indice there is the measured absorbance.
    @param dataSet: dataset of absorbance values for for every extracted file
    @param wavelength: array of extracted wavelength
    @param dataSetProperties: dictionary with dataSet-changing values
    @param path: path to the current working directory
    """
    if not os.path.exists(path+'/Origin/'):
        os.mkdir(path+'/Origin/')
    originfile = open(path+'/Origin/Spectra.txt', 'w')
    text = 'wavelength'
    safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, dataSetProperties)
                
    for dArray in safedataSet:
        text += ' '+str(dArray.getNumber())
    text += '\n'
    for j in range(len(safewavelength)):
        text += str(safewavelength[j])
        for dArray in safedataSet:
            text +=' '+str(dArray[j])
        text += '\n'
        
    originfile.write(text)
    originfile.close()
    print "INFO: saved to ORIGIN/spectra.txt"
    print "      wavelength: ("+str(dataSetProperties['wl'][0])+','+str(dataSetProperties['wl'][1])+')'
    print "      files: ("+str(dataSetProperties['file'][0])+','+str(dataSetProperties['file'][1])+')'
    
def Igor_safefile(dataSet, wavelength, dataSetProperties, name, path):
    """
    Safes data to IGOR compatible files.
    Each file contains wavelength on first column and TR on second column.
    It considers the file limits, they can be set with
    >> setlimits file [min, max]
    """
    
    if not os.path.exists(path+'/IGOR/'):
        os.mkdir(path+'/IGOR/')
        
    safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, dataSetProperties)
        
    pzf = int(np.ceil(np.log10(dataSetProperties['file'][1])))+1
    for dArray in safedataSet:
        nr = dArray.getNumber()
        prozero = ((nr<10)*(pzf>1)*('0')+
                    (nr<100)*(pzf>2)*('0')+
                    (nr<1000)*(pzf>3)*('0')+
                    (nr<10000)*(pzf>4)*('0')+
                    (nr<100000)*(pzf>5)*('0')
                    )
        newfile = open(path+'/IGOR/'+name+"_"+prozero+str(nr)+"_IGOR.txt", 'w')
        newfile.write('wavelength TR\n')
        for j, wl in enumerate(safewavelength):
            newfile.write(str(wl)+' '+str(dArray[j])+'\n')
        newfile.close()
    
    print "INFO: saved to IGOR/NAME_IGOR.txt"
    print "      wavelength: ("+str(dataSetProperties['wl'][0])+','+str(dataSetProperties['wl'][1])+')'
    print "      files: ("+str(dataSetProperties['file'][0])+','+str(dataSetProperties['file'][1])+')'


""" Visualisation """
#####################
def plot2d(dataSet, wavelength, plotProperties, dataSetProperties):
    
    safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, dataSetProperties)
    
    fig, ax = plt.subplots()

    cax = ax.imshow(safedataSet)
    ax.set_title(plotProperties['name'])
    ax.set_xlabel(plotProperties['2d_axinfo'][0])
    ax.set_ylabel(plotProperties['2d_axinfo'][1])
    
    cbar = fig.colorbar(cax)
    cbar.set_label(plotProperties['2d_axinfo'][2])
    fig.canvas.set_window_title('UVvis Reader 2D ColorPlot')

    plt.show()
    
def multipleLines2dPlot(dataSet, wavelength, plotProperties, dataSetProperties):
    safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, dataSetProperties)
    for dArray in safedataSet:
        plt.plot(safewavelength, dArray, color=np.random.rand(3,1), label=dArray.getNumber())
    plt.xlabel(plotProperties['n_axinfo'][0])
    plt.ylabel(plotProperties['n_axinfo'][1])
    plt.title(plotProperties['name'])
    plt.legend()
    plt.show()

def plot3d(dataSet, wavelength, plotProperties, dataSetProperties):
    safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, dataSetProperties)
    X = []
    for dArray in safedataSet:
        X.append(dArray.getNumber())
    X = np.array(X)
    Y = np.array(safewavelength)
    Y, X = np.meshgrid(Y, X)
    Z = list(safedataSet[:])
    fig = plt.figure()
    fig.canvas.set_window_title('UVvis Reader 3D Plot')
    ax = Axes3D(fig)
    ax.set_xlabel(plotProperties['3d_axinfo'][0])
    ax.set_ylabel(plotProperties['3d_axinfo'][1])
    ax.set_zlabel(plotProperties['3d_axinfo'][2])
    ax.set_title(plotProperties['name'])
    ax.plot_surface(X, Y, Z, rstride = 5, cstride = 5, cmap=plt.get_cmap('nipy_spectral'), linewidth=0, antialiased=True)
    plt.show()


""" Command Line """
######################################################

def convertCommand(inputcommand):
    if type(inputcommand) != str:
        return ['']
    stripinput = inputcommand.strip()
    return stripinput.split(' ')

def getValuesFromStringArray(strArray):
    elements = strArray[1:-1].split(',')
    for i in elements:
        i = i.strip()
    return elements

def handleCommands(inputArray, dataSet_raw, wavelength_raw, dataSet, wavelength,
                   dataSetProperties_raw, dataSetProperties, plotProperties, path):
    
    
    if inputArray == ['']:
        return False
    
    
    elif inputArray[0] == 'help':
        helpCommand()
        
        
    elif inputArray[0] == 'exit':
        return True
    
    
    elif inputArray[0] == 'plot':
        if dataSet== []:
            print "ERROR: no data. Try 'load file.txt'"
            return False
        elif len(inputArray) == 1:
            print "ERROR: No valid plot object. [2d/3d/n]"
            return False
        elif inputArray[1] == '2d':
            plot2d(dataSet, wavelength, plotProperties, dataSetProperties)
        elif inputArray[1] == '3d':
            plot3d(dataSet, wavelength, plotProperties, dataSetProperties)
        elif inputArray[1] == 'n':
            multipleLines2dPlot(dataSet, wavelength, plotProperties, dataSetProperties)
        else:
            print "ERROR: No valid plot object. [2d/3d/n]"
            
            
    elif inputArray[0] == 'cutoff':
        if dataSet== []:
            print "ERROR: no data. Try 'load file.txt'"
            return False
        if len(inputArray) == 1:
            print "ERROR: No valid cutoff object. [files/wavelength]"
            return False 
        if inputArray[1] == 'wavelength':
            cutoffWavelengthCommand(inputArray[2], dataSet_raw, dataSet, dataSetProperties_raw, dataSetProperties)
        elif inputArray[1] == 'files':
            cutoffFilenumbersCommand(inputArray[2], dataSet_raw, dataSet, dataSetProperties_raw, dataSetProperties)
        else:
            print "ERROR: No valid cutoff object. [files/wavelength]"
        
        
    elif inputArray[0] == 'load':
        if len(inputArray) == 1:
            print 'ERROR: please enter a file name. Example: load test001.txt'
            return False
        ds, wl, filename = extractFiles(inputArray[1], path)
        if ds == [] or wl == []:
            print 'ERROR: data could not be loaded. Got path and filename correct?'
            return False
        #copy raw data to working data
        dataSet[:] = ds[:]
        wavelength[:] = wl[:]
        dataSet_raw[:] = ds[:]
        wavelength_raw[:] = wl[:]
        #set file name for plot properties
        plotProperties['name'] = filename
        #set dataSet properties raw  and working file
        dataSetProperties_raw['wl'][0] = wavelength_raw[0]
        dataSetProperties_raw['wl'][1] = wavelength_raw[-1]
        dataSetProperties['wl'][0] = wavelength_raw[0]
        dataSetProperties['wl'][1] = wavelength_raw[-1]
        dataSetProperties_raw['file'][0] = dataSet_raw[0].getNumber()
        dataSetProperties_raw['file'][1] = dataSet_raw[-1].getNumber()
        dataSetProperties['file'][0] = dataSet_raw[0].getNumber()
        dataSetProperties['file'][1] = dataSet_raw[-1].getNumber()
        
        
    elif inputArray[0] == 'show':
        showCommands(inputArray)
        
    
    elif inputArray[0] == 'safe':
        if dataSet== []:
            print "ERROR: no data. Try 'load file.txt'"
            return False
        if len(inputArray) == 1:
            print 'ERROR: which format shall be saved? [origin/igor]'
            return False
        if inputArray[1] == 'igor':
            Igor_safefile(dataSet, wavelength, dataSetProperties, plotProperties['name'], path)
        elif inputArray[1] == 'origin':
            origin_safefile(dataSet, wavelength, dataSetProperties, path)
        else:
            print 'ERROR: No valid safe format. Try [origin/igor]'
            return False
        
    elif inputArray[0] == 'axis':
        axisCommand(inputArray, plotProperties)
        
        
        
        
        
    elif inputArray[0] == 'selector':
        if len(inputArray) == 1:
            print 'ERROR: enter how many files to be skipped.'
            return False
        try:
            int(inputArray[1])
            dataSetProperties['selector'] = int(inputArray[1])
        except:
            print 'ERROR: selector must be of type integer.'
            return False
        
    elif inputArray[0] == 'smooth':
        if len(inputArray) == 1:
            print 'ERROR: enter smoothing factor.'
            return False
        smoothCommand(inputArray[1], dataSet, dataSet_raw, wavelength, dataSetProperties)
    
    
    else:
        print "ERROR: No valid command. Try 'help' for list of commands" 
            
                

""" Commands """
######################################################

def helpCommand():
    print '########################################'
    print "UVvisReader commands"
    print ' '
    print "help"
    print "** this list **"
    print ' '
    print "load FILENAME.txt"
    print "** loads files with same filename, beginning with 1,"
    print "   ending with last found filenumber **"
    print ' '
    print "exit"
    print "** close UVvisReader **"
    print ' '
    print "plot [2d/3d/n]"
    print "** plot different figures. Example: UVvisReader>>plot 3d **"
    print ' '
    print "cutoff [wavelength/files] (start,stop)"
    print "** Set the limits for wavelength or filenumber."
    print "   WARNING: No spaces in (start,stop) allowed. **"
    print ' '
    print "axis [2d/3d/n] [x/y/z] VALUE"
    print "** Change axis description **"
    print ' '
    print "safe [origin/igor]"
    print "** safe modified dataSet as origin/igor importable. **"
    print ' '
    print 'selector VALUE'
    print "** Set how many files shall be skipped **"
    print ' '
    print "smooth VALUE"
    print "** Set smoothing value **"
    print ' '
    print "show [dataSet/dataSet_raw/wavelength/wavelength_raw] [NONE/len)"
    print "** DEBUG function: Show the actual values of intern variables **"
    print ' '
    print "show [dataSetproperties/dataSetProperties_raw/plotProperties]"
    print "** DEBUG function: Show the actual values of intern variables **"
    print '########################################'
    print ' '
    
def cutoffWavelengthCommand(min_max, dataSet_raw, dataSet, dataSetProperties_raw, dataSetProperties):
    mymin, mymax = getValuesFromStringArray(min_max)
    mymin = float(mymin)
    mymax = float(mymax)
    if mymax <= mymin:
        return
    if mymax > dataSetProperties_raw['wl'][1]:
        dataSetProperties['wl'][1] = dataSetProperties_raw['wl'][1]
        print 'INFO: max value out of bond. Setting to highest measured wavelength'
    else:
        dataSetProperties['wl'][1] = mymax
        
    if mymin < dataSetProperties_raw['wl'][0]:
        dataSetProperties['wl'][0] = dataSetProperties_raw['wl'][0]
        print 'INFO: min value out of bond. Setting to lowest measured wavelength'
    else:
        dataSetProperties['wl'][0] = mymin

def cutoffFilenumbersCommand(min_max, dataSet_raw, dataSet, dataSetProperties_raw, dataSetProperties):
    mymin, mymax = getValuesFromStringArray(min_max)
    mymin = int(mymin)
    mymax = int(mymax)
    if mymax <= mymin:
        return
    if mymax > dataSetProperties_raw['file'][1]:
        dataSetProperties['file'][1] = dataSetProperties_raw['file'][1]
        print 'INFO: max value out of bond. Setting to highest file number'
    else:
        dataSetProperties['file'][1] = mymax
        
    if mymin < dataSetProperties_raw['file'][0]:
        dataSetProperties['file'][0] = dataSetProperties_raw['file'][0]
        print 'INFO: min value out of bond. Setting to lowest file number'
    else:
        dataSetProperties['file'][0] = mymin
        
        
def smoothCommand(smoo, dataSet, dataSet_raw, wavelength, dataSetProperties):
    smoo = int(smoo)
    if smoo <=1:
        dataSet[:] = dataSet_raw[:]
        dataSetProperties['smooth'] = 1
    else:
        ds, wl = smooth_files(dataSet_raw, wavelength_raw, smoo)
        dataSet[:] = ds[:]
        wavelength[:] = wl[:]
        dataSetProperties['smooth'] = smoo
        
        
def axisCommand(inputArray, plotProperties):
    if len(inputArray) < 4:
        print "ERROR: wrong format. Example: axis 2d x wavelength (nm)"
        return False
    myvalue = inputArray[3]
    if len(inputArray) > 4:
        for i in range(4,len(inputArray)):
            myvalue += ' '+ inputArray[i]
    if inputArray[1] in ['2d', '3d', 'n']:
        inputArray[1] += '_axinfo'
    else:
        print 'ERROR: No valid plot for axis change'
        return False
        
    if inputArray[2] not in ['x', 'y', 'z']:
        print 'ERROR: No valid coordinate. Try [x/y/z].'
        return False
    for i, val in enumerate(['x', 'y', 'z']):
        if val == inputArray[2]:
            plotProperties[inputArray[1]][i] = myvalue
        
        

def showCommands(inputArray):
    if inputArray[1] == 'dataSet':
        if len(inputArray)==2:
            print dataSet
        elif inputArray[2] == 'len':
            print len(dataSet)
    elif inputArray[1] == 'dataSet_raw':
        if len(inputArray)==2:
            print dataSet_raw
        elif inputArray[2] == 'len':
            print len(dataSet_raw)
    elif inputArray[1] == 'wavelength':
        if len(inputArray)==2:
            print wavelength
        elif inputArray[2] == 'len':
            print len(wavelength)
    elif inputArray[1] == 'wavelength_raw':
        if len(inputArray)==2:
            print wavelength_raw
        elif inputArray[2] == 'len':
            print len(wavelength_raw)
    elif inputArray[1] == 'plotProperties':
        print plotProperties
    elif inputArray[1] == 'dataSetProperties':
        print dataSetProperties
    elif inputArray[1] == 'dataSetProperties_raw':
        print dataSetProperties_raw
    else:
        print '######################################'
        print 'Invalid value to be shown. Values are:'
        print 'dataSet, dataSet_raw'
        print 'wavelength, wavelength_raw'
        print 'dataSetProperties, dataSetProperties_raw'
        print 'plotProperties'
        print '######################################'
        print ' '
    
    
    
    
    
    
##################################
#########     CODE     ###########
##################################  

path = os.getcwd()
 
dataSet = []
dataSet_raw = []
wavelength = []
wavelength_raw = []
dataSetProperties = {'wl': [None, None], 'file': [None, None], 'selector':1, 'smooth':1}
dataSetProperties_raw = {'wl': [None, None], 'file': [None, None], 'selector':1, 'smooth':1}
plotProperties = {'name':'',
                  'n_axinfo': ['wavelength (nm)', 'Transmission/Reflection (a.u.)'],
                  '2d_axinfo': ['file number', 'wavelength (nm)', 'Transmission/Reflection (a.u.)'],
                  '3d_axinfo': ['file number', 'wavelength (nm)', 'Transmission/Reflection (a.u.)']
                  }

print "################ UVvisReader ###################"
print "# by Jannik Woehnert and Matthias Schwartzkopf #"
print "# Version 0.2 beta                             #"
print "# try 'help' for more information              #"
print "################################################" 
myexit = False
while myexit != True:
    commandLine = raw_input('UVvisReader >>')
    myexit = handleCommands(convertCommand(commandLine),
                            dataSet_raw = dataSet_raw, wavelength_raw = wavelength_raw,
                            dataSet = dataSet, wavelength=wavelength,
                            dataSetProperties = dataSetProperties, dataSetProperties_raw = dataSetProperties_raw,
                            plotProperties = plotProperties, path = path
                            )

