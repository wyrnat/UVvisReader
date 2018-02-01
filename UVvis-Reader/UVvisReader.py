'''
Created on 14.09.2016

@author: Jannik Woehnert
@version: 0.4 beta
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
        self.fileName = None
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
        
    def setFileName(self, name):
        self.fileName = name
    
    def getFileName(self):
        return self.fileName
             


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
        if numberlen == len(nameList):
            numberlen = None
            break
        if nameList[numberlen] == '0':
            break
        else:
            numberlen +=1
    if numberlen == None:
        filenumber = ''
        filename = ''.join(nameList)  
    else:   
        filenumber = ''.join(nameList[numberlen:])
        filename = ''.join(nameList[:numberlen])
    filetype = nameAndType[-1]
    if len(pathArray) == 1:
        additionalPath = ''
    else:
        additionalPath = '/'+pathArray[0]
    
    return additionalPath, filename, filenumber, filetype



def extractFiles(fileNameString, filepath, start = 1, end = None, header=88, Cwl = 1, Ctr = 8):
    """
    Multifunctional file handling method. Adaptable to different UVvis safefiles.
    @param filenamestring: the path to the UVvis files fromthe view of UVvisReader
    @param filepath: gloabal path of the UVvisReader
    @param start: which file to import first. default 1.
    @param end: last file to import. default None.
    @param header: How many lines in file to skip before data. default 88.
    @param Cwl: column of wavelength. default 1.
    @param Ctr: column of Transmission/Reflection. default 8.
    @return:  dataArray, wavelength array, path from UVvisReader to files,
    @return:  raw filename, filetype, number of digits in filename
    """
    waveLength = []
    dataBib = {'TR': []}
    
    additionalPath, filename, filenumber, filetype = extractFileName(fileNameString)
    pzf = len(filenumber)   #prozero factor                      
    
    k = start
    while True:
        partArray = dataArray()
        partArray.setNumber(k)
        prozero = ((k<10)*(pzf>1)*('0')+
                   (k<100)*(pzf>2)*('0')+
                   (k<1000)*(pzf>3)*('0')+
                   (k<10000)*(pzf>4)*('0')+
                   (k<100000)*(pzf>5)*('0')
                   )
        if (not os.path.isfile(filepath+additionalPath+'/'+filename+prozero+str(k)+"."+filetype)
            or ( end != None and k == end+1)):
            print 'INFO: Finished Extraction at:'
            print filepath+additionalPath+'/'+filename+prozero+str(k-1)+"."+filetype
            break
        myfile = open(filepath+additionalPath+'/'+filename+prozero+str(k)+"."+filetype,'r')
        
        waveLength, partArray = extraction(partArray, waveLength, myfile, k, start, header, Cwl, Ctr)
        
        k+=1
        
        if partArray.boolZeroMeasure() == False:
            dataBib['TR'].append(partArray)
        
        myfile.close()
    return dataBib['TR'], waveLength, additionalPath, filename, filetype, pzf

def extraction(partArray, waveLength, myfile, k, start, header, Cwl, Ctr):
    """
    Help function for extract files
    @param partArray: dataArray with TR values
    @param wavelength: wavelength Array
    @param myfile: stream to data file
    @param k: file number of UVvis data
    @param start: first file number
    @param header: header lines in data file
    @param Cwl, Ctr: column lines of wavelength and TR in data file
    @return: filled wavelength and dataArray
    """
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
                if k == start:
                    waveLength.append(float(newline[Cwl]))
                partArray.append(float(newline[Ctr]))
    
        j+=1
    return waveLength, partArray

def fillVariables(ds, wl, additionalPath, filename, filetype, pzf,
                  dataSet_raw, dataSet, wavelength_raw, wavelength,
                  dataSetProperties_raw, dataSetProperties, fileProperties):
    
    #copy raw data to working data
    dataSet[:] = ds[:]
    wavelength[:] = wl[:]
    dataSet_raw[:] = ds[:]
    wavelength_raw[:] = wl[:]
    #set file properties
    fileProperties['name'] = filename
    fileProperties['additionalPath'] = additionalPath
    fileProperties['type'] = filetype
    #set dataSet properties raw  and working file
    dataSetProperties_raw['wl'][0] = wavelength_raw[0]
    dataSetProperties_raw['wl'][1] = wavelength_raw[-1]
    dataSetProperties['wl'][0] = wavelength_raw[0]
    dataSetProperties['wl'][1] = wavelength_raw[-1]
    dataSetProperties_raw['file'][0] = dataSet_raw[0].getNumber()
    dataSetProperties_raw['file'][1] = dataSet_raw[-1].getNumber()
    dataSetProperties['file'][0] = dataSet_raw[0].getNumber()
    dataSetProperties['file'][1] = dataSet_raw[-1].getNumber()
    
    
def extractSubstrat(fileNameString, filepath):
    additionalPath, filename, filenumber, filetype = extractFileName(fileNameString)
    completefilepath = filepath+additionalPath+'/'+filename+filenumber+"."+filetype
    if os.path.isfile(completefilepath):  
        myfile = open(completefilepath, 'r')
        wavelength, partArray = extraction([], [], myfile, 0, 88, 1, 8)
        print "INFO: extracted substrate at:"
        print completefilepath
        return partArray
    else:
        print "ERROR: substrate file not found."
    


""" Handle Data"""
#######################################################################

def smooth(dArray, box_pts):
    """
    Smoothes data with numpy.convole package
    @param y: data set
    @param box_pts: number of measure points in both directions to smooth the value
    @return: y_smooth: smoothed data set
    """
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(dArray, box, mode='valid')
    return y_smooth


def wlExtrem(wavelength, dataSetProperties):
    
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
        
    return wlMinIndex, wlMaxIndex

def fileExtrem(dataSet, dataSetProperties):
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
    
    return fileMinIndex, fileMaxIndex

def prepareDataSet(dataSet, wavelength, substrate, dataSetProperties, fileNameProperties):
    """
    Cut off files, wavelength or skip files.
    This function provides the choosen dataSet for plotting and saving
    @param dataSet: list of dataArray objects of the measure
    @param wavelength: list of loaded wavelength
    @param dataSetProperties: dictionary with dataSet-changing values
    @return: prepared dataSet and wavelength with values from dataSetProperties
    """
    nameFunc = lambda dArray: fileNameProperties['factor']*dArray.getNumber() + fileNameProperties['offset']

    #extract extrema of wavelength and file number
    wlMinIndex, wlMaxIndex =  wlExtrem(wavelength, dataSetProperties)
    fileMinIndex, fileMaxIndex = fileExtrem(dataSet, dataSetProperties)
    slct = dataSetProperties['selector']
    #cut dataSet
    cutFileArray = dataSet[fileMinIndex:fileMaxIndex:slct]
    
    # set smoothed and cutted wavelength
    preparedwavelength = wavelength[wlMinIndex:wlMaxIndex]
    preparedwavelength = smooth(preparedwavelength, dataSetProperties['smooth'])
    
    #set smoothed, cutted and substrate reduced dataSet
    
    if dataSetProperties['substrate']:
        preparedsubstrate = substrate[wlMinIndex:wlMaxIndex]
    preparedArray = []
    for dArray in cutFileArray:
        preArray = dArray[wlMinIndex:wlMaxIndex]
        if dataSetProperties['substrate']:
            newArray = dataArray(np.array(preArray) * np.array(preparedsubstrate)/100.)
        else:
            newArray = preArray
        newArray = dataArray(smooth(newArray, dataSetProperties['smooth']))   
        newArray.setNumber(dArray.getNumber())
        newArray.setFileName( str(nameFunc(dArray)) + ' ' + fileNameProperties['unit'] )
        preparedArray.append(newArray)
        
    return preparedArray, preparedwavelength



def findExtremum(dataSet, wavelength, fileNameProperties, start, stop, minmax):
    nameFunc = lambda dArray: fileNameProperties['factor']*dArray.getNumber() + fileNameProperties['offset']
    exList = []
    nameList =[]
    wlMinIndex, wlMaxIndex = wlExtrem(wavelength, {'wl':[start, stop]})
    for dArray in dataSet:
        ds = dArray[wlMinIndex:wlMaxIndex]
        wl = wavelength[wlMinIndex:wlMaxIndex]
        if minmax == 'min':
            i = ds.index(min(ds))
        elif minmax == 'max':
            i = ds.index(max(ds))
        exList.append(wl[i])
        nameList.append(nameFunc(dArray))
        
    return nameList, exList
        
        


""" Safe Files """
####################

def origin_safefile(dataSet, wavelength, substrate, dataSetProperties, fileProperties):
    """
    Creates a .txt File with the header wavelength and the extracted file indices.
    In the wavelength column there are all extracted wavelength.
    Under each extracted file indice there is the measured absorbance.
    @param dataSet: dataset of absorbance values for for every extracted file
    @param wavelength: array of extracted wavelength
    @param dataSetProperties: dictionary with dataSet-changing values
    @param filepath: path to the current working directory
    """
    filepath = fileProperties['path']
    additionalPath = fileProperties['additionalPath']
    if not os.path.exists(filepath+additionalPath+'/ORIGIN/'):
        os.mkdir(filepath+additionalPath+'/ORIGIN/')
    originfile = open(filepath+additionalPath+'/ORIGIN/Spectra.txt', 'w')
    text = 'wavelength'
    safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, substrate, dataSetProperties, fileNameProperties)
                
    for dArray in safedataSet:
        text += ' '+dArray.getFileName()
    text += '\n'
    for j in range(len(safewavelength)):
        text += str(safewavelength[j])
        for dArray in safedataSet:
            text +=' '+str(dArray[j])
        text += '\n'
        
    originfile.write(text)
    originfile.close()
    print "INFO: saved to "+additionalPath+ "/ORIGIN/spectra.txt"
    print "      wavelength: ("+str(dataSetProperties['wl'][0])+','+str(dataSetProperties['wl'][1])+')'
    print "      files: ("+str(dataSetProperties['file'][0])+','+str(dataSetProperties['file'][1])+')'
    print "      smooth: "+str(dataSetProperties['smooth'])
    print "      selector: "+ str(dataSetProperties['selector'])
    print "      substrate substracted: "+ ("Yes" if dataSetProperties['substrate'] else "No")
    
def Igor_safefile(dataSet, wavelength, substrate, dataSetProperties, fileProperties):
    """
    Safes data to IGOR compatible files.
    Each file contains wavelength on first column and TR on second column.
    It considers the file limits, they can be set with
    >> setlimits file [min, max]
    """
    filepath = fileProperties['path']
    additionalPath = fileProperties['additionalPath']
    pzf = fileProperties['pzf']
    name = fileProperties['name']
    type = fileProperties['type']
    
    if not os.path.exists(filepath+additionalPath+'/IGOR/'):
        os.mkdir(filepath+additionalPath+'/IGOR/')
        
    safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, substrate, dataSetProperties, fileNameProperties)
        
    for dArray in safedataSet:
        nr = dArray.getNumber()
        prozero = ((nr<10)*(pzf>1)*('0')+
                    (nr<100)*(pzf>2)*('0')+
                    (nr<1000)*(pzf>3)*('0')+
                    (nr<10000)*(pzf>4)*('0')+
                    (nr<100000)*(pzf>5)*('0')
                    )
        newfile = open(filepath+additionalPath+'/IGOR/'+name+"_"+prozero+str(nr)+"_IGOR."+type, 'w')
        newfile.write('wavelength TR\n')
        for j, wl in enumerate(safewavelength):
            newfile.write(str(wl)+' '+str(dArray[j])+'\n')
        newfile.close()
    
    print "INFO: saved to "+additionalPath+"/IGOR/NAME_IGOR.txt"
    print "      wavelength: ("+str(dataSetProperties['wl'][0])+','+str(dataSetProperties['wl'][1])+')'
    print "      files: ("+str(dataSetProperties['file'][0])+','+str(dataSetProperties['file'][1])+')'
    print "      smooth: "+str(dataSetProperties['smooth'])
    print "      selector: "+ str(dataSetProperties['selector'])
    print "      substrate substracted: "+ ("Yes" if dataSetProperties['substrate'] else "No")


def special_safefile(specialData, name, fileProperties):
    filepath = fileProperties['path']
    additionalPath = fileProperties['additionalPath']
    if not os.path.exists(filepath+additionalPath+'/SPECIAL/'):
        os.mkdir(filepath+additionalPath+'/SPECIAL/')
    
    filenameList = specialData[name][0]
    extremumList = specialData[name][1]
    print specialData[name][0]
        
    myfile = open(filepath+additionalPath+'/SPECIAL/'+name+".txt", 'w')
    myfile.write('file wavelength\n')
    for i, filename in enumerate(filenameList):
        myfile.write(str(filename) + ' ' + str(extremumList[i])+'\n')
    myfile.close()

""" Visualisation """
#####################
def plot2d(dataSet, wavelength, substrate, plotProperties, dataSetProperties, fileProperties):
    
    safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, substrate, dataSetProperties, fileNameProperties)
    
    fig, ax = plt.subplots()
    
    cax = ax.imshow(safedataSet)
    ax.set_title(fileProperties['name'])
    ax.set_xlabel(plotProperties['2d_axinfo'][0])
    ax.set_ylabel(plotProperties['2d_axinfo'][1])
    
    cbar = fig.colorbar(cax)
    cbar.set_label(plotProperties['2d_axinfo'][2])
    fig.canvas.set_window_title('UVvis Reader 2D ColorPlot')

    plt.show()
    
def multipleLines2dPlot(dataSet, wavelength, substrate, plotProperties, dataSetProperties, fileProperties):
    safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, substrate, dataSetProperties, fileNameProperties)
    mycm = plt.get_cmap('nipy_spectral')
    N = float(len(safedataSet))
    for i, dArray in enumerate(safedataSet):
        plt.plot(safewavelength, dArray, color=mycm(i/N), label=dArray.getFileName())
    plt.xlabel(plotProperties['n_axinfo'][0])
    plt.ylabel(plotProperties['n_axinfo'][1])
    plt.title(fileProperties['name'])
    plt.legend()
    plt.show()

def plot3d(dataSet, wavelength, substrate, plotProperties, dataSetProperties, fileProperties):
    safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, substrate, dataSetProperties, fileNameProperties)
    X = []
    #func = lambda dArray: fileNameProperties['factor']*dArray.getNumber() + fileNameProperties['offset']
    for dArray in safedataSet:
        #X.append(func(dArray))
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
    ax.set_title(fileProperties['name'])
    ax.plot_surface(X, Y, Z, rstride = plotProperties['3d_precision'][0], cstride = plotProperties['3d_precision'][1],
                    cmap=plt.get_cmap('nipy_spectral'), linewidth=0, antialiased=True
                    )
    plt.show()
    
def plotSpecial(name, specialData):
    nameList = specialData[name][0]
    exList = specialData[name][1]
    plt.plot(nameList, exList, color='r')
    plt.xlabel('filename')
    plt.ylabel("wavelength [nm]")
    plt.title("Special plot: "+name)
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
                   dataSetProperties_raw, substrate, dataSetProperties,
                   plotProperties, fileNameProperties, fileProperties):
    
    
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
            plot2d(dataSet, wavelength, substrate, plotProperties, dataSetProperties, fileProperties)
        elif inputArray[1] == '3d':
            plot3d(dataSet, wavelength, substrate, plotProperties, dataSetProperties, fileProperties)
        elif inputArray[1] == 'n':
            multipleLines2dPlot(dataSet, wavelength, substrate, plotProperties, dataSetProperties, fileProperties)
        elif inputArray[1] == 'special':
            if len(inputArray) == 3:
                if inputArray[2] in specialData:
                    plotSpecial(inputArray[2], specialData)
                else:
                    print "ERROR: special data not found."
            else:
                print "ERROR: invalid format. Example: 'plot special firstmin'."
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
        loadCommand(inputArray, dataSet_raw, dataSet,
                wavelength_raw, wavelength, dataSetProperties_raw,
                dataSetProperties, fileProperties)
        
        
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
            Igor_safefile(dataSet, wavelength, substrate, dataSetProperties, fileProperties)
        elif inputArray[1] == 'origin':
            origin_safefile(dataSet, wavelength, substrate, dataSetProperties, fileProperties)
        elif inputArray[1] == 'special':
            if (len(inputArray) >= 3 and inputArray[2] in specialData):
                special_safefile(specialData, inputArray[2], fileProperties)
            else:
                print "ERROR: special plot format invalid or data not found."
        else:
            print 'ERROR: No valid safe format. Try [origin/igor]'
            return False
        
    elif inputArray[0] == 'axis':
        axisCommand(inputArray, plotProperties)
        
        
    elif inputArray[0] == 'filename':
        if len(inputArray) < 3:
            print "ERROR: wrong syntax. Try 'filename [factor/offset/unit] VALUE"
            return False
        fileNameCommand(inputArray, dataSet, fileNameProperties)
        
        
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
        smoothCommand(inputArray[1], dataSetProperties)
    
    elif inputArray[0] == '3d':
        if inputArray[1] == 'precision':
            try:
                x, y = getValuesFromStringArray(inputArray[2])
                x = int(x)
                y = int(y)
                plotProperties['3d_precision'] = [x,y]
            except:
                print "ERROR: could not read Integer Values. Example: (4,10)"
                return False
        else:
            print "ERROR: no valid format. Try '3d precision (x,y)'."
            return False
    
    elif inputArray[0] == 'substrate':
        if len(inputArray) != 2:
            print "ERROR: wrong syntax. Example: 'substrate 1'"
        if inputArray[1] == '0':
            dataSetProperties['substrate'] = False
        elif inputArray[1] == '1':
            dataSetProperties['substrate'] = True
        else:
            print "ERROR: substrate value either '0' or '1'."
            
    elif inputArray[0] == 'find':
        if len(inputArray) != 4:
            print "ERROR: wrong format. Example: 'find min (250,350) firstmin'."
            return False
        findCommand(inputArray, dataSet, wavelength, substrate, dataSetProperties, fileNameProperties)
    
    else:
        print "ERROR: No valid command. Try 'help' for list of commands" 
            
                

""" Commands """
######################################################

def loadCommand(inputArray, dataSet_raw, dataSet,
                wavelength_raw, wavelength, dataSetProperties_raw,
                dataSetProperties, fileProperties):
    if len(inputArray) < 3:
        print 'ERROR: please enter a file name and specfify what to load.'
        print '       Example: load data test001.txt'
        return False
    if inputArray[1] == 'substrate':
        sb = extractSubstrat(inputArray[2], fileProperties['path'])
        substrate[:] = sb[:]
        dataSetProperties['substrate'] = True
    elif inputArray[1] == 'data':
        ds, wl, additionalPath, filename, filetype, pzf = extractFiles(inputArray[2], fileProperties['path'])
        if ds == [] or wl == []:
            print 'ERROR: data could not be loaded. Got path and filename correct?'
            return False
        fillVariables(ds, wl, additionalPath, filename, filetype, pzf,
                      dataSet_raw, dataSet, wavelength_raw, wavelength,
                      dataSetProperties_raw, dataSetProperties, fileProperties
                      )

def helpCommand():
    print '################################################'
    print "UVvisReader commands"
    print ' '
    print "3d precision (x,y)"
    print "** Set precision of 3d plot for x/y axis **"
    print ' '
    print "axis [2d/3d/n] [x/y/z] VALUE"
    print "** Change axis description **"
    print ' '
    print "cutoff [wavelength/files] (start,stop)"
    print "** Set the limits for wavelength or filenumber."
    print "   WARNING: No spaces in (start,stop) allowed. **"
    print ' '
    print "exit"
    print "** close UVvisReader **"
    print ' '
    print "filename [factor/offset/unit] VALUE"
    print "** Chance the filename for 'plot n' and 'safe origin'"
    print "   function: (factor*FILENUMBER + offset) unit **"
    print ' '
    print "help"
    print "** this list **"
    print ' '
    print "load [data/substrate] FILENAME.txt"
    print "** loads UVvis data, that means all files it can found in numerical order,"
    print "   as well as the substrate measure to gain the real UVvis values. **"
    print ' '
    print "plot [2d/3d/n/special] {special:NAME}"
    print "** plot different figures. Example: 'plot 3d'"
    print "   pLot special data (find min/max). Example: 'plot special firstmin' **"
    print ' '
    print "safe [origin/igor/special] {special:NAME}"
    print "** safe modified dataSet as origin/igor importable."
    print "   safe special data. Example: 'safe special firstmin' **"
    print ' '
    print 'selector VALUE'
    print "** Set how many files shall be skipped **"
    print ' '
    print "show [dataSet/dataSet_raw/wavelength/wavelength_raw] [NONE/len)"
    print "** DEBUG function: Show the actual values of intern variables **"
    print ' '
    print "show [dataSetproperties/dataSetProperties_raw/plotProperties]"
    print "** DEBUG function: Show the actual values of intern variables **"
    print ' '
    print "smooth VALUE"
    print "** Set smoothing value **"
    print '################################################'
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
        
        
def smoothCommand(smoo, dataSetProperties):
    smoo = int(smoo)
    if smoo <=1:
        dataSetProperties['smooth'] = 1
    else:
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
    elif inputArray[1] == 'fileNameProperties':
        print fileNameProperties
    elif inputArray[1] == 'fileProperties':
        print fileProperties
    elif inputArray[1] == 'substrate':
        if len(inputArray) == 2:
            print substrate
        elif inputArray[2] == 'len':
            print len(substrate)
    elif inputArray[1] == 'specialData':
        if len(inputArray) == 2:
            print specialData.keys()
        else:
            if inputArray[2] in specialData:
                print specialData[inputArray[2]]
    else:
        print '######################################'
        print 'Invalid value to be shown. Values are:'
        print 'dataSet, dataSet_raw'
        print 'wavelength, wavelength_raw'
        print 'dataSetProperties, dataSetProperties_raw'
        print 'plotProperties, fileNameProperties'
        print 'substrate, specialData'
        print '######################################'
    
    
def fileNameCommand(inputArray, dataSet, fileNameProperties):
    if inputArray[1] == 'unit':
        unitstr = inputArray[2]
        if len(inputArray) > 3:
            for word in inputArray[3:]:
                unitstr += ' '+word
        fileNameProperties['unit'] = unitstr
    elif inputArray[1] in ['factor', 'offset']:
        try:
            valuefloat = float(inputArray[2])
            fileNameProperties[inputArray[1]] = valuefloat
        except:
            print 'ERROR: VALUE not float.'
            return False
            
    

def findCommand(inputArray, dataSet, wavelength, substrate, dataSetProperties, fileNameProperties):
    if inputArray[1] in ['min', 'max']:
        safedataSet, safewavelength = prepareDataSet(dataSet, wavelength, substrate,
                                                     dataSetProperties, fileNameProperties)
        try:
            minwl, maxwl = getValuesFromStringArray(inputArray[2])
            minwl = float(minwl)
            maxwl = float(maxwl)
        except:
            print "ERROR: wavelength area not float"
            return False
        nameList, exList = findExtremum(safedataSet, safewavelength, fileNameProperties, minwl, maxwl, inputArray[1])
        specialData[inputArray[3]] = [nameList, exList]
            
            
    
##################################
#########     CODE     ###########
##################################  
 
dataSet = []
dataSet_raw = []
wavelength = []
wavelength_raw = []
specialData = {}
substrate = []
dataSetProperties = {'wl': [None, None], 'file': [None, None], 'selector':1, 'smooth':1, 'substrate':False}
dataSetProperties_raw = {'wl': [None, None], 'file': [None, None], 'selector':1, 'smooth':1, 'substrate': False}
plotProperties = {'name':'',
                  'n_axinfo': ['wavelength (nm)', 'Transmission/Reflection (a.u.)'],
                  '2d_axinfo': ['wavelength (nm)', 'file number', 'Transmission/Reflection (a.u.)'],
                  '3d_axinfo': ['file number', 'wavelength (nm)', 'Transmission/Reflection (a.u.)'],
                  '3d_precision': [10,10]
                  }
fileProperties = {'path': os.getcwd(), 'additionalPath': '', 'name': '', 'pzf': 0, 'type': '', 'safefile':''}
fileNameProperties = {'factor':1, 'offset':0, 'unit': ''}

print "################ UVvisReader ###################"
print "# by Jannik Woehnert                           #"
print "# thanks to Matthias Schwartzkopf, Marc Gensch #"
print "# Version 0.4 beta                             #"
print "# try 'help' for more information              #"
print "################################################" 
myexit = False
while myexit != True:
    commandLine = raw_input('UVvisReader >>')
    myexit = handleCommands(convertCommand(commandLine),
                            dataSet_raw = dataSet_raw, wavelength_raw = wavelength_raw,
                            dataSet = dataSet, wavelength=wavelength,
                            dataSetProperties = dataSetProperties, dataSetProperties_raw = dataSetProperties_raw,
                            substrate = substrate, plotProperties = plotProperties,
                            fileNameProperties =fileNameProperties, fileProperties = fileProperties
                            )

