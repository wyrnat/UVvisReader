import os
path = os.getcwd()
if not os.path.exists(path+'/NEW/'):
    os.mkdir(path+'/NEW/')
def convertFiles(name, start, stop):
    for i in range(start,stop+1):
        prozero = (i==0)*('0')+(i<10)*('0')+(i<100)*('0')+(i<1000)*('0')
        myfile = open(path+'/'+name+"_"+prozero+str(i)+".txt",'r')
        newfile = open(path+'/NEW/'+name+"_"+prozero+str(i)+"_NEW.txt", 'w')

        line = "marc"
        j = 0
        while j < 2137:
              line = myfile.readline()
              if (j > 88):
                  newline = line.split(';')
                  for i, element in enumerate(newline):
                      newline[i] = element.replace('   ','0')
                  writtenline =''
                  for element in newline:  
                      writtenline += (element+' ')
                  newfile.write(writtenline)
              j+=1
        newfile.close()
        myfile.close()
          
