'''
Created on Sep 30, 2016
@author: Maximilian Pfeiffer
@argv: filename, -seperator: file to be opened and unicode string seperating data entries, respectively.
@argv: quantileWidth, baselineWidth: size of buckets to be scanned for the baseline / quantiles, respectively.
TODO: clean up of code! Output should be generated in one call
    general process outline:
    
    1. finish deconvolution
    2. generation of plots / transient property files into seperate functions, as well as INPUT
    3. All normalization / correction methods shouldn't take more input than necessary and have defaults
    4. 2 & 3 should collapse main.py to an overseeable format! -> if not: rinse & repeat
    5. COMPARE PERFORMANCES
        -how fast is it compared to alternatives, OPTIMIZE
        -how close is the result? OPTIMIZE
    6. Put methods in a library-esque format
    7. Documentation! 
        - README for the done script
        - explanation for methods
TIMELINE:
    1: by (including) 24th of April
    2: 26.4
    3: 28.4
    4: 1.5
    5: 7.5
    6. 9.5
    7. 11.5
'''
import time, argparse, cal_neuroIm
from numpy import mean,ceil
#'/home/maximilian/unistuff/paris_ens/cal_neuroim/driftData/a5 KO anaesthetized mice data/2014.02.04 a5 cage 18/TSeries-02042014-1019-4562/Analysis/Results.xls'
start_time = time.time() # get sys time to calculate runtime
 
descriptionString = "Takes a list of itensity matrices.\nNormalizes the input (for -q != 0), defines events in the data and converts those to spike trains.\nReturns spike trains and transient/event properties in respective files."
inputParser = argparse.ArgumentParser(description=descriptionString)
inputParser.add_argument('inputFile',nargs="+", help="Input matrices containing mean Intensities",type=str)        
inputParser.add_argument('outputFilePath',help="Full file path to save the output to", type=str)
inputParser.add_argument("-q",'-quantileWidth',"-qWidth", help="Width of the window to be considered for quantile normalization. Set to 0 to skip quantile Norm", nargs='?', type=int,default=150)        
inputParser.add_argument("-b",'-baselineWidth',"-bWidth",help="Width of the window to be considered for baseline normalization", nargs='?',type=int,default=150)                
inputParser.add_argument('-seperator',"-sep",help="Seperator used for entries in the input matrix", nargs='?',default="\t")        
inputParser.add_argument('-csv',help="Is the file in .csv format? (ImageJ Results.xls files are NOT .xls but csv, despite the file name!)",default=False,action='store_true')
args = inputParser.parse_args()

for iterationNumber, filename in enumerate(args.inputFile):
    
    loadObject = cal_neuroIm.importMatrix(filename,args.seperator,args.csv)
    rawMatrix = loadObject[0]
    numOfVals = loadObject[1]
    print("on file:" + str(filename) + " with " + str(numOfVals) + " ROIs")
    
    baseObject = cal_neuroIm.pushToBaseline(rawMatrix.copy(),args.b)
    baselineMatrix = baseObject[0]
    baselineCoordinates = baseObject[1]
    eventDetectObject = cal_neuroIm.eventDetect(baselineMatrix.copy(), args.q,len(rawMatrix)/150)
    quantileMatrix = eventDetectObject[1]
    transientMatrix = eventDetectObject[0]
    meanKernel = cal_neuroIm.createMeanKernel(transientMatrix)
    #TODO: array containing only the transient data -> check negative transients first and filter!
    #ALSO: transient kernel is what excatly, where do we get it from? -> Alpha function which was fit over single peak transients
       
    #for k in reversed(delList):
    #    rawMatrix = delete(rawMatrix, k, axis=1)
    #    quantileMatrix = delete(quantileMatrix,k,axis=1)
    #    correctionMatrix = delete(correctionMatrix,k,axis=1)
    #numOfVals -= len(delList)
    #transientList = cal_neuroIm.generateOutput(rawMatrix, baseObject, tempData, filename, numOfVals)

    for i in range(numOfVals):
        import matplotlib.pyplot as plt
        f,(axarr0,axarr1,axarr2) = plt.subplots(3, sharex=True)
        eventCoordinates = [transientMatrix[i][j].getCoordinates() for j in range(len(transientMatrix[i]))] # ugly hack that defeats the purpose. REFACTOR
        spikeTrain = (cal_neuroIm.deconvolve(transientMatrix[i], meanKernel))
        axarr2.plot(ceil(spikeTrain))
        axarr1.plot(baselineMatrix[:,i])
        axarr0.plot(rawMatrix[:,i],"r")
        minVal = float(min(quantileMatrix[:,i]))
        thisMean2 = mean(baselineMatrix[:,i])
        #minVal2 = float(min(baselineMatrix[:,i]))
        #maxVal2 = float(max(baselineMatrix[:,i]))
        if(transientMatrix[i]): # is there a transient in this time series?
            for j in transientMatrix[i]: # yes, iterate over transient objects
                coords = j.getCoordinates()
                axarr0.plot(coords,[minVal,minVal],'r-', lw=2)
                if ((coords[1]) != 0):
                    thisSlice = (j.getData())
                    #transientList.append(cal_neuroIm.Transient(thisSlice,j[0],j[1])) # call to transient kills 
                    #axarr1.plot([coords[0],coords[0]],[minVal2,maxVal2],'r-', lw=1)
        axarr0.plot(baselineCoordinates,[thisMean2,thisMean2],'k.',lw=1)
        #axarr0.plot([baselineCoordinates[0],baselineCoordinates[0]],[minVal2,maxVal2],'k.', lw=1)
        #axarr0.plot([baselineCoordinates[1],baselineCoordinates[1]],[minVal2,maxVal2],'k.', lw=1)
        plt.savefig(str(filename) + str(i)+".png")
        plt.close()
    '''try:    
        with open(str(args.outputFilePath) + str(filename).split('/')[-1] + ".csv", 'w') as outFile:
            outFile.write("Amplitude\tRiseTime\tDecayTime\tmeanIntensity\ttotalLength\tnumOfPeaks\tstartTime\tendTime\n")    
            for num,t in enumerate(transientList):
                outFile.write(str(t.getAmplitude()) + "\t" + str(t.getRiseTime())+ "\t" + str(t.getDecayTime()) + "\t" + str(t.getMeanIntensity()) + "\t" + str(t.getTotalTime()) +"\t" + str(t.getNumOfPeaks())+"\t" + str(t.getStartTime()) + "\t"+ str(t.getEndTime()) + "\n")
                if(t.getIsLast()):
                    outFile.write("__________________________ROI" + str(num) + "__________________________\n")
    except IOError:
        print(str(args.outputFilePath) +"output" +  str(filename).split('\\')[-1] + " is not a valid file path for output!")
        exit()
        '''
    #cal_neuroIm.writeOut(baselineMatrix, str(filename) + "base")
    #cal_neuroIm.writeOut(quantileMatrix, str(filename) + "quant")

    
print("--- %s seconds ---" % (time.time() - start_time))

