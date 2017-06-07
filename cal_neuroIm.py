'''
Created on Oct 3, 2016
script contains methods used for main.py
@author: maximilian
'''
from numpy import percentile,zeros,var, mean, arange, concatenate, savetxt, sign, ceil,\
matrix, std, around, fft, divide, abs, ones, resize, exp, nditer,trapz, argmax,\
    append, sum, argmin
from sys import maxint
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import norm
from scipy.optimize import curve_fit
import pandas, xlrd, time
from xlrd.biffh import XLRDError
from scipy.stats.stats import mode

def importMatrix(filename,valSeperator):
    '''
    imports csv or xls formatted file, returns ndArray 
    @input: filename: input path; valSeperator: seperator between entries, treats file as xlrd if none
    @return: numOfVals: number of collumns in file (cells); dataMatrix: ndArray containing file data (dtype float)
    '''
    dataMatrix = []
    try:
        if(valSeperator):
            frame = pandas.read_csv(filename,sep = valSeperator)
            for x in frame.values:
                dataMatrix.append(x)
            numOfVals = x.size
        else:
            try:
                book = xlrd.open_workbook(filename)
                sh = book.sheet_by_index(0)
                numOfVals = sh.ncols
                for rx in range(sh.nrows):
                    dataMatrix.append(getValues(sh.row(rx)))
            except XLRDError:
                print("Wrong format! Try adding -csv to the invocation")
                exit()
    except IOError:
        print(str(filename) + " is not a valid file path !")
        exit()
    from numpy import matrix as m
    dataMatrix = m(dataMatrix).A
    return(dataMatrix, numOfVals)

'''
@param numOfVals: number of Values per row (e.g. ROIs)
@param row: current row of csvfile
@return: Vector of length 'numOfVals', containing the entry values row-vise
actually yields one entry less than numOfVals - ImageJ result.xls file headers start with a whitespace (value seperator)
TODO: this looks extremely un-pythonic and stupid ... what am I using it for again?'''
def getValues (row):
    thisRow = []
    for i in row:
        thisRow.append(float(i.value))
    return thisRow;

def meanSlope(inArray):
    '''
    Given a time series, return the slope (which only depends on the start and end point
    '''
    iterNum = len(inArray)
    avgSlope = (inArray[-1] - inArray[0])/iterNum
    return (float(avgSlope))

def eventDetect (dataMatrix, quantileWidth,slopeWidth):
    '''
    TODO: clean up and turn into ndArray ? -> check if it is the bottleneck first!
    
    @param dataMatrix: NxM matrix of the data containing transients and drift(possibly)
    @param windowSize: size of the window considered for the quantile normalization (none if q = 0)
    @param slopeWidth: size of the window considered for the calculation of slopes, critical parameter for event detection
    @return: dataMatrix, eventEndCoordinates, delList, correctionMatrix, baseLineArray
            dataMatrix:  quantile normalized (if q>0) NxM matrix
            eventEndCoordinates: tuples of on- and offset coordinates of transients
            delList: list of ROIs (denominated via their vertical position integer in the matrix) which were not considered and should be dropped
            correctionMatrix: NxM Matrix containing the subtracted quantiles for visualization of the normalization
            baseLineArray: 1xM array containing the start coordinates of the sections considered as the baseline 
    ''' 
    transients = [] # global list for lists of transients
    numOfVals = dataMatrix.shape[1] # number of ROIs in file
    numOfEntries = len(dataMatrix) # number of entries in file
    startNumOfHits = 5 # number of slope-hits to be made before an event is declared as beginning
    minEventLength = 30 # minimum length of events
    thresholdList = []
    slopeDistributions = []
    correctionMatrix = zeros([numOfEntries,numOfVals])
    
    # here we generate a global noise / slope threshold for all ROIs of the file
    for collumn in nditer(dataMatrix, order='F',flags=['external_loop'],op_flags=['readonly']): 
        slopeList = zeros(numOfEntries - slopeWidth)
        #slopeThreshold += abs(meanSlope(detectPeakline(collumn, slopeWidth)[0])) # threshold to be passed for start of significance
        for i in range(len(slopeList)):
            slopeList[i]= (meanSlope(collumn[i:i+slopeWidth]))  # threshold to be passed for start of significance
        # best fit of data according to the log-likelihood maximization estimate
        thisMode = mode(around(slopeList,3),axis=None)[0]
        slopeList2 = []
        for i in slopeList:
            if i < thisMode:
                slopeList2.append(i)
        
        mu, sigma = norm.fit(slopeList2)
        # the histogram of the data
        thresholdList.append(sigma*2) 
        import matplotlib.pyplot as plt
        import matplotlib.mlab as mlab
        n, bins, patches = plt.hist(slopeList, 100, normed=1, facecolor='green', alpha=0.75)
        y = mlab.normpdf( bins, mu, sigma)
        slopeDistributions.append((bins,y,slopeList,sigma*2))
        plt.close()
        
        #plt.plot(bins, y, 'r--', linewidth=2)
        #plt.plot([sigma*2,sigma*2],[0,5],'k--')
        #plot
        '''
        plt.title(r'$\mathrm{test:}\ \mu=%.3f,\ \sigma=%.3f,\ thisMode=%.3f$' %(mu, sigma,thisMode))
        plt.grid(True)
        plt.savefig("Slope_distrib_" + str(num)+".png")
        plt.close()
        '''
        #noiseThreshold += (var(detectBaseline(collumn, baselineWidth)[0])) # threshold to be passed for end of significance
    #slopeThreshold = slopeThreshold/numOfVals
    # re-iterate over the data ROI-wise for determination of event positions - in reverse order because we delete skipped rows immediately
    for horizontalPosition, collumn in enumerate(dataMatrix.T): # generate local noise / slope threshold and consider if we should skip
        
        collumn = collumn.T # looks stupid but makes the general code more pythonic (instead of referencing collumns via slicing)
        thisSlopeThreshold = thresholdList[horizontalPosition]# threshold to be passed for start of significance 
        isEvent = False # bool to remember wether we are in event-mode
        eventStart = 0
        eventEnd = 0
        threshHit = 0
        #theseEventEndCoordinates = [] # local list for event coordinates, to be appended to the bigger version
        theseTransients = [] # local list for transients
        if(quantileWidth != 0):
            shiftValue = percentile(collumn,8)#shift of data to avoid negative values and the corresponding inversion
        else:
            shiftValue = 0
        for verticalPosition in range(numOfEntries): # iterate over number of pictures in file, e.g. row-wise
                
            if(verticalPosition+40 <= numOfEntries): # get correct slice of data matrix
                thatSlice = dataMatrix[verticalPosition:verticalPosition+slopeWidth,horizontalPosition]
            else:
                thatSlice = dataMatrix[numOfEntries-slopeWidth:,horizontalPosition]
            thisSlope = abs(meanSlope(thatSlice))   
            if (quantileWidth != 0):
                if (verticalPosition+quantileWidth<=numOfEntries):
                    correctionVal = percentile(collumn[verticalPosition:verticalPosition+quantileWidth],q=8)
                else:
                    correctionVal = correctionMatrix[verticalPosition-1,horizontalPosition]
                collumn[verticalPosition] -= correctionVal
                correctionMatrix[verticalPosition,horizontalPosition] = correctionVal
                
            if (not isEvent): 
                # we're not in an event ... so check if one is starting
                if (thisSlope > thisSlopeThreshold):
                    # it is! count hit
                    threshHit += 1
                    if(threshHit == startNumOfHits):# and dataMatrix[verticalPosition,horizontalPosition] > 0):
                        eventStart = verticalPosition - startNumOfHits +1
                        eventEnd = verticalPosition
                        isEvent = True
    
            else: 
                # we are in an event!
                eventEnd += 1
                # does it end?
                if (collumn[eventStart] > collumn[verticalPosition]): # check if it ends
                    isEvent = False
                    threshHit = 0
                    if(eventEnd-eventStart > minEventLength):
                        theseTransients.append(Transient(collumn[eventStart:eventEnd],eventStart,eventEnd,numOfEntries))
            
        if(isEvent):
            # event lasted until end of the data, correct the last bit
            #theseEventEndCoordinates.append([eventStart,eventEnd])
            theseTransients.append(Transient(collumn[eventStart:eventEnd],eventStart,eventEnd,numOfEntries))
        '''
        if (len(theseEventEndCoordinates)== 0):
            delList.append(horizontalPosition)
        else:
            eventEndCoordinates.append(theseEventEndCoordinates[:])
        '''
        #eventEndCoordinates.append(theseEventEndCoordinates[:])
        transients.append(theseTransients)
        collumn += shiftValue # shift of data to avoid negatives
        #if (theseEventEndCoordinates):axarr2
            #collumn = discardNonEvent(collumn, theseEventEndCoordinates,baseLineArray[horizontalPosition])
    return (transients, correctionMatrix,slopeDistributions);

def thresholdEventDetect(dataMatrix, quantileWidth, emptyplaceholder):
    '''
    alternative to eventDetect, uses a simple threshold (baseline mean + 5 standard deviations of baseline) for transient detection
    '''
    numOfEntries = dataMatrix.shape[1]
    eventEndCoordinates = []
    numOfVals = len(dataMatrix)
    correctionMatrix = zeros([numOfVals,numOfEntries])
    minEventlength = 50
    for position1,collumn in enumerate(dataMatrix.T):
        theseEventEndCoordinates = []
        eventBool = True
        collumn = collumn.T
        base = detectBaseline(collumn,200)[0]
        baseMean = mean(base)
        baseStd = std(base)
        threshold = baseStd*5 + baseMean
        eventStart = 0
        eventEnd = 0
        for position2,value in enumerate(collumn):
            if (quantileWidth > 0):
                if(position2+quantileWidth < numOfVals):
                    correctionValue = percentile(collumn[position2:position2+quantileWidth],8)
                else:
                    correctionValue = percentile(collumn[numOfVals-quantileWidth:],8)
                    value -= correctionValue
                    correctionMatrix[position1,position2] = correctionValue
            if eventBool: # we are not in an event, check if one is starting
                if threshold < value:# yes: below, no: continue
                    eventBool = False
                    eventStart = position2
            else:
                if value < baseMean+ baseStd and minEventlength < position2 - eventStart:
                    eventBool = True
                    eventEnd = position2
                    theseEventEndCoordinates.append([eventStart,eventEnd])
        if(not eventBool):
            theseEventEndCoordinates.append([eventStart,numOfVals-1])
        eventEndCoordinates.append(theseEventEndCoordinates)
    return(dataMatrix,eventEndCoordinates,emptyplaceholder, correctionMatrix)

def quantileCorrect (inputArray, eventCoordinates, windowSize):
    ''' 
    DEPRECATED - non-transient data is ignored anyway, quantile correction is handled by eventDetect
    takes input array and coordinates of transients, quantile corrects non-transient data 
    '''
    correctionArray = zeros(len(inputArray))
    if(len(eventCoordinates) == 0):
        return(inputArray)
    dataSize = len(inputArray)
    priorEnd = 0
    iterationStart = 0
    iterationEnd = len(eventCoordinates)     
    if (eventCoordinates[0][0] == 0): 
        priorEnd = eventCoordinates[0][1]
        iterationStart += 1 
        
    for j in range(iterationStart,iterationEnd): # iterates over num of events
        posteriorStart = eventCoordinates[j][0] 
        posteriorEnd = eventCoordinates[j][1] 
        
        for k in range(priorEnd,posteriorStart): # from last event end to next event start    
            if (k + windowSize >= posteriorStart): # check if there is an event in the k + windowSize range (which we have to wrap around)

                if(posteriorEnd + windowSize > dataSize): # check if we have a conflict with the end of the data
                    correctionValue = percentile(inputArray[k:posteriorStart], 8)
                    inputArray[k] -= correctionValue
                    correctionArray[k] = correctionValue                       
                else: # we have no conflict with end of data, wrap around event normally
                    concatSlice = concatenate([inputArray[k:posteriorStart],inputArray[posteriorEnd:posteriorEnd+windowSize - posteriorStart + k]])
                    correctionValue = percentile(concatSlice, 8)
                    inputArray[k] -=correctionValue
                    correctionArray[k] = correctionValue
                    
            else: # no event in range
                if(k+windowSize > dataSize): # are we at the end of the data?
                    correctionValue = percentile(inputArray[dataSize-k:], 8)
                    inputArray[k] -= correctionValue
                    correctionArray[k] = correctionValue
                    
                else: # no, just point correct
                    correctionValue = percentile(inputArray[k:k+windowSize],8)
                    inputArray[k] -= correctionValue
                    correctionArray[k] = correctionValue
        priorEnd = posteriorEnd
        
    if (eventCoordinates[-1][1] != dataSize-1): # no event covering the end, need to norm last portion of data
        for k in range(eventCoordinates[-1][1],dataSize):
            correctionValue = percentile(inputArray[k:k+windowSize], q=8)
            inputArray[k] -= correctionValue
            correctionArray[k] = correctionValue         
    return (inputArray,correctionArray);

def quantileNorm (inputArray, eventCoordinates, windowSize):
    ''' Alternative version of quantileCorrect which does not skip over events- not in use, eventDetect does this anyway
    '''
    if(len(eventCoordinates) == 0):
        return(inputArray)
    position = 0
    correctionArray = zeros(len(inputArray))
    
    for event in eventCoordinates:
        while(position < event[0]):
            correctionValue = percentile(inputArray[position:position+windowSize], 8)
            inputArray[position] -= correctionValue
            correctionArray[position] = (correctionValue)
            position += 1
        position = event[1]
    if (eventCoordinates[-1][1] != len(inputArray)-1):
        #thisSlice = inputArray[len(inputArray)-windowSize:len(inputArray-1)]
        for i in range(eventCoordinates[-1][1],len(inputArray)):
            correctionValue = percentile(inputArray[i:i+windowSize], 8)
            inputArray[i] -= correctionValue
            correctionArray[position] = (correctionValue)
    return(inputArray,correctionArray);


def writeOut(data, fileName):
    '''
    Small method to print matrix with "meanN"- header, as found in the .cls files
    '''
    numOfCols = (data.shape[1])
    headerString = ''
    for i in range(numOfCols):
        headerString += 'mean' + str(i) + '\t'
    savetxt(fileName, data, fmt='%10.9f',delimiter="\t",header=headerString)


def pushToBaseline (dataMatrix, bucketSize):
    '''
    @param dataMatrix: matrix of the data whose baseline we're looking for
    @param bucketSize: size of the baseline bins
    @return: baseline corrected version of the data, currently also expressed data in terms of sigma
    '''
    coordList = []
    for roi in nditer(dataMatrix,order='F',flags=['external_loop'],op_flags=['readwrite']): # iterate over number of ROIs in file, e.g. collumn-wise
        baseLineArray, coordinates= detectBaseline(roi, bucketSize)
        meanVal = abs(mean(baseLineArray))
        coordList.append(coordinates)
        roi[...] = roi/meanVal -1 # baseline correction 
            #dataMatrix[:,i] = dataMatrix[:,i]/(std(baseLineArray))
    return (dataMatrix,coordList);

def pushToBaseLine2(inputMatrix, baselineArray):
    numOfVals = inputMatrix.shape[1]
    for i in range(0,numOfVals): # iterate over number of ROIs in file, e.g. collumn-wise
        meanVal = abs(mean(inputMatrix[baselineArray[i][0]:baselineArray[i][1]])) # calculate mean over baseline area, slicing in python treats arrays as gaulois-fields, so index out of bounds doesn't occur
        try:
            inputMatrix[:,i] = (inputMatrix[:,i]-meanVal)/meanVal # baseline correction
        except ZeroDivisionError:
            print(str(meanVal) + " at position: " + str(i))
            exit()
    return(inputMatrix)

def detectBaseline (data, bucketSize):
    '''
    @param data: the array out of which the region with the lowest noise is to be identified
    @param bucketSize: size of the bins to be checked
    @return: bin with the lowest noise and its starting coordinate, in a tuple
    '''
    numOfEntries = len(data)
    lowestSigma = maxint # for size comparasion
    baselineArray = zeros(bucketSize)
    coordinate = []
    for j in range(0,int(numOfEntries-bucketSize),int(numOfEntries/(bucketSize*2))): # iterate over 1 out of 10 bucket positions
        thisStd = std(data[j:j+bucketSize])#(axisStd(data[j:j+bucketSize])) # current deviation
        if (thisStd < lowestSigma): # new min deviation found
            lowestSigma = thisStd
            coordinate = (j,j+bucketSize)
        baselineArray = data[coordinate[0]:coordinate[1]]
    return(baselineArray, coordinate)

def detectPeakline (data, bucketSize):
    '''
    @param data: the array out of which the region with the highest noise is to be identified
    @param bucketSize: size of the bins to be checked
    @return: bin with the highest noise and its starting coordinate, in a tuple
    '''
    coordinate = 0
    numOfEntries = len(data)
    highestSigma = 0 # for size comparasion
    peaklineArray = zeros(bucketSize) # stores portion of data which we consider the peak line
    # print(str(data) + " " + str(len(data)))
    for j in range(0,int(numOfEntries-bucketSize)): # iterate over all possible bucket positions
        thisStd = (var(data[j:j+bucketSize])) # current deviation
        if (thisStd > highestSigma): # new min deviation found  
            highestSigma = thisStd 
            peaklineArray = data[j:j+bucketSize]
            coordinate = j
    return(peaklineArray,coordinate)

def maxAmp(inputData):
    max_index = argmax(inputData)
    return(inputData[max_index], max_index);

def detectNumOfPeaks(data):
    if(len(data) <2):
        return(1)
    numOfPeaks = 0
    data = gaussian_filter(data, 9)
    oldSign = sign(meanSlope([data[0],data[1]]))
    
    for i in range(len(data)-1):
        newSign = sign(meanSlope([data[i],data[i+1]]))
        if (newSign != oldSign):
            numOfPeaks += 1
        oldSign = newSign
        
    numOfPeaks = max([int(ceil(numOfPeaks/2)),1])
    return(numOfPeaks);

def getSinglePeakTransient(transients):
    lowestAmplitude = maxint
    lowestPeakNum = maxint
    peakProperties = []
    
    for t in transients:
        if (t.getNumOfPeaks() > lowestPeakNum):
            continue
        if (t.getAmplitude() < lowestAmplitude):
            lowestAmplitude = t.getAmplitude
            peakProperties = t
    return(peakProperties)

def createMeanKernel(transientMatrix):
    #Function is called for one ROI, iterates over the transients. Problem: requires information from ALL ROI's and transients for computation.

    singlePeakArray = [i for j in transientMatrix for i in j if i.getNumOfPeaks()==1 ]
    risetimeArray = [i.getRiseTime() for i in singlePeakArray]
    decaytimeArray = [i.getDecayTime() for i in singlePeakArray]
    amplitudeArray = [i.getAmplitude() for i in singlePeakArray]
    
    risetimeStd = std(risetimeArray)
    risetimeMean = mean(risetimeArray)
    decaytimeStd = std(decaytimeArray)
    decaytimeMean = mean(decaytimeArray)
    amplitudeStd = std(amplitudeArray)
    amplitudeMean = mean(amplitudeArray)
    risetimeThreshUp = risetimeStd*2 + risetimeMean
    decaytimeThreshUp = decaytimeStd*2 + decaytimeMean 
    amplitudeThreshUp = amplitudeStd*2 + amplitudeMean
    risetimeThreshDown =  risetimeMean - risetimeStd
    decaytimeThreshDown = decaytimeMean - decaytimeStd  
    amplitudeThreshDown =  amplitudeMean -amplitudeStd*2 
    print ("amp: %f\t%f\nRT: %f\t%f \nDT: %f\t%f" % (amplitudeThreshDown,amplitudeThreshUp,risetimeThreshDown,risetimeThreshUp,decaytimeThreshDown,decaytimeThreshUp))
    # block below: mean transients falling within 2 STD of 1 peak transient mean
    import matplotlib.pyplot as plt
    meanArray = zeros(0,dtype=float)
    divideArray = zeros(transientMatrix[0][0].getNumOfFrames(),dtype=int)
    constrainRelaxVariable = 1
    while (not meanArray.any()):
        for i in singlePeakArray: # go over ROIs
            if (i.getNumOfPeaks() == 1 and amplitudeThreshDown/constrainRelaxVariable < i.getAmplitude() < amplitudeThreshUp*constrainRelaxVariable and 
                    decaytimeThreshDown/constrainRelaxVariable < i.getDecayTime() < decaytimeThreshUp*constrainRelaxVariable and 
                    risetimeThreshDown/constrainRelaxVariable < i.getRiseTime() < risetimeThreshUp*constrainRelaxVariable): # we already know which ones are the single peak transients, but re-checking is cheaper than a new array
                    
                    currentTransient = i.getData()
                    plt.plot(currentTransient,"grey",alpha=0.2)
                    if(len(currentTransient)>len(meanArray)):
                        meanArray.resize(currentTransient.shape)
                        tempArray = ones(currentTransient.shape,dtype=int)
                        tempArray.resize(divideArray.shape,refcheck=False)
                        divideArray += tempArray  
                        meanArray += currentTransient
                    else:
                        tempArray = ones(currentTransient.shape,dtype=int)
                        tempArray.resize(divideArray.shape)
                        divideArray += tempArray
                        zerosArray = zeros(len(meanArray)-len(currentTransient))
                        currentTransient = concatenate((currentTransient,zerosArray))
                        meanArray += currentTransient
        if (not meanArray.any()):
            print "No single peak transients within 2 amplitude/ rise time/ decay time standard deviations of each other found!\nRelaxing constraints by 50% for deconvolution"
            constrainRelaxVariable += 0.5
    # loop up top: finished adding single peak transients, need to divide by our divideArray (after cropping it accordingly)
    
    divideArray.resize(len(meanArray))
    meanArray = meanArray / (divideArray)
    tVariable = arange(0,len(meanArray),1.)
    #tVariable.resize(meanArray.shape)
    popt, pcov = curve_fit(alphaKernel, tVariable,meanArray, p0=[80,100,50],maxfev=10000)
    plt.plot(tVariable, meanArray, label="meanArray")
    plt.plot(alphaKernel(tVariable, *popt), 'r-', label='fit')
    plt.savefig('alphaKernel.png')
    plt.close()
    return(alphaKernel(arange(i.getNumOfFrames()),*popt))

def alphaKernel (t, A,t_A,t_B):
    return A*(exp(-t/t_A)-exp(-t/t_B))

def deconvolve(transients, kernel):
    '''
    IDEA: using kernel (e.g. exp(-t/tau)), model corrected fluorescence observations as exp. decay with
    point-wise increase for each spike
    '''
    if (not transients):
        return(zeros(kernel.size))
    spikeSignal = zeros(transients[0].getNumOfFrames())
    for i in (transients):
        transientData = i.getData()
        
        if transientData.size > kernel.size:
            thisKernel = append(kernel,zeros(transientData.size-kernel.size))
        else:
            thisKernel = resize(kernel,transientData.shape)
        #spikeTrain = fft.ifft(fft.fft(transientData)/ fft.fft(transientKernel))# with ifft = inverse fast fourier transform and fft = fast fourier transform
        spikeSignal[i.getStartTime():i.getEndTime()] = (fft.ifft(divide(fft.fft(transientData), fft.fft(thisKernel))))
        spikeSignal[i.getStartTime():i.getEndTime()] = generateSpiketrainFromSignal(spikeSignal[i.getStartTime():i.getEndTime()])
    return(spikeSignal)

def generateSpiketrainFromSignal(spikeSignal):
    #TODO: change method so it assignes ceil(sum(spikeSignal)) amount of APs to peaks of deconvo signal
    
    spikeTrain = zeros(spikeSignal.size)
    numOfSpikes = ceil(sum(spikeSignal))
    for i in arange(numOfSpikes):
        maxVal = argmax(spikeSignal)
        spikeTrain[maxVal] = 1
        spikeSignal[maxVal] = 0
        
    return(spikeTrain)

def negativeTransientsCorrect(transients):
    #TODO: this.
    numOfNeg = 0
    meanTrans = ()
    for i in transients:
        meanTrans.append(mean(i.getData))
        if mean(i.getData)< 0:
            numOfNeg += 1
    for j in range(numOfNeg):
        transients.remove(argmin(meanTrans))
    return(transients)

def aucCorrect(transients, kernel):
    '''
    naive approach, less accurate for overlapping activity
    '''
    kernelAUC = trapz(kernel)
    for i in (transients):
        transientAUC = trapz(i.getData())
        #spikeTrain = fft.ifft(fft.fft(transientData)/ fft.fft(transientKernel))# with ifft = inverse fast fourier transform and fft = fast fourier transform
        if (transientAUC/kernelAUC < 0.1):
            transients.remove(i)
    return(transients)

class Transient:
    '''
    class to store transient properties, such as rise time, decay time, amplitude and duration
    '''
    def __init__(self, data, startTime,endTime, numOfFrames):
        
        maxVal = maxAmp(data)
        self.isLast = False
        self.startTime = startTime
        self.endTime = endTime
        self.std = std(data)
        self.amplitude = float(maxVal[0]) 
        self.riseTime = maxVal[1]
        self.decayTime = len(data) + 1 - self.riseTime
        self.totalTime = len(data)
        self.numOfPeaks = detectNumOfPeaks(data)
        self.data = data
        self.numOfFrames = numOfFrames
    
    def getData(self):
        return(self.data)
    def getNumOfFrames(self):
        return(self.numOfFrames)
    def getCoordinates(self):
        return[self.startTime,self.endTime]
    
    def getStartTime(self):
        return(self.startTime)
    def getEndTime(self):
        return(self.endTime)
    def getRiseTime(self):
        return self.riseTime;
        
    def getDecayTime(self):
        return self.decayTime;
    def getNumOfPeaks(self):
        return(self.numOfPeaks)
    def getStd(self):
        return(self.std);
    
    def getIsLast(self):
        return self.isLast;
    def setIsLast(self):
        self.isLast = True
        
    def getAmplitude(self):
        return self.amplitude;
    
    def getTotalTime(self):
        return (self.totalTime);

import numpy as np
import scipy.signal 
import scipy.linalg

from warnings import warn    
try:
    from cvxopt import matrix, spmatrix, spdiag, solvers
    import picos
except ImportError:
    raise ImportError('Constrained Foopsi requires cvxopt and picos packages.')


#%%
def constrained_foopsi(fluor, 
                     b = None, 
                     c1 = None, 
                     g = None, 
                     sn = None, 
                     p= 2, 
                     method = 'cvx', 
                     bas_nonneg = True, 
                     noise_range = [.25,.5],
                     noise_method = 'logmexp',
                     lags = 5, 
                     resparse = 0,
                     fudge_factor = 1, 
                     verbosity = False):

    """
    Infer the most likely discretized spike train underlying a fluorescence
    trace, using a noise constrained deconvolution approach
    Inputs
    ----------
    fluor   : nparray
        One dimensional array containing the fluorescence intensities with
        one entry per time-bin.
    b       : float, optional
        Fluorescence baseline balue. If no value is given, then b is estimated 
        from the data
    c1      : 
    g       : float, optional
        Parameters of the AR process that models the fluorescence impulse response.
        Estimated from the data if no value is given
    sn      : float, optional
        Standard deviation of the noise distribution.  If no value is given, 
        then sn is estimated from the data.        
    options : dictionary
        list of user selected options (see more below)
    'p'             :         2, # AR order 
    'method'        :     'cvx', # solution method (no other currently supported)
    'bas_nonneg'    :      True, # bseline strictly non-negative
    'noise_range'   :  [.25,.5], # frequency range for averaging noise PSD
    'noise_method'  : 'logmexp', # method of averaging noise PSD
    'lags'          :         5, # number of lags for estimating time constants
    'resparse'      :         0, # times to resparse original solution (not supported)
    'fudge_factor'  :         1, # fudge factor for reducing time constant bias
    'verbosity'     :     False, # display optimization details
    
    Returns
    -------
    c            : ndarray of float
        The inferred denoised fluorescence signal at each time-bin.
    b, c1, g, sn : As explained above
    sp           : ndarray of float
        Discretized deconvolved neural activity (spikes)
    
    References
    ----------
    * Pnevmatikakis et al. 2015. Submitted (arXiv:1409.2903).
    * Machado et al. 2015. Cell 162(2):338-350
    """

    
    if g is None or sn is None:        
        g,sn = estimate_parameters(fluor, p=p, sn=sn, g = g, range_ff=noise_range, method=noise_method, lags=lags, fudge_factor=fudge_factor)

  
    
    T = len(fluor)
    # construct deconvolution matrix  (sp = G*c) 
    G = spmatrix(1.,range(T),range(T),(T,T))

    for i in range(p):
        G = G + spmatrix(-g[i],np.arange(i+1,T),np.arange(T-i-1),(T,T))
        
    gr = np.roots(np.concatenate([np.array([1]),-g.flatten()])) 
    gd_vec = np.max(gr)**np.arange(T)  # decay vector for initial fluorescence
    gen_vec = G * matrix(np.ones(fluor.size))  
    
    # Initialize variables in our problem
    prob = picos.Problem()
    
    # Define variables
    calcium_fit = prob.add_variable('calcium_fit', fluor.size)    
    cnt = 0
    if b is None:
        flag_b = True
        cnt += 1
        b = prob.add_variable('b', 1)
        if bas_nonneg:
            b_lb = 0
        else:
            b_lb = np.min(fluor)
            
        prob.add_constraint(b >= b_lb)
    else:
        flag_b = False

    if c1 is None:
        flag_c1 = True
        cnt += 1
        c1 = prob.add_variable('c1', 1)
        prob.add_constraint(c1 >= 0)
    else:
        flag_c1 = False
    
    # Add constraints    
    prob.add_constraint(G * calcium_fit >= 0)
    #this line takes FOREVER. what's going on here?
    res = abs(matrix(fluor.astype(float)) - calcium_fit - b*matrix(np.ones(fluor.size)) - matrix(gd_vec) * c1)
    prob.add_constraint(res < sn * np.sqrt(fluor.size))
    prob.set_objective('min', calcium_fit.T * gen_vec)
    
    # solve problem
    try:
        prob.solve(solver='mosek', verbose=verbosity)
        sel_solver = 'mosek'
#        prob.solve(solver='gurobi', verbose=verbosity)
#        sel_solver = 'gurobi'
    except ImportError:
        warn('MOSEK is not installed. Spike inference may be VERY slow!')
        sel_solver = []
        prob.solver_selection()
        prob.solve(verbose=verbosity)
        
    # if problem in infeasible due to low noise value then project onto the cone of linear constraints with cvxopt
    if prob.status == 'prim_infeas_cer' or prob.status == 'dual_infeas_cer' or prob.status == 'primal infeasible':
        warn('Original problem infeasible. Adjusting noise level and re-solving')   
        # setup quadratic problem with cvxopt        
        solvers.options['show_progress'] = verbosity
        ind_rows = range(T)
        ind_cols = range(T)
        vals = np.ones(T)
        if flag_b:
            ind_rows = ind_rows + range(T) 
            ind_cols = ind_cols + [T]*T
            vals = np.concatenate((vals,np.ones(T)))
        if flag_c1:
            ind_rows = ind_rows + range(T)
            ind_cols = ind_cols + [T+cnt-1]*T
            vals = np.concatenate((vals,np.squeeze(gd_vec)))            
        P = spmatrix(vals,ind_rows,ind_cols,(T,T+cnt))
        H = P.T*P
        Py = P.T*matrix(fluor.astype(float))
        sol = solvers.qp(H,-Py,spdiag([-G,-spmatrix(1.,range(cnt),range(cnt))]),matrix(0.,(T+cnt,1)))
        xx = sol['x']
        c = np.array(xx[:T])
        sp = np.array(G*matrix(c))
        c = np.squeeze(c)
        if flag_b:
            b = np.array(xx[T+1]) + b_lb
        if flag_c1:
            c1 = np.array(xx[-1])
        sn = np.linalg.norm(fluor-c-c1*gd_vec-b)/np.sqrt(T)   
    else: # readout picos solution
        c = np.squeeze(calcium_fit.value)
        sp = np.squeeze(np.asarray(G*calcium_fit.value))        
        if flag_b:    
            b = np.squeeze(b.value)        
        if flag_c1:    
            c1 = np.squeeze(c1.value)                    

    return c,b,c1,g,sn,sp


def estimate_parameters(fluor, p = 2, sn = None, g = None, range_ff = [0.25,0.5], method = 'logmexp', lags = 5, fudge_factor = 1):
    """
    Estimate noise standard deviation and AR coefficients if they are not present
    """
    
    if sn is None:
        sn = GetSn(fluor,range_ff,method)
        
    if g is None:
        g = estimate_time_constant(fluor,p,sn,lags,fudge_factor)

    return g,sn

def estimate_time_constant(fluor, p = 2, sn = None, lags = 5, fudge_factor = 1):
    """    
    Estimate AR model parameters through the autocovariance function    
    Inputs
    ----------
    fluor        : nparray
        One dimensional array containing the fluorescence intensities with
        one entry per time-bin.
    p            : positive integer
        order of AR system  
    sn           : float
        noise standard deviation, estimated if not provided.
    lags         : positive integer
        number of additional lags where he autocovariance is computed
    fudge_factor : float (0< fudge_factor <= 1)
        shrinkage factor to reduce bias
        
    Return
    -----------
    g       : estimated coefficients of the AR process
    """    
    

    if sn is None:
        sn = GetSn(fluor)
        
    lags += p
    xc = axcov(fluor,lags)        
    xc = xc[:,np.newaxis]
    
    A = scipy.linalg.toeplitz(xc[lags+np.arange(lags)],xc[lags+np.arange(p)]) - sn**2*np.eye(lags,p)
    g = np.linalg.lstsq(A,xc[lags+1:])[0]
    if fudge_factor < 1:
        gr = fudge_factor*np.roots(np.concatenate([np.array([1]),-g.flatten()]))
        gr = (gr+gr.conjugate())/2
        gr[gr>1] = 0.95
        gr[gr<0] = 0.15
        g = np.poly(gr)
        g = -g[1:]        
        
    return g.flatten()
    
def GetSn(fluor, range_ff = [0.25,0.5], method = 'logmexp'):
    """    
    Estimate noise power through the power spectral density over the range of large frequencies    
    Inputs
    ----------
    fluor    : nparray
        One dimensional array containing the fluorescence intensities with
        one entry per time-bin.
    range_ff : (1,2) array, nonnegative, max value <= 0.5
        range of frequency (x Nyquist rate) over which the spectrum is averaged  
    method   : string
        method of averaging: Mean, median, exponentiated mean of logvalues (default)
        
    Return
    -----------
    sn       : noise standard deviation
    """
    

    ff, Pxx = scipy.signal.welch(fluor)
    ind1 = ff > range_ff[0]
    ind2 = ff < range_ff[1]
    ind = np.logical_and(ind1,ind2)
    Pxx_ind = Pxx[ind]
    sn = {
        'mean': lambda Pxx_ind: np.sqrt(np.mean(Pxx_ind/2)),
        'median': lambda Pxx_ind: np.sqrt(np.median(Pxx_ind/2)),
        'logmexp': lambda Pxx_ind: np.sqrt(np.exp(np.mean(np.log(Pxx_ind/2))))
    }[method](Pxx_ind)

    return sn

def axcov(data, maxlag=5):
    """
    Compute the autocovariance of data at lag = -maxlag:0:maxlag
    Parameters
    ----------
    data : array
        Array containing fluorescence data
    maxlag : int
        Number of lags to use in autocovariance calculation
    Returns
    -------
    axcov : array
        Autocovariances computed from -maxlag:0:maxlag
    """
    
    data = data - np.mean(data)
    T = len(data)
    bins = np.size(data)
    xcov = np.fft.fft(data, np.power(2, nextpow2(2 * bins - 1)))
    xcov = np.fft.ifft(np.square(np.abs(xcov)))
    xcov = np.concatenate([xcov[np.arange(xcov.size - maxlag, xcov.size)],
                           xcov[np.arange(0, maxlag + 1)]])
    #xcov = xcov/np.concatenate([np.arange(T-maxlag,T+1),np.arange(T-1,T-maxlag-1,-1)])
    return np.real(xcov/T)
    
def nextpow2(value):
    """
    Find exponent such that 2^exponent is equal to or greater than abs(value).
    Parameters
    ----------
    value : int
    Returns
    -------
    exponent : int
    """
    
    exponent = 0
    avalue = np.abs(value)
    while avalue > np.power(2, exponent):
        exponent += 1
    return exponent 