'''
Created on Oct 3, 2016
script contains methods used for main.py
@author: maximilian
'''
from numpy import percentile,zeros,var, mean, arange, concatenate, savetxt, sign, ceil,\
 std, around, fft, divide, abs, ones, resize, exp, nditer,trapz, argmax,\
    append, sum, argmin, matrix, real
from sys import maxint
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import norm
from scipy.optimize import curve_fit
import pandas, xlrd
from xlrd.biffh import XLRDError
from scipy.stats.stats import mode

def _importMatrix(filename,valSeparator):
    '''
    imports csv or xls formatted file, returns ndArray 
    @input: filename: input path; valSeperator: seperator between entries, treats file as xlrd if none
    @return: numOfVals: number of collumns in file (cells); dataMatrix: ndArray containing file data (dtype float)
    '''
    try:
        if(valSeparator):
            frame = pandas.read_csv(filename,sep = valSeparator,engine='python',skipinitialspace=True)
            numOfRois = frame.shape[1]
            dataMatrix = frame.as_matrix()
            dataMatrix = matrix(dataMatrix,dtype=float).A
        else:
            try:
                dataMatrix = []
                book = xlrd.open_workbook(filename)
                sh = book.sheet_by_index(0)
                numOfRois = sh.ncols
                for rx in range(sh.nrows):
                    dataMatrix.append(map(float,sh.row(rx)))
                dataMatrix = matrix(dataMatrix).A
            except XLRDError:
                print("Wrong format! The file might contain comma-seperated values. Try using the -sep argument (E.G. -sep '\\t'")
                exit()
    except IOError:
        print(str(filename) + " is not a valid file path !")
        exit()
    
    return(dataMatrix, numOfRois)

'''
@param row: current row of csvfile
@return: List of length 'numOfVals', containing the entry values row-vise
actually yields one entry less than numOfVals - ImageJ result.xls file headers start with a whitespace (value seperator)
'''
def _getValues (row):
    thisRow = []
    for i in row:
        thisRow.append(float(i.value))
    return thisRow;

def _meanSlope(inArray):
    '''
    Given a time series, return the slope (which only depends on the start and end point
    '''
    iterNum = len(inArray)
    avgSlope = (inArray[-1] - inArray[0])/iterNum
    return (float(avgSlope))

def eventDetect (dataMatrix, quantileWidth,slopeWidth):
    '''
    
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
    transients = [] # list of lists of transients
    numOfVals = dataMatrix.shape[1] # number of ROIs in file
    numOfEntries = len(dataMatrix) # number of entries in file
    
    startNumOfHits = 4 # number of hits necessary for an event start to be declared 
    minEventLength = 15 # number of data points necessary above start level for an event to be kept as such
    
    thresholdList = []
    slopeDistributions = []
    correctionMatrix = zeros([numOfEntries,numOfVals]) # placeholder for correction visualization - doesn't do anything
    
    # here we generate a global noise / slope threshold for all ROIs of the file
    for collumn in nditer(dataMatrix, order='F',flags=['external_loop'],op_flags=['readonly']): 
        slopeList = zeros(numOfEntries - slopeWidth)
        for i in range(len(slopeList)):
            slopeList[i]= (_meanSlope(collumn[i:i+slopeWidth]))  # threshold to be passed for start of significance
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
    # re-iterate over the data ROI-wise for determination of event positions - in reverse order because we delete skipped rows immediately
    for horizontalPosition, collumn in enumerate(dataMatrix.T): # generate local noise / slope threshold and consider if we should skip
        baselineMean = mean(_detectBaseline(collumn, 500)[0]) # we use this to end transients earlier - and have a new one start if activity continues
        collumn = collumn.T # looks stupid but makes the general code more pythonic (instead of referencing collumns via slicing)
        thisSlopeThreshold = thresholdList[horizontalPosition]# threshold to be passed for start of significance 
        isEvent = False # bool to remember wether we are in event-mode
        eventStart = 0
        eventEnd = 0
        threshHit = 0
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
            thisSlope = abs(_meanSlope(thatSlice))   
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
                if (collumn[eventStart] > collumn[verticalPosition] or collumn[verticalPosition] < baselineMean): # check if it ends
                    isEvent = False
                    threshHit = 0
                    if(eventEnd-eventStart > minEventLength):
                        theseTransients.append(_Transient(collumn[eventStart:eventEnd],eventStart,eventEnd,numOfEntries))
            
        if(isEvent):
            # event lasted until end of the data, correct the last bit
            #theseEventEndCoordinates.append([eventStart,eventEnd])
            theseTransients.append(_Transient(collumn[eventStart:eventEnd],eventStart,eventEnd,numOfEntries))

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
        base = _detectBaseline(collumn,200)[0]
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

def _quantileCorrect (inputArray, eventCoordinates, windowSize):
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

def _quantileNorm (inputArray, eventCoordinates, windowSize):
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


def _writeOut(data, fileName):
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
    for roi in nditer(dataMatrix,order='F',flags=['external_loop','refs_ok'],op_flags=['readwrite']): # iterate over number of ROIs in file, e.g. collumn-wise
        baseLineArray, coordinates= _detectBaseline(roi, bucketSize)
        meanVal = abs(mean(baseLineArray))
        coordList.append(coordinates)
        roi[...] = roi/meanVal -1 # baseline correction 
            #dataMatrix[:,i] = dataMatrix[:,i]/(std(baseLineArray))
    return (dataMatrix,coordList);

def _detectBaseline (data, bucketSize):
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

def _detectPeakline (data, bucketSize):
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

def _maxAmp(inputData):
    max_index = argmax(inputData)
    return(inputData[max_index], max_index);

def _detectNumOfPeaks(data):
    if(len(data) <10):
        return(1)
    numOfPeaks = 0
    steps = len(data)/10
    data2 = [sum(data[i:i+steps-1])/steps for i in range(len(data)-steps)]
    data = gaussian_filter(data2, 9)    
    prevPoint = 0
    oldSign = 1
    for i in data:
        newSign = sign(i- prevPoint)
        if (newSign != oldSign):
            numOfPeaks += 1
        oldSign = newSign
        prevPoint = i
    return(numOfPeaks);

def createMeanKernel(transientMatrix):
    '''
    @transientMatrix: list (or list of lists) of transients
    @return: alpha kernel of length n
    The transients fed into this function should come from comparable traces,
     i.e. same cell type / experimenting run.
    It selects single-peaked transients and then further retains only those 
    within the [0.2,0.8]-percentile of rise time, decay time,
    amplitude and area under the curve.
    The resulting mean-1-peak-transient is used as a template to fit an 
    alpha distribution to, which is suitable for deconvolution.
    '''
    
    if (any(isinstance(el, list) for el in transientMatrix)):
        singlePeakArray = [i for j in transientMatrix for i in j if i.numOfPeaks==1]
    else:
        singlePeakArray = [i for i in transientMatrix if i.numOfPeaks ==1]
        
    risetimeArray = [i.riseTime for i in singlePeakArray]
    decaytimeArray = [i.decayTime for i in singlePeakArray]
    amplitudeArray = [i.amplitude for i in singlePeakArray]
    aucArray = [i.auc for i in singlePeakArray]
    
    aucThreshUp = percentile(aucArray, 80)
    risetimeThreshUp = percentile(risetimeArray, 80)
    decaytimeThreshUp = percentile(decaytimeArray, 80)
    amplitudeThreshUp = percentile(amplitudeArray, 80)
    aucThresDown = percentile(aucArray, 20)
    risetimeThreshDown =  percentile(risetimeArray, 20)
    decaytimeThreshDown = percentile(decaytimeArray, 20)
    amplitudeThreshDown =  percentile(amplitudeArray, 20)
    print ("amp: %f\t%f\nRT: %f\t%f \nDT: %f\t%f" % (amplitudeThreshDown,amplitudeThreshUp,risetimeThreshDown,risetimeThreshUp,decaytimeThreshDown,decaytimeThreshUp))
   
    import matplotlib.pyplot as plt
    meanArray = zeros(0,dtype=float)
    divideArray = zeros(transientMatrix[0][0].numOfFrames,dtype=int)
    for i in singlePeakArray: # go over ROIs
        if (amplitudeThreshDown < i.amplitude < amplitudeThreshUp and 
                decaytimeThreshDown < i.decayTime < decaytimeThreshUp and 
                risetimeThreshDown < i.riseTime < risetimeThreshUp and
                aucThresDown < i.auc < aucThreshUp): # we already know which ones are the single peak transients, but re-checking is cheaper than a new array
            #TODO: should I keep the AUC check or not? more testing -> get the full data set in here
            currentTransient = i.data
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
    # loop up top: finished adding single peak transients, need to divide by our divideArray (after cropping it accordingly)
    
    divideArray.resize(len(meanArray))
    meanArray = meanArray / (divideArray)
    tVariable = arange(0,len(meanArray),1.)
    #tVariable.resize(meanArray.shape)
    popt, pcov = curve_fit(_alphaKernel, tVariable,meanArray, p0=[80,100,50],maxfev=100000)
    plt.plot(tVariable, meanArray, label="meanArray")
    plt.plot(_alphaKernel(tVariable, *popt), 'r-', label='fit')
    plt.savefig('/home/maximilian/unistuff/paris_ens/cal_neuroim/testData/alphaKernel.png')
    plt.close()
    return(_alphaKernel(arange(i.numOfFrames),*popt))

def _alphaKernel (t, A,t_A,t_B):
    return A*(exp(-t/t_A)-exp(-t/t_B))

def deconvolve(transients, kernel):
    '''
    @transients: list or list of lists of transient-objects, like returned by eventDetect and thresholdEventDetect
    @kernel: kernel to be used for deconvolution, i.e. as returned by createMeanKernel  
    @return: 1xn or dxn array (corresponding to the input) of integers, 
    with n= number of frames, d= number of lists.
    Integer values correspond to the number of action potentials in the respective frame.
    '''
    if (not transients): # no transients found - the input is empty.
        return(zeros(kernel.size))
    
    if (any(isinstance(el, list) for el in transients)): #this is true, we have a list of lists -> iterate over lists and transients
        spikeSignal = zeros((len(transients),transients[0][0].numOfFrames))# problem: crashes if first ROI doesn't contain a transient
        for i,transientList in enumerate(transients):
            for transient in transientList:
                
                if transient.data.size > kernel.size:
                    thisKernel = append(kernel,zeros(transient.data.size-kernel.size))
                else:
                    thisKernel = resize(kernel,transient.data.shape)
                spikeSignal[i][transient.startTime:transient.endTime] = (fft.ifft(divide(fft.fft(transient.data), fft.fft(thisKernel))))
                spikeSignal[i][transient.startTime:transient.endTime] = _generateSpiketrainFromSignal(spikeSignal[i][transient.startTime:transient.endTime])
        
    else: # the above clause is false, we have a list of transients and can iterate directly
        spikeSignal = zeros(transients[0].numOfFrames)
        for transient in (transients):
            
            if transient.data.size > kernel.size:
                thisKernel = append(kernel,zeros(transient.data.size-kernel.size))
            else:
                thisKernel = resize(kernel,transient.data.shape)
            spikeSignal[transient.startTime:transient.endTime] = real(fft.ifft(divide(fft.fft(transient.data), fft.fft(thisKernel))))
            spikeSignal[transient.startTime:transient.endTime] = _generateSpiketrainFromSignal(spikeSignal[transient.startTime:transient.endTime])
        
    return(spikeSignal)

def _generateSpiketrainFromSignal(spikeSignal):
    
    spikeTrain = zeros(spikeSignal.size)
    numOfSpikes = int(ceil(sum(spikeSignal)))
    
    for i in range(numOfSpikes):
        maxVal = argmax(spikeSignal)
        spikeTrain[maxVal] = 1
        spikeSignal[maxVal] /= 2 #TODO: this doesn't work! need to flatten stuff more extremely
    return(spikeTrain)

def negativeTransientsCorrect(transients):
    #kick out fraction of positive going transients relative to number of negative going ones
    numOfNeg = 0
    for i in transients:
        meanVal = mean(i.data)
        if meanVal< 0:
            transients.remove(i)
            numOfNeg += 1
            
    for j in range(numOfNeg):
        del transients[argmin(meanVal)]
    return(transients)

def aucCorrect(transients, kernel):
    '''
    naive approach, less accurate for overlapping activity
    kicks out transients which lack surface area
    '''
    kernelAUC = trapz(kernel)
    for i in (transients):
        transientAUC = trapz(i.data)
        if (transientAUC/kernelAUC < 0.1):
            transients.remove(i)
    return(transients)

class _Transient:
    '''
    class to store transient properties, such as rise time, decay time, amplitude and duration
    '''
    def __init__(self, data, startTime,endTime, numOfFrames):
        
        self.maxVal = _maxAmp(data)
        self.startTime = startTime
        self.endTime = endTime
        self.amplitude = float(self.maxVal[0]) 
        self.riseTime = int(self.maxVal[1])
        self.decayTime = int(len(data) + 1 - self.riseTime)
        self.totalTime = len(data)
        self.numOfPeaks = _detectNumOfPeaks(data)
        self.data = data
        self.numOfFrames = numOfFrames
        self.coordinates = (self.startTime, self.endTime)
        self.auc = trapz(data)