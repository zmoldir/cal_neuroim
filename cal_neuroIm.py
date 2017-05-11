'''
Created on Oct 3, 2016
script contains methods used for main.py
@author: maximilian
'''
from numpy import percentile,zeros,var, mean, arange, concatenate, savetxt, sign, ceil,\
matrix, std, subtract, repeat, around, fft, divide, abs, ones, resize, exp, nditer
from sys import maxint
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import norm
from scipy.optimize import curve_fit
import pandas, xlrd, time
from xlrd.biffh import XLRDError
from math import sqrt
from scipy.stats.stats import mode
def importMatrix(filename,valSeperator,csvBool):
    dataMatrix = [] # TODO: it's busted. Turn this into a proper pandas routine instead of some csv-reader BS
    try:
        if(csvBool):
            frame = pandas.read_csv(filename)
            
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
    rawMatrix = matrix(dataMatrix).A
    del(dataMatrix)
    return(rawMatrix, numOfVals)

'''
@param numOfVals: number of Values per row (e.g. ROIs)
@param row: current row of csvfile
@return: Vector of length 'numOfVals', containing the entry values row-vise
actually yields one entry less than numOfVals - ImageJ result.xls file headers start with a whitespace (value seperator)
TODO: this looks extremely un-pythonic and stupid ... what am I using it for again?'''
def getValues (row):
    thisRow = []
    for i in range(0,len(row)):
        thisRow.append(float(row[i].value))
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
    startNumOfHits = 2 # number of slope-hits to be made before an event is declared as beginning
    stopNumOfHits = 2
    minEventLength = 30 # minimum length of events
    delList = []
    correctionMatrix = zeros([numOfEntries,numOfVals])
    # here we generate a global noise / slope threshold for all ROIs of the file
    for num,collumn in enumerate(dataMatrix.T): 
        slopeList = zeros(numOfEntries - slopeWidth)
        collumn = collumn.T
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
        delList.append(sigma*2) #TODO don't store this in delList
        '''
        import matplotlib.pyplot as plt
        import matplotlib.mlab as mlab
        n, bins, patches = plt.hist(slopeList, 100, normed=1, facecolor='green', alpha=0.75)
        y = mlab.normpdf( bins, mu, sigma)
        plt.plot(bins, y, 'r--', linewidth=2)
        plt.plot([sigma*2,sigma*2],[0,5],'k--')
        #plot
        plt.title(r'$\mathrm{test:}\ \mu=%.3f,\ \sigma=%.3f,\ thisMode=%.3f$' %(mu, sigma,thisMode))
        plt.grid(True)
        plt.savefig("Slope_distrib_" + str(num)+".png")
        plt.close()'''
        #noiseThreshold += (var(detectBaseline(collumn, baselineWidth)[0])) # threshold to be passed for end of significance
    #slopeThreshold = slopeThreshold/numOfVals
    # re-iterate over the data ROI-wise for determination of event positions - in reverse order because we delete skipped rows immediately
    for horizontalPosition, collumn in enumerate(dataMatrix.T): # generate local noise / slope threshold and consider if we should skip
        
        collumn = collumn.T # looks stupid but makes the general code more pythonic (instead of referencing collumns via slicing)
        thisSlopeThreshold = delList[horizontalPosition]# threshold to be passed for start of significance 
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
                        noiseHit = 0

            else: 
                # we are in an event!
                eventEnd += 1
                # does it end?
                if (collumn[eventStart] > collumn[verticalPosition]): # check if it ends
                    noiseHit += 1
                    if(noiseHit >= stopNumOfHits):
                        isEvent = False
                        threshHit = 0
                        if(eventEnd-eventStart > minEventLength):
                            #theseEventEndCoordinates.append([eventStart,eventEnd])
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
        #if (theseEventEndCoordinates):
            #collumn = discardNonEvent(collumn, theseEventEndCoordinates,baseLineArray[horizontalPosition])
    return (transients, correctionMatrix);

def thresholdEventDetect(dataMatrix, quantileWidth, emptyplaceholder):
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
    ''' Alternative version of quantileCorrect which does not skip over events
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



def arraySlopeCorrect4 (inputArray, coordinates,quantileWidth):
    ''' takes an array of values and a list of start-end coordinates of events located on said array, which are then corrected using the base line drift
    takes full transient quantile ... why?
    @param inputArray: array with events to be corrected
    @param coordinates: list of start / end coordinates of the events
    @param correctionArray: place holder to keep track of performed corrections
    @param baselineSlope: slope of the baseline used for corrections
    @return inputArray with corrections, correctionArray of the same size to visualize corrections
    '''
    #returnArray = zeros(len(inputArray))
    correctionArray = zeros(len(inputArray))
    for x in (coordinates):
        for i in range(x[1]-x[0]):
            if (x[0] + i <= len(inputArray)):
                correctionVal = percentile(inputArray[x[0]+i:x[0]+i+quantileWidth],q=8)
                inputArray[x[0]+i] = inputArray[x[0]+i] - correctionVal
                correctionArray[x[0]+i] = correctionVal
            else:
                correctionVal = percentile(inputArray[len(inputArray+1)-500:],q=8)
            inputArray[x[0]+i] = inputArray[x[0]+i] - correctionVal
            correctionArray[x[0]+i] = correctionVal
    return(inputArray, correctionArray)

def arraySlopecorrect (inputArray, coordinates, eventCorrectionValues, correctionArray):
    '''
    Modifies values by subtracting a line going from the starting value to the end value of the respective event, 
    thus keeping the shape (unlike e.g. quantile normalization)
    @param: inputArray array containing events to be corrected
    @param: coordinates 2xN list containing the start and end coordinates of said events
    @param: eventCorrectionValues events which overlap with the end or start of the array have to be corrected seperately using the values found
            at their start (end of array) or end (start of array), respectively.
    @return: array with corrected events
    '''
    for x in coordinates:
        counter = 0
        eventlen = x[1] - x[0]
        thisTransient = inputArray[x[0]:x[1]]
        slope = (thisTransient[-1]- thisTransient[0])/eventlen
        diffArray = arange(thisTransient[0],thisTransient[-1],slope)
        if(x[0] == 0):
            inputArray[x[0]+1:x[1]] -=  eventCorrectionValues[counter]
            correctionArray[x[0]+1:x[1]] = eventCorrectionValues[counter]
            counter += 1
            continue
        if(x[1] == len(inputArray) -1):
            inputArray[x[0]:x[1]+1] -=  eventCorrectionValues[counter]
            correctionArray[x[0]:x[1]+1] = eventCorrectionValues[counter]
            counter += 1
            continue
        for i in range(eventlen):
            inputArray[x[0] + i] -= diffArray[i]
            correctionArray[x[0] + i] = diffArray[i]
    return(inputArray, correctionArray)

def writeOut(data, fileName):
    '''
    Small method to print matrix with "meanN"- header, as found in the .cls files
    '''
    numOfCols = (data.shape[1])
    headerString = ''
    for i in range(numOfCols):
        headerString += 'mean' + str(i) + '\t'
    savetxt(fileName, data, fmt='%10.9f',delimiter="\t",header=headerString)

def meanShiftCorrect(inputArray, eventCoordinates, windowSize):
    '''
    @param inputArray: array of data with events to be mean-shifted
    @param eventCoordinates: array containing the starting points of the events to be corrected
    @param windowSize size of the baseline whose mean will be used to calculate the size of the shift
    @return array with the baseline-mean minus event min-shifted events 
    '''
    baseLineCoords = detectBaseline(inputArray, windowSize)[1]
    baseMean = mean(inputArray[baseLineCoords:baseLineCoords+windowSize])
    for x in eventCoordinates:
        event = inputArray[x[0]:x[1]]
        correctionValue =  abs(min(event)) - abs(baseMean)
        inputArray[x[0]:x[1]] -= correctionValue#+= correctionValue
    return(inputArray)

def pushToBaseline (dataMatrix, bucketSize):
    '''
    @param dataMatrix: matrix of the data whose baseline we're looking for
    @param bucketSize: size of the baseline bins
    @return: baseline corrected version of the data, currently also expressed data in terms of sigma
    '''
    numOfVals = dataMatrix.shape[1]
    coordinates = zeros(2)
    for i in range(0,numOfVals): # iterate over number of ROIs in file, e.g. collumn-wise
        baselineObject = detectBaseline(dataMatrix[:,i], bucketSize)
        coordinates[0] = baselineObject[1][0]
        coordinates[1] = baselineObject[1][1]
        baseLineArray = baselineObject[0]
        meanVal = abs(mean(baseLineArray))
        try:
            dataMatrix[:,i] = (dataMatrix[:,i]-meanVal)/meanVal # baseline correction
            dataMatrix[:,i] = dataMatrix[:,i]/(std(baseLineArray))
        except ZeroDivisionError:
            print(str(meanVal) + " at position: " + str(i))
            exit()

    return (dataMatrix,coordinates);

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
    for j in range(0,int(numOfEntries-bucketSize)): # iterate over all possible bucket positions
        thisStd = (axisStd(data[j:j+bucketSize])) # current deviation
        if (thisStd < lowestSigma): # new min deviation found
            lowestSigma = thisStd
            coordinate = (j,j+bucketSize)
        baselineArray = data[coordinate[0]:coordinate[1]]
    return(baselineArray, coordinate)

def axisStd(inArray):
    slopeArray = repeat(inArray[0],len(inArray))
    thisStd = sqrt(mean(sum(subtract(inArray,slopeArray)**2)))
    return (thisStd)

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
    maxVal = -maxint
    maxPosition = 0
    for position, value in enumerate(inputData):
        if (value > maxVal):
            maxVal = float(value)
            maxPosition = position 
    return(maxVal, maxPosition);

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
    risetimeArray = []
    decaytimeArray = []
    amplitudeArray = []
    
    for transients in transientMatrix: # go over cells
        for i in transients: # go over ROIs
            if (i.getNumOfPeaks() == 1): # select 1-peak entries
                risetimeArray.append(i.getRiseTime())
                decaytimeArray.append(i.getDecayTime())
                amplitudeArray.append(i.getAmplitude())
                            
    risetimeStd = std(risetimeArray)
    risetimeMean = mean(risetimeArray)
    decaytimeStd = std(decaytimeArray)
    decaytimeMean = mean(decaytimeArray)
    amplitudeStd = std(amplitudeArray)
    amplitudeMean = mean(amplitudeArray)
    risetimeThreshUp = risetimeStd*2 + risetimeMean
    decaytimeThreshUp = decaytimeStd*2 + decaytimeMean 
    amplitudeThreshUp = amplitudeStd*2 + amplitudeMean
    risetimeThreshDown = risetimeStd*2 - risetimeMean
    decaytimeThreshDown = decaytimeStd*2 - decaytimeMean 
    amplitudeThreshDown = amplitudeStd*2 - amplitudeMean
    print ("amp: %f\t%f\nRT: %f\t%f \nDT: %f\t%f" % (amplitudeThreshDown,amplitudeThreshUp,risetimeThreshDown,risetimeThreshUp,decaytimeThreshDown,decaytimeThreshUp))
    # block below: mean transients falling within 2 STD of 1 peak transient mean
    import matplotlib.pyplot as plt
    meanArray = zeros(0,dtype=float)
    divideArray = zeros(transientMatrix[0][0].getNumOfFrames(),dtype=int)
    constrainRelaxVariable = 1
    while (not meanArray.any()):
        for transients in transientMatrix: # go over cells
            for i in transients: # go over ROIs
                if(i.getNumOfPeaks() == 1):
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
    fitTime = time.time()
    popt, pcov = curve_fit(alphaKernel, tVariable,meanArray, p0=[80,100,50],maxfev=10000)
    plt.plot(tVariable, meanArray, label="meanArray")
    plt.plot(divideArray)
    plt.plot(alphaKernel(tVariable, *popt), 'r-', label='fit')
    plt.savefig('alphaKernel.png')
    plt.close()
    return(alphaKernel(tVariable,*popt))

def alphaKernel (t, A,t_A,t_B):
    return A*(exp(-t/t_A)-exp(-t/t_B))

def deconvolve(transients, kernel):
    '''TODO: this.
    IDEA: using kernel (e.g. exp(-t/tau)), model corrected fluorescence observations as exp. decay with
    point-wise increase for each spike
    '''
    spikeTrain = zeros(transients[0].getNumOfFrames())
    for i in (transients):
        transientData = i.getData()
        thisKernel = resize(kernel, transientData.shape)
        #spikeTrain = fft.ifft(fft.fft(transientData)/ fft.fft(transientKernel))# with ifft = inverse fast fourier transform and fft = fast fourier transform
        spikeTrain[i.getStartTime():i.getEndTime()] = abs(fft.ifft(divide(fft.fft(transientData), fft.fft(thisKernel))))
    return(spikeTrain)

def generateOutput(rawMatrix, baseObject,tempData, filename, numOfVals):
    #ugly AF placeholder, refactor this so its flexible with input
    transientData = []
    baselineMatrix = baseObject[0]
    baselineCoordinates = baseObject[1]
    eventCoordinates = tempData[1]
    quantileMatrix = matrix(tempData[0])
    correctionMatrix = tempData[3]
    
    import matplotlib.pyplot as plt
    #for j in range(numOfVals):
    #    plt.hist(delList[j], bins=100)
    #    plt.show()
    #    plt.close()
    #TODO: extract this crap into a function
    for i in range(numOfVals):
        f,(axarr0,axarr1,axarr2) = plt.subplots(3, sharex=True)
        axarr1.plot(quantileMatrix[:,i])
        axarr1.plot(correctionMatrix[:,i])
        axarr0.plot(rawMatrix[:,i],"r")
        thisMean = mean(quantileMatrix[:,i])
        minVal = float(min(quantileMatrix[:,i]))
        minVal2 = float(min(baselineMatrix[:,i]))
        maxVal2 = float(max(baselineMatrix[:,i]))
        if(eventCoordinates[i]):
            for j in eventCoordinates[i]:
                axarr1.plot(j,[minVal,minVal],'r-', lw=1)
                if ((j[1]) != 0):
                    #transientList.append(cal_neuroIm.Transient(thisSlice,j[0],j[1])) # call to transient kills 
                    axarr1.plot([j[0],j[0]],[minVal2,maxVal2],'r-', lw=1)
                    axarr2.plot([j[1],j[1]],[minVal2,maxVal2],'k-', lw=1)
                    #axarr2.annotate(str(transientList[-1].getStartTime()), xy=(j[0],transientList[-1].getNumOfPeaks()))
            #transientList[-1].setIsLast()
        #cal_neuroIm.getSinglePeakTransient(transientList)
        axarr2.plot(baselineMatrix[:,i])
        axarr1.plot([baselineCoordinates[i][0],baselineCoordinates[i][1]],[thisMean,thisMean],'k-',lw=1)
        axarr1.plot([baselineCoordinates[i][0],baselineCoordinates[i][0]],[minVal2,maxVal2],'k-', lw=1)
        axarr1.plot([baselineCoordinates[i][1],baselineCoordinates[i][1]],[minVal2,maxVal2],'k-', lw=1)
        plt.savefig(str(filename) + str(i)+".png")
        plt.close()
    try:    
        with open(str(filename) + str(filename).split('/')[-1] + ".csv", 'w') as outFile:
            outFile.write("Amplitude\tRiseTime\tDecayTime\tmeanIntensity\ttotalLength\tnumOfPeaks\tstartTime\tendTime\n")    
            for num,t in enumerate(transientData):
                outFile.write(str(t.getAmplitude()) + "\t" + str(t.getRiseTime())+ "\t" + str(t.getDecayTime()) + "\t" + str(t.getMeanIntensity()) + "\t" + str(t.getTotalTime()) +"\t" + str(t.getNumOfPeaks())+"\t" + str(t.getStartTime()) + "\t"+ str(t.getEndTime()) + "\n")
                if(t.getIsLast()):
                    outFile.write("__________________________ROI" + str(num) + "__________________________\n")
    except IOError:
        print(str(filename) +"output" +  str(filename).split('\\')[-1] + " is not a valid file path for output!")
        exit()
    #cal_neuroIm.writeOut(baselineMatrix, str(filename) + "base")
    #cal_neuroIm.writeOut(quantileMatrix, str(filename) + "quant")
    return(transientData)

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
    