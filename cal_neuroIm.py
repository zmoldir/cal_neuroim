'''
Created on Oct 3, 2016
script contains methods used for main.py
@author: maximilian
'''
import numpy
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import norm
from scipy.optimize import curve_fit
import pandas, xlrd, sys
from xlrd.biffh import XLRDError
from scipy.stats.stats import mode

def _importMatrix(filename,valSeparator):
    '''
    imports csv or xls formatted file, returns ndArray 
    @input: filename: input path; valSeperator: seperator between entries, treats file as xlrd if none
    @return: numOfVals: number of collumns in file (ROIs); dataMatrix: ndArray containing file data (dtype float)
    '''
    if(valSeparator):
        try:
            frame = pandas.read_csv(filename,sep = valSeparator,engine='python',skipinitialspace=True)    
        except (IOError,pandas.errors.ParserError) as err:
            print(err)
            print("See error message above! Either wrong file path or wrong separator.")
            exit()
        numOfVals = frame.shape[1]
        dataMatrix = frame.as_matrix()
        dataMatrix = numpy.nan_to_num(dataMatrix)
        numOfRois = len(dataMatrix)
        if(numOfVals < len(dataMatrix)):
            numOfRois = numOfVals
            numOfVals = len(dataMatrix)
        else:
            dataMatrix = dataMatrix.transpose()
    else:
        try:
            dataMatrix = []
            book = xlrd.open_workbook(filename)
            sh = book.sheet_by_index(0)
            numOfRois = sh.ncols
            for rx in range(sh.nrows):
                dataMatrix.append(map(float,sh.row(rx)))
            dataMatrix = numpy.matrix(dataMatrix).A
        except XLRDError:
            print("Wrong format! The file might contain comma-seperated values. Try using the -sep argument (E.G. -sep '\\t'")
            exit()
    if(not dataMatrix.dtype == "float64"):
        if(dataMatrix.dtype == "int64"):
            return(dataMatrix,numOfRois)
        print ("Can't convert file input to float! Wrong separator argument?")
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
    Given a time series, return the slope (which only depends on the start and end point)
    '''
    avgSlope = (inArray[-1] - inArray[0])/(len(inArray)-1)
    return (float(avgSlope))

def eventDetect (dataMatrix, quantileWidth=None,slopeWidth=None,cutoff=None,startNum=None,minEventLength=None,sigma=None):
    '''
    @dataMatrix: NxM matrix of the data containing transients and drift(possibly)
    @windowSize: size of the window considered for the quantile normalization (none if q = 0)
    @slopeWidth: size of the window considered for the calculation of slopes, critical parameter for event detection
    @cutoff: distance from slope distribution mean in standard deviance, where the value is taken as the transient threshold
    @startNum: number of threshold hits necessary for an event start to be declared
    @minEventLength: minimum number of frames above baseline for an event to be accepted
    @sigma: sigma of the gaussian filter for number of peak determination in the transients
    @return: transients: list (or list of lists) of objects of the type transient, one list per ROI (i.e. row in dataMatrix)
            slopeDistribution: distribution of slopes within the data, not used for further calculations but nice for visualizations
    ''' 
    #if- block below: default parameter values if none are supplied - overwritten if you hand the function anything
    if sigma == None: sigma = 7
    if quantileWidth == None: quantileWidth = 0
    if slopeWidth == None: slopeWidth = 3
    if cutoff == None: cutoff = 2
    if startNum == None: startNum = 3
    if minEventLength == None:  minEventLength = 5 # number of data points necessary above start level for an event to be kept as such
    
    transients = [] # list of lists of transients
    numOfEntries = len(dataMatrix) # number of entries in file
    startNumOfHits = startNum # number of hits necessary for an event start to be declared 
    
    thresholdList = []
    slopeDistributions = []
    # here we generate a global noise / slope threshold for all ROIs of the file
    for collumn in dataMatrix.T: 
        collumn = numpy.trim_zeros(collumn,'b')
        collumn = pandas.DataFrame(collumn)
        slopeList = collumn.rolling(slopeWidth,min_periods=1).apply(lambda x : _meanSlope(x))
        slopeList = slopeList.rolling(slopeWidth,min_periods=1).mean().as_matrix()
        slopeList = numpy.nan_to_num(slopeList)
        
        # best fit of data according to the log-likelihood maximization estimate
        thisMode = mode(numpy.around(slopeList,3),axis=None)[0]
        slopeList2 = [i for i in slopeList if i > thisMode]
        
        mu, deviance = norm.fit(slopeList2)
        # the histogram of the data
        thresholdList.append(deviance*cutoff) 
        
        import matplotlib.pyplot as plt
        import matplotlib.mlab as mlab

        #TODO: crash here is "max must be larger than min in range parameter"...
        n, bins, patches = plt.hist(slopeList, 100, normed=1, facecolor='green', alpha=0.75)
        y = mlab.normpdf( bins, mu, sigma)
        slopeDistributions.append((bins,y,slopeList,sigma*cutoff))
        plt.close()
    # re-iterate over the data ROI-wise for determination of event positions - in reverse order because we delete skipped rows immediately
    for horizontalPosition, collumn in enumerate(dataMatrix.T): # generate local noise / slope threshold and consider if we should skip
        base = _detectBaseline(collumn,400)[0]
        baseMean = numpy.mean(base)
        collumn = numpy.trim_zeros(collumn,'b').T # looks stupid but makes the general code more pythonic (instead of referencing collumns via slicing)
        numOfEntries = len(collumn)
        thisSlopeThreshold = thresholdList[horizontalPosition]# threshold to be passed for start of significance 
        isEvent = False # bool to remember wether we are in event-mode
        eventStart = 0
        eventEnd = 0
        threshHit = 0
        theseTransients = [] # local list for transients
        if(quantileWidth != 0):
            shiftValue = numpy.percentile(collumn,8)#shift of data to avoid negative values and the corresponding inversion
        else:
            shiftValue = 0
        
        for verticalPosition,val in enumerate(collumn): # iterate over number of pictures in file, e.g. row-wise
                
            thisSlope = slopeDistributions[horizontalPosition][2][verticalPosition]
            
            if (quantileWidth != 0):
                if (verticalPosition+quantileWidth<=numOfEntries):
                    correctionVal = numpy.percentile(collumn[verticalPosition:verticalPosition+quantileWidth],q=8)
                else:
                    correctionVal = numpy.percentile(collumn[numOfEntries-quantileWidth:],q=8)
                dataMatrix[verticalPosition,horizontalPosition] -= correctionVal
                
            if (not isEvent): 
                # we're not in an event ... so check if one is starting
                if (thisSlope > thisSlopeThreshold):
                    # it is! count hit
                    threshHit += 1
                    if(threshHit == startNumOfHits):
                        eventStart = verticalPosition - startNumOfHits +1
                        isEvent = True
            else: 
                # we are in an event, does it end?
                if (collumn[eventStart] > collumn[verticalPosition] or collumn[eventStart] < baseMean 
                    or collumn[eventStart < 1]): # check if it ends
                    isEvent = False
                    eventEnd = verticalPosition
                    threshHit = 0
                    if(eventEnd-eventStart > minEventLength):
                        
                        eventStart -= _getOnsetCoordinate(collumn[max(eventStart-100,0):eventStart],baseMean)
                        theseTransients.append(_Transient(collumn[eventStart:eventEnd],eventStart,eventEnd,numOfEntries,sigma))
                        
            if (verticalPosition >= len(collumn)-slopeWidth-1):
                break
        if(isEvent):
            # event lasted until end of the data, correct the last bit
            #theseEventEndCoordinates.append([eventStart,eventEnd])
            eventStart -= _getOnsetCoordinate(collumn[max(eventStart-100,0):eventStart],baseMean)
            theseTransients.append(_Transient(collumn[eventStart:],eventStart,numOfEntries,numOfEntries,sigma))

        #eventEndCoordinates.append(theseEventEndCoordinates[:])
        transients.append(theseTransients)
        dataMatrix[:,horizontalPosition] += shiftValue # shift of data to avoid negatives if we quantile corrected
        #if (theseEventEndCoordinates):axarr2
            #collumn = discardNonEvent(collumn, theseEventEndCoordinates,baseLineArray[horizontalPosition])
    return (transients,slopeDistributions);

def thresholdEventDetect(dataMatrix, quantileWidth=None, minEventLength= None, thresholdDeviance= None,sigma=None):
    '''
    alternative to eventDetect, uses a simple threshold (baseline mean + 5 standard deviations of baseline) for transient detection
    '''
    if sigma == None: sigma = 7
    if quantileWidth == None: quantileWidth = 0
    if minEventLength == None: minEventLength = 15
    if thresholdDeviance == None: thresholdDeviance = 2.5
    
    transients = []
    numOfVals =  len(dataMatrix)
    for collumn in (dataMatrix.T):
        theseTransients = []
        eventBool = True
        collumn = numpy.trim_zeros(collumn,'b').T
        eventStart = 0
        eventEnd = 0
        numOfEntries = len(collumn)
        
        if(quantileWidth > 0):    
            
            for position2,value in enumerate(collumn):
                if (quantileWidth > 0):
                    if(position2+quantileWidth < numOfVals):
                        correctionValue = numpy.percentile(collumn[position2:position2+quantileWidth],8)
                    else:
                        correctionValue = numpy.percentile(collumn[numOfVals-quantileWidth:],8)
                    value[...] -= correctionValue
        base = _detectBaseline(collumn,400)[0]
        baseMean = numpy.mean(base)
        deviance = numpy.std(collumn)
        threshold =  baseMean + deviance * thresholdDeviance
        for position2,value in enumerate(collumn):
            if eventBool: # we are not in an event, check if one is starting
                if threshold < value:# yes: below, no: continue
                    #print("event started at val: %f vs thres %f" % (value, threshold))
                    eventBool = False
                    onset = _getOnsetCoordinate(collumn[max(position2-100,0):position2],baseMean+ deviance)
                    eventStart = position2 - onset
            else:
                if baseMean + deviance > value and minEventLength <  position2 - eventStart:
                    eventBool = True
                    eventEnd = position2
                    theseTransients.append(_Transient(collumn[eventStart:eventEnd]-baseMean,eventStart,eventEnd,numOfEntries,sigma))
        if(not eventBool):
            theseTransients.append(_Transient(collumn[eventStart:]-1,eventStart,len(collumn),numOfEntries,sigma))

        transients.append(theseTransients)
    return(transients)

def _getOnsetCoordinate(roi,baseMean):
    '''
    short code snippet to find onset coordinate in front of the threshold-passing data point
    '''
    for onset,val in enumerate(reversed(roi)):
        if val <= baseMean:
            return onset
    return False

def _quantileNorm (inputArray, windowSize):
    '''DEPRECATED Alternative version of quantileCorrect which does not skip over events- not in use, eventDetect does this anyway
    '''
    correctionArray = numpy.zeros(inputArray.shape)
    for pos,i in enumerate(inputArray.T):
        if(isinstance(i,numpy.ndarray)):
            shiftValue = numpy.percentile(inputArray[:,pos],92)
            for j,val in enumerate(i):
                correctionVal = numpy.percentile(inputArray[j:j+windowSize,pos],5)
                inputArray[j,pos] -= correctionVal
                correctionArray[j,pos] = correctionVal
            i[...] += shiftValue
        else:
            shiftValue = numpy.percentile(inputArray,92)
            correctionVal = numpy.percentile(inputArray[pos:pos+windowSize],92)
            i[...] -= correctionVal
            correctionArray[pos] = correctionVal
    return(inputArray,correctionArray);


def _writeOut(data, fileName):
    '''
    Small method to print matrix with "meanN"- header, as found in .cls files
    '''
    numOfCols = (data.shape[1])
    headerString = ''
    for i in range(numOfCols):
        headerString += 'mean' + str(i) + '\t'
    numpy.savetxt(fileName, data, fmt='%10.9f',delimiter="\t",header=headerString)


def pushToBaseline (dataMatrix, bucketSize=None):
    '''
    @param dataMatrix: matrix of the data whose baseline we're looking for
    @param bucketSize: size of the baseline bins
    @return: baseline corrected version of the data: 
            coordList: tuple with start / end coordinate of the region determined as baseline activity, for visualization purposes
    '''
    if bucketSize == None: bucketSize = 300
    coordList = []
    for roi in (dataMatrix.T): # iterate over number of ROIs in file, e.g. collumn-wise
        baseLineArray, coordinates= _detectBaseline(roi, bucketSize)
        meanVal = abs(numpy.mean(baseLineArray))
        coordList.append(coordinates)
        roi[...] = roi/meanVal - numpy.mean(baseLineArray/meanVal)# baseline correction 
        
    return (dataMatrix,coordList);

def _detectBaseline (data, bucketSize):
    '''
    @param data: the array out of which the region with the lowest noise is to be identified
    @param bucketSize: size of the bins to be checked
    @return: bin with the lowest noise and its starting coordinate, in a tuple
    '''
    data = numpy.trim_zeros(data, 'b')
    numOfEntries = len(data)
    lowestSigma = sys.maxsize # for size comparasion
    #baselineArray = numpy.zeros(bucketSize)
    coordinate = []
    for j in range(0,int(numOfEntries-bucketSize),int(numOfEntries/(bucketSize*2))): 
        thisStd = numpy.std(data[j:j+bucketSize])#(axisStd(data[j:j+bucketSize])) # current deviation
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
    peaklineArray = numpy.zeros(bucketSize) # stores portion of data which we consider the peak line
    # print(str(data) + " " + str(len(data)))
    for j in range(0,int(numOfEntries-bucketSize)): # iterate over all possible bucket positions
        thisStd = (numpy.var(data[j:j+bucketSize])) # current deviation
        if (thisStd > highestSigma): # new min deviation found  
            highestSigma = thisStd 
            peaklineArray = data[j:j+bucketSize]
            coordinate = j
    return(peaklineArray,coordinate)

def _maxAmp(inputData):
    max_index = numpy.argmax(inputData)
    return(inputData[max_index], max_index);

def _detectNumOfPeaks(data,sigma):
    if(len(data) <10):
        return(1)
    numOfPeaks = 0
    steps = int(len(data)/10)

    data = pandas.DataFrame(data)
    data.rolling(steps,min_periods=min(steps,3)).mean().as_matrix()
    data = gaussian_filter(data, sigma)    
    prevPoint = 0
    oldSign = 1
    for i in data:
        newSign = numpy.sign(i- prevPoint)
        if (newSign != oldSign):
            numOfPeaks += 1
        oldSign = newSign
        prevPoint = i
    
    return(numpy.ceil(numOfPeaks/2.));

def createMeanKernel(transientMatrix):
    '''
    @transientMatrix: list (or list of lists) of transients
    @return: alpha kernel of length n
    The transients fed into this function should come from comparable traces,
     i.e. same cell type / experimenting run.
    It selects single-peaked transients and then further retains only those 
    within the [0.2,0.8]-numpy.percentile of rise time, decay time,
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
    
    aucThreshUp = numpy.percentile(aucArray, 80)
    risetimeThreshUp = numpy.percentile(risetimeArray, 80)
    decaytimeThreshUp = numpy.percentile(decaytimeArray, 80)
    amplitudeThreshUp = numpy.percentile(amplitudeArray, 80)
    aucThresDown = numpy.percentile(aucArray, 20)
    risetimeThreshDown =  numpy.percentile(risetimeArray, 20)
    decaytimeThreshDown = numpy.percentile(decaytimeArray, 20)
    amplitudeThreshDown =  numpy.percentile(amplitudeArray, 20)
    
    import matplotlib.pyplot as plt
    meanArray = numpy.zeros(0,dtype=float)
    divideArray = numpy.zeros(transientMatrix[0][0].numOfFrames,dtype=int)
    for i in singlePeakArray: # go over ROIs
        if (amplitudeThreshDown < i.amplitude < amplitudeThreshUp and 
                decaytimeThreshDown < i.decayTime < decaytimeThreshUp and 
                risetimeThreshDown < i.riseTime < risetimeThreshUp and
                aucThresDown < i.auc < aucThreshUp): # we already know which ones are the single peak transients, but re-checking is cheaper than a new array
            currentTransient = i.data
            plt.plot(currentTransient,"grey",alpha=0.2)
            if(len(currentTransient)>len(meanArray)):
                meanArray.resize(currentTransient.size)
                tempArray = numpy.ones(currentTransient.size,dtype=int)
                tempArray.resize(divideArray.shape,refcheck=False)
                divideArray += tempArray  
                meanArray += currentTransient
            else:
                tempArray = numpy.ones(currentTransient.size,dtype=int)
                tempArray.resize(divideArray.shape)
                divideArray += tempArray
                zerosArray = numpy.zeros(len(meanArray)-len(currentTransient))
                currentTransient = numpy.concatenate((currentTransient, zerosArray))
                meanArray += currentTransient
    # loop up top: finished adding single peak transients, need to divide by our divideArray (after cropping it accordingly)
    divideArray.resize(len(meanArray))
    meanArray = meanArray / (divideArray)
    tVariable = numpy.arange(0,len(meanArray),1.)
    popt, pcov = curve_fit(_alphaKernel, tVariable,meanArray, p0=[80,100,50],maxfev=100000)
    # print popt
    #plt.plot(tVariable, meanArray, label="meanArray")
    #plt.plot(_alphaKernel(tVariable*2, *popt), 'r-', label='fit')
    #plt.savefig('/home/maximilian/unistuff/paris_ens/cal_neuroim/alphaKernel.png')
    #plt.close()
    #newKernel = numpy.tile(_alphaKernel(tVariable*2, *popt), i.numOfFrames/len(meanArray))
    return(_alphaKernel(numpy.arange(0,i.numOfFrames,1), *popt))

def _alphaKernel (t, A,t_A,t_B):
    return A*(numpy.exp(-t/t_A)-numpy.exp(-t/t_B))

def deconvolve(transients, kernel):
    '''
    @transients: list or list of lists of transient-objects, like returned by eventDetect and thresholdEventDetect
    @kernel: kernel to be used for deconvolution, i.e. as returned by createMeanKernel  
    @return: 1xn or dxn array (corresponding to the input) of integers, 
    with n= number of frames, d= number of lists.
    Integer values correspond to the number of action potentials in the respective frame.
    '''
    interpolSteps = 1
    if (not transients): # no transients found - the input is empty.
        return(numpy.zeros(kernel.size))
    if (any(isinstance(el, list) for el in transients)): #this is true, we have a list of lists -> iterate over lists and transients
        spikeSignal = numpy.zeros((len(transients),transients[0][0].numOfFrames))# problem: crashes if first ROI doesn't contain a transient
        for i,transientList in enumerate(transients):
            for transient in transientList:
                #problem we check for below: can we interpolate 10-point-wise, or is the transient too big?
                '''if (transient.data.size*10 < transient.numOfFrames): 
                    newTransient = numpy.interp(numpy.arange(0,len(transient.data),interpolSteps), range(len(transient.data)), transient.data)
                    thisKernel1 = numpy.resize(kernel,len(newTransient))#transient.data.shape)
                    deconvolved2 = numpy.real(numpy.fft.ifft(numpy.divide(numpy.fft.fft(newTransient), numpy.fft.fft(thisKernel1))))
                    spikeSignal[i][transient.startTime:transient.endTime] = _generateSpiketrainFromSignal(deconvolved2,transient.data.size)
                else:'''
                thisKernel = numpy.resize(kernel,transient.data.shape)
                deconvolved = numpy.real(numpy.fft.ifft(numpy.divide(numpy.fft.fft(transient.data), numpy.fft.fft(thisKernel))))
                spikeSignal[i][transient.startTime:transient.endTime] = _generateSpiketrainFromSignal(deconvolved)   
    else: # the above clause is false, we have a list of transients and can iterate directly
        spikeSignal = numpy.zeros(transients[0].numOfFrames)
        for transient in (transients):
            '''if (transient.data.size*10 < transient.numOfFrames): 
                newTransient = numpy.interp(numpy.arange(0,len(transient.data),interpolSteps), range(len(transient.data)), transient.data)
                thisKernel1 = numpy.resize(kernel,len(newTransient))#transient.data.shape)
                deconvolved2 = numpy.real(numpy.fft.ifft(numpy.divide(numpy.fft.fft(newTransient), numpy.fft.fft(thisKernel1))))
                spikeSignal[transient.startTime:transient.endTime] = _generateSpiketrainFromSignal(deconvolved2,transient.data.size)

            else:'''
            thisKernel = numpy.resize(kernel,transient.data.shape)
            deconvolved = numpy.real(numpy.fft.ifft(numpy.divide(numpy.fft.fft(transient.data), numpy.fft.fft(thisKernel))))
            spikeSignal[transient.startTime:transient.endTime] = _generateSpiketrainFromSignal(deconvolved)
    return(spikeSignal)

def _generateSpiketrainFromSignal(spikeSignal,originalLength = None):
    interpolationFactor = 1
    if originalLength == None: 
        originalLength = spikeSignal.size
        interpolationFactor = 1
    if (numpy.isnan(spikeSignal).any()):
        return 0
    
    spikeTrain = numpy.zeros(originalLength)
    numOfSpikes = int(numpy.ceil(sum(spikeSignal)))
    spikeTrain[0] = 1
    for i in range(1,numOfSpikes):
        maxVal = numpy.argmax(spikeSignal)
        spikeTrain[maxVal/interpolationFactor] = 1
        spikeSignal[maxVal] /= 2 
    return(spikeTrain)

def negativeTransientsCorrect(transients):
    '''
    @transients: list of transients to be checked/corrected
    @return: same list, with the smallest transients removed in 
    proportion to the number of negative going transients
    '''
    numOfNeg = 0
    for i in transients:
        meanVal = numpy.mean(i.data)
        if meanVal< 0:
            transients.remove(i)
            numOfNeg += 1
            
    for j in range(numOfNeg):
        del transients[numpy.argmin(meanVal)]
    return(transients)

def aucCorrect(transients, kernel):
    '''
    @transients: list of transients to be checked/corrected
    @kernel: kernel used for deconvolution
    @return: same transients as before, minus those which have 
    less than 10% of the kernel surface area
    -naive approach, less accurate for overlapping activity
    '''
    
    kernelAUC = numpy.trapz(kernel)
    if (any(isinstance(el, list) for el in transients)):
        for j in transients:
            for i in j:
                transientAUC = numpy.trapz(i.data)
                if (transientAUC/kernelAUC < 0.2):
                    j.remove(i)
    else:
        for i in (transients):
            transientAUC = numpy.trapz(i.data)
            if (transientAUC/kernelAUC < 0.2):
                transients.remove(i)

    return(transients)

class _Transient:
    '''
    class to store transient properties, such as rise time, decay time, amplitude and duration
    '''
    def __init__(self, data, startTime,endTime, numOfFrames,sigma):
        
        self.maxVal = _maxAmp(data)
        self.startTime = startTime
        self.endTime = endTime
        self.amplitude = float(self.maxVal[0]) 
        self.riseTime = int(self.maxVal[1])
        self.decayTime = int(len(data) + 1 - self.riseTime)
        self.totalTime = len(data)
        self.numOfPeaks = _detectNumOfPeaks(data,sigma)
        self.data = data
        self.numOfFrames = numOfFrames
        self.coordinates = (self.startTime, self.endTime)
        self.auc = numpy.trapz(data) 