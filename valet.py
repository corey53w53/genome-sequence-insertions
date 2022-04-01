###########################
## valet.py
##
## Module that contains helper functions useful for valet
###########################
import math

def iqr(valsin, weightsin=None, thresh=1.5, scale=True):
    """
    NAME: iqr()

    PURPOSE:
        Computers inter-quartile range for a list of numbers, potentially
        with an associated set of weights.

    :param valsin: values for which the IQR is computed
    :type valsin: list
    :param weightsin: weights for the values
    :type weightsin: list
    :param thresh: IQR threshold
    :type thresh: float
    :param scale: Scale by coverage
    :type scale: bool
    :return: tuple containing the 25% and 75% quantiles of the set of numbers
    :rtype: tuple
    """

    ## Note, this implementation is based on sorting - there are more efficient ways of doing this
    if weightsin and len(valsin) != len(weightsin) :
        raise ValueError("iqr: Input list and weights must have the same number of elements")

    avgcvg = 0
    totsize = 0

    tuples = []
    if weightsin :
        for i in range(len(valsin)):
            tuples.append([valsin[i], weightsin[i]])
            avgcvg += valsin[i] * weightsin[i]
            totsize += weightsin[i]
    else:
        for i in range(len(valsin)):
            tuples.append([valsin[i]])
            avgcvg += valsin[i]
            totsize += 1

    avgcvg /= totsize

    tuples.sort(key=lambda x: x[0]) # sort by value

#    print(tuples)

    lowval = None
    highval = None
    if weightsin :   # weighted option
        totval = 0
        for v in weightsin:
            totval += v

        lower = totval * 0.25
        higher = totval * 0.75

#        print("lower = {l}, higher={h}".format(l=lower, h=higher))

        totv = 0
        for i in range(len(tuples)):
            totv += tuples[i][1]
            if not lowval and totv >= lower:
                lowval = tuples[i-1][0]
            if not highval and totv >= higher:
                highval = tuples[i][0] # value that pushed us over the edge
    else: # just plain
        lower = int(len(valsin) * 0.25)
        higher = int(len(valsin) * 0.75)
        lowval = tuples[lower][0]
        highval = tuples[higher][0]

    if not lowval or not highval:
        print("Could not find quartiles?")
        exit(1)

    if lowval > highval:
        print("Low {l} is higher than high {h}".format(l=lowval, h=highval))
        exit(1)

    adj = thresh * (highval - lowval)
    if scale:
        adj /= avgcvg

#    print("lo = {l}, hi = {h}, adj = {a}".format(l=lowval, h=highval, a=adj))
    lowval -= adj
    if lowval < 0:
        lowval = 0

    highval += adj
    if highval > max(valsin):
        highval = max(valsin)

    return int(lowval), int(highval)  # nominally the 1.5 IQR calculation

#################################
def poisswin(listin, totlen, winsize=300, mtesting=True, pthresh=0.05):
    """
    NAME:
        poisswin()

    PURPOSE:
        Take in a list of the positions along a sequence of
        length totlen and returns the coordinates of windows
        that contain more events than expected according to poisson
        statistics

    :param listin: The list of coordinates. This is assumed
                to be in sorted order
    :type listin: list
    :param totlen: The total length of the sequence
    :type totlen: int
    :param winsize: The size of the window (default: 300)
    :type winsize: int
    :param mtesting: Account for multiple testing (multiply p-value by # of windows) (default: True)
    :type mtesting: bool
    :param pthresh: P-value threshold (default: 0.05)
    :type pthresh: float

    :return:
        A list of tuples that contain information about each window:
        (start, end, best_start, best_end, pvalue)
        start and end are locations in the list rather than the actual coordinates
    :rtype: list
    """

# Code starts here
    rt = len(listin) / totlen    # rate parameter for poisson statistic

    outlist = []  # this is what we'll return
    best = {}  # current best window
    last = {} # last significant window
    i = 0
    while i < len(listin) : # start a window at each "event"
        j = i  # this will be the end of the window
        while j < len(listin) and (listin[j] - listin[i]) <= winsize:
            if j > 0 and listin[j] < listin[j - 1]: # list is not in order
                raise ValueError("Input list is not sorted")
            j += 1

        # At this point j is one past the end of the window
        # number of events in window is j - i
        n = j - i
        w = winsize
        if j == len(listin) :  #dealing with a partial window, must adjust window size
            w = totlen - listin[i]

        ## TODO fix the math in case we have large n (i.e., where math.factorial will fail)
        pval = 0
        try:
            pval = (rt * w)**n * math.exp(-rt * w) / math.factorial(n)
        except OverflowError:
            if n > rt * w * 4: # at least four times the rate...yes, it's a hack
                pval = 0
            else:
                print("N too large {n} {r}".format(n=n, r=rt*w))
                exit(1)

        if mtesting :
            pval *= len(listin)  # adjust by number of tests we have made

        if pval <= pthresh : # significant interval
            # check if current window overlaps with best window
            if not best : # best list is empty
                best = {'s':i, 'e':j, 'p':pval}
            elif best['p'] > pval :  # update best if necessary
                best = {'s':i, 'e':j, 'p':pval}
                if last : # need to update the best window stored in last
                    last['bs'] = best['s']
                    last['be'] = best['e']
                    last['bp'] = best['p']

            if not last : # if there is no last significant window
                last = {'s': i, 'e': j, 'bs': best['s'], 'be': best['e'], 'bp': best['p']}
            elif i == last['e'] : # if we are at the next window
                last['e'] = j # update last
            else: # need to start new best and last
#                best = [i, j, pval]
                outlist.append(last)
                last = {'s': i, 'e': j, 'bs': best['s'], 'be': best['e'], 'bp': best['p']}
            i = j # skip over the whole window and start again
        else: # not a significant window
            if last:
                outlist.append(last)
            last = {}
            best = {}
            i += 1

    # now we must output the final best
    if last :
        outlist.append(last)

    return outlist

##############################

## coverage code here...
## TODO add window "weight" score - sum of depth times width
def flagCoverage(inlist, totlen, findmax=True, findall=False) :
    """
    NAME:
        flagCoverage()
        
    PURPOSE:
        Computes depth of coverage along a chromosome, and outputs
        regions that are either deeper or shallower than expected.
        
        It either reports the interval with a maximal/minimal coverage
        within a window that is unusually deep/shallow, or the entire
        window.
        
    
        
    :param inlist: list of lists, each representing a different interval
                   the individual lists must have the first and second
                   elements represent the left/right (in order) ends of the
                   interval. The other elements are ignored.
    :type inlist: list
    :param totlen: total length of the sequence within which intervals are located
    :type totlen: int
    :param findmax: determines if it will flag regions with unusually high (findmax=True)
                    or unusually low (findmax=False) depth of coverage. (default: True)
    :type findmax: bool
    :param findall:    determines if to return all the regions with unusual coverage (all=True)
                   or just the window with the highest coverage per region (default)
    :return:       returns a list of tuples, each containing the coordinates of a flagged
                   regions as well as the depth of coverage
    :rtype: list
    """


    # some variables that have to be adjusted based on whether we're looking for
    # min or max
    adj = 1  # multiplier for adjusting looking for max (1) or min (-1)
#    zthresh = 3 # z value threshold
    if not findmax: # looking for minimum
        adj = -1
#        zthresh = 2

    # Starts and ends of intervals
    starts = []
    ends = []
    for i in range(len(inlist)) :
        starts.append(inlist[i][0])
        ends.append(inlist[i][1])

    starts.sort()
    ends.sort()

    s = 0
    e = 0
    cvg = 0
    coverages = []
    last = 0
    while s < len(starts) and e < len(ends) : # the usual merge sort happening here
        if starts[s] <= ends[e] : # found a start
            coverages.append([last, starts[s], cvg])
            last = starts[s]
            cvg += 1
            s += 1
        else : # found an end
            coverages.append([last, ends[e], cvg])
            last = ends[e]
            cvg -= 1
            e += 1

    # At this point we should only have ends left
    if s < len(starts) :
        print("Weird - starts ended before ends")
        exit(1)

    while e < len(ends):
        coverages.append([last, ends[e], cvg])
        last = ends[e]
        cvg -= 1
        e += 1

    coverages.append([last, totlen, cvg])

#    print(coverages)
    vals = []
    wghts = []
    for i in range(len(coverages)) :
        vals.append(coverages[i][2]) # coverage
        wghts.append(coverages[i][1]-coverages[i][0]) # length of interval

#    (mean, stdev) = wmom(vals, wghts, calcerr = True)
    (low, high) = iqr(vals, wghts)
#   print(wmom(vals, wghts, calcerr=True))
#    print("Mean = {m}, SD = {s}".format(m=mean, s=stdev))
#    print("Low = {l}, High = {h}".format(l=low, h=high))

    # now print all regions with > 3 stdev from mean
    outlist = []
    curr_window = {} # contains fields: s -start, e - end, bs - best start, be- best end, c -best coverage
                     # w - window weight (coverage times width)
    windowWeight = 0

    for i in range(len(coverages)) :
#        if adj*(coverages[i][2] - mean) > zthresh*stdev :
        if adj > 0 and coverages[i][2] > high or adj < 0 and coverages[i][2] < low:  # outlier
            if curr_window: # a window is currently open
                if findall : # simply extend it
                    if adj*(curr_window['c'] - coverages[i][2]) < 0:
                        curr_window['c'] = coverages[i][2]
                        curr_window['bs'] = coverages[i][0]
                        curr_window['be'] = coverages[i][1]

                    windowWeight += coverages[i][2]*(coverages[i][1] - coverages[i][0])
                    curr_window['e'] = coverages[i][1]
                    curr_window['w'] = windowWeight
                elif adj*(curr_window['c'] - coverages[i][2]) < 0: # replace window
                    windowWeight = coverages[i][2]*(coverages[i][1] - coverages[i][0])
                    curr_window = {'s': coverages[i][0], 'e': coverages[i][1], 'c': coverages[i][2],
                                   'bs': coverages[i][0], 'be': coverages[i][1], 'w': windowWeight}
            else:
                windowWeight = coverages[i][2] * (coverages[i][1] - coverages[i][0])
                curr_window = {'s': coverages[i][0], 'e': coverages[i][1], 'c': coverages[i][2],
                               'bs': coverages[i][0], 'be': coverages[i][1], 'w': windowWeight}
        else:  # not an outlier
            if curr_window: # output the window and reset the current one
                outlist.append(curr_window)
                curr_window = {}

    if curr_window:  # if last window is outlier
        outlist.append(curr_window)

    return outlist

############################