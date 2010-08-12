"""
GeneTrack Lite
a.k.a., genetrack_lite
2/27/2010 (c) MIT License (http://www.opensource.org/licenses/mit-license.php)
Shinichiro Wachi

=============================================================================
The MIT License

Copyright (c) 2010 Shinichiro Wachi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
=============================================================================

Why?
Created mainly to avoid going through the gene-track web interface
by allowing command-line scripting of peak calling.

How?
Ripped out from old genetrack code (MIT Open Source License), this module produces the 
"gaussian kernel" smoothed values and picks peaks, just like the web-interface counterpart.

For citation, please refer to:
Albert I, Wachi S, Jiang C, Pugh BF, Genetrack--a genomic data processing and visualization framework.
Bioinformatics 2008 May 15;24(10):1305-6
PMID: 18388141
"""
import numpy
import operator
from math import log, exp
from itertools import izip, count, islice, imap, ifilter, tee, repeat

# ####### ######## ######## ######## ######## ######## ######## ######## ######## ######## 
# ####### ######## ######## ######## ######## ######## ######## ######## ######## ######## 
# ####### ######## ######## ######## ######## ######## ######## ######## ######## ######## 
# Modified functions from GeneTrack by I.Albert follow:

### modified from genetrack: smooth.py
def normal_function( sigma, width ):
    """
    Defaulf fitting function, it returns values from a normal distribution
    """
    log2    = log(2)
    sigma2  = float(sigma)**2
    lo, hi  = width, width+1
    
    def normal_func(value, index):
        return value * exp( -index*index/sigma2 * log2 )    
    
    values = [ normal_func(1, x) for x in range(-lo, hi) ]
    values = numpy.array( values )

    return lo, hi, values

def safe_add(fit, ix, lo, hi, values):
    """
    Mutates fit data, adds values to the slice
    """
    actual = len ( fit[ix-lo:ix+hi] )
    fit[ix-lo:ix+hi] += values[:actual]
    
### Modified from genetrack: peack_detection.py
def make_list( size, default=0, pad=0, ntype=numpy.int8):
    if default==0:
        return numpy.zeros( size+pad, ntype )
    else:
        return numpy.ones( size+pad, ntype )
    
def detect_peaks( data ):
    """
    Detects peaks (local maxima) from an iterator that generates 
    float values. Input data is a list of tuples containing 
    the index and the value at the index.
    
    Returns a generator yielding tuples that corresponds 
    to the peak index and peak value.
    
    >>> values  = [ 0.0, 1.0, 2.5, 1.0, 3.5, 1.0, 0.0, 0.0, 10.5, 2.0, 1.0, 0.0 ]
    >>> data  = izip( count(), values )
    >>> peaks = detect_peaks( data )
    >>> peaks = list(peaks)
    >>> peaks
    [(2, 2.5), (4, 3.5), (8, 10.5)]
    >>> select_peaks( peaks, width=1)
    [(2, 2.5), (4, 3.5), (8, 10.5)]
    >>> select_peaks( peaks, width=2)
    [(4, 3.5), (8, 10.5)]
    """
    # a more functional approach using less memory
    def is_local_maxima( row ):
        "Selects local maxima"
        return row[0][1] < row[1][1] >= row[2][1]

    teedata   = tee(data, 3)
    mapped    = imap( islice, teedata , count(), repeat(None) )
    zipped    = izip( *mapped )
    filtered  = ifilter( is_local_maxima, zipped )
    for left, mid, right in filtered:
        yield mid

def select_peaks( peaks, width, level=0):
    """
    Selects successive non-overlapping peaks 
    with a given exclusion zone.
    Takes as input a list of (index, value) tuples.
    Returns a list of tuples = ( midpoint, maxima )
    """
    # order by peak height
    work  = [ (y, x) for x, y in peaks if y >= level and x > width ]
    work.sort()
    work.reverse()

    # sanity check
    if not work:
        return []

    #largest index        
    xmax = max(work, key=operator.itemgetter(1) )[1]
    
    # selected peaks
    selected = []

    # keeps track of empty regions
    empty = make_list( size=xmax, pad=width+1, default=1 )
    
    half = width/2
    for peaky, peakx in work:
        if empty[peakx]:
            left  = peakx - width
            right = peakx + width
           
            # store into output data
            selected.append( ( peakx, peaky ) )
            
            # block the region
            empty[left:right] = numpy.zeros (right - left)

    selected.sort()
    return selected 


def smooth_data(indices_, tag_counts_, sigma=20, width=120, size=None):
    """
    #sigma = 20
    #width = 120
    #epsilon = 0 #some minimum peak height
    #size, epsilon = data_size, minimum_peak_size
    """
    if size is None:
        size = max(indices_)
    
    lo, hi, normal = normal_function(sigma, width)
    
    pad = hi + 1
    fit = numpy.zeros( size + pad )  # watson fit
    
    for idx, _tag_count in zip(indices_, tag_counts_):
        if _tag_count > 0.0:
            safe_add( fit=fit, ix=idx, lo=lo, hi=hi, values=normal*_tag_count )
    
    return fit

# ####### ######## ######## ######## ######## ######## ######## ######## ######## ######## 
# ####### ######## ######## ######## ######## ######## ######## ######## ######## ######## 
# ####### ######## ######## ######## ######## ######## ######## ######## ######## ######## 
# Original functions for GeneTrack Lite follow:

def example_usage_for_peak_detection():
    values  = [ 0.0, 1.0, 2.5, 1.0, 3.5, 1.0, 0.0, 0.0, 10.5, 2.0, 1.0, 0.0 ]
    data  = izip( count(), values )
    peaks = detect_peaks( data )
    peaks = list(peaks)
    
    print select_peaks( peaks, width=2)
    
def get_tag_count_data(infilename, chromosome='chr01'):
    """
    Reads in tag count of a particular chromosome.
    Input format is the same as that of GeneTrack/LionDB3
    """
    k_idx = []
    k_fwd = []
    k_rev = []
    for i in  open(infilename):
        isplit = i.split()
        if isplit[0]==chromosome:
            k_idx.append( int(isplit[1]))
            k_fwd.append( int(isplit[2]))
            k_rev.append( int(isplit[3]))
    return k_idx, k_fwd, k_rev

def get_peaks_iter(smoothed_tag_counts):
    return detect_peaks( izip(count(), smoothed_tag_counts))

def get_peaks_list(smoothed_tag_counts):
    """
    e.g.: 
    peaks = list(detect_peaks( izip(count(), wxfit)))
    selected_peaks = select_peaks(peaks, width=400, level =0)
    """
    return list(get_peaks_iter(smoothed_tag_counts))
def get_selected_peaks(peak_iter, width=400, level =0):
    """
    selected_peaks = select_peaks(peak_iter, width=400, level =0)
    """
    return select_peaks(peak_iter, width=400, level =0)

def get_selected_peaks_from_fitted_data(fit_data, width=400, level =0):
    peak_iter = get_peaks_iter(fit_data)
    return select_peaks(peak_iter, width=width, level=level)

def get_peaks(idx, tag_counts, sigma, width, level=0):
    xfit = smooth_data(idx, tag_counts, sigma, width)
    peak_iter = get_peaks_iter(xfit)
    selected_peaks = select_peaks(peak_iter, width=width, level=level) # format: (index, height)
    return selected_peaks

def sum_indexed_tag_counts(dict_list):
    """
    Usage:
    infilename = index_file_root%files[0]
    k_idx, k_fwd, k_rev = genetrack_like.get_tag_count_data(infilename, chromosome)
    dict1 = dict(zip(k_idx, zip(k_fwd, k_rev)))
    
    infilename = index_file_root%files[1]
    k_idx, k_fwd, k_rev = genetrack_like.get_tag_count_data(infilename, chromosome)
    dict2 = dict(zip(k_idx, zip(k_fwd, k_rev)))
    
    infilename = index_file_root%files[2]
    k_idx, k_fwd, k_rev = genetrack_like.get_tag_count_data(infilename, chromosome)
    dict3 = dict(zip(k_idx, zip(k_fwd, k_rev)))
    
    k_idx, k_fwd, k_rev = sum_indexed_tag_counts([dict1,dict2,dict3])
    """ 
    def sum_tuples(pair1, pair2):
        """
        >>> pair1 = (1,2)
        >>> pair2 = (3,4)
        >>> z12 = zip(pair1, pair2)
        >>> z12
        [(1, 3), (2, 4)]
        >>> [x + y for x, y in z12]
        [4, 6]
        """
        k = zip(pair1, pair2)
        return [x + y for x, y in k]
    
    def get_a_pair(dictX, idx):
        """
        inits a pair (0,0) if index doesn't exist, return value if it does.
        """
        if idx in dictX.keys():
            return dictX[idx]
        else:
            return (0,0)
    
    def sum_a_dict(dict1, dict2):
        """
        """
        summed_dict = {}
        for idx in list(set(dict1.keys()).union(set(dict2.keys()))):
            pair1 = get_a_pair(dict1, idx)
            pair2 = get_a_pair(dict2, idx)
            summed_dict[idx]  = sum_tuples(pair1, pair2)
        return summed_dict
    
    summed_dict = reduce(sum_a_dict, dict_list)#[dict1,dict2,dict3])
    summed_dict_idx = summed_dict.keys()
    summed_dict_idx.sort()
    summed_dict_wx, summed_dict_cx = zip(*[(summed_dict[x][0],summed_dict[x][1])  for x in summed_dict_idx])

    return summed_dict_idx, summed_dict_wx, summed_dict_cx

def write_wcfit_file_compact(outfilename, wxfit, cxfit):
    """
    e.g.:
    # output wxfit
    outfilename =  "/Users/shin/workspace/main_Solexa17_100108_HWUSI-EAS610_0004_100125/output/fit_data/test.txt"
    genetrack_lite.write_wcfit_file_compact(outfilename, wxfit, cxfit)
    """
    outfile = open(outfilename,'w')
    counter = 1
    precision_limit = 0.000001
    for i, j in zip(wxfit, cxfit):
        if i>precision_limit and j>precision_limit:
            outfile.write("%d\t%f\t%f\n"%(counter, i,j))
        if i>precision_limit and j<=precision_limit:
            outfile.write("%d\t%f\t0\n"%(counter, i))
        if j>precision_limit and i<=precision_limit:
            outfile.write("%d\t0\t%f\n"%(counter, j))
        counter+=1
    outfile.close()

def get_idx_fwd_rev_dict(infilename, chromosome):
    """
    returns indexed dictionary full of fwd and reverse tag counts
    """
    k_idx, k_fwd, k_rev = get_tag_count_data(infilename, chromosome)
    return dict(zip(k_idx, zip(k_fwd, k_rev)))

def get_peak_pairs(peaks_w, peaks_c, width, variance):
    """
    usage:
    for w_peak_location, c_peak_location, w_peak_height, c_peak_height in genetrack_lite.get_peak_pairs(selected_peaks, selected_peaks_cx, width, variance):
           row = [chromosome]
           row.extend( [ ("%r"%x).replace("'","") for x in [w_peak_location, w_peak_height, c_peak_location, w_peak_height ]])
           outfile_peaks.write("%s\n"%('\t'.join(row)))
   
    """
    #need to sort 
    peaks_w.sort(key=operator.itemgetter(0), reverse=False) #sort by ascending index
    peaks_c.sort(key=operator.itemgetter(0), reverse=False) #sort by ascending index
    
    wx_generator =      (x for x in peaks_w)
    wx_generator_next = (x for x in peaks_w)
    cx_generator =      (x for x in peaks_c)
    
    # determine C-W distance between peaks:
    # rules:
    # given a W tag index, look for downstream C tag.  In order for a C tag to be eligible, it needs tobe
    # a) within 2xWidth
    # b) if over 1xWidth, not overstepping another idx
    idx_pair=[]
    wx_generator_next.next()
    for idx,val in wx_generator:    
        try:
            idx_nxt,val_nxt = wx_generator_next.next()
            #if val_nxt > 
        except: # should only happen a few times
            idx_nxt = None
            val_nxt = None
        for idx_cx,val_cx in cx_generator:
            #if idx_nxt is not None and idx_cx  > idx_nxt and idx_cx-idx > width:
            #    break
            if idx_cx-idx > (width+variance):
                break
            if idx_cx > idx:
                idx_pair.append((idx,idx_cx, val, val_cx))
                break
    
    return idx_pair # w_peak_location, c_peak_location, w_peak_height, c_peak_height
#################################################################################################

def example1():
    """
    Usage example.
    """
    import module_100227_analysis_solexa_peak_detection_module as genetrack_lite
    import module_r_interface_module as rmodule
    import rpy2.robjects as robjects
    
    chromosomes = 'chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 chr13 chr14 chr15 chr16'.split()
    width = 192
    variance = 100 # eyeballed here, make it like 2stddev later.
    threshold = 0
    
    sigma = 20
    
    distance_list = []
    
    
    chromosome = chromosomes[0]
    fileset = [index_file_root%x for x in files[0:3]] # list of filenames
    
    dict_list = [genetrack_lite.get_idx_fwd_rev_dict(infilename, chromosome) for infilename in fileset[0:2]]#[0:3]] # leave out 6
    
    # combine indexed tag-count files into one
    k_idx, k_fwd, k_rev = genetrack_lite.sum_indexed_tag_counts(dict_list) 
    
    wxfit = genetrack_lite.smooth_data(k_idx, k_fwd, sigma, width)
    cxfit = genetrack_lite.smooth_data(k_idx, k_rev, sigma, width)
    
    peak_iter = genetrack_lite.get_peaks_iter(wxfit)
    selected_peaks = genetrack_lite.select_peaks(peak_iter, width=width, level =0) # format: (index, height)
    
    # get top peaks
    # selected_peaks.sort(key=operator.itemgetter(1), reverse=True) #sort by peak height (hi2lo)
    rmodule.set_float_vector_in_R("peak_height",[x[1] for x in selected_peaks])
    
    
    print "now processing rev reads (cxfit)"
    peak_iter_cx = genetrack_lite.get_peaks_iter(cxfit)
    selected_peaks_cx = genetrack_lite.select_peaks(peak_iter_cx, width=width, level=0)
    rmodule.set_float_vector_in_R("peak_height_rev",[x[1] for x in selected_peaks_cx])
    
    
    rmodule.set_x11_window_params(width=16,height=10)
    rmodule.open_x11_window()
    
    rmodule.run_rcode("""
    par(mfrow=c(2,2))
    hist(peak_height,100)
    hist(peak_height_rev,100)
    """)
    k = rmodule.get_array_from_R("dev.list()")
    
    
    selected_peaks.sort(key=operator.itemgetter(0), reverse=False) #sort by ascending index
    selected_peaks_cx.sort(key=operator.itemgetter(0), reverse=False) #sort by ascending index
    wx_generator = (x for x in selected_peaks)
    wx_generator_next = (x for x in selected_peaks)
    cx_generator = (x for x in selected_peaks_cx)
    
    # determine C-W distance between peaks:
    # rules:
    # given a W tag index, look for downstream C tag.  In order for a C tag to be eligible, it needs tobe
    # a) within 2xWidth
    # b) if over 1xWidth, not overstepping another idx
    idx_pair=[]
    wx_generator_next.next()
    for idx,val in wx_generator:    
        try:
            idx_nxt = wx_generator_next.next()
        except:
            idx_nxt = None
        for idx_cx,val_cx in cx_generator:
            #if idx_nxt is not None and idx_cx  > idx_nxt and idx_cx-idx > width:
            #    break
            if idx_cx-idx > (width+variance):
                break
            if idx_cx > idx:
                idx_pair.append((idx,idx_cx, val, val_cx))
                break

    print 'k'
    
    
    rmodule.set_float_vector_in_R("idx_paired_wx", [x[0] for x in idx_pair])
    rmodule.set_float_vector_in_R("idx_paired_cx", [x[1] for x in idx_pair])
    rmodule.set_float_vector_in_R("val_paired_wx", [x[2] for x in idx_pair])
    rmodule.set_float_vector_in_R("val_paired_cx", [x[3] for x in idx_pair])
    
    
    rmodule.set_x11_window_params(width=8,height=6)
    rmodule.open_x11_window()
    rmodule.run_rcode("""
    #idx_paired <- cbind(idx_paired_wx, idx_paired_cx)
    x <- idx_paired_cx-idx_paired_wx
    hist(x, 50, xlab="C-W distance (bp)", ylab="Count", main="Criteria 1: C-W")
    curve(dnorm(x), add=T)
    plot(x, val_paired_wx )
    points(x, val_paired_cx, col='red')
    """)


def example2():
    """
    """
    
    index_file_root = ""
    
    files = """""".split();
    outfilename =  ""
    
    
    sigma = 20
    width = 120
    #epsilon = 0 #some minimum peak height (never used)
    
    lo, hi, normal = normal_function(sigma, width)
    
    
    infilename = index_file_root%files[0]
    
    data = open(infilename).readlines()
    
    k = []
    k_idx = []
    k_fwd = []
    k_rev = []
    for i in data:
        isplit = i.split()
        if isplit[0] == 'chr02':
            k.append(i)
            k_idx.append( int(isplit[1]))
            k_fwd.append( int(isplit[2]))
            k_rev.append( int(isplit[3]))
    
    size = max(k_idx)
    pad = hi + 1
    wxfit = numpy.zeros( size + pad )  # watson fit
    
    #size, epsilon = data_size, minimum_peak_size
    for ix, wx in zip(k_idx, k_fwd):
        if wx > 0.0:
            safe_add( fit=wxfit, ix=ix, lo=lo, hi=hi, values=normal*wx )
    
    
    cxfit = numpy.zeros( size + pad )  # crick fit
    for ix, cx in zip(k_idx, k_rev):
        if cx > 0.0:
            safe_add( fit=cxfit, ix=ix, lo=lo, hi=hi, values=normal*cx )
    
    
    outfile = open(outfilename,'w')
    for i, j in zip(wxfit, cxfit):
        outfile.write("%f\t%f\n"%(i,j))
    outfile.close()
    
    print 'k'
    
    peaks = list(detect_peaks( izip(count(), wxfit)))
    
    selected_peaks = select_peaks(peaks, width=400, level =0)
    
    
    print 'k'

if __name__ == "__main__":
    print "This version is intended to be used as a module.  Command-line options will be added later."