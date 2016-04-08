# CEll Growth ANAlyzer (CEGANA) 1.2
# Copyright (c) 2015 Edouard A. Harris.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this program, included as the two text files
# "COPYING" and "COPYING.LESSER".  If not, see
# <http://www.gnu.org/licenses/>.

# For inquiries or bugs, please contact the corresponding author,
# David R. McMillen, at david.mcmillen@utoronto.ca.


# Applies Bayesian inference to the cell counting process, using the toolkit
# built for BioLogica.

# New in Version 1.2
# - Fixed a bug that caused CEGANA not to plot 2-d marginals when only two
# parameters are to be plotted
# - Fixed a bug that caused CEGANA not to support contour plots of 2-d
# marginals when only two parameters are to be plotted

import warnings
warnings.filterwarnings('ignore',category = DeprecationWarning)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import itertools as it
import time
import pickle as pk
import collections as cl
import random as rd
import scipy.optimize as opt
import sys
import scipy.stats as sts

scriptpath = os.path.abspath(sys.argv[0])
folderpath = os.path.dirname(scriptpath)
parentfolder = os.path.abspath(folderpath)

datafolder = parentfolder + '/Datasets/'
picklefolder = parentfolder + '/Saved Distributions/'

# Main method of CEGANA.
# The data file must contain the count data in five columns, with the first
# row giving the title of each column, and subsequent rows giving the data
# values.  Column 1 is the time ti at which each count was taken; Column 2 is
# the number of cells Mi*yi seeded for the population; Column 3 is the volume
# fraction yi of cells in the parent population that were seeded; Column 4 is the
# number of cells ni counted in the population sample; and Column 5 is the volume
# fraction xi of the population that was sampled for counting.  "Exact" biological
# replicates are handled by assigning a new row with the same ti values.
# Version 1.2 will only be able to handle a flat prior.  The user must input
# his or her estimated value range and grid point density in such a way as to find
# the most probable values efficiently.

def cegana_main(data_file,gridpoints,rmin,rmax,Kmin,Kmax,alphamin,alphamax):
    global datafolder
    start = time.time()

    if not os.path.isfile(datafolder + data_file):
        print "WARNING: the specified data file could not be found."
    
    args = [rmin,rmax,Kmin,Kmax,alphamin,alphamax]

    i = 0
    while i < len(args):
        args[i] = checkfloat(args[i])
        args[i + 1] = checkfloat(args[i + 1])

        if None in [args[i],args[i + 1]]:
            print "WARNING: your domain boundaries should be ordered as follows: minimum growth rate, maximum growth rate, minimum carrying capacity, maximum carrying capacity, minimum extrinsic error, maximum extrinsic error."
        elif args[i + 1] <= args[i]:
            print "WARNING: your domain boundaries must comprise the minimum values, followed by the maximum values, in that order.  The minimum and maximum values cannot be identical."

        i = i + 2

    inputdata = np.genfromtxt(datafolder + data_file,dtype = float,skip_header = 1)

    parameterlist = []
    parameternames = ['Growth rate, r','Carrying capacity, K','Extrinsic error, alpha']

    i = 0
    while i < len(args):
        params = np.linspace(args[i],args[i + 1],gridpoints)
        parameterlist = parameterlist + [params]
        
        i = i + 2

# We define the shape in this generalized way so that in the future we can
# easily accommodate parameter searches with different numbers of gridpoints
# along each dimension (that is, for each parameter).

    inputshape = []
    i = 0
    while i < len(parameterlist):
        inputshape = inputshape + [len(parameterlist[i])]
        i = i + 1

# Total number of calculations that must be done to populate the search space.
# Used to display progress.

    operations = 1
    numbases = []
    for ops in inputshape:
        operations = operations*ops
        numbases = numbases + [operations]

    inputshape = tuple(inputshape)

    logallprobs = np.zeros(operations)

    print "Calculating inference.."

    i = 0
    for parameters in it.product(*parameterlist):
        percentcompletion = round((float(i)/float(operations))*100,1)
        string = str(percentcompletion) + '% complete.'
        sys.stdout.write('\r' + string)
        sys.stdout.flush()

        inputs = list(parameters)

        logtotalprob = 0.

        k = 0

        while k < len(inputdata):

            instanceinputs = list(inputdata[k])
            allinputs = tuple(instanceinputs + inputs + [i,k])#DELETE + [i,k]
            logpreprob = log_likelihood_cellcount_with_counterror(*allinputs)
            logtotalprob = logtotalprob + logpreprob

            k = k + 1

        logallprobs[i] = logtotalprob

        i = i + 1

    logallprobs = np.reshape(logallprobs,inputshape)

    string = '100.0% complete.'
    sys.stdout.write('\r' + string)
    sys.stdout.flush()
    sys.stdout.write('\n')

    maxlogprob = np.max(logallprobs)

    if maxlogprob == -1e309:
        maxlogprob = 0.
        print "Maximum log probability = -inf."
    else:
        print "Maximum log probability = " + str(maxlogprob) + "."

    normlogallprobs = logallprobs - maxlogprob

    normallprobs = np.exp(normlogallprobs)
    trueallprobs = normallprobs/np.sum(normallprobs)

    print ""
    end = time.time()
    elapsed = end - start
    print "Time elapsed: " + str(elapsed) + " s."

    return [parameterlist,parameternames,trueallprobs,normlogallprobs,logallprobs]

# Plots a best-fit line specifically for the cell counting situation.
# Given a data file in the format specified above, data points are
# plotted as the number of cells counted divided by the volume fraction
# of the sample, normalized to seeding number Mi*yi in each case.  The
# model line is a given set of parameters r, K, plotted with a
# seeding number of <Miyi>, plugged into the logistic growth equation,
# and divided by the mean value <Miyi>.

def plot_fit(data_file,r,K):
    global datafolder
    
    data = np.genfromtxt(datafolder + data_file,dtype = float,skip_header = 1)
    data = data.T

    ti = data[0]
    Miyi = data[1]
    ni = data[3]
    xi = data[4]

    Ni = ni/xi

    xdata = ti
    ydata = Ni/Miyi
        
    xmodel = np.linspace(np.min(xdata),np.max(xdata),1000)
    avseed = np.mean(Miyi)
    ymodel = logistic(xmodel,avseed,K,r)/avseed

    plt.plot(xmodel,ymodel,'b-',xdata,ydata,'g.',markersize = 15)
    plt.show()

    return [xdata,ydata,xmodel,ymodel]

# Plots every model tested in distribution against the data contained
# in data_file.  Model curves are plotted according to their
# posterior probability; the higher the probability, the darker the
# curve.  A least squares fit may be optionally plotted, in red, for
# comparison with the model curves.  To do this, set show_leasetsq to
# True.

def plot_all_fits(data_file,distribution,show_leastsq = False):
    global datafolder

    resolution = 100

    data = np.genfromtxt(datafolder + data_file,dtype = float,skip_header = 1)
    data = data.T

    trueallprobs = distribution[2]
    parameters = distribution[0]
    r = parameters[0]
    K = parameters[1]
    
    probabilities = marginalize2d(trueallprobs,0,1)
    probabilities = probabilities/np.max(probabilities)

    ti = data[0]
    Miyi = data[1]
    ni = data[3]
    xi = data[4]

    Ni = ni/xi

    xdata = ti
    ydata = Ni/Miyi

    xmodels = np.linspace(np.min(xdata),np.max(xdata),resolution)
    avseed = np.mean(Miyi)

    n = 0

    colorstrings = []#
    ymodels = np.empty((len(r)*len(K),resolution))
    colornumbers = np.array([],dtype = int)
    while n < len(r):
        m = 0
        while m < len(K):
            ymodel = logistic(xmodels,avseed,K[m],r[n])/avseed
            colornumber = int(probabilities[n][m]*256.)
            if colornumber == 256:
                colornumber = 255
            elif colornumber == 0:
                colornumber = 1

            colornumber = 256 - colornumber

            ymodels[n*len(K) + m] = ymodel
            colornumbers = np.append(colornumbers,colornumber)

            m = m + 1

        n = n + 1

    sortindices = np.argsort(colornumbers)[::-1]
    colornumbers = colornumbers[sortindices]
    ymodels = ymodels[sortindices]

    n = 0
    while n < len(colornumbers):

        colornumber = hex(colornumbers[n])[2:-1]
        
        if colornumbers[n] < 16:
            colorstring = '#0' + str(colornumber) + '0' + str(colornumber) + '0' + str(colornumber)
        else:
            colorstring = '#' + str(colornumber) + str(colornumber) + str(colornumber)

        colorstrings = colorstrings + [colorstring]
        plt.plot(xmodels,ymodels[n],'-',color = colorstring,linewidth = 4)

        n = n + 1

    if show_leastsq:
        def fit_func(x,r,K):
            return logistic(x,avseed,K,r)/avseed

# New in Version 1.1: better initial guesses for r and K.
        rguess = sts.linregress(xdata,np.log(Ni))[0]
        Kguess = np.max(Ni)
        
        rsq,Ksq = opt.curve_fit(fit_func,xdata,ydata,p0 = np.array([rguess,Kguess]))[0]
    
        ylstsq = logistic(xmodels,avseed,Ksq,rsq)/avseed
        plt.plot(xmodels,ylstsq,'r-',linewidth = 2)

        print "Least squares values:"
        print "K = " + str(Ksq)
        print "r = " + str(rsq)
        
    plt.plot(xdata,ydata,'g.',markersize = 15)
    plt.show()

# Saves a posterior distribution in the Pickle Jar under a given
# file_name. This is useful so that full probability distributions
# can be looked at later without having to recalculate the
# inference.

def save_distribution(distribution,file_name):
    global picklefolder

    file_name = str(file_name)

    output = open(picklefolder + file_name + '.pk1','w')
    pk.dump(distribution,output)
    output.close()

# Lists the names of all files stored in the Pickle Jar folder.
# Gives the names, less the .pk1 file extensions.  Note that
# we delete the default Mac file .DS_Store from the list of
# files since it is irrelevant to what we are trying to list.

def list_saved():
    global picklefolder

    filelist = os.listdir(picklefolder)
    
    if '.DS_Store' in filelist:
        backindex = filelist.index('.DS_Store')
        del filelist[backindex]
    
    if len(filelist) == 0:
        print "There are no saved datasets at the moment."
    else:
        print "Saved datasets:"
    
        for name in filelist:
            print name[0:-4]

        print ""

# Loads a distribution with a given name, or emits a warning message
# if no such distribution exists.

def load_distribution(file_name):
    global picklefolder

    file_name = str(file_name)

    if not os.path.isfile(picklefolder + file_name + '.pk1'):
        print "WARNING: The specified file could not be found."
    else:
        file1 = open(picklefolder + file_name + '.pk1')
        distribution = pk.load(file1)
        file1.close()

# New in Version 1.1: return distribution is moved into the else
# statement.
        return distribution

# Given a distribution of the format generated by cellcrunch_main(),
# plots each parameter's probability density against all of the others
# as a 3D surface plot.  plotr, plotK, and plotalpha are Booleans:
# if plotr = True and the others are false, the marginal posterior for
# r is plotted in a 2d plot; if two are true, their joint marginal
# posterior is plotted in a 3d plot; if all three are true, the joint
# marginal posteriors of all 3 possible pairs are plotted in 3d
# plots.  countour set to False makes a surface plot, contour set to
# True generates a contour plot.

def visualize(distribution,plotr,plotK,plotalpha,contour = False):
    
    plotr = checkbool(plotr)
    plotK = checkbool(plotK)
    plotalpha = checkbool(plotalpha)

    whatplots = np.array([plotr,plotK,plotalpha])
    whereplots = np.where(whatplots)[0]

    if None in whatplots:
        print "WARNING: plotr, plotK, and plotalpha should all be Boolean variables, with values of either True or False."
        
    parameterlist = distribution[0]
    parameternames = distribution[1]
    trueallprobs = distribution[2]
    normlogallprobs = distribution[3]
    logallprobs = distribution[4]

# Plots each parameter against all others.  For n parameters, produces
# (1/2)*n*(n - 1) two-parameter plots.

    if len(whereplots) == 1:
        dist1d = marginalize1d(trueallprobs,whereplots[0])
        
        plt.plot(parameterlist[whereplots[0]],dist1d,'b-')
        plt.xlabel(parameternames[whereplots[0]])
        plt.ylabel('Probability density')
        plt.show()
        
    elif len(whereplots) == 2:
        dist2d = marginalize2d(trueallprobs,whereplots[0],whereplots[1])

        fig = plt.figure(1)
        X = parameterlist[whereplots[1]]
        Y = parameterlist[whereplots[0]]
        X, Y = np.meshgrid(X,Y)
        Z = dist2d

# New in Version 1.2: added support for contour plots to the
# len(whereplots) == 2 case.
        if contour:
            cs = plt.contour(X,Y,Z,colors = 'k')
            plt.xlabel(parameternames[j])
            plt.ylabel(parameternames[i])
            plt.clabel(cs)

        else:
            ax = fig.gca(projection = '3d')
            ax.set_xlabel(parameternames[whereplots[1]])
            ax.set_ylabel(parameternames[whereplots[0]])
            ax.set_zlabel('Probability density')
            
            surf = ax.plot_surface(X,Y,Z,rstride = 1,cstride = 1,linewidth = 0,antialiased = False)
# New in Version 1.2: added plt.show().
        plt.show()

    elif len(whereplots) == 3:
        i = 0
        k = 1
        while i < len(parameterlist):
            j = i + 1
            while j < len(parameterlist):
                dist2d = marginalize2d(trueallprobs,i,j)

                fig = plt.figure(k)
                X = parameterlist[j]
                Y = parameterlist[i]
                X, Y = np.meshgrid(X,Y)
                Z = dist2d
                
                if contour:
                    cs = plt.contour(X,Y,Z,colors = 'k')
                    plt.xlabel(parameternames[j])
                    plt.ylabel(parameternames[i])
                    plt.clabel(cs)
                else:
                    ax = fig.gca(projection = '3d')
                    ax.set_xlabel(parameternames[j])
                    ax.set_ylabel(parameternames[i])
                    ax.set_zlabel('Probability')
                    
                    surf = ax.plot_surface(X,Y,Z,rstride = 1,cstride = 1,linewidth = 0,antialiased = False)
                
                k = k + 1
                j = j + 1

            i = i + 1

        plt.show()

# Produces a one-dimensional marginalization of an n >= 1-
# dimensional array.  We add up the values of the array along
# every axis except the one whose index is given in the inputs.
# Used chiefly for the visualize() function, but may have a
# more general utility.

def marginalize1d(input_array,index):
    n = len(np.shape(input_array))
    
    if index >= n:
        print "WARNING: index is too large for the dimensionality of your input."
    elif len(np.shape(input_array)) == 0:
        print "WARNING: the input array must have at least one dimension."

    if len(np.shape(input_array)) == 1:
        output_array = input_array
    else:
        input_array = np.swapaxes(input_array,0,index)

        while len(np.shape(input_array)) > 1:
            input_array = np.sum(input_array,1)

        output_array = input_array.T

    return output_array


# Produces a two-dimensional marginalization of an n >= 2-
# dimensional array.  We add up the values of the array along
# every axis except the ones whose indices are given in the
# inputs. Used chiefly for the visualize() function, but may
# have a more general utility.

def marginalize2d(input_array,index1,index2):
    n = len(np.shape(input_array))
    
    if index1 >= n or index2 >= n:
        print "WARNING: one or both of your indices is too large for the dimensionality of your input."
    elif len(np.shape(input_array)) < 2:
        print "WARNING: the input array must have at least two dimensions."

    if len(np.shape(input_array)) == 2:
        output_array = input_array.T
    else:
        input_array = np.swapaxes(input_array,0,index1)

        if index2 == 0:
            input_array = np.swapaxes(input_array,1,index1)
        else:
            input_array = np.swapaxes(input_array,1,index2)

        while len(np.shape(input_array)) > 2:
            input_array = np.sum(input_array,2)

        output_array = input_array

    return output_array

def log_likelihood_cellcount_with_counterror(ti,Miyi,yi,ni,xi,r,K,alpha,i,k):#DELETE i,k
    randomthresh = 0.00#

    randomcheck = (rd.random() < randomthresh)#
    
    def gfuncN(Ni):
        return ((Ni*xi)**2.)*((alpha**2.) + (1. - xi)/(Ni*xi))

    def gfuncM():
        return ((Mi*yi)**2.)*((alpha**2.) + (1. - yi)/(Mi*yi))

    def ffuncN(Ni):
        return (K*Ni*np.exp(-r*ti))/(K + Ni*(np.exp(-r*ti) - 1.))

    def log_prob(Ni):
# New in Version 1.1: see just below.
        scalar = False
        
        try:
            len(Ni)
        except TypeError:
            Ni = np.array([Ni])
# New in Version 1.1: detect the scalar nature of an Ni input.
            scalar = True

        ans = np.empty(len(Ni))
        ans.fill(-1e309)
            
# New in Version 1.1: handles the cases in which alpha == 0 (no extrinsic
# error) and either xi == 1 (the sample is the entire population), or
# yi == 1 (the entire parent flask is seeded), or both.  These are the
# next three if and elif statements.  The else statement will be the
# most common, nominal case.
        if alpha == 0. and xi == 1. and yi != 1.:
# New in Version 1.1: if xi == 1 and alpha == 0, then P(ni | Ni, xi, alpha)
# is just a Kroenecker delta(ni, Ni).  There is only one nonzero term
# in the sum, corresponding to Ni == ni.
            gfM = gfuncM()
            ffN = ffuncN(ni)

# If ni == 0, ffuncN(ni)/ni is undefined, since it's zero divided by zero.
# We therefore calculate the limit manually.
            if ni == 0.:
                nonzeroterm = -0.5*np.log(gfM) + 2.*np.log(K*np.exp(-r*ti)/K) + r*ti - ((Mi*yi)**2.)/(2.*gfM)
            else:
                nonzeroterm = -0.5*np.log(gfM) + 2.*np.log(ffN/ni) + r*ti - ((ffN - Mi*yi)**2.)/(2.*gfM)

            ans[Ni == ni] = nonzeroterm

# New in Version 1.1: if yi == 1 and alpha == 0, then P(Ni | ...) is just
# a Kroenecker delta(Ni, f(Mi)).  Note that we are just setting this
# probability as a Kroenecker, without including the extra multiplicative
# terms that would result from taking the explicit limit.  This is because
# in this special case we are treating Ni as an integer, rather than as
# being drawn from a continuum (the change of variables factor would apply
# in this latter case).
        elif alpha == 0. and xi != 1. and yi == 1.:
            N = float(int((K*Mi*np.exp(r*ti))/(K + Mi*(np.exp(r*ti) - 1.))))
            gfN = gfuncN(N)

            nonzeroterm = -0.5*np.log(gfN) - ((ni - N*xi)**2.)/(2.*gfN)

            ans[Ni == N] = nonzeroterm

# New in Version 1.1: if alpha == 0, xi == 0 and yi == 0, then there is
# perfect seeding and sampling.  Therefore, the log probability is zero
# (so the probability is unity) if and only if the number sampled is
# exactly equal to the logistic growth from seeding number Mi, and
# zero in every other case.
        elif alpha == 0. and xi == 1. and yi == 1.:
            if ni == float(int((K*Mi*np.exp(r*ti))/(K + Mi*(np.exp(r*ti) - 1.)))):
                ans[Ni == ni] = 0.

        else:
            leftlimit = ni
            if K == 0.:
                rightlimit = 0.
            elif r == 0.:
                rightlimit = 1e309
            else:
                rightlimit = float(int(K/(1. - np.exp(-r*ti))))
                
# Note that we subtract 1 from rightlimit because the log_prob function
# breaks if we evaluate it at the edge of the domain and it happens
# to be an exact integer by incredible coincidence (but has happened!).
# Note also that we add the constraint Ni > 0 in Version 1.1, so that
# we can handle the ni == 0 case separately.
            finiteindex = np.where((Ni > 0.) & (Ni >= leftlimit) & (Ni < rightlimit) & (Ni != 1e309))[0]

            gfN = gfuncN(Ni[finiteindex])
            gfM = gfuncM()
            ffN = ffuncN(Ni[finiteindex])

            ans[finiteindex] = -0.5*np.log(gfN*gfM) + 2.*np.log(ffN/Ni[finiteindex]) + r*ti - ((ni - Ni[finiteindex]*xi)**2.)/(2.*gfN) - ((ffN - Mi*yi)**2.)/(2.*gfM)

# New in Version 1.1.  If ni == 0, P(ni | Ni, xi, alpha) is forcefully
# set to 1 for Ni == 0.  If this were not the case, gfuncN() would give
# a divide by 0 error.
            if ni == 0.:
                zeroindex = np.where(Ni == 0.)[0]
                ans[zeroindex] = -0.5*np.log(gfM) - ((ffuncN(Ni[zeroindex]) - Mi*yi)**2.)/(2.*gfM) + r*ti + 2.*np.log(K*np.exp(-r*ti)/(K + Ni[zeroindex]*(np.exp(-r*ti) - 1.)))

# New in Version 1.1: revert to scalar iff the input was a scalar.
        if scalar:
            ans = ans[0]
        
        return ans


# salamibite is the maximum allowed size of an array.  Around
# 5 million entries are comfortably handled by numpy at any
# one time, so this is the default.
# sparsemultiplier is inversely proportional to the sparse
# evaluation resolution.  A sparsemultiplier of 100. means
# that 1% of points are sampled in the sparse evaluation;
# a sparsemultiplier of 200. means that 0.5% of points are
# sampled, etc.
    sparsemultiplier = 5000
    logthreshold = 20.
    salamibite = 5e6
    deltax = 1.
    maxerror = 0.001
    
    Mi = float(Miyi)/float(yi)

    maximum1 = ni
    maximum2 = float(int(ni/xi))
    
# Regardless of what r is, K = 0 means that all cells just die.
# We force this expicitly because, if r = 0 as well, the formulae for
# maximum3 and maximum4 are undefined until we state the order of our
# limits.
    if K == 0.:
        maximum3 = 0.
        maximum4 = 0.
    else:
        maximum3 = float(int((K*Mi*yi*np.exp(r*ti))/(K + Mi*yi*(np.exp(r*ti) - 1.))))
# If r == 0, then there is no change in the population at all: it
# should remain at the seeding number, which is unknown and
# unbounded.  Hence, the absolute maximum cell number is infinity.
        if r == 0:
            maximum4 = 2.*max(maximum2,maximum3)
        else:
# This is the limit N0 -> inf.
            maximum4 = float(int(K/(1. - np.exp(-r*ti))))

    if maximum1 > maximum4:
        return -1e309

# Begin with a sparse evaluation.
    Ni = np.arange(maximum1,maximum4 + sparsemultiplier*deltax,sparsemultiplier*deltax)
    evalvalues = log_prob(Ni)

# Identify the sparse global maximum.
    Nimax = np.argmax(evalvalues)

# Carry out a dense evaluation of points on either side of the sparse maximum.
# New in Version 1.1: handle the exceptional case in which Ni has a length of
# unity.
    if len(Ni) == 1:
        highpoints = Ni
    elif Nimax == 0:
        highpoints = np.arange(Ni[Nimax],Ni[Nimax + 1],deltax)
    elif Nimax == len(Ni) - 1:
        highpoints = np.arange(Ni[Nimax - 1],Ni[Nimax],deltax)
    else:
        highpoints = np.arange(Ni[Nimax - 1],Ni[Nimax + 1],deltax)

    highvals = log_prob(highpoints)

# Identify the dense global maximum.  This is the true maximum.
    overallmaxvalue = np.max(highvals)

# Use the dense maximum to find the floor above which the integral will be
# densely evaluated.
    evalfloor = overallmaxvalue - logthreshold

# Identify all sparse points that lie above and just adjacent to the evalfloor.
# We do this by finding the indices of all points that are either above the
# evalfloor, or have immediate neighbours that lie above the evalfloor.
    padleft = np.concatenate((np.array([-1e309,-1e309]),evalvalues))
    padleft[-1] = -1e309
    padright = np.concatenate((evalvalues,np.array([-1e309,-1e309])))
    padright[0] = -1e309
    padboth = np.concatenate((np.array([-1e309]),evalvalues,np.array([-1e309])))

    sparseabove = np.where((padleft > evalfloor) | (padright > evalfloor) | (padboth > evalfloor))[0]
# To compensate for the padding of infinities.
    sparseabove = sparseabove - 1

# If there are no sparse points above threshold (which can occur in a region
# of rapid change, if the global max picked out by the highvals is way higher
# than the maximum of the sparse points), then we just look at the highvals
# and highpoints.
    if len(sparseabove) == 0:
        fullNi = highpoints
        
        fullevalvalues = highvals

# New in Version 1.1: if the maximum likelihood on the Ni domain is zero, then
# naturally the summed likelihood over the Ni domain must be nil as well.
        if overallmaxvalue == -1e309:
            normlikelihood = 0.
        else:
            fullevalvalues = fullevalvalues - overallmaxvalue

            normlikelihood = np.sum(np.exp(fullevalvalues))*deltax

    else:
# Find differences between successive indices.
        jumplist = np.ediff1d(sparseabove)

# Identify big jumps; these correspond to hops between domains above the evalfloor.
        bigjumps = np.where(jumplist > 1)[0]

        if evalvalues[0] > evalfloor:
            lowerlimit = Ni[0]
            lowerindex = 0
        else:
            lowerlimit = Ni[sparseabove[0]]
            lowerindex = sparseabove[0]

        if evalvalues[-1] > evalfloor:
            upperlimit = Ni[-1]
            upperindex = -1
        else:
            upperlimit = Ni[sparseabove[-1]]
            upperindex = sparseabove[-1]

        if len(bigjumps) == 0:
            limits = [lowerlimit,upperlimit]
            indices = [lowerindex,upperindex]

        else:
            limits = [lowerlimit,Ni[sparseabove[bigjumps[0]]]]
            indices = [lowerindex,sparseabove[bigjumps[0]]]
            
            j = 1
            while j < len(bigjumps):
                limits = limits + [Ni[sparseabove[bigjumps[j - 1] + 1]],Ni[sparseabove[bigjumps[j]]]]
                indices = indices + [sparseabove[bigjumps[j - 1] + 1],sparseabove[bigjumps[j]]]
                j = j + 1

            limits = limits + [Ni[sparseabove[bigjumps[j - 1] + 1]],upperlimit]
            indices = indices + [sparseabove[bigjumps[j - 1] + 1],upperindex]
        
        normlikelihood = 0.
        j = 0
        while j < len(limits):
            sparselocalvals = evalvalues[indices[j]:indices[j + 1] + 1]
            sparselocalNi = Ni[indices[j]:indices[j + 1] + 1]

            fullNi = dynamic_approximator(sparselocalNi,sparselocalvals,deltax,maxerror,log_prob)#delete log_prob

            allmultipliers = np.ediff1d(fullNi)
            
            fullevalvalues = log_prob(fullNi)
            fullevalvalues = fullevalvalues - overallmaxvalue

# We take the average of successive fullevalvalues to get a centered
# Riemann sum.
            fullexpvalues = np.exp(fullevalvalues)
            localnormlikelihood = np.sum(0.5*(fullexpvalues[:-1] + fullexpvalues[1:])*allmultipliers)*deltax
            normlikelihood = normlikelihood + localnormlikelihood
            
            j = j + 2

# New in Version 1.1: if normlikelihood is zero, then naturally the logarithm
# of the likelihood must be negative infinity.
    if normlikelihood == 0.:
        loglikelihood = -1e309
    else:
        lognormlikelihood = np.log(normlikelihood)
        loglikelihood = lognormlikelihood + overallmaxvalue

## DELETE BELOW
    if randomcheck:
        print ""
        print ""
        print "i = " + str(i)
        print "k = " + str(k)
        print str((ti,Miyi,yi,ni,xi,r,K,alpha,i,k))
        print "domain size = " + str(maximum4 - maximum1)
        
        Ni2 = arange_salami(maximum1, maximum4 + deltax,deltax,salamibite)
        
        evalvalues2 = salami_function(Ni2,log_prob)

        overallmaxvalue2 = salami_max(evalvalues2)
        evalvalues2 = salami_subtract(evalvalues2,overallmaxvalue2)
        normlikelihood2 = salami_sum(salami_function(evalvalues2,np.exp))*deltax
        lognormlikelihood2 = np.log(normlikelihood2)
        loglikelihood2 = lognormlikelihood2 + overallmaxvalue2

        logerror = (loglikelihood2 - loglikelihood)/abs(loglikelihood)
        error = (np.exp(loglikelihood2) - np.exp(loglikelihood))/np.exp(loglikelihood)

        print "log error = " + str(logerror)
        print "lin error = " + str(error)
        print ""
## DELETE ABOVE
            
    return loglikelihood

# The next four methods are designed to check whether the given input
# has the indicated type, or can be converted to the indicated type.
# If yes, return the converted variable; else, return None.

def checkstr(x):
    try:
        out = str(x)
        return out
    except ValueError:
        pass

def checkint(x):
    try:
        out = int(x)
        return out
    except ValueError:
        pass

def checkfloat(x):
    try:
        out = float(x)
        return out
    except ValueError:
        pass

def checkbool(x):
    try:
        out = bool(x)
        return out
    except ValueError:
        pass

# Chops up a very large array into bite sized chunks that can
# then be fed into numpy for greater memory efficiency.
def array_salami(array,bitesize):
    bitesize = int(bitesize)
    salami = []
    i = 0
    while i < len(array)/bitesize:
        salami = salami + [array[i*bitesize:(i + 1)*bitesize]]
        i = i + 1

    if len(array)%bitesize != 0:
        salami = salami + [array[i*bitesize:len(array)]]

    return salami

def arange_salami(start,stop,step,bitesize):
    bitesize = int(bitesize)
    step = float(step)
    salami = []
    i = 0
    numsteps = int((stop - start)/step)
    while i < abs(numsteps/bitesize):
        salami = salami + [np.arange(start + i*step*bitesize,start + (i + 1)*step*bitesize,step)]
        i = i + 1

    if abs(numsteps%bitesize) != 0:
        salami = salami + [np.arange(start + i*step*bitesize,stop,step)]

    return salami

def salami_function(array_salami,function):
    salami = []
    for entry in array_salami:
        salami = salami + [function(entry)]

    return salami

def salami_max(array_salami):
    maxes = []
    for entry in array_salami:
        maxes = maxes + [np.max(entry)]

    overallmax = np.max(np.array(maxes))

    return overallmax

# Subtracts number from all entries of array_salami.
def salami_subtract(array_salami,number):
    salami = []
    for entry in array_salami:
        salami = salami + [entry - number]

    return salami

def salami_sum(array_salami):
    overallsum = 0.
    for entry in array_salami:
        overallsum = overallsum + np.sum(entry)

    return overallsum

# Returns a single argument, that labels the element of the
# salami array as though it were a regular array.
def salami_argmax(array_salami):
    maxes = np.zeros(len(array_salami))
    maxargs = []
    i = 0
    while i < len(array_salami):
        maxargs = maxargs + [np.argmax(array_salami[i])]
        maxes[i] = array_salami[i][maxargs[i]]
        i = i + 1

    outerarg = np.argmax(maxes)
    fullarg = outerarg*len(array_salami[0]) + maxargs[outerarg]

    return fullarg

# Similar to np.where(), condition() is a function that evaluates
# elements of the array_salami to a Boolean value.
def salami_where(array_salami,condition):
    indices = []
    i = 0
    while i < len(array_salami):
        indices = indices + [i*len(array_salami) + np.where(condition(array_salami))]
        i = i + 1

    return indices

# Pulls the nth element from the salami array with n numbered
# as though it were a normal array.  Array overflow error is
# handled naturally by Python itself.
def salami_element(array_salami,n):
    outerarg = n/len(array_salami[0])
    innerarg = n%len(array_salami[0])

    return array_salami[outerarg][innerarg]

def salami_ediff1d(array_salami):
    diffsalami = []
    i = 0
    while i < len(array_salami):
        slicei = np.ediff1d(array_salami[i])
        
        if i < len(array_salami) - 1:
            lastentry = array_salami[i + 1][0] - array_salami[i][-1]
            slicei = np.append(slicei,lastentry)

        diffsalami = diffsalami + [slicei]

        i = i + 1

    return diffsalami

# Builds a new approximation that can be used to integrate function
# over the domain given by Ni.  It assigns the highest density of
# points to regions of function that have the highest absolute value
# of the derivative (and so are changing the fastest and need to be
# sampled the most), given that function(Ni) = evalvalues.  There is
# a maximum limit to the density of points given by deltax, the
# minimum distance between successive points.  maxerror is the maximum
# acceptable amount of error, NORMALIZED TO THE LENGTH OF THE LINE SEGMENT
# BEFORE IT IS PASSED, which is calculated by looking at the
# upside Riemann difference; often it will actually be less than len(Ni).
# NOTE: it is assumed that the spacings between the Ni are the same,
# i.e., that Ni[i] - Ni[i - 1] == Ni[j] - Ni[j - 1] for all i, j.
# Hence, the relative derivatives can be evaluated as the differences
# between successive evalvalues, without dividing by Ni, which is
# computationally cheaper to do.
def dynamic_approximator(Ni,evalvalues,deltax,maxerror,function):

    slices = 10
    
    deltasparse = Ni[1] - Ni[0]
    maxerror = maxerror*deltasparse

    derivatives = np.abs(np.ediff1d(evalvalues)/deltasparse)
    infinds = np.where(derivatives == 1e309)[0]
    finiteinds = np.where(derivatives != 1e309)[0]
    derivfinite = derivatives[finiteinds]

    rawnumpoints = np.zeros(len(derivatives))
    rawnumpoints[infinds] = deltasparse
    rawnumpoints[finiteinds] = (maxerror/(derivfinite*deltax))*(np.sqrt(1. + ((deltasparse*derivfinite*deltax/maxerror)**2.)) - 1.)

    denseareas = Ni[rawnumpoints >= 1.]
    densenumpoints = rawnumpoints[rawnumpoints >= 1.]
    densemultipliers = deltasparse/(densenumpoints*deltax)
    densemultipliers = densemultipliers.astype(int)

    if len(densenumpoints) == 0:
        denseNi = np.array([])
    else:
        if len(denseareas) <= slices:
            maxsize = int(np.max(densenumpoints)) + 1
            denseNi = np.empty((len(densenumpoints),maxsize - 1))
            denseNi[:] = np.arange(1.,maxsize,1.)
            denseNi = denseNi.T

            denseNi = denseNi*densemultipliers*deltax
            denseNi[denseNi > deltasparse] = np.nan
            denseNi = denseNi + denseareas

            denseNi = np.reshape(denseNi,np.size(denseNi))
            denseNi = denseNi[np.isnan(denseNi) == False]

        else:
            slicesize = len(denseareas)/slices

            sortindices = np.argsort(densenumpoints)
            densenumpoints = densenumpoints[sortindices]
            densemultipliers = densemultipliers[sortindices]
            denseareas = denseareas[sortindices]
            
            denseNi = np.array([])

            i = 0
            while i < slices:
                if i < slices - 1:
                    densenumpointslice = densenumpoints[i*slicesize:(i + 1)*slicesize]
                    densemultiplierslice = densemultipliers[i*slicesize:(i + 1)*slicesize]
                    denseareaslice = denseareas[i*slicesize:(i + 1)*slicesize]
                else:
                    densenumpointslice = densenumpoints[i*slicesize:]
                    densemultiplierslice = densemultipliers[i*slicesize:]
                    denseareaslice = denseareas[i*slicesize:]
                
                maxsizeslice = int(np.max(densenumpointslice)) + 1
                denseNislice = np.empty((len(densenumpointslice),maxsizeslice - 1))
                denseNislice[:] = np.arange(1.,maxsizeslice,1.)
                denseNislice = denseNislice.T

                denseNislice = denseNislice*densemultiplierslice*deltax
                denseNislice[denseNislice > deltasparse] = np.nan
                denseNislice = denseNislice + denseareaslice

                denseNislice = np.reshape(denseNislice,np.size(denseNislice))
                denseNislice = denseNislice[np.isnan(denseNislice) == False]

                denseNi = np.concatenate((denseNi,denseNislice))

                i = i + 1

    sparseareas = Ni[rawnumpoints < 1.]

    if len(sparseareas) == 0:
        sparseNi = np.array([])
    else:
        sparsemean = np.mean(rawnumpoints[rawnumpoints < 1.])
        
        if sparsemean == 0.:
            sparseNi = np.array([])
        else:
            sparsemultiplier = int(deltasparse/(sparsemean*deltax))

            sparseNi = np.arange(Ni[0],Ni[-1] + sparsemultiplier*deltax,sparsemultiplier*deltax)

# We make sure to include also the first and last points explicitly,
# because when we take the average of the exponentials, we need to make
# sure that the low values on either side of the high values are
# included, otherwise the sum will "miss" part of the high value.
    fullNi = np.concatenate((denseNi,sparseNi,np.array([Ni[0],Ni[-1]])))
    fullNi = np.unique(fullNi)

    return fullNi
    
# The logistic growth function.  Takes arrays or floats as inputs.

def logistic(input_data,n0,K,r):
    output = K*n0*np.exp(r*input_data)/(K + n0*(np.exp(r*input_data) - 1))

    return output
