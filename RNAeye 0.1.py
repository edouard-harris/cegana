# RNAeye v 0.1; Created by Ed Harris.

from pylab import *
from datetime import datetime
import os
import urllib
import time
import random
from scipy.optimize import curve_fit

# If you are opening RNAeye for the first time, please enter a homefolder
# to specify the file path to RNAeye on your hard disk; this will activate
# the program.  Failure to specify this file path correctly will
# result in an error message and may prevent RNAeye from processing some
# commands correctly.
# Enter the file path in single quotes, as in the following example:
# homefolder = '/Users/edharris/Documents/RNAeye_0.1/'

homefolder = '/Users/edharris/Documents/PhD_Project/Data/RNAeye_0.1/'

# The following are default parameter values for threshold alignment scores
# and threshold energies of the miRanda algorithm used by RNAeye to compute
# interference hits.  These are the defaults recommended for local
# inplementations of miRanda; however, they may be changed at the convenience
# of the user.

thresholdenergy = -7
thresholdalign = 150














mirandalocation = 'miRanda-3.3a/src/miranda'

if os.path.isfile(homefolder + 'User/RNAeye 0.1.py'):
    print "Welcome to RNAeye."
else:
    print "Error: RNAeye's home folder is corrupted or the file path was incorrectly specified. Some functions may not work as indicated."

# List of functions currently present in active memory.
activefunctions = []
activedatasets = []

#################### SOME PARAMETERS BELOW ###############################

# Generalized fudge factor for floating-point comparisons.
fudge = 0.0000001

# Gibbs free energies of nearest-neighbors for TD stability calculations.
# From Wiki: Nucleic Acid Thermodynamics; primary source in footnote [1].
AAE = -4.26
ATE = -3.67
TAE = -2.5
CAE = -6.12
GTE = -6.09
CTE = -5.4
GAE = -5.51
CGE = -9.07
GCE = -9.36
GGE = -7.66

# end: the ends of the target are compared for thermodynamic stability. The
# variable 'end' is the number of nucleotides that are included in the
# definition of the ends of the target (same number of nucleotides are used for
# 3' and 5' ends). Should be an integer between 1 and 19 (more genenerally,
# between 1 and the length of the target).  Nominally (according to Schwarz et
# al., 2003), end should be around 4, since the siRNA duplex is unwound by
# a helicase that can 'see' about 4 bases at a time.
end = 4

# stab: the weighting factor for the relative stability of the 3' end of a
# target as compared to its 5' end. The more stable the 5' end is compared
# to the 3' end, the better. Should be a floating-point number >=0.
stab = 0.01

# ind18: the weighting factor for the presence of an A or a T at position 18 of
# the target sequence. A is better than T; this factor determines how much better.
# Should be a floating-point number >=0.
ind18 = 1000.

# ind7: the weighting factor for the presence of a T residue at position 7 of the
# target sequence. T is better than anything else; this factor determines how
# much better. Should be a floating-point number >=0.
ind7 = 1000.

# flank: the weighting factor for the presence of A or T residues at the 4
# nucleotides flanking the 3' end of a target sequence. At each of those
# positions, A or T are better than anything else; this factor determines how
# much better. Should be a floating-point number >=0.
flank = 1.

# This subroutine downloads promoter sequences of genes whose refseq IDs are
# listed in the file homefolder + 'USER Gene IDs.txt'. It saves them in txt
# format in the file folder homefolder + 'Gene_Sequences'. The download is the
# entire raw sequence 10 kb upstream of the start codon of interest, downloaded
# from mybioinfo.info.

#########################################################################

def DownloadPromoters(geneids,homefolder):
    i = 0
    outflag = True
    while i < len(geneids):
        if os.path.isfile(homefolder + 'Gene_Sequences/' + geneids[i] + 'pro.fasta') == False:
            print "Downloading promoter for sequence " + geneids[i] + "."

        
            f = urllib.urlopen('http://www.mybioinfo.info/display.php?search_type=dna_refseq&search_key=' + geneids[i] + '&B3=Search')
            query = f.read()

            Downloadstart = query.find('This gene is on') + 42
            Downloadend = query.find('">',Downloadstart)
            if Downloadstart == 41:
                print "RNAeye was unable to find the promoter for target sequence " + geneids[i] + "."
                outflag = False

            else:
                Download = query[Downloadstart:Downloadend]
                f = urllib.urlopen(Download)
                rawseq = f.read()

                returns = []
                lastreturn = rawseq.rfind('\n')
                returnpos = 0

                while returnpos < lastreturn:
                    returnpos = rawseq.find('\n', returnpos + 1)
                    returns = returns + [returnpos]

                intro = '>' + geneids[i] + ' promoter.\n'
                j = 0
                finalseq = ''
                while j < (len(returns) - 1):
                    finalseq = finalseq + rawseq[returns[j] + 1:returns[j + 1]]
                    j = j + 1

                finalseq = finalseq[-1200:-200]

                seqfile = open(homefolder + 'Gene_Sequences/' + geneids[i] + 'pro.fasta', 'w')
                seqfile.write(intro)
                seqfile.write(finalseq)
                seqfile.close()

        i = i + 1

    return(outflag)

# This subroutine downloads transcripts from the refseq database using the same
# technique as the promoter downloader above. (Also from mybioinfo.info.)

def DownloadTranscripts(geneids,homefolder):
    i = 0
    outflag = True
    while i < len(geneids):
        if os.path.isfile(homefolder + 'Gene_Sequences/' + geneids[i] + 'trans.fasta') == False:
            print "Downloading transcript for sequence " + geneids[i] + "."

            f = urllib.urlopen('http://www.mybioinfo.info/display.php?search_type=dna_refseq&search_key=' + geneids[i] + '&B3=Search')
            query = f.read()

            marker = query.find(geneids[i])
            reducedsoup = query[0:marker]
            Downloadstart = reducedsoup.rfind('<a href="') + 9
            Downloadend = query.find('">',marker)

            if marker == -1 or Downloadstart == 8 or query.find('No result found for') != -1:
                 print "RNAeye was unable to find the transcript for target sequence " + geneids[i] + "."
                 outflag = False

            else:
                Download = query[Downloadstart:Downloadend]
                f = urllib.urlopen(Download)
                rawseq = f.read()

                returns = []
                lastreturn = rawseq.rfind('\n')
                returnpos = 0

                while returnpos < lastreturn:
                    returnpos = rawseq.find('\n', returnpos + 1)
                    returns = returns + [returnpos]

                j = 0
                finalseq = '>' + geneids[i] + ' transcript.\n'
                while j < (len(returns) - 1):
                    finalseq = finalseq + rawseq[returns[j] + 1:returns[j + 1]]
                    j = j + 1

                finalseq = finalseq.replace('T','U')
                
                seqfile = open(homefolder + 'Gene_Sequences/' + geneids[i] + 'trans.fasta', 'w')
                seqfile.write(finalseq)
                seqfile.close()

        i = i + 1

    return(outflag)

# This subroutine checks whether the lengths of the Gene IDs.txt and
# microRNA IDs.txt files are the same; this must be done because the variations
# of miRandatuRbo run by the script associate a unique miRNA to a gene
# promoter or transcript, so this subroutine is required to catch index out of
# bounds errors before they occur.

def lengthcheck(microids,geneids):
    flag = True
    
    if len(geneids) != len(microids):
        flag = False

    return(flag)

# This subroutine downloads miRNA sequences from miRBase.  It begins by
# searching the site directly, then goes to an miRNA's home page if that page
# is accessible from the search page for a particular miRNA.  If the identifier
# -3p or -5p is given to a miRNA, the program will search for that sequence
# specifically; if not, it will search for the sequence of the strand that is
# thought to be the likeliest candidate for the guide strand (the strand that
# is actually complementary to the 3' UTR of interest).  This strand is
# identified as that without a star (*) in its previous ID, as this symbol is
# used to denote the passenger strand.  The format of search queries in miRBase
# is highly inconsistent, so this program incorporates many different checks
# at each stage of the search.

def DownloadMicros(microids,homefolder):
    i = 0
    identifiers = []
    indices = []
    while i < len(microids):
        if os.path.isfile(homefolder + 'microRNA_Sequences/' + microids[i] + '.fasta') == False:
            indices = indices + [i]
            print "Downloading sequence for microRNA " + microids[i] + "."

            f = urllib.urlopen('http://www.mirbase.org/cgi-bin/query.pl?terms=' + microids[i])
            query = f.read()

            initsearch = query.find('Mature sequence')
            flag = True
            link = 'http://www.mirbase.org/cgi-bin/query.pl?terms=' + microids[i]
    
            if initsearch == -1:
                query = query.lower()
                reducedstart = query.find('restore to the original order')
                basicmicroid = microids[i].lower()
        
                if microids[i].find('3p') != -1 or microids[i].find('5p') != -1:
                    basicmicroid = basicmicroid[0:len(basicmicroid) - 3]
                    locator = query.find('>' + basicmicroid + '<', reducedstart)
                    if locator == -1:
                        locator = query.find('>' + basicmicroid + ' <', reducedstart)
                    if locator == -1:
                        locator = query.find('>' + basicmicroid + '*<', reducedstart)

                    if locator != -1:
                        reducedsoup = query[reducedstart:locator + 1]
                        linkstart = reducedsoup.rfind('<a href="') + 9
                        linkend = reducedsoup.rfind('">')
                        link = reducedsoup[linkstart:linkend]

                    elif query.find(basicmicroid + '-', reducedstart) != -1:
                        locator = query.find(basicmicroid + '-', reducedstart)
                        reducedsoup = query[reducedstart:locator]
                        linkstart = reducedsoup.rfind('<a href="') + 9
                        linkend = reducedsoup.rfind('">')
                        link = reducedsoup[linkstart:linkend]

                    elif query.find(basicmicroid + 'a', reducedstart) != -1:
                        locator = query.find(basicmicroid + 'a', reducedstart)
                        reducedsoup = query[reducedstart:locator]
                        linkstart = reducedsoup.rfind('<a href="') + 9
                        linkend = reducedsoup.rfind('">')
                        link = reducedsoup[linkstart:linkend]

                    else:
                        flag = False
                        identifiers = identifiers + ['Not found']
                        print "RNAeye was unable to find the microRNA " + microids[i] + " in the search page."
                
                else:
                    locator = query.find('>' + basicmicroid + '<', reducedstart)
                    if locator == -1:
                        locator = query.find('>' + basicmicroid + ' <', reducedstart)
                    if locator == -1:
                        locator = query.find('>' + basicmicroid + '*<', reducedstart)

                    if locator != -1:
                        reducedsoup = query[reducedstart:locator + 1]
                        linkstart = reducedsoup.rfind('<a href="') + 9
                        linkend = reducedsoup.rfind('">')
                        link = reducedsoup[linkstart: linkend]

                    elif query.find(basicmicroid + '-', reducedstart) != -1:
                        locator = query.find(basicmicroid + '-', reducedstart)
                        reducedsoup = query[reducedstart:locator]
                        linkstart = reducedsoup.rfind('<a href="') + 9
                        linkend = reducedsoup.rfind('">')
                        link = reducedsoup[linkstart:linkend]

                    elif query.find(basicmicroid + 'a', reducedstart) != -1:
                        locator = query.find(basicmicroid + 'a', reducedstart)
                        reducedsoup = query[reducedstart:locator]
                        linkstart = reducedsoup.rfind('<a href="') + 9
                        linkend = reducedsoup.rfind('">')
                        link = reducedsoup[linkstart:linkend]

                    else:
                        flag = False
                        identifiers = identifiers + ['Not found']
                        print "RNAeye was unable to find the microRNA " + microids[i] + " in the search page."

            if flag:
                if initsearch == -1:
                    link = 'http://www.mirbase.org' + link
        
                f = urllib.urlopen(link)
                query = f.read()
                secstart = query.find('Mature sequence')

                if microids[i].find('3p') != -1 or microids[i].find('5p') != -1:
                    reducedstart = query.find(microids[i], secstart)

                    namestart = secstart + 16
                    nameend = query[namestart:namestart + 100].find('h2') + namestart - 3
                    truename = query[namestart:nameend]

                    reducedend = query.find('>Get sequence', reducedstart) + 2
                    reducedsoup = query[reducedstart:reducedend]

                    downloadstart = reducedsoup.rfind('<a href="') + 9
                    downloadend = reducedsoup.rfind('">')
                    download = reducedsoup[downloadstart:downloadend]

                    if reducedstart == -1 or namestart == 0 or nameend == -1 or reducedend == 1 or downloadstart == 8 or downloadend == -1:
                        flag = False
                        identifiers = identifiers + ['Not found']
                        print "RNAeye was unable to find the microRNA " + microids[i] + " on the main page."

                else:
                    reducedstart = query.find(microids[i] + '-5p', secstart)

                    namestart = secstart + 16
                    nameend = query[namestart:namestart + 100].find('h2') + namestart - 3
                    truename = query[namestart:nameend]

                    reducedend = query.find('>Get sequence', reducedstart) + 2
                    reducedsoup = query[reducedstart:reducedend]

                    downloadstart = reducedsoup.rfind('<a href="') + 9
                    downloadend = reducedsoup.rfind('">')
                    download = reducedsoup[downloadstart:downloadend]

                    previous = reducedsoup.rfind('Previous IDs')
                    previousstart = reducedsoup.find('">', previous) + 2
                    previousend = reducedsoup.find('<', previousstart)
                    previousname = reducedsoup[previousstart:previousend]
                    star = previousname.find('*')

                    if reducedstart == -1 or namestart == 0 or nameend == -1 or reducedend == 1 or downloadstart == 8 or downloadend == -1 or star != -1:
                        reducedstart = query.find(microids[i] + '-3p', secstart)

                        namestart = secstart + 16
                        nameend = query[namestart:namestart + 100].find('h2') + namestart - 3
                        truename = query[namestart:nameend]

                        reducedend = query.find('>Get sequence', reducedstart) + 2
                        reducedsoup = query[reducedstart:reducedend]

                        downloadstart = reducedsoup.rfind('<a href="') + 9
                        downloadend = reducedsoup.rfind('">')
                        download = reducedsoup[downloadstart:downloadend]

                        previous = reducedsoup.rfind('Previous IDs')
                        previousstart = reducedsoup.find('">', previous) + 2
                        previousend = reducedsoup.find('<', previousstart)
                        previousname = reducedsoup[previousstart:previousend]
                        star = previousname.find('*')

                        if reducedstart == -1 or namestart == 0 or nameend == -1 or reducedend == 1 or downloadstart == 8 or downloadend == -1 or star != -1:

                            starflag = True
                            while starflag and secstart != 19:
                                reducedstart = query.find(microids[i], secstart)
                                secstart = query.find('Mature sequence', secstart) + 20

                                namestart = secstart + 16
                                nameend = query[namestart:namestart + 100].find('h2') + namestart - 3
                                truename = query[namestart:nameend]

                                reducedend = query.find('>Get sequence', reducedstart) + 2
                                reducedsoup = query[reducedstart:reducedend]

                                downloadstart = reducedsoup.rfind('<a href="') + 9
                                downloadend = reducedsoup.rfind('">')
                                download = reducedsoup[downloadstart:downloadend]

                                previous = reducedsoup.rfind('Previous IDs')
                                previousstart = reducedsoup.find('">', previous) + 2
                                previousend = reducedsoup.find('<', previousstart)
                                previousname = reducedsoup[previousstart:previousend]
                                star = previousname.find('*')

                                if reducedstart != -1 and namestart != 0 and nameend != -1 and reducedend != 1 and downloadstart != 8 and downloadend != -1 and star == -1:
                                    starflag = False

                            if starflag:
                                flag = False
                                identifiers = identifiers + ['Not found']
                                print "RNAeye was unable to find the microRNA " + microids[i] + " on the main page."

                if flag:            
                    f = urllib.urlopen('http://www.mirbase.org' + download)
                    seqpage = f.read()

                    firstreturn = seqpage.find('\n')
                    seqstart = seqpage.find('\n', firstreturn + 1) + 1
                    seqend = seqpage.find('\n', seqstart + 1)
                    sequence = seqpage[seqstart:seqend]

                    if seqstart == -1 or seqend == -1:
                        identifiers = identifiers + ['Not found']
                        print "An error occurred when RNAeye tried to download the sequence for miRNA " + microids[i] + "."

                    else:
                        identifiers = identifiers + [truename]

                        seqfile = open(homefolder + 'microRNA_Sequences/' + microids[i] + '.fasta', 'w')
                        seqfile.write('>microRNA ' + microids[i] + '; found under the name ' + truename + ' on miRBase.\n')
                        seqfile.write(sequence)
                        seqfile.close()

        i = i + 1

    i = 1

    while os.path.isfile(homefolder + 'User Archives/Sequence Downloads/' + str(datetime.now().day) + '-'+ str(datetime.now().month) + '-' + str(datetime.now().day) + '_' + str(i) + '.txt'):
        i = i + 1

    idfile = open(homefolder + 'User Archives/Sequence Downloads/' + str(datetime.now().day) + '-' + str(datetime.now().month) + '-' + str(datetime.now().year) + '_' + str(i) + '.txt', 'w')
    idfile.write('miRNA name (as given)\tmiRNA name (as found).\n')
    i = 0

    while i < len(indices):
        idfile.write(microids[indices[i]] + '\t' + identifiers[i] + '\n')
        i = i + 1

    idfile.close()

    return(identifiers)

# This subroutine runs miRandatuRbo, scanning all miRNAs in the list against
# the associated promoters & transcripts.  The geneids and microids lists are
# scanned against one another in order.  The output files list only hybrids
# with an alignment score higher than thresholdalign, and with a binding
# energy less than thresholdenergy.  It also uses the 'strict' option of
# miRanda, which considers only matches for which the seed sequence of the
# miRNA matches exactly (or very closely) the complementary sequence of the
# target.  Default values: thresholdalign = 150; thresholdenergy = -7, as
# typically recommended for local implementations of miRanda.

def miRandatuRboint(mirandalocation,thresholdenergy,thresholdalign,microids,geneids,homefolder):    
    i = 0
    while i < len(microids):
        if os.path.isfile(homefolder + 'microRNA_Sequences/' + microids[i] + '.fasta') and os.path.isfile(homefolder + 'Gene_Sequences/' + geneids[i] + 'trans.fasta'):
            os.system(homefolder + mirandalocation + ' ' + homefolder + 'microRNA_Sequences/' + microids[i] + '.fasta ' + homefolder + 'Gene_Sequences/' + geneids[i] + 'trans.fasta -out ' + homefolder + 'miRanda_Matches/miRNA_' + microids[i] + '_gene_' + geneids[i] + 'trans.txt -en ' + str(thresholdenergy) + ' -sc ' + str(thresholdalign) + ' -strict')
            #print "Now scanning miRNA " + microids[i] + " against gene " + geneids[i] + "."
        else:
            print "Could not find microRNA " + microids[i] + " with gene " + geneids[i] + "."
        i = i + 1

# This subroutine reads the output from miRanda and extracts the energy of
# the hybrid having the lowest binding energy, and the alignments score of the
# hybrid having the highest score.  These data are then ported to a text file
# in the User Archives folder, named according to the date on which the data
# was saved (format: yyyy-mm-dd_n where n is a number; that is, this is the
# nth dataset of this type saved on this date).

def miRandascoResint(microids,geneids,homefolder):
    maxscores = []
    maxenergies = []

    i = 0
    while i < len(microids):        
        if os.path.isfile(homefolder + 'miRanda_Matches/miRNA_' + microids[i] + '_gene_' + geneids[i] + 'trans.txt'):
            f = open(homefolder + 'miRanda_Matches/miRNA_' + microids[i] + '_gene_' + geneids[i] + 'trans.txt')
            alignments = f.read()
            f.close()

            if 'No Hits Found above Threshold' in alignments:
                maxscores = maxscores + ['No hits']
                maxenergies = maxenergies + ['No hits']

            else:
                start = alignments.rfind('>>')
                end = alignments.rfind('Scan Complete')
                reducedfile = alignments[start:end]

                scorestart = 0
                j = 0
                while j < 4:
                    scorestart = reducedfile.find('\t', scorestart + 1)
                    j = j + 1

                scoreend = reducedfile.find('\t', scorestart + 1)
                maxscore = reducedfile[scorestart + 1:scoreend]

                enstart = scoreend
                enend = reducedfile.find('\t', enstart + 1)
                maxenergy = reducedfile[enstart:enend]

                maxscores = maxscores + [maxscore]
                maxenergies = maxenergies + [maxenergy]

        else:
            maxscores = maxscores + ['Not found']
            maxenergies = maxenergies + ['Not found']

        i = i + 1

    i = 1
    while os.path.isfile(homefolder + 'User Archives/Interference Scores/' + str(datetime.now().day) + '-' + str(datetime.now().month) + '-' + str(datetime.now().year) + '_' + str(i) + '.txt'):
        i = i + 1

    f = open(homefolder + 'User/External Data.txt')
    data = f.read()
    data = data.strip()

    if data == 'External_Data' or len(shape(data)) == 2:
        data = []
    else:
        data = genfromtxt(homefolder + 'User/External Data.txt', dtype = str, skip_header = 1, usecols = 0)
    
    outfile = open(homefolder + 'User Archives/Interference Scores/' + str(datetime.now().day) + '-' + str(datetime.now().month) + '-' + str(datetime.now().year) + '_' + str(i) + '.txt', 'w')
    outfile.write('microRNA Interference Scores and Energies (maximal; as calculated by miRanda)\n\n')

    if len(data) != len(microids):
        if len(data) != 0:
            print "Note: the length of your list of external data was different in length from your list of microRNAs, so your external data was not included in the compilation of results in the User Archives."

        string = '%8s\t%8s\t%8s\t%8s'%('miRNA','Gene','Alignment Score','Binding Energy')
        outfile.write(string + '\n')
    
        j = 0
        while j < len(microids):
            string = '%8s\t%8s\t%8s\t%8s'%(microids[j],geneids[j],maxscores[j],maxenergies[j])
            outfile.write(string + '\n')
            j = j + 1

    else:
        string = '%8s\t%8s\t%8s\t%8s\t%8s' %('miRNA','Gene','Alignment Score','Binding Energy','External Data')
        outfile.write(string + '\n')

        j = 0
        while j < len(microids):
            string = '%8s\t%8s\t%8s\t%8s\t%8s'%(microids[j],geneids[j],maxscores[j],maxenergies[j],data[j])
            outfile.write(string + '\n')
            j = j + 1

    outfile.close()

    return(maxscores + maxenergies)

# This subroutine receives the list of microRNA names (IDs) and
# produces a list of files in the Temporary_Files folder that contain
# the reversed (NOT reverse-complement!) sequences for all these miRNAs.
# This is necessary to use as input for miRanda for the activating RNA
# scan, since activating RNAs should possess a greater thermodynamic
# stability at their 5' ends than at their 3' ends (that is, their genomic
# targets ought to be more stable at the 3' end than at the 5' end), which
# is the reverse of the pattern found in successful RNA interference
# targets.  Since miRanda does not include an explicit feature allowing
# it to reverse its preferred thermodynamically stable end, the program
# must be fed pre-processed sequence input of this type (reversed, but
# not complemented).  The program does the same thing to the sequences
# of genes that these RNAs are to be scanned against (also supplied by
# the user).

def MakeReverseSet(microids,geneids,homefolder):
    i = 0
    while i < len(microids):
        if os.path.isfile(homefolder + 'microRNA_Sequences/' + microids[i] + '.fasta'):
            f = open(homefolder + 'microRNA_Sequences/' + microids[i] + '.fasta')
            query = f.read()

            seqstart = query.rfind('\n') + 1
            fwdsequence = query[seqstart:len(query)]
            revsequence = fwdsequence[::-1]

            output = open(homefolder + 'Temporary_Files/REV_' + microids[i] + '.fasta', 'w')
            output.write('>microRNA ' + microids[i] + '.\n')
            output.write(revsequence)
            output.close()

        i = i + 1

    i = 0
        
    while i < len(microids):
        if os.path.isfile(homefolder + 'Gene_Sequences/' + geneids[i] + 'pro.fasta'):
            f = open(homefolder + 'Gene_Sequences/' + geneids[i] + 'pro.fasta')
            query = f.read()

            seqstart = query.rfind('\n') + 1
            fwdsequence = query[seqstart:len(query)]
            revsequence = fwdsequence[::-1]

            output = open(homefolder + 'Temporary_Files/REV_' + geneids[i] + 'pro.fasta', 'w')
            output.write('>' + geneids[i] + ' promoter.\n')
            output.write(revsequence)
            output.close()

        i = i + 1

# This subroutine runs miRandatuRbo, scanning all miRNAs in the list against
# the associated promoters & transcripts.  The geneids and microids lists are
# scanned against one another in order.  The output files list only hybrids
# with an alignment score higher than thresholdalign, and with a binding
# energy less than thresholdenergy.  All values used are miRanda's own internal
# defaults.  This subroutine is intended to take as input the reversed sequences
# generated by the MakeReverseSet subroutine, and its output should be post-
# processed to re-order the sequences of the output file correctly.

def miRandatuRboact(mirandalocation,microids,geneids,homefolder):    
    i = 0
    while i < len(microids):
        if os.path.isfile(homefolder + 'Temporary_Files/REV_' + microids[i] + '.fasta') and os.path.isfile(homefolder + 'Temporary_Files/REV_' + geneids[i] + 'pro.fasta'):
            os.system(homefolder + mirandalocation + ' ' + homefolder + 'Temporary_Files/REV_' + microids[i] + '.fasta ' + homefolder + 'Temporary_Files/REV_' + geneids[i] + 'pro.fasta -out ' + homefolder + 'Temporary_Files/REV_miRNA_' + microids[i] + '_gene_' + geneids[i] + 'pro.txt -sc 120 -en -5')
            #print "Now scanning microRNA " + microids[i] + " against gene " + geneids[i] + "."
        else:
            print "Could not find microRNA " + microids[i] + " with gene " + geneids[i] + "."
        i = i + 1

# This subroutine unscrambles the output from the miRanda program, under the
# assumption that miRanda has been fed reversed miRNA and target sequences
# initially.  Briefly, the subroutine changes the starting and ending positions
# of both the query and the reference sequences directly in the output of the
# miRanda program, and then sends the modified output to a new file folder for
# viewing by the user.  It also reverses the order of the query and target,
# allowing the "real" alignment to be read off easily.  Furthermore, the
# routine also picks up the relevant data (sequences of the query and reference,
# position of the match along the target, number of "skips" in the query and
# reference sequences), and returns it in an array that combines all these
# data together; the information can be extracted afterwards by the main
# method.

def UnscrambleOutput(microids,geneids,homefolder):
    i = 0
    targetstarts = []
    targetends = []
    micros = []
    targets = []
    bindings = []
    energies = []
    alignments = []
    
    while i < len(microids):
        if os.path.isfile(homefolder + 'Temporary_Files/REV_miRNA_' + microids[i] + '_gene_' + geneids[i] + 'pro.txt'):
            f = open(homefolder + 'Temporary_Files/REV_miRNA_' + microids[i] + '_gene_' + geneids[i] + 'pro.txt')
            query = f.read()
            f.close()

            #print "Now unscrambling " + microids[i] + " and " + geneids[i] + "."

            f = open(homefolder + 'Temporary_Files/REV_' + microids[i] + '.fasta')
            micro = f.read()
            f.close()

            microseq = micro[micro.rfind('\n') + 2:len(micro)]
            microlen = len(microseq)

            f = open(homefolder + 'Temporary_Files/REV_' + geneids[i] + 'pro.fasta')
            gene = f.read()
            f.close()

            geneseq = gene[gene.rfind('\n') + 2:len(gene)]
            genelen = len(geneseq)
            
            output = open(homefolder + 'miRanda_Matches/miRNA_' + microids[i] + '_gene_' + geneids[i] + 'pro.txt', 'w')

            cutoff = query.rfind('=-=-=-=-=-=') + 10
            output.write(query[0:cutoff])
            
            docend = query.rfind('>>') + 1
            fulldata = query[cutoff:docend]

            targetstarts = targetstarts + [[]]
            targetends = targetends + [[]]
            micros = micros + [[]]
            targets = targets + [[]]
            bindings = bindings + [[]]
            energies = energies + [[]]
            alignments = alignments + [[]]
            
            startsec = 0
            marker = fulldata.find('>')
            endsec = fulldata.find('\n\n', marker)
            lastmarker = fulldata.rfind('>')
            k = 0
            
            while marker != lastmarker:
                
                diagram = fulldata[startsec:marker]
                
                queryfirststart = diagram.find('Q:') + 2
                queryfirstend = diagram.find(' ', queryfirststart)
                queryfirst = diagram[queryfirststart:queryfirstend]
                queryfirst = int(queryfirst)
                newqueryfirst = microlen - queryfirst + 1

                querylaststart = queryfirstend + 4
                querylastend = diagram.find(' ', querylaststart + 1)
                querylast = diagram[querylaststart:querylastend]
                querylast = int(querylast)
                newquerylast = microlen - querylast + 1

                referencefirststart = diagram.find('R:') + 2
                referencefirstend = diagram.find(' ', referencefirststart)
                referencefirst = diagram[referencefirststart:referencefirstend]
                referencefirst = int(referencefirst)
                newreferencefirst = genelen - referencefirst + 1

                referencelaststart = referencefirstend + 4
                referencelastend = diagram.find(' ', referencelaststart + 1)
                referencelast = diagram[referencelaststart:referencelastend]
                referencelast = int(referencelast)
                newreferencelast = genelen - referencelast + 1

                queryseqstart = diagram.find("3'") + 3
                queryseqend = diagram.find("5'", queryseqstart) - 1
                queryseq = diagram[queryseqstart:queryseqend]
                newqueryseq = queryseq[::-1]

                bindingstart = queryseqend + 21
                bindingend = diagram.find('\n', bindingstart) - 1
                binding = diagram[bindingstart:bindingend]
                newbinding = binding[::-1]

                refseqstart = diagram.find("5'", bindingend) + 3
                refseqend = diagram.find("3'", refseqstart) - 1
                refseq = diagram[refseqstart:refseqend]
                newrefseq = refseq[::-1]

                energystart = diagram.find('Energy:') + 9
                energyend = diagram.find(' ', energystart)
                energy = float(diagram[energystart:energyend])

                alignmentstart = diagram.find('Score:') + 7
                alignmentend = diagram.find(' ', alignmentstart)
                alignment = float(diagram[alignmentstart:alignmentend])
                
                scores = fulldata[marker:endsec]
                firstposition = scores.find('\t')
                j = 0
                while j < 3:
                    firstposition = scores.find('\t', firstposition + 1)
                    j = j + 1

                lastposition = firstposition
                j = 0
                while j < 2:
                    lastposition = scores.find('\t', lastposition + 1)
                    j = j + 1

                lastposition = lastposition + 2

                output.write(diagram[0:queryfirststart] + str(newquerylast) + ' to ' + str(newqueryfirst) + '  R:' + str(newreferencelast) + ' to ' + str(newreferencefirst))
                output.write(diagram[referencelastend:queryseqstart] + newqueryseq + diagram[queryseqend:bindingstart] + newbinding + diagram[bindingend:refseqstart] + newrefseq)
                output.write(diagram[refseqend:len(diagram) - 1] + '\n')
                output.write(scores[0:firstposition] + '\t' + str(newquerylast) + ' ' + str(newqueryfirst) + '\t' + str(newreferencelast) + ' ' + str(newreferencefirst) + '\t')
                output.write(scores[lastposition:len(scores)])

                targetstarts[len(targetstarts) - 1] = targetstarts[len(targetstarts) - 1] + [newreferencelast]
                targetends[len(targetends) - 1] = targetends[len(targetends) - 1] + [newreferencefirst + 1]
                micros[len(micros) - 1] = micros[len(micros) - 1] + [newqueryseq]
                bindings[len(bindings) - 1] = bindings[len(bindings) - 1] + [newbinding]
                targets[len(targets) - 1] = targets[len(targets) - 1] + [newrefseq]
                energies[len(energies) - 1] = energies[len(energies) - 1] + [energy]
                alignments[len(alignments) - 1] = alignments[len(alignments) - 1] + [alignment]

                startsec = endsec
                marker = fulldata.find('>', marker + 1)
                endsec = fulldata.find('\n\n', marker)
                k = k + 1

            output.close()

        else:
            targetstarts = targetstarts + [['Not found.']]
            targetends = targetends + [['Not found.']]
            micros = micros + [['Not found.']]
            targets = targets + [['Not found.']]
            bindings = bindings + [['Not found.']]
            energies = energies + [['Not found.']]
            alignments = alignments + [['Not found.']]
        i = i + 1

    alldata = targetstarts + targetends + micros + targets + bindings + energies + alignments

    return(alldata)
        
# Valid subroutine. Finds CpG islands in a submitted sequence, using roughly
# the same algorithm as used by D. Takai and P. A. Jones. Briefly, a CpG
# island is any stretch of DNA >=200 nucleotides in length, in which the GC
# content is >50%, the CpG content (number of CG (ordered) pairs) is greater
# than 0.6 times the expected value of this content (expected value is
# calculated as (number of C)*(number of G)/(length of sequence)) and that
# contains at least 7 CpG sites (to avoid mathematical artefacts). The
# output of this subroutine differs slightly from that of the Takai/Jones
# program, since this implementation includes no minimum distance between
# consecutive islands.

def CpGislanddetector(sequence):
    if len(sequence) > 200:

        CpGverdicts = []
        GC = (sequence.count('G', 0, 200) + sequence.count('C', 0, 200))/200.
        CpG = sequence.count('CG', 0, 200)

        CpGexp = (sequence.count('C', 0, 200)*sequence.count('G', 0, 200))/200.

        if GC > 0.5 and CpG > 0.6*CpGexp and CpG > 7:
            CpGverdicts = CpGverdicts + [True]
        else:
            CpGverdicts = CpGverdicts + [False]

        i = 1
        
        while i < len(sequence) - 199:
            GCchange = 0
            if sequence[i - 1] == 'G' or sequence[i - 1] == 'C':
                GCchange = GCchange - 1
            if sequence[i + 199] == 'G' or sequence[i + 199] == 'C':
                GCchange = GCchange + 1
            GC = ((GC*200.) + GCchange)/200.
            CpGexp = (sequence.count('C',i, i + 200)*sequence.count('G',i, i + 200))/200.

            CpGchange = 0
            if sequence[i - 1] == 'C' and sequence[i] == 'G':
                CpGchange = CpGchange - 1
            if sequence[i + 198] == 'C' and sequence[i + 199] == 'G':
                CpGchange = CpGchange + 1
            CpG = CpG + CpGchange

            if GC > 0.5 and CpG > 0.6*CpGexp and CpG > 7:
                CpGverdicts = CpGverdicts + [True]
            else:
                CpGverdicts = CpGverdicts + [False]

            i = i + 1

        rawislandboundaries = []
        i = 0
        j = 0
        while i < len(CpGverdicts):
            if CpGverdicts[i]:
                rawislandboundaries = rawislandboundaries + [[]]
                rawislandboundaries[j] = rawislandboundaries[j] + [i, i + 199]
                j = j + 1
            i = i + 1

        if rawislandboundaries == []:
            islandboundaries = []
        else:    
            islandboundaries = [rawislandboundaries[0]]
            i = 0
            j = 0
            while i < len(rawislandboundaries):
                if rawislandboundaries[i][0] < islandboundaries[j][1]:
# Note: Python has a strange feature whereby assigning a list to a variable in a list (as in the case of
# [rawislandboundaries[0]] being assigned to islandboundaries above), the two lists are forever linked;
# that is, any variable assignment to one list is transitively ported to the other list as well. This
# means that the first entry corresponding to each "real" CpG island in the list rawislandboundaries will
# always be the full length of the real CpG island, whereas all the other entries in this list will
# correspond, as expected, to mini-CpG island windows of 200 bp.
                    islandboundaries[j][1] = rawislandboundaries[i][1]

                else:
                    j = j + 1
                    islandboundaries = islandboundaries + [rawislandboundaries[i]]
                i = i + 1

    else:
        print ""
        print "Warning: the submitted sequence is less than 200 bp long. This makes positive CpG island identification impossible."

    return islandboundaries

# Valid subroutine. Checks whether a target, whose starting position is
# submitted, lies within the boundaries of the CpG islands that have also been
# submitted to the subroutine. If any part of the target lies within any of the
# boundary sets in the list islandboundaries, return True; else, return False.
# Note that the target sequence may overlap with the CpG island at ONE base
# (either the first or last, because of the strict inequalities in the two
# 'and' clauses), but no more (i.e., it may lie on the edge of the island).

def inisland(targetstart,targetend,islandboundaries):
    i = 0
    flag = False
    while i < len(islandboundaries) and flag == False:
        if (targetstart >= islandboundaries[i][0] and targetstart < islandboundaries[i][1]) or (targetend > islandboundaries[i][0] and targetend <= islandboundaries[i][1]):
            flag = True
        i = i + 1

    return flag
                
# Valid subroutine. Given a sequence, returns a list with the starting positions
# of all CpG sites in the sequence. Note that two adjacent CpG sites are
# considered as separate, and the position number returned is the position of
# the C residue in the submitted sequence.

def CpGsitedetector(sequence):
    i = 1
    CpGs = []
    while i < len(sequence):
        if sequence[i - 1:i + 1] == 'CG':
            CpGs = CpGs + [i - 1]
        i = i + 1

    return CpGs

# Valid subroutine. Given the starting position of a target and a list of
# the C positions of CpG sites, returns True if the target has any overlap with
# any of the 2-nucleotide CpG sites, and false otherwise.

def insite(targetstart,targetend,CpGs):
    i = 0
    flag = False
    while i < len(CpGs) and flag == False:
        if (CpGs[i] >= targetstart and CpGs[i] <= targetend) or (CpGs[i] + 1 >= targetstart and CpGs[i] + 1 <= targetend):
            flag = True
        i = i + 1

    return flag

# Given a sequence, this subroutine determines whether it contains more than four
# consecutive characters (read: nucleotides), regardless of case.  NOTE: this means
# in practice that when this method is applied on miRanda output, it will return
# False in the case of sequences with more than four dash '-' characters!! This
# subroutine is not intended to be used this way, but may be, if so desired.

def searchconsecutive(sequence):
    sequence = sequence.upper()
    i = 0
    flag = False

    while i < len(sequence) - 4 and flag == False:
        if sequence[i] == sequence[i + 1] and sequence[i] == sequence[i + 2] and sequence[i] == sequence[i + 3] and sequence[i] == sequence[i + 4]:
            flag = True

        i = i + 1

    return(flag)

# This subroutine preprocesses the data generated from the miRanda results on
# the activating sequences.  It performs the first round of eliminations, seeking
# out sequences whose targets are in CpG islands, include CpG sites, or include more
# than four consecutive nucleotides (or whose targets include more than four
# consecutive nucleotides).  It then returns an array of alignments scores; however,
# if a particular alignment is to be eliminated due to presence in a CpG island or
# CpG site or to the presence of more than four consecutive nucleotides, the
# subroutine replaces the alignment score with a negative number.  -1 corresponds
# to a CpG site, -2 to a CpG island, and -3 to the presence of more than four
# consecutive nucleotides on either the microRNA or the target sequence.

def ActivPreAnalyzer(targetstarts,targetends,micros,targets,alignments,energies,microids,geneids,homefolder):
    modalignments = alignments
    modenergies = energies
    i = 0

    while i < len(microids):
        #print "Analyzing sequence pair " + microids[i] + " and " + geneids[i] + "."
        
        if modalignments[i] != ['Not found.'] and modalignments[i] != []:
            f = open(homefolder + 'Gene_Sequences/' + geneids[i] + 'pro.fasta')
            query = f.read()

            seqstart = query.find('\n') + 1
            sequence = query[seqstart:len(query)]

            CpGs = CpGsitedetector(sequence)
            islandboundaries = CpGislanddetector(sequence)

            j = 0

            while j < len(targetstarts[i]):
                augmentedtarget = sequence[targetstarts[i][j] - 3:targetends[i][j] + 3]

                conmicroflag = searchconsecutive(micros[i][j])
                contargetflag = searchconsecutive(augmentedtarget)
                siteflag = insite(targetstarts[i][j],targetends[i][j],CpGs)
                islandflag = inisland(targetstarts[i][j],targetends[i][j],islandboundaries)
                if siteflag == True:
                    modalignments[i][j] = -1
                    modenergies[i][j] = -1
                elif islandflag == True:
                    modalignments[i][j] = -2
                    modenergies[i][j] = -2
                elif conmicroflag == True or contargetflag == True:
                    modalignments[i][j] = -3
                    modenergies[i][j] = -3

                j = j + 1

        i = i + 1

    modparameters = modalignments + modenergies
    
    return(modparameters)

# This subroutine reads the output from the preprocessor and extracts the energy of
# the hybrid having the lowest (negative) binding energy, and the alignments score
# of the hybrid having the highest score, subject to the constraints that these matches
# should not be in CpG sites, or CpG islands, and should not include more than four
# consecutive nucleotides (this information is passed to the routine in the
# modparameters list).  These data are then ported to a text file
# in the User Archives folder, named according to the date on which the data
# was saved (format: yyyy-mm-dd_n where n is a number; that is, this is the
# nth dataset of this type saved on this date).

def miRandascoResact(modalignments,modenergies,microids,geneids,homefolder):
    maxscores = []
    maxenergies = []
    
    i = 0
    
    while i < len(microids):
        if modalignments[i] == ['Not found.']:
            maxscores = maxscores + ['Not found']
            maxenergies = maxenergies + ['Not found']

        elif modalignments[i] == []:
            maxscores = maxscores + ['No hits']
            maxenergies = maxenergies + ['No hits']

        elif max(modalignments[i]) == -1:
            maxscores = maxscores + ['CpG site']
            maxenergies = maxenergies + ['CpG site']

        elif max(modalignments[i]) == -2:
            maxscores = maxscores + ['CpG island']
            maxenergies = maxenergies + ['CpG island']

        elif max(modalignments[i]) == -3:
            maxscores = maxscores + ['5+ consecutive']
            maxenergies = maxenergies + ['5+ consecutive']

        else:
            maxscores = maxscores + [max(modalignments[i])]
            maxenergies = maxenergies + [min(modenergies[i])]

        i = i + 1

    i = 1

    while os.path.isfile(homefolder + 'User Archives/Activation Scores/' + str(datetime.now().day) + '-' + str(datetime.now().month) + '-' + str(datetime.now().year) + '_' + str(i) + '.txt'):
        i = i + 1

    f = open(homefolder + 'User/External Data.txt')
    data = f.read()
    data = data.strip()

    if data == 'External_Data' or len(shape(data)) == 2:
        data = []
    else:
        data = genfromtxt(homefolder + 'User/External Data.txt', dtype = str, skip_header = 1, usecols = 0)

    outfile = open(homefolder + 'User Archives/Activation Scores/' + str(datetime.now().day) + '-' + str(datetime.now().month) + '-' + str(datetime.now().year) + '_' + str(i) + '.txt', 'w')
    outfile.write('microRNA Activation Scores and Energies (maximal; as calculated by a modification of miRanda)\n\n')

    if len(data) != len(microids):
        if len(data) != 0:
            print "Note: The length of your list of external data was different in length from your list of microRNAs, so your external data was not included in the compilation of results in the User Archives."

        string = '%8s\t%8s\t%8s\t%8s' %('miRNA','Gene','Alignment Score','Binding Energy')
        outfile.write(string + '\n')

        j = 0
        while j < len(microids):
            string = '%8s\t%8s\t%8s\t%8s'%(microids[j],geneids[j],maxscores[j],maxenergies[j])
            outfile.write(string + '\n')
            j = j + 1

    else:
        string = '%8s\t%8s\t%8s\t%8s\t%8s' %('miRNA','Gene','Alignment Score','Binding Energy','External Data')
        outfile.write(string + '\n')

        j = 0
        while j < len(microids):
            string = '%8s\t%8s\t%8s\t%8s\t%8s'%(microids[j],geneids[j],maxscores[j],maxenergies[j],data[j])
            outfile.write(string + '\n')
            j = j + 1

    outfile.close()

    return(maxscores + maxenergies)

# This subroutine adds microRNA sequences to the user's library, with the sequences submitted
# by the user, rather than found and downloaded online.  The sequences are saved in FASTA format
# and the names of the files correspond to the microRNA names given by the user.  The routine
# returns a warning message if the sequences submitted by the user is already in the user's
# library; in this case, the file is not overwritten.

def addmicros(microids,microseqs,homefolder):
    i = 0
    while i < len(microids):
        if os.path.isfile(homefolder + 'microRNA_Sequences/' + microids[i] + '.fasta'):
            print "Warning: a file by the name " + microids[i] + ".fasta already exists in your library."

        else:
            microseqs[i] = microseqs[i].upper()
            microseqs[i] = microseqs[i].replace('T','U')
            
            j = 0
            flag = True
            while j < len(microseqs[i]) and flag:
                if microseqs[i][j] != 'A' and microseqs[i][j] != 'G' and microseqs[i][j] != 'C' and microseqs[i][j] != 'U':
                    flag = False

                j = j + 1
                
            if flag == False:
                print "Warning: the given sequence for microRNA " + microids[i] + " is invalid." 

            else:
                revcomp = RCgenRNA(microseqs[i])
                outfile = open(homefolder + 'microRNA_Sequences/' + microids[i] + '.fasta', 'w')
                outfile.write('>microRNA ' + microids[i] + '; sequence was user-submitted.\n')
                outfile.write(revcomp)

        i = i + 1

# This subroutine adds promoter and transcript sequences submitted by the user to the user's
# library.  As above, the sequences are user-submitted and saved in FASTA format.  The subroutine
# also returns an error message if the user's library already contains a file by the name
# submitted; in this case, the file is not overwritten.

def addgenes(geneids,geneproseqs,genetransseqs,homefolder):
    i = 0
    while i < len(geneids):
        if os.path.isfile(homefolder + 'Gene_Sequences/' + geneids[i] + 'pro.fasta'):
            print "Warning: a file by the name " + geneids[i] + "pro.fasta already exists in your library."

        else:
            geneproseqs[i] = geneproseqs[i].upper()
            geneproseqs[i] = geneproseqs[i].replace('U','T')

            j = 0
            flag = True
            while j < len(geneproseqs[i]) and flag:
                if geneproseqs[i][j] != 'A' and geneproseqs[i][j] != 'G' and geneproseqs[i][j] != 'C' and geneproseqs[i][j] != 'T':
                    flag = False

                j = j + 1

            if flag == False:
                print "Warning: the given promoter sequence for gene " + geneids[i] + " is invalid."

            else:
                outfile = open(homefolder + 'Gene_Sequences/' + geneids[i] + 'pro.fasta', 'w')
                outfile.write('>' + geneids[i] + ' promoter; sequence was user-submitted.\n')
                outfile.write(geneproseqs[i])

        if os.path.isfile(homefolder + 'Gene_Sequences/' + geneids[i] + 'trans.fasta'):
            print "Warning: a file by the name " + geneids[i] + "trans.fasta already exists in your library."

        else:
            genetransseqs[i] = genetransseqs[i].upper()
            genetransseqs[i] = genetransseqs[i].replace('T','U')

            j = 0
            flag = True
            while j < len(genetransseqs[i]) and flag:
                if genetransseqs[i][j] != 'A' and genetransseqs[i][j] != 'G' and genetransseqs[i][j] != 'C' and genetransseqs[i][j] != 'U':
                    flag = False

                j = j + 1

            if flag == False:
                print "Warning: the given transcript sequence for gene " + geneids[i] + " is invalid."
                
            else:
                outfile = open(homefolder + 'Gene_Sequences/' + geneids[i] + 'trans.fasta', 'w')
                outfile.write('>' + geneids[i] + ' transcript; sequence was user-submitted.\n')
                outfile.write(genetransseqs[i])

        i = i + 1

# Given a gene's refseq ID, this subroutine determines the name of the gene
# and its acronym, and prints the information in the shell.  It also attempts to
# give reasons to the user in case of failure, and returns a flag of False
# if a search fails, and True if the search is successful.

def QueryGene(name):
    flag = True
    
    f = urllib.urlopen('http://www.mybioinfo.info/display.php?search_type=dna_refseq&search_key=' + name + '&B3=Search')
    query = f.read()

    end = query.find(name)
    reducedsoup = query[0:end]
    end = reducedsoup.rfind('</td>')
    reducedsoup = reducedsoup[0:end]
    end = reducedsoup.rfind('</td>')
    reducedsoup = reducedsoup[0:end]

    descriptstart = reducedsoup.rfind('">') + 2
    description = reducedsoup[descriptstart:end]

    end = reducedsoup.rfind('</td>')
    reducedsoup = reducedsoup[0:end]

    acrostart = reducedsoup.rfind('">') + 2
    acronym = reducedsoup[acrostart:end]

    if descriptstart == 1 or acrostart == 1 or end == -1 or 'This gene is on' not in query:
        flag = False

        if "No result found for" in query:
            print "The queried refseq ID could not be found in the mybioinfo database."

        else:
            print "The queried refseq ID could not be found for an unknown reason."

        print "The gene may be searched manually at http://www.mybioinfo.info/DNA_refseq_search.php."

    else:
        print "The queried gene refseq ID corresponds to:"
        print description + " (" + acronym + ")."

    return(flag)

# Given a microRNA's standard ID, this subroutine determines the true name of the
# microRNA (as specified on miRBase) and its full sequence, and prints this
# information in the shell.  In case of search failure, the subroutine attempts
# to identify the reason for the failure, and conveys this information to the user
# through the shell.  It returns a flag of True if the search was successful,
# and False if not.

def QueryMicro(name):
    link = 'http://www.mirbase.org/cgi-bin/query.pl?terms=' + name
    
    f = urllib.urlopen(link)
    queryi = f.read()

    initsearch = queryi.find('Mature sequence')
    flag = True

    if initsearch == -1:
        query = queryi.lower()
        reducedstart = query.find('restore to the original order')
        basic = name.lower()

        if '3p' in name or '5p' in name:
            basic = basic[0:len(basic) - 3]
            locator = query.find('>' + basic + '<', reducedstart) + 1
            if locator == 0:
                locator = query.find('>' + basic + ' <', reducedstart) + 1
            if locator == 0:
                locator = query.find('>' + basic + '*<', reducedstart) + 1
            if locator == 0:
                locator = query.find(basic + '-', reducedstart)
            if locator == -1:
                locator = query.find(basic + 'a', reducedstart)
                    
            if locator != -1:
                reducedsoup = query[reducedstart:locator]
                linkstart = reducedsoup.rfind('<a href="') + 9
                linkend = reducedsoup.rfind('">')
                link = reducedsoup[linkstart:linkend]

        else:
            locator = query.find('>' + basic + '<', reducedstart) + 1
            if locator == 0:
                locator = query.find('>' + basic + ' <', reducedstart) + 1
            if locator == 0:
                locator = query.find('>' + basic + '*<', reducedstart) + 1
            if locator == 0:
                locator = query.find(basic + '-', reducedstart)
            if locator == -1:
                locator = query.find(basic + 'a', reducedstart)

            if locator != -1:
                reducedsoup = query[reducedstart:locator]
                linkstart = reducedsoup.rfind('<a href="') + 9
                linkend = reducedsoup.rfind('">')
                link = reducedsoup[linkstart:linkend]
                
        if locator == -1:
            flag = False
            print "RNAeye was unable to find this microRNA."

            if 'dead mirna entry' in query:
                deadfind = query.find('dead mirna entry')
                marker = query.find('comment:',deadfind)
                commentstart = query.find('\n\t', marker) + 2
                commentend = query.find('\n</dd>',commentstart)
                comment = queryi[commentstart:commentend]

                print "Reason: " + comment

            elif 'here are a few tips on how to get the most' in query:
                print "Reason: The microRNA " + name + " is not registered on miRBase."

            else:
                print "The reason for the search failure is unknown."

    if flag:
        if initsearch == -1:
            link = 'http://www.mirbase.org' + link

        f = urllib.urlopen(link)
        query = f.read()
        secstart = query.find('Mature sequence')

        if '3p' in name or '5p' in name:
            reducedstart = query.find(name, secstart)
            namestart = secstart + 16
            nameend = query[namestart:namestart + 100].find('h2') + namestart - 3
            truename = query[namestart:nameend]

            reducedend = query.find('>Get sequence', reducedstart) + 2
            reducedsoup = query[reducedstart:reducedend]

            downloadstart = reducedsoup.rfind('<a href="') + 9
            downloadend = reducedsoup.rfind('">')
            download = reducedsoup[downloadstart:downloadend]

            if reducedstart == -1 or namestart == 0 or nameend == -1 or reducedend == 1 or downloadstart == 8 or downloadend == -1:
                flag = False
                print "RNAeye was unable to find the specified microRNA on miRBase; however, you may be able to download its sequence by visiting " + link + "."

        else:
            reducedstart = query.find(name + '-5p', secstart)
            namestart = secstart + 16
            nameend = query[namestart:namestart + 100].find('h2') + namestart - 3
            truename = query[namestart:nameend]

            reducedend = query.find('>Get sequence', reducedstart) + 2
            reducedsoup = query[reducedstart:reducedend]

            downloadstart = reducedsoup.rfind('<a href="') + 9
            downloadend = reducedsoup.rfind('">')
            download = reducedsoup[downloadstart:downloadend]

            previous = reducedsoup.rfind('Previous IDs')
            previousstart = reducedsoup.find('">', previous) + 2
            previousend = reducedsoup.find('<', previousstart)
            previousnames = reducedsoup[previousstart:previousend]
            star = previousnames.find('*')

            if reducedstart == -1 or namestart == 0 or nameend == -1 or reducedend == 1 or downloadstart == 8 or downloadend == -1 or star != -1:
                reducedstart = query.find(name + '-3p', secstart)
                namestart = secstart + 16
                nameend = query[namestart:namestart + 100].find('h2') + namestart - 3
                truename = query[namestart:nameend]

                reducedend = query.find('>Get sequence', reducedstart) + 2
                reducedsoup = query[reducedstart:reducedend]

                downloadstart = reducedsoup.rfind('<a href="') + 9
                downloadend = reducedsoup.rfind('">')
                download = reducedsoup[downloadstart:downloadend]

                previous = reducedsoup.rfind('<a href="') + 9
                previousstart = reducedsoup.find('">', previous) + 2
                previousend = reducedsoup.find('<', previousstart)
                previousnames = reducedsoup[previousstart:previousend]
                star = previousnames.find('*')

                if reducedstart == -1 or namestart == 0 or nameend == -1 or reducedend == 1 or downloadstart == 8 or downloadend == -1 or star != -1:
                    starflag = True
                    while starflag and secstart != 19:
                        reducedstart = query.find(name, secstart)
                        secstart = query.find('Mature sequence', secstart)

                        namestart = secstart + 16
                        nameend = query[namestart:namestart + 100].find('h2') + namestart - 3
                        truename = query[namestart:nameend]

                        reducedend = query.find('>Get sequence', reducedstart) + 2
                        reducedsoup = query[reducedstart:reducedend]

                        downloadstart = reducedsoup.rfind('<a href="') + 9
                        downloadend = reducedsoup.rfind('">')
                        download = reducedsoup[downloadstart:downloadend]

                        previous = reducedsoup.rfind('Previous IDs')
                        previousstart = reducedsoup.find('">', previous) + 2
                        previousend = reducedsoup.find('<', previousstart)
                        previousnames = reducedsoup[previousstart:previousend]
                        star = previousnames.find('*')

                        if reducedstart != -1 and namestart != 0 and nameend != -1 and reducedend != -1 and downloadstart != 8 and downloadend != -1 and star == -1:
                            starflag = False

                    if starflag:
                        flag = False
                        print "RNAeye was unable to find the specified microRNA on miRBase; however, you may be able to download its sequence by visiting " + link + "."

    if flag:
        f = urllib.urlopen('http://www.mirbase.org' + download)
        seqpage = f.read()

        firstreturn = seqpage.find('\n')
        seqstart = seqpage.find('\n', firstreturn + 1) + 1
        seqend = seqpage.find('\n', seqstart + 1)
        sequence = seqpage[seqstart:seqend]

        if seqstart == -1 or seqend == -1:
            print "An unexpected error occurred when RNAeye tried to download the sequence of this microRNA.  You may run the method again, or go to http://www.mirbase.org" + download + " to download the sequence manually."

        else:
            print "microRNA " + name + " was successfully located on miRBase under the name " + truename + "."
            print "The sequence for this microRNA is"
            print sequence

    return(flag)

# This program automatically scans a list of DNA sequences in .txt format and
# BLASTs each sequence successively in BLASTn, with a 3-second delay between
# sequences to comply with NCBI bandwidth usage guidelines.  This program is
# specifically designed to search only short sequences, < 30 nucleotides.
# Longer sequences will require a modification of the search parameters.
# NCBI servers are pinged for results once every 60 seconds, also to comply
# with NCBI guidelines.  Thus, high-throughput BLAST searches take approximately
# 63 seconds per sequence to complete (100 sequences will take 1 h 45 min).

def BLASTER(seqlist,shortflag,goodflag,homefolder):
    RID = []
    flags = []

# The following loop issues BLASTn search orders to the NCBI mainframe, with
# a delay of 3 seconds between them, and stores the consequent Request ID
# numbers issued to it by the mainframe.  This is the "Put" part of the program.
    i = 0
    if shortflag:
        print ""
        if goodflag:
            while i < len(seqlist):
                print "RNAeye is submitting the sequence for microRNA " + seqlist[i] + " to the NCBI server."

                if os.path.isfile(homefolder + 'microRNA_Sequences/' + seqlist[i] + '.fasta'):
                    f = open(homefolder + 'microRNA_Sequences/' + seqlist[i] + '.fasta')
                    rawdoc = f.read()

                    seqstart = rawdoc.find('\n') + 1
                    seqend = len(rawdoc)
                    sequence = rawdoc[seqstart:seqend]
                    sequence = RCgenRNA(sequence)

# Explanation of the URL below:
# CMD = Put means that a command is issued to the NCBI server to submit a
# sequence to BLAST.
# QUERY = sequence is the actual sequence that is to be submitted to the
# server.
# PROGRAM = blastn means that this sequence will be searched by the nucleotide
# BLAST program.
# DATABASE = refseq_rna means that the database to be BLASTed against is that
# of all reference RNA sequences; in other words, we are looking for targeted
# transcripts.
# ENTREZ_QUERY = human%5Borgn%5D is URL-friendly script for "human[orgn]", and
# ensures that BLAST narrows its search of the refseq_rna database it is
# searching to only human entries.
# BLAST_PROGRAM = MegaBlast means that blastn will search using MegaBlast, which
# is optimized to find highly similar sequences.
# SHOW_OVERVIEW = false means that BLAST will not return a graphical overview of
# the hits.  This feature should be turned-off during all script-based searches.
                    try:
                        f = urllib.urlopen('http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&QUERY=' + sequence + '&PROGRAM=blastn&DATABASE=refseq_rna&ENTREZ_QUERY=human%5Borgn%5D&BLAST_PROGRAM=MegaBlast&SHOW_OVERVIEW=false')
                    except IOError:
                        print "RNAeye was unable to detect an Internet connection; the NCBI server could not be accessed."
                        
                    putpage = f.read()
                    RIDposition = putpage.find('RID%3D') + 6

                    if RIDposition == 5:
                        print "An error occurred when RNAeye tried to access the Request ID for microRNA sequence " + seqlist[i] + "."
                        flags = flags + [False]
                        RID = RID + ['']

                    else:
                        flags = flags + [True]
                        j = 0
                        RID = RID + [putpage[RIDposition:RIDposition + 11]]

                else:
                    flags = flags + [False]
                    RID = RID + ['']

                i = i + 1
                time.sleep(3)

        else:
            print "RNAeye could not submit your microRNA sequences to NCBI."

    else:
        print ""
        if goodflag:
            while i < len(seqlist):
                print "RNAeye is submitting the transcript sequence for gene " + seqlist[i] + " to the NCBI server."

                if os.path.isfile(homefolder + 'Gene_Sequences/' + seqlist[i] + 'trans.fasta'):
                    f = open(homefolder + 'Gene_Sequences/' + seqlist[i] + 'trans.fasta')
                    rawdoc = f.read()

                    seqstart = rawdoc.find('\n') + 1
                    seqend = len(rawdoc)
                    sequence = rawdoc[seqstart:seqend]

# Explanation of the URL below:
# CMD = Put means that a command is issued to the NCBI server to submit a
# sequence to BLAST.
# QUERY = sequence is the actual sequence that is to be submitted to the
# server.
# PROGRAM = blastn means that this sequence will be searched by the nucleotide
# BLAST program.
# DATABASE = refseq_rna means that the database to be BLASTed against is that
# of all reference RNA sequences; in other words, we are looking for targeted
# transcripts.
# ENTREZ_QUERY = human%5Borgn%5D is URL-friendly script for "human[orgn]", and
# ensures that BLAST narrows its search of the refseq_rna database it is
# searching to only human entries.
# SHOW_OVERVIEW = false means that BLAST will not return a graphical overview of
# the hits.  This feature should be turned-off during all script-based searches.

                    f = urllib.urlopen('http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&QUERY=' + sequence + '&PROGRAM=blastn&DATABASE=refseq_rna&ENTREZ_QUERY=human%5Borgn%5D&SHOW_OVERVIEW=false')
                    putpage = f.read()
                    RIDposition = putpage.find('RID%3D') + 6

                    if RIDposition == 5:
                        print "An error occurred when RNAeye tried to access the Request ID for gene sequence " + seqlist[i] + "."
                        flags = flags + [False]
                        RID = RID + ['']

                    else:
                        flags = flags + [True]
                        RID = RID + [putpage[RIDposition:RIDposition + 11]]

                else:
                    flags = flags + [False]
                    RID = RID + ['']

                i = i + 1
                time.sleep(3)

        else:
            print "RNAeye could not submit your gene transcript sequences to NCBI."

    if goodflag:
        print ""
        i = 0
        print "Accessing BLAST results now..."

# The following loop retrieves the results of the searches issued in the
# previous loop, using the RIDs collected therefrom.  This is the "Get" part of
# the program.

        while i < len(seqlist):
            print "Retrieving BLAST results for sequence " + seqlist[i] + "."

            if flags[i] == False:
                print "The results for the sequence " + seqlist[i] + " could not be retrieved."

            else:
#Explanation of the URL below:
# CMD = Get means that a command is issued to the NCBI server to retrieve the
# results of a BLAST search.
# FORMAT_OBJECT = Alignment sets the object to be formatted by the FORMAT_TYPE
# variable to be the alignment of sequences found.
# FORMAT_TYPE = Text specifies the format to be .txt, so that Python will have
# an easier time reading it than if it were in HTML.
# NOHEADER = true removes the header from the results page, simplifying it
# considerably and easing reading by the Python script.
# RID = RID[i] sets the Request ID; that is, it tells BLAST to retrieve the
# results of the query that has this particular RID.
# SHOW_OVERVIEW = false means that the graphical overview is turned off, as it
# always should be for script-assisted searches.
                time.sleep(60)
                f = urllib.urlopen("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=Alignment&FORMAT_TYPE=Text&NOHEADER=true&RID=" + RID[i] + "&SHOW_OVERVIEW=false")
                getpage = f.read()

                if shortflag == True:
                    outfile = open(homefolder + 'BLAST_Results/microRNA_' + seqlist[i] + '.txt', 'w')

                else:
                    outfile = open(homefolder + 'BLAST_Results/Gene_' + seqlist[i] + '.txt', 'w')

                outfile.write(getpage)
                outfile.close()
                    
            i = i + 1

# When given no arguments, this subroutine adds to your library all microRNAs currently listed
# in your document "miRNA IDs.txt.  This method can detect whether the sequences are
# user-submitted (in which case it passes the sequences directly to its library) or absent (in which
# case it searches for the sequences online based on the identifiers provided by the user).

def addlibmicro():
    outflag = True
    f = open(homefolder + 'User/miRNA IDs.txt')
    microdoc = f.read()
    microdoc = microdoc.strip()
    
    if microdoc == 'miRNA_ID\tmiRNA_Sequence_(Optional)':
        microinput = []
    else:
        microinput = loadtxt(homefolder + 'User/miRNA IDs.txt', dtype = str, skiprows = 1)

    if len(shape(microinput)) == 1:
        microids = microinput

        if shape(microinput)[0] != 0:            
            print "Downloading microRNA sequences.."
            DownloadMicros(microids,homefolder)
            print ""
            print "The submitted sequences have been added to your library."

        else:
            print "Your list of microRNA sequences is blank."
            outflag = False

    elif len(shape(microinput)) == 2:
        if shape(microinput)[1] == 2:
            microids = microinput[:,0]
            microseqs = microinput[:,1]

            flag = lengthcheck(microids,microseqs)
            if flag == False:
                print "Error: your list of microRNA names is different in length from your list of microRNA sequences."
                outflag = False
            else:
                print "Adding microRNA sequences to library..."
                addmicros(microids,microseqs,homefolder)

    else:
        print "Error: your microRNA IDs list contains an incorrect number of columns."
        outflag = False

    return(outflag)

# When given no arguments, this subroutine adds to your library all genes currently listed
# in the file 'Gene IDs.txt.'  This method can detect whether the sequences are
# user-submitted (in which case it passes the sequences directly to its library) or absent (in which
# case it searches for the sequences online based on the identifiers provided by the user).

def addlibgene():
    outflag = True
    f = open(homefolder + 'User/Gene IDs.txt')
    genedoc = f.read()
    genedoc = genedoc.strip()
    
    if genedoc == 'GeneID\t\tPromoter_Sequence_(Optional)\t\tTranscript_Sequence_(Optional)':
        geneinput = []
    else:
        geneinput = loadtxt(homefolder + 'User/Gene IDs.txt', dtype = str, skiprows = 1)

    if len(shape(geneinput)) == 1:
        geneids = geneinput

        if shape(geneinput)[0] != 0:
            print "Downloading promoter sequences..."
            DownloadPromoters(geneids,homefolder)
            print ""
            print "Downloading transcript sequences..."
            DownloadTranscripts(geneids,homefolder)
            print ""

        else:
            print "Your list of gene IDs is blank."
            outflag = False

    elif len(shape(geneinput)) == 2:
        if shape(geneinput)[1] == 3:
            geneids = geneinput[:,0]
            geneproseqs = geneinput[:,1]
            genetransseqs = geneinput[:,2]

            flag = lengthcheck(geneids,geneproseqs)
            if flag == False:
                print "Error: your list of gene names is different in length from your list of gene promoter sequences."
                outflag = False
            else:
                flag = lengthcheck(geneids,genetransseqs)
                if flag == False:
                    print "Error: your list of gene names is different in length from your list of gene transcript sequences."
                    outflag = False
                else:
                    print "Adding gene sequences to library..."
                    addgenes(geneids,geneproseqs,genetransseqs,homefolder)
                    print "Adding microRNA sequences to library..."
                    addmicros(microids,microseqs,homefolder)

    else:
        print "Error: your Gene IDs list contains an incorrect number of columns."
        outflag = False

    return(outflag)

# This subroutine modifies the scores list so that it can be correlated numerically.
# Specifically, it eliminates all entries that read "Not found" and sets all entries
# whose scores are nullified by lack of hits, CpG presence, or 5 or more consecutive
# nucleotides to zero.

def reducescores(scores,indices):
    newscores = []
    i = 0
    j = 0
    while i < len(indices):
        if i == indices[j]:
            if scores[i] == 'No hits' or scores[i] == 'CpG site' or scores[i] == 'CpG island' or scores[i] == '5+ consecutive':
                newscores = newscores + [0]
            else:
                newscores = newscores + [float(scores[i])]

            j = j + 1

        i = i + 1

    return(newscores)

# This subroutine returns all the indices of valid score lists; that is, an entry
# of 'Not found' will not have its index included in the output, but any other
# entry in the list will have its index included.

def indexscores(scores):
    indices = []
    i = 0
    while i < len(scores):
        if scores[i] != 'Not found':
            indices = indices + [i]

        i = i + 1

    return(indices)

# This subroutine generates the reverse complement of an input sequence of RNA.
# It includes no checks on the validity of the sequence (i.e. absence of non-
# RNA characters), so it should only be called by subroutines that already have
# such checks in place.

def RCgenRNA(sequence):
    revcomp = ''
    i = len(sequence) - 1
    while i > -1:
        if sequence[i] == 'A':
            revcomp = revcomp + 'U'
        elif sequence[i] == 'G':
            revcomp = revcomp + 'C'
        elif sequence[i] == 'C':
            revcomp = revcomp + 'G'
        else:
            revcomp = revcomp + 'A'

        i = i - 1

    return(revcomp)

# This subroutine generates the reverse complement of an input sequence of DNA.
# It includes no checks on the validity of the sequence (i.e. absence of non-
# DNA characters), so it should only be called by subroutines that already have
# such checks in place.

def RCgenDNA(sequence):
    revcomp = ''
    i = len(sequence) - 1
    while i > -1:
        if sequence[i] == 'A':
            revcomp = revcomp + 'T'
        elif sequence[i] == 'G':
            revcomp = revcomp + 'C'
        elif sequence[i] == 'C':
            revcomp = revcomp + 'G'
        else:
            revcomp = revcomp + 'A'

        i = i - 1

    return(revcomp)

def calculate(newactivalign,newactivenergies,newinteralign,newinterenergies):
    finalscores = []
    i = 0
    while i < len(newactivalign):
        
# Scoring formula; to be optimized:

        score = -newactivalign[i]*newactivenergies[i] + newinteralign[i]*newinterenergies[i]
        
        finalscores = finalscores + [score]

        i = i + 1

    return(finalscores)

def calladdtolib(sequence_name):
    flag = True
    if 'XM_' in sequence_name or 'NM_' in sequence_name:
        if os.path.isfile(homefolder + 'Gene_Sequences/' + sequence_name + 'pro.fasta') and os.path.isfile(homefolder + 'Gene_Sequences/' + sequence_name + 'trans.fasta'):
            print "The gene " + sequence_name + " is already in your library."
            
        else:
            proflag = DownloadPromoters([sequence_name],homefolder)       
            transflag = DownloadTranscripts([sequence_name],homefolder)
            print ""
        
            if proflag == True and transflag == True:
                print "The gene " + sequence_name + " has been added to your library."

            else:
                print "The download of gene " + sequence_name + " could not be completed successfully."
                flag = False

    else:
        if os.path.isfile(homefolder + 'microRNA_Sequences/' + sequence_name + '.fasta'):
            print "The microRNA " + sequence_name + " is already in your library."

        else:
            indicator = DownloadMicros([sequence_name],homefolder)
            print ""

            if indicator != ['Not found']:
                print "The microRNA " + sequence_name + " has been added to your library; found as " + indicator[0] + " on miRBase."

            else:
                print "The download of microRNA " + sequence_name + " could not be completed successfully."
                flag = False

    return(flag)

# This subroutine creates a hairpin, given the desired passenger strand (or targeted
# sequence) of that hairpin.  The generated haripin includes the loop region and the
# restriction sites needed for cloning into AgeI/EcoRI restriction sites.  In essence,
# the sequences are ready-to-order.  Note that the passenger strand is identical to
# the genomic target sequence; and the guide strand is the reverse complement to this
# genomic sequence.

def callmaketargethp(targetseq):
    guideseq = RCgenDNA(targetseq)
    sense = 'CCGGT' + targetseq + 'CTCGAG' + guideseq + 'TTTTTG'
    asense = 'AATTCAAAAA' + targetseq + 'CTCGAG' + guideseq + 'A'
    hairpin = [sense, asense]
    return hairpin

# This subroutine takes a sequence and converts it to all caps, and changes
# all Ts to Us.  If the sequence includes characters other than a, g, c, t or u,
# the subroutine returns the string, 'failure'.

def checkRNAsequence(sequence):
    sequence = sequence.upper()
    
    i = 0
    flag = True
    while i < len(sequence) and flag:
        if sequence[i] != 'A' and sequence[i] != 'G' and sequence[i] != 'C' and sequence[i] != 'T' and sequence[i] != 'U':
            flag = False
        i = i + 1

    if flag:
        sequence = sequence.replace('T','U')
    else:
        sequence = 'failure'
        
    return sequence

# This subroutine takes a sequence and converts it to all caps, and changes
# all Us to Ts.  If the sequence includes characters other than a, g, c, t or u,
# the subroutine returns the string, 'failure'.

def checkDNAsequence(sequence):
    sequence = sequence.upper()
    
    i = 0
    flag = True
    while i < len(sequence) and flag:
        if sequence[i] != 'A' and sequence[i] != 'G' and sequence[i] != 'C' and sequence[i] != 'T' and sequence[i] != 'U':
            flag = False
        i = i + 1

    if flag:
        sequence = sequence.replace('U','T')
    else:
        sequence = 'failure'
        
    return sequence

# This subroutine takes a sequence in the correct format and generates a list
# of all n-base-long subsequences within that target. If the sequence length is
# m, then it generates m - n target subsequences, each with n nucleotides.
# Here, n = targlength.

def maketargets(sequence,targlength):
    i = 0
    targets = []
    while i < len(sequence) - targlength - 1:
        targets = targets + [sequence[i:i + targlength]]
        i = i + 1
    return targets

# This subroutine tests whether a particular target primer has EITHER a G
# or a C at its last two 3' positions, AND either a T, a U or an A at its third-
# to-last 3' position.  If all of these are true, it returns True, else it
# returns False.

def threeprimetest(primer):
    flag = False
    if (primer[-1] == 'G' or primer[-1] == 'C') and (primer[-2] == 'G' or primer[-2] == 'C') and (primer[-3] == 'A' or primer[-3] == 'T' or primer[-3] == 'U'):
        flag = True
    return flag

# This subroutine tests whether a particular primer sequence has 4 or more
# consecutive nucleotides in it.  If yes, return True; if no, return False.

def checkconsecutive(sequence):
    i = 3
    flag = False
    while i < len(sequence) and flag == False:
        if sequence[i - 3] == sequence[i - 2] and sequence[i - 2] == sequence[i - 1] and sequence[i - 1] == sequence[i]:
            flag = True

        i = i + 1

    return flag

# This subroutine takes a target sequence (of any length) and calculates the
# proportion of nucleotides that are either G or C, then returns that
# number.

def countGC(target):
    i = 0
    GC = 0
    while i < len(target):
        if target[i] == 'G' or target[i] == 'C':
            GC = GC + 1
        i = i + 1

    GC = float(GC)
    GC = GC/len(target)
    return GC

# Calculates melting temperature of an oligo; valid for oligos < 14 bases.
# Assumes the reaction is carried out in the presence of 50 mM monovalent
# cations.

def smalloligoTm(sequence):
    GC = countGC(sequence)*len(sequence)
    AT = (1. - countGC(sequence))*(len(sequence))

    Tm = 4.*GC + 2.*AT

    return Tm

# Calculates melting temperature of an oligo; valid for oligos >= 14 bases.

def largeoligoTm(sequence):
    GC = countGC(sequence)*len(sequence)

    Tm = 64.9 + 41.*(GC - 16.4)/len(sequence)

    return Tm

# Valid subroutine. Given a target sequence submitted 5' to 3', calculates the
# difference in the thermodynamic stabilities of the 3' and 5' ends using the
# nearest-neighbor method. Variable parameter: end (see top), which defines
# the number of nucleotides that make up the 3' end and the 5' end. The
# difference is 3' end energy - 5' end energy. All energies are negative values
# because the system is bound. It is desirable for the 5' end to be more stable
# than the 3' end. If the 5' end is more stable (energy more negative) than the
# 3' end, then the difference as defined by the program will be positive; this
# is the desirable condition.

def endstability(target):
    i = 0
    E5prime = 0.
    while i < end - 1:
        if target[i:i + 2] == 'AA' or target[i:i + 2] == 'UU':
            energy = AAE
        elif target[i:i + 2] == 'AU':
            energy = ATE
        elif target[i:i + 2] == 'UA':
            energy = TAE
        elif target[i:i + 2] == 'CA' or target[i:i + 2] == 'AC':
            energy = CAE
        elif target[i:i + 2] == 'GU' or target[i:i + 2] == 'UG':
            energy = GTE
        elif target[i:i + 2] == 'CU' or target[i:i + 2] == 'UC':
            energy = CTE
        elif target[i:i + 2] == 'GA' or target[i:i + 2] == 'AG':
            energy = GAE
        elif target[i:i + 2] == 'CG':
            energy = CGE
        elif target[i:i + 2] == 'GC':
            energy = GCE
        elif target[i:i + 2] == 'GG' or target[i:i + 2] == 'CC':
            energy = GGE
        else:
            print ""
            print "Error: none of the recognized nearest-neighbor pairings were detected."
            
        E5prime = E5prime + energy
        i = i + 1

    i = len(target) - end
    E3prime = 0.
    while i < len(target) - 1:
        if target[i:i + 2] == 'AA' or target[i:i + 2] == 'UU':
            energy = AAE
        elif target[i:i + 2] == 'AU':
            energy = ATE
        elif target[i:i + 2] == 'UA':
            energy = TAE
        elif target[i:i + 2] == 'CA' or target[i:i + 2] == 'AC':
            energy = CAE
        elif target[i:i + 2] == 'GU' or target[i:i + 2] == 'UG':
            energy = GTE
        elif target[i:i + 2] == 'CU' or target[i:i + 2] == 'UC':
            energy = CTE
        elif target[i:i + 2] == 'GA' or target[i:i + 2] == 'AG':
            energy = GAE
        elif target[i:i + 2] == 'CG':
            energy = CGE
        elif target[i:i + 2] == 'GC':
            energy = GCE
        elif target[i:i + 2] == 'GG' or target[i:i + 2] == 'CC':
            energy = GGE
        else:
            print ""
            print "Error: none of the recognized nearest-neighbor pairings were detected."

        E3prime = E3prime + energy 
        i = i + 1

# Note: we want the 3' end to have a lower TD stability than the 5' end. Since
# all the free energies added here are negative, that means that a POSITIVE
# value of stabilitydiff is what we want; the higher the positive value, the
# more stable the 5' end is than the 3' end.
    stabilitydiff = E3prime - E5prime
    return stabilitydiff

# Valid subroutine. Checks whether the given target has an 'A' at its 19th
# position. If it does, return True, else return False.

def checkpos19(target):
    flag = False
    if target[18] == 'A':
        flag = True
    return flag

# Valid subroutine. Checks whether the given target has an 'A' or a 'U' at its
# 18th position. If it has 'A', return 2; if 'U', return 1; if neither, return
# 0.

def checkpos18(target):
    indicator = 0
    if target[17] == 'A':
        indicator = 2
    elif target[17] == 'U':
        indicator = 1
    return indicator

# Valid subroutine. Checks whether the given target has a 'U' at its 7th
# position. If it does, return 1; else, return 0.

def checkpos7(target):
    indicator = 0
    if target[6] == 'U':
        indicator = 1
    return indicator

# Valid subrountine.  Given a target sequence submitted 5' to 3', determines
# whether the 3' or 5' end base is more energetically stable.  If the 5' end
# is more stable, returns True (since this means that the 5' end of the guide
# strand is less stable); otherwise, returns False.

def endbase(target):
    if target[0] == 'A' or target[0] == 'U':
        energy = AAE/2
    elif target[0] == 'G' or target[0] == 'C':
        energy = GGE/2
    else:
        print ""
        print "Error: none of the standard bases were detected."

    E5prime = energy
    #print E5prime#

    if target[-1] == 'A' or target[-1] == 'U':
        energy = AAE/2
    elif target[-1] == 'G' or target[-1] == 'C':
        energy = GGE/2
    else:
        print ""
        print "Error: none of the stardard bases were detected."

    E3prime = energy
    #print E3prime#

    if E3prime - E5prime > fudge:
        recog = True
    else:
        recog = False

    return recog

# Valid subroutine. Checks whether the four nucleotides immediately 3' of the
# end of a target sequence (beginning at nucleotide number 'position') are
# 'A' or 'U'; returns a score equal to the number of 'A' or 'U' instances
# found among these nucleotides.

def checkflanks(sequence,position):
    i = position + 19
    score = 0
    while i < position + 23:
        if sequence[i] == 'A' or sequence[i] == 'U':
            score = score + 1
        i = i + 1
    return score

# Generates a list of random sequences, DNA or RNA, with a length of the user's
# choice.

def generate_sequences(number, length = 19, RNA = True):
    if RNA:
        string = 'AGCU'
    else:
        string = 'AGCT'

    sequences = []
    i = 0
    while i < number:
        sequence = ''
        j = 0
        while j < length:
            sequence = sequence + random.choice(string)
            j = j + 1
        sequences = sequences + [sequence]
        i = i + 1

    return sequences

# Validates hits according to Place et al.'s saRNA detection heuristic.

def validate(sequence,targetstart,targetend,checkislands = True):
    targetseq = sequence[targetstart:targetend]
    
    CpGs = CpGsitedetector(sequence)
    
    flag = True

    if checkconsecutive(targetseq):
        flag = False

    if endstability(targetseq) < 0:
        flag = False

    if countGC(targetseq) < 0.4 or countGC(targetseq) > 0.65:
        flag = False

    if checkpos19(targetseq) == False:
        flag = False

    if checkpos18(targetseq) == 0:
        flag = False

    if checkpos7(targetseq) == 0:
        flag = False

    if endbase(targetseq) == False:
        flag = False

    if insite(targetstart,targetend,CpGs):
        flag = False

    if checkislands:
        islandboundaries = CpGislanddetector(sequence)
        
        if inisland(targetstart,targetend,islandboundaries):
            flag = False

    return flag

# Incorporates rules for elimination from the Qiagen website
# http://www.qiagen.com/products/genesilencing/hponguardsirnadesign.aspx
# In particular, seed region is nucleotides 2-7 on the antisense strand
# meaning it's 13-18 (inclusive) on the target site.  This, plus 10 more
# matches anywhere else on the siRNA, constitutes sufficient evidence for
# cross-targeting to eliminate the match.

def seedfinder(BLAST_result_19_bp):
    badmatch = False
    page = BLAST_result_19_bp

    if page.find('No significant similarity found.') == -1:
        resultstart = page.find('ALIGNMENTS')
        resultpage = page[resultstart:-1]
        resultstart = resultpage.find('>ref',2)
        resultpage = resultpage[resultstart:-1]
        
        while resultpage.find('>ref') != -1 and badmatch == False:
            querylinestart = resultpage.find('Query')
            querylineend = resultpage.find('\n',querylinestart)
            queryline = resultpage[querylinestart:querylineend]

            space = queryline.rfind(' ')
            endbp = queryline[space + 1:len(queryline)]
            endbp = int(endbp)

            bondlinestart = querylineend + 2
            bondlineend = resultpage.find('\n',bondlinestart)
            bondline = resultpage[bondlinestart:bondlineend]

            if bondline.count('|') > 15:
                if endbp == 19 and bondline[-7:-1] == '||||||':
                    badmatch = True
                elif endbp == 18 and bondline[-6:len(bondline)] == '||||||':
                    badmatch = True

            resultstart = resultpage.find('>ref',2)
            resultpage = resultpage[resultstart:-1]

    return badmatch

# Incorporates rules for elimination from the Qiagen website
# http://www.qiagen.com/products/genesilencing/hponguardsirnadesign.aspx
# In particular, seed region is nucleotides 2-7 on the antisense strand
# meaning it's 13-18 (inclusive) on the target site.  This, plus 10 more
# matches anywhere else on the siRNA, constitutes sufficient evidence for
# cross-targeting to eliminate the match.  Returns every gene flagged as
# an off-target hit for the input sequence.

def seedidentifier(BLAST_result_19_bp):
    badmatches = []
    page = BLAST_result_19_bp

    if page.find('No significant similarity found.') == -1:
        resultstart = page.find('ALIGNMENTS')
        resultpage = page[resultstart:-1]
        resultstart = resultpage.find('>ref',2)
        resultpage = resultpage[resultstart:-1]
        
        while resultpage.find('>ref') != -1:
            querylinestart = resultpage.find('Query')
            querylineend = resultpage.find('\n',querylinestart)
            queryline = resultpage[querylinestart:querylineend]

            space = queryline.rfind(' ')
            endbp = queryline[space + 1:len(queryline)]
            endbp = int(endbp)

            bondlinestart = querylineend + 2
            bondlineend = resultpage.find('\n',bondlinestart)
            bondline = resultpage[bondlinestart:bondlineend]

            if bondline.count('|') > 15:
                if endbp == 19 and bondline[-7:-1] == '||||||':
                    genestart = resultpage.find('(')
                    geneend = resultpage.find(')')
                    geneid = resultpage[genestart + 1:geneend]
                    badmatches = badmatches + [geneid]
                elif endbp == 18 and bondline[-6:len(bondline)] == '||||||':
                    genestart = resultpage.find('(')
                    geneend = resultpage.find(')')
                    geneid = resultpage[genestart + 1:geneend]
                    badmatches = badmatches + [geneid]

            resultstart = resultpage.find('>ref',2)
            resultpage = resultpage[resultstart:-1]

    return badmatches
                

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

# The following subroutines are designed to be used directly from the IDLE shell by the user.

# addtolib: This adds a sequence, whose name you provide, to your library.  addlib can determine whether
# the identifier passed to it is that of a gene or a micro, downloads the sequence from its reference
# databases, and returns an error message if the sequence cannot be found.

def addtolib(sequence_name):
    flag = calladdtolib(sequence_name)

# addlib: When given no arguments, addlib adds to your library all microRNAs and genes currently listed
# in the documents within the User file folder.  This method can detect whether the sequences are
# user-submitted (in which case it passes the sequences directly to its library) or absent (in which
# case it searches for the sequences online based on the identifiers provided by the user).
# NOTE: addlib assumes that user-submitted miRNA sequences refer to the target region, NOT to the
# guide strand of the miRNA itself; as a result, these sequences are reverse-complemented, then
# stored in your library (so you library will contain the REVERSE COMPLEMENTS of the sequences
# submitted to RNAeye; that is, the guide strands).  miRNAs that have been downloaded online, however,
# are kept as-is, because online sequences refer to the actual guide strand.

def addlib():
    f = open(homefolder + 'User/miRNA IDs.txt')
    microdoc = f.read()
    microdoc = microdoc.strip()
    
    if microdoc == 'miRNA_ID\tmiRNA_Sequence_(Optional)':
        microinput = []
    else:
        microinput = loadtxt(homefolder + 'User/miRNA IDs.txt', dtype = str, skiprows = 1)

    f = open(homefolder + 'User/Gene IDs.txt')
    genedoc = f.read()
    genedoc = genedoc.strip()
    
    if genedoc == 'GeneID\t\tPromoter_Sequence_(Optional)\t\tTranscript_Sequence_(Optional)':
        geneinput = []
    else:
        geneinput = loadtxt(homefolder + 'User/Gene IDs.txt', dtype = str, skiprows = 1)

    if len(shape(microinput)) == 1 and len(shape(geneinput)) == 1:
        microids = microinput
        geneids = geneinput

        if shape(geneinput)[0] != 0:
            print "Downloading promoter sequences..."
            DownloadPromoters(geneids,homefolder)
            print ""
            print "Downloading transcript sequences..."
            DownloadTranscripts(geneids,homefolder)
            print ""

        else:
            print "Your list of gene IDs is blank."

        if shape(microinput)[0] != 0:            
            print "Downloading microRNA sequences.."
            DownloadMicros(microids,homefolder)
            print ""
            print "The submitted sequences have been added to your library."

        else:
            print "Your list of microRNAs is blank."

    elif len(shape(microinput)) == 2 and len(shape(geneinput)) == 2:
        if shape(microinput)[1] == 2 and shape(geneinput)[1] == 3:
            geneids = geneinput[:,0]
            geneproseqs = geneinput[:,1]
            genetransseqs = geneinput[:,2]
            microids = microinput[:,0]
            microseqs = microinput[:,1]

            flag = lengthcheck(geneids,geneproseqs)
            if flag == False:
                print "Error: your list of gene names is different in length from your list of gene promoter sequences."
            else:
                flag = lengthcheck(geneids,genetransseqs)
                if flag == False:
                    print "Error: your list of gene names is different in length from your list of gene transcript sequences."
                else:
                    flag = lengthcheck(microids,microseqs)
                    if flag == False:
                        print "Error: your list of microRNA names is different in length from your list of microRNA sequences."
                    else:
                        flag = lengthcheck(geneids,microids)
                        if flag == False:
                            print "Error: your list of microRNA names is different in length from your list of gene names."
                        else:
                            print "Adding gene sequences to library..."
                            addgenes(geneids,geneproseqs,genetransseqs,homefolder)
                            print "Adding microRNA sequences to library..."
                            addmicros(microids,microseqs,homefolder)

        else:
            print "Error: one or both of your sequence lists contain an incorrect number of columns."

    elif len(shape(microinput)) == 2 and len(shape(geneinput)) == 1:        
        if shape(microinput)[1] == 2:
            microids = microinput[:,0]
            microseqs = microinput[:,1]
            geneids = geneinput

            flag = lengthcheck(microids,microseqs)
            if flag == False:
                print "Error: your list of microRNA names is different in length from your list of microRNA sequences."
            
            else:
                if shape(geneinput)[0] != 0:
                    print "Downloading promoter sequences..."
                    DownloadPromoters(geneids,homefolder)
                    print ""
                    print "Downloading transcript sequences..."
                    DownloadTranscripts(geneids,homefolder)
                    print ""
                    
                print "Adding microRNA seqeuences to library..."
                addmicros(microids,microseqs,homefolder)
                print ""
                print "The submitted sequences have been added to your library."

        else:
            print "Error: your list of miRNA IDs list contains an incorrect number of columns."

    else:
        if len(shape(microinput)) == 1 and shape(geneinput)[1] == 3:
            microids = microinput
            geneids = geneinput[:,0]
            geneproseqs = geneinput[:,1]
            genetransseqs = geneinput[:,2]
        
            flag = lengthcheck(geneids,geneproseqs)
            if flag == False:
                print "Error: your list of gene names is different in length from your list of gene promoter sequences."
            else:
                flag = lengthcheck(geneids,genetransseqs)
                if flag == False:
                    print "Error: your list of gene names is different in length from your list of gene transcript sequences."

                else:
                    print "Adding gene sequences to library..."
                    addgenes(geneids,geneproseqs,genetransseqs,homefolder)
                    print ""

                    if shape(microinput) != 0:
                        print "Downloading microRNA sequences..."
                        DownloadMicros(microids,homefolder)
                        print ""
                        
                    print "The submitted sequences have been added to your library."

        else:
            print "Error: your Gene IDs list contains an incorrect number of columns."

# fullscan: This is the method responsible for downloading all necessary sequences into
# your library (or recording sequences it has been given), analyzing all sequences
# for activation and interference potential, and returning the alignment scores and
# maximal binding energies to the Activation Scores and Interference Scores folders;
# essentially, this is a catch-all method.

def fullscan():
    addlib()

    microids = genfromtxt(homefolder + 'User/miRNA IDs.txt', dtype = str, skip_header = 1, usecols = 0)
    geneids = genfromtxt(homefolder + 'User/Gene IDs.txt', dtype = str, skip_header = 1, usecols = 0)

    print ""
    print "Running miRanda on interfering sequences and transcripts..."
    miRandatuRboint(mirandalocation,thresholdenergy,thresholdalign,microids,geneids,homefolder)
    print "Compiling miRanda interference scores..."
    interscores = miRandascoResint(microids,geneids,homefolder)

    interalign = interscores[0:len(microids)]
    interenergies = interscores[len(microids):2*len(microids)]

    print ""
    print "Running miRanda on activating sequences and promoters..."
    MakeReverseSet(microids,geneids,homefolder)
    miRandatuRboact(mirandalocation,microids,geneids,homefolder)
    alldata = UnscrambleOutput(microids,geneids,homefolder)

    targetstarts = alldata[0:len(microids)]
    targetends = alldata[len(microids):2*len(microids)]
    micros = alldata[2*len(microids):3*len(microids)]
    targets = alldata[3*len(microids):4*len(microids)]
    alignments = alldata[6*len(microids):7*len(microids)]
    energies = alldata[5*len(microids):6*len(microids)]

    print "Compiling miRanda activation scores..."
    modparameters = ActivPreAnalyzer(targetstarts,targetends,micros,targets,alignments,energies,microids,geneids,homefolder)
    modalignments = modparameters[0:len(microids)]
    modenergies = modparameters[len(microids):2*len(microids)]

    activscores = miRandascoResact(modalignments,modenergies,microids,geneids,homefolder)

    activalign = activscores[0:len(microids)]
    activenergies = activscores[len(microids):2*len(microids)]
    print ""
    print "RNAeye has finished scoring the sequences provided."
    
"""
    data = genfromtxt(homefolder + 'User/External Data.txt', dtype = str, skip_header = 1, usecols = 0)

    activindices = indexscores(activalign)
    interindices = indexscores(interalign)

    newinteralign = reducescores(interalign, interindices)
    newinterenergies = reducescores(interenergies, interindices)
    newactivalign = reducescores(activalign, activindices)
    newactivenergies = reducescores(activenergies, activindices)
    newdata = reducescores(data, activindices)

    finalscores = calculate(newactivalign,newactivenergies,newinteralign,newinterenergies)
    
    figure(1)
    plot(newinteralign,newdata,'b.')
    xlabel('Maximal Alignment Scores')
    ylabel('External Data')
    title('Alignment Scores vs. External Data')

    figure(2)
    plot(newinterenergies,newdata,'b.')
    xlabel('Maximal Binding Energies')
    ylabel('External Data')
    title('Binding Energies vs. External Data')

    show()
    
    print ""
    print "RNAeye has finished scoring the sequences provided."
"""
# alignall: This method has all the functionality of fullscan, except that
# it does not download the sequences given, but merely aligns the input
# sequences that are already present in your library.

def alignall():
    f = open(homefolder + 'User/miRNA IDs.txt')
    microdoc = f.read()
    microdoc = microdoc.strip()

    f = open(homefolder + 'User/Gene IDs.txt')
    genedoc = f.read()
    genedoc = genedoc.strip()

    if microdoc == 'miRNA_ID\tmiRNA_Sequence_(Optional)' or genedoc == 'GeneID\t\tPromoter_Sequence_(Optional)\t\tTranscript_Sequence_(Optional)':
        print "Error: either your list of microRNA names, your list of gene names, or both, are blank."

    else:
        microids = genfromtxt(homefolder + 'User/miRNA IDs.txt', dtype = str, skip_header = 1, usecols = 0)
        geneids = genfromtxt(homefolder + 'User/Gene IDs.txt', dtype = str, skip_header = 1, usecols = 0)

        if lengthcheck(microids,geneids) == False:
            print "Error: your list of microRNA names is different in length from your list of gene names."

        else:
            print ""
            print "Running miRanda on interfering sequences and transcripts..."
            miRandatuRboint(mirandalocation,thresholdenergy,thresholdalign,microids,geneids,homefolder)
            print "Compiling miRanda interference scores..."
            interscores = miRandascoResint(microids,geneids,homefolder)

            interalign = interscores[0:len(microids)]
            interenergies = interscores[len(microids):2*len(microids)]

            print ""
            print "Running miRanda on activating sequences and promoters..."
            MakeReverseSet(microids,geneids,homefolder)
            miRandatuRboact(mirandalocation,microids,geneids,homefolder)
            alldata = UnscrambleOutput(microids,geneids,homefolder)

            targetstarts = alldata[0:len(microids)]
            targetends = alldata[len(microids):2*len(microids)]
            micros = alldata[2*len(microids):3*len(microids)]
            targets = alldata[3*len(microids):4*len(microids)]
            alignments = alldata[6*len(microids):7*len(microids)]
            energies = alldata[5*len(microids):6*len(microids)]

            print "Compiling miRanda activation scores..."
            modparameters = ActivPreAnalyzer(targetstarts,targetends,micros,targets,alignments,energies,microids,geneids,homefolder)
            modalignments = modparameters[0:len(microids)]
            modenergies = modparameters[len(microids):2*len(microids)]

            activscores = miRandascoResact(modalignments,modenergies,microids,geneids,homefolder)

            activalign = activscores[0:len(microids)]
            activenergies = activscores[len(microids):2*len(microids)]

"""
            data = genfromtxt(homefolder + 'User/External Data.txt', dtype = str, skip_header = 1, usecols = 0)

            activindices = indexscores(activalign)
            interindices = indexscores(interalign)

            newinteralign = reducescores(interalign, interindices)
            newinterenergies = reducescores(interenergies, interindices)
            newactivalign = reducescores(activalign, activindices)
            newactivenergies = reducescores(activenergies, activindices)
            newdata = reducescores(data, activindices)

            finalscores = calculate(newactivalign,newactivenergies,newinteralign,newinterenergies)
            
            figure(1)
            plot(newinteralign,newdata,'b.')
            xlabel('Maximal Alignment Scores')
            ylabel('External Data')
            title('Alignment Scores vs. External Data')

            figure(2)
            plot(newinterenergies,newdata,'b.')
            xlabel('Maximal Binding Energies')
            ylabel('External Data')
            titlee('Binding Energies vs. External Data')

            show()
            
            print ""
            print "RNAeye has finished scoring the sequences provided."
"""
  
# query: this searches online for a specific sequence (gene or micro) and
# downloads it if it can be found online and if it is not already in the
# user's library.  Querying a sequence for which an online search failed
# yields information on the reason for the search failure and in some cases
# generates a url that the user can access to potentially download the
# desired sequence manually.

def query(sequence_name):
    if 'XM_' in sequence_name or 'NM_' in sequence_name:
        print "Querying gene " + sequence_name + "..."
        print ""
        flag = QueryGene(sequence_name)
        
    else:
        print "Querying microRNA " + sequence_name + "..."
        print ""
        flag = QueryMicro(sequence_name)

    if flag:
        print ""
        addtolib(sequence_name)

    print ""
    print "Query completed."

# microblaster: This method BLASTs all sequences in the miRNA IDs document, downloading them
# as needed, and stores the results in the BLAST_Results folder.  This method, as well
# as all other methods that call the BLASTER subroutine, is in strict compliance with
# NCBI usage guidelines regarding script-assisted searches; see
# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node2.html
# for details.

def microblaster():    
    goodflag = addlibmicro()
    microids = genfromtxt(homefolder + 'User/miRNA IDs.txt', dtype = str, skip_header = 1, usecols = 0)
    decision = 'y'

    if datetime.now().hour > 4 and datetime.now().hour < 21 and len(microids) > 10:
        decision = raw_input("Warning: NCBI guidelines concerning bandwidth usage require that high-throughput BLAST scripts only be activated during off-peak hours; that is, between 9 pm and 5 am Eastern Standard Time.  Failure to comply with these guidelines may result in NCBI suspending service to your IP address.  Do you wish to continue? (y/n) ")

        while decision != 'n' and decision != 'y':
            decision = raw_input("Invalid response; please try again: ")

    if decision == 'y':
    
        BLASTER(microids,True,goodflag,homefolder)

        print ""
        print "RNAeye has finished BLASTing your microRNA sequences.  They can be retrieved in the folder 'BLAST_Results'."

# geneblaster: This method BLASTs all sequences in the Gene IDs document, downloading them
# as needed, and stores the results in the BLAST_Results folder.  This method, as well
# as all other methods that call the BLASTER subroutine, is in strict compliance with
# NCBI usage guidelines regarding script-assisted searches; see
# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node2.html
# for details.

def geneblaster():
    goodflag = addlibgene()
    geneids = genfromtxt(homefolder + 'User/Gene IDs.txt', dtype = str, skip_header = 1, usecols = 0)
    decision = 'y'

    if datetime.now().hour > 4 and datetime.now().hour < 21 and len(geneids) > 10:
        decision = raw_input("Warning: NCBI guidelines concerning bandwidth usage require that high-throughput BLAST scripts only be activated during off-peak hours; that is, between 9 pm and 5 am Eastern Standard Time.  Failure to comply with these guidelines may result in NCBI suspending service to your IP address.  Do you wish to continue? (y/n) ")

        while decision != 'n' and decision != 'y':
            decision = raw_input("Invalid response; please try again: ")

    if decision == 'y':
        
        BLASTER(geneids,False,goodflag,homefolder)

        print ""
        print "RNAeye has finished BLASTing your gene transcript sequences.  They can be retrieved in the folder 'BLAST_Results'."

# blaster: This method BLASTs a sequence submitted by the user, downloading it
# as needed, and stores the result in the BLAST_Results folder.  This method, as well
# as all other methods that call the BLASTER subroutine, is in strict compliance with
# NCBI usage guidelines regarding script-assisted searches; see
# http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node2.html
# for details.

def blaster(sequence_name):
    if 'NM_' in sequence_name or 'XM_' in sequence_name:
        flag = False
    else:
        flag = True

    goodflag = calladdtolib(sequence_name)
    
    BLASTER([sequence_name],flag,goodflag,homefolder)

    if goodflag:
        print ""
        print "Your sequence results are available in the folder 'BLAST_Results'."

# This subroutine manufactures a hairpin, given the hairpin's targeted sequence in the genome.
# Both strands of the hairpin are output in the shell, as ready-to-order sequences.

def makehairpin(sequence):
    sequence = checkDNAsequence(sequence)
    if sequence == 'failure':
        print "The submitted sequence included non-standard characters.  A hairpin could not be produced."
    else:
        hairpin = callmaketargethp(sequence)
        reverse = hairpin[1][::-1]
        
        bars = ""
        i = 0
        while i < len(sequence):
            bars = bars + "|"
            i = i + 1

        spaces = ""
        i = 0
        while i < len(sequence):
            spaces = spaces + " "
            i = i + 1
            
        print "Sense strand:"
        print hairpin[0]
        print ""
        print "Antisense strand:"
        print hairpin[1]
        print ""
        print ""
        print "Hairpin (genomic):"
        print ""
        print "5'  " + hairpin[0] + "      3'"
        print "        |" + bars + "||||||" + bars + "||||||"
        print "3'      " + reverse + "  5'"
        print ""
        print ""
        print "Hairpin (transcript):"
        print spaces + "    C"
        print "5'  " + sequence + " T"
        print "    " + bars + "  C"
        print "3'  " + RCgenDNA(sequence)[::-1] + "  G"
        print spaces + "    GA"

# This subroutine converts a DNA sequence to RNA.

def toRNA(sequence):
    sequence = checkRNAsequence(sequence)
    if sequence == 'failure':
        print "The submitted sequence included non-standard characters and could not be converted to RNA."
    else:
        return sequence

# This subroutine converts an RNA sequence to DNA.

def toDNA(sequence):
    sequence = checkDNAsequence(sequence)
    if sequence == 'failure':
        print "The submitted sequence included non-standard characters and could not be converted to DNA."
    else:
        return sequence

# This subroutine delivers the reverse complement of a submitted sequence of RNA or DNA,
# and is capable of detecting errors such as inconsistencies in bases ('U' vs 'T'), as well
# as the absence of these ambiguous bases.  Its reverse complement is either DNA or RNA,
# depending upon the sequence submitted.

def revcomp(sequence):
    sequence = sequence.upper()
    if sequence.find('T') == -1 and sequence.find('U') == -1:
        sequence = checkDNAsequence(sequence)
        if sequence == 'failure':
            print "The submitted sequence includes non-standard characters and cannot be processed."
        else:
            revcomp = RCgenDNA(sequence)
            print revcomp
            print ""
            print "Warning: No 'T' or 'U' was detected, so revcomp defaulted to a DNA complement.  If an RNA complement is specifically desired, please use the revcompR method."
    elif sequence.find('T') != -1 and sequence.find('U') != -1:
        sequence = checkDNAsequence(sequence)
        if sequence == 'failure':
            print "The submitted sequence includes non-standard characters and cannot be processed."
        else:
            print "The submitted sequence includes both 'T' and 'U' bases, so RNAeye could not determine whether a DNA or RNA complement is desired.  Please use the toRNA or toDNA methods to preprocess your sequence."
    elif sequence.find('T') == -1 and sequence.find('U') != -1:
        sequence = checkRNAsequence(sequence)
        if sequence == 'failure':
            print "The submitted sequence includes non-standard characters and cannot be processed."
        else:
            revcomp = RCgenRNA(sequence)
            print revcomp
    else:
        sequence = checkDNAsequence(sequence)
        if sequence == 'failure':
            print "The submitted sequence includes non-standard characters and cannot be processed."
        else:
            revcomp = RCgenDNA(sequence)
            print revcomp

# This subroutine delivers the reverse complement of a submitted sequence of RNA.  Any 'T'
# base is force-converted to a 'U' base.

def revcompR(sequence):
    sequence = checkRNAsequence(sequence)
    if sequence == 'failure':
        print "The submitted sequence includes non-standard characters and cannot be processed."
    else:
        revcomp = RCgenRNA(sequence)
        print revcomp

# This subroutine delivers the reverse complement of a submitted sequence of DNA.  Any 'U'
# base is force-converted to a 'T' base.

def revcompD(sequence):
    sequence = checkDNAsequence(sequence)
    if sequence == 'failure':
        print "The submitted sequence includes non-standard characters and cannot be processed."
    else:
        revcomp = RCgenDNA(sequence)
        print revcomp

# This subroutine manufactures hairpins from the miRNA IDs.txt file.  The entries
# in this file are expected to be the SENSE strand targets of the hairpins,
# NOT the guide strands!! The resulting hairpins are saved in Hairpin_List_ + file_name.txt.

def makehairpins(file_name):
    addlibmicro()
    print ""
    
    targetids = genfromtxt(homefolder + 'User/miRNA IDs.txt', dtype = str, skip_header = 1, usecols = 0)

    i = 0
    hairpinlist = []
    while i < len(targetids):
        f = open(homefolder + 'microRNA_Sequences/' + targetids[i] + '.fasta')
        rawdoc = f.read()

        seqstart = rawdoc.find('\n') + 1
        seqend = len(rawdoc)
        sequence = rawdoc[seqstart:seqend]
        sequence = checkDNAsequence(sequence)
        sequence = RCgenDNA(sequence)
        
        if sequence == 'failure':
            hairpinlist = hairpinlist + [[]]
            print "Could not process sequence " + targetids[i] + "."
        else:
            hairpin = callmaketargethp(sequence)
            hairpinlist = hairpinlist + [hairpin]

        i = i + 1

    output = open(homefolder + 'User Archives/Hairpin_List_' + file_name + '.txt', 'w')

    i = 0
    while i < len(hairpinlist):
        seqlength = (len(hairpinlist[i][0]) - 17)/2
        j = 0
        bars = ''
        spaces = ''
        while j < seqlength:
            bars = bars + '|'
            spaces = spaces + ' '
            j = j + 1
        
        if hairpinlist[i] != []:
            seq = hairpinlist[i][0][5:seqlength + 5]
            antiseq = RCgenDNA(hairpinlist[i][0][5:seqlength + 5])
            
            output.write('----------------------------------------------------------------------\n')
            output.write(targetids[i] + '\n\n')
            output.write('Sense strand:\n')
            output.write(hairpinlist[i][0] + '\n\n')
            output.write('Antisense strand:\n')
            output.write(hairpinlist[i][1] + '\n\n')
            output.write('Hairpin (genomic):\n\n')
            output.write('5\'  ' + hairpinlist[i][0] + '      3\'\n')
            output.write('        |' + bars + '||||||' + bars + '||||||\n')
            output.write('3\'      ' + hairpinlist[i][1][::-1] + '  5\'\n\n')
            output.write('Hairpin (transcript):\n')
            output.write(spaces + '    C\n')
            output.write('5\'  ' + seq + ' T\n')
            output.write('    ' + bars + '  C\n')
            output.write('3\'  ' + antiseq[::-1] + '  G\n')
            output.write(spaces + '    GA\n\n')

        i = i + 1

    output.close()
        
    print "Hairpins constructed.  Sequences may be found in the file 'Hairpin_List_" + file_name + ".txt'."

# This subroutine deletes the sequence (whether gene or miRNA) named
# sequence_name from your library.  Any other files relating to the sequence,
# such as BLAST result logs, binding energy analysis, or miRanda results, are
# preserved.

def wipe(sequence_name):
    propath = homefolder + 'Gene_Sequences/' + sequence_name + 'pro.fasta'
    transpath = homefolder + 'Gene_Sequences/' + sequence_name + 'trans.fasta'
    mirnapath = homefolder + 'microRNA_Sequences/' + sequence_name + '.fasta'
    if os.path.isfile(mirnapath):
        if os.path.isfile(propath) or os.path.isfile(transpath):
            print "Warning: " + sequence_name + " is both a gene and a miRNA.  Please specify which you would like to delete using either the wipeG or the wipeM method."
        else:
            os.remove(mirnapath)
            print "File " + sequence_name + ".fasta has been deleted from your library."
    else:
        if os.path.isfile(propath):
            os.remove(propath)
            print "File " + sequence_name + "pro.fasta has been deleted from your library."
        if os.path.isfile(transpath):
            os.remove(transpath)
            print "File " + sequence_name + "trans.fasta has been deleted from your library."
        if os.path.isfile(propath) == False and os.path.isfile(transpath) == False:
            print "No file corresponding to sequence " + sequence_name + " could be found in your library."

# This subroutine delets the gene gene_name from your library, in the same way
# as the wipe() method; the difference is that, if there is a miRNA in your
# library of the same name, it will not be deleted (it only searches for genes).

def wipeG(gene_name):
    propath = homefolder + 'Gene_Sequences/' + gene_name + 'pro.fasta'
    transpath = homefolder + 'Gene_Sequences/' + gene_name + 'trans.fasta'
    if os.path.isfile(propath):
        os.remove(propath)
        print "File " + gene_name + "pro.fasta has been deleted from your library."
    if os.path.isfile(transpath):
        os.remove(transpath)
        print "File " + gene_name + "trans.fasta has been deleted from your library."
    if os.path.isfile(propath) == False and os.path.isfile(transpath) == False:
        print "No file corresponding to gene " + gene_name + " could be found in your library."

# This subroutine deletes the miRNA microRNA_name from your library, in the same
# way that wipeG() deletes a gene from your library.

def wipeM(microRNA_name):
    mirnapath = homefolder + 'microRNA_Sequences/' + microRNA_name + '.fasta'
    if os.path.isfile(mirnapath):
        os.remove(mirnapath)
        print "File " + microRNA_name + ".fasta has been deleted from your library."
    else:
        print "No file corresponding to miRNA " + microRNA_name + " could be found in your library."

# This subroutine deletes from your library all of the genes in Gene IDs.txt that
# are currently present; this includes any BLAST results about these genes that
# have been accumulated, so this method should be used with some caution.

def wipegenes():
    geneids = genfromtxt(homefolder + 'User/Gene IDs.txt', dtype = str, skip_header = 1, usecols = 0)

    i = 0
    while i < len(geneids):
        propath = homefolder + 'Gene_Sequences/' + geneids[i] + 'pro.fasta'
        transpath = homefolder + 'Gene_Sequences/' + geneids[i] + 'trans.fasta'
        blastpath = homefolder + 'BLAST_Results/Gene_' + geneids[i] + '.txt'
        if os.path.isfile(propath):
            os.remove(propath)
        if os.path.isfile(transpath):
            os.remove(transpath)
        if os.path.isfile(blastpath):
            os.remove(blastpath)
        i = i + 1

    print "All genes in Gene IDs.txt have been deleted from your library."

# This subroutine deletes from your library all of the microRNAs in miRNA IDs.txt that
# are currently present; this includes any BLAST results about these miRNAs that
# have been accumulated, so this method should be used with some caution.

def wipemirnas():
    mirnaids = genfromtxt(homefolder + 'User/miRNA IDs.txt', dtype = str, skip_header = 1, usecols = 0)

    i = 0
    while i < len(mirnaids):
        mirpath = homefolder + 'microRNA_Sequences/' + mirnaids[i] + '.fasta'
        if os.path.isfile(mirpath):
            os.remove(mirpath)
        blastpath = homefolder + 'BLAST_Results/microRNA_' + mirnaids[i] + '.txt'
        if os.path.isfile(blastpath):
            os.remove(blastpath)
        i = i + 1

    print "All miRNAs in miRNA IDs.txt have been deleted from your library."

# This subroutine deletes from your library all of the genes in Gene IDs.txt
# and all of the microRNAs in miRNA IDs.txt that are currently present; this
# includes any BLAST results about these genes that have been accumulated,
# so this method should be used with some caution.

def wipeall():
    wipegenes()
    wipemirnas()

# This subroutine clears the entire contents of your library, subject to
# final confirmation by the user.  All data will be wiped, and is non-
# recoverable, so this method should be used only when clearing out
# a library to port the build to another machine.

def clearlib():
    decision = raw_input('Warning: this command will delete the entire contents of your library.  Proceed? (y/n):')

    while decision != 'y' and decision != 'n':
        decision = raw_input('Invalid response; please try again:')

    if decision == 'y':
        genepath = homefolder + 'Gene_Sequences/'
        mirnapath = homefolder + 'microRNA_Sequences/'
        blastpath = homefolder + 'BLAST_Results/'

        for targetfile in os.listdir(genepath):
            filepath = os.path.join(genepath, targetfile)
            os.remove(filepath)

        for targetfile in os.listdir(mirnapath):
            filepath = os.path.join(mirnapath, targetfile)
            os.remove(filepath)

        for targetfile in os.listdir(blastpath):
            filepath = os.path.join(blastpath, targetfile)
            os.remove(filepath)

        print "Your library has been completely cleared of data."

# This subroutine is designed to manufacture sequencing primers for
# a submitted DNA sequence, with spaces between the primers as determined
# by the user.

def sequencingprimers(primer_length, target_sequence, distance_between_primers = 600):
    global seqsurvivors
    global seqpositions
    
    if primer_length > 30:
        print "Your desired primer length is too long."
    elif primer_length < 18:
        print "Your desired primer length is too short."
    else:
        if len(target_sequence) < 100:
            targets = maketargets(target_sequence,primer_length)
        else:
            targets = maketargets(target_sequence[0:100], primer_length)

        i = 0
        survivors = [[]]
        positions = [[]]
        while i < len(targets):
            flag = True
                
            if threeprimetest(targets[i]) == False:
                flag = False
            if checkconsecutive(targets[i]):
                flag = False
            GC = countGC(targets[i])
            if GC < 0.4 or GC > 0.6:
                flag = False

            if flag:
                survivors[0] = survivors[0] + [targets[i]]
                positions[0] = positions[0] + [i + 1]

            i = i + 1

        runs = len(target_sequence)/distance_between_primers
            
        if runs != 0:
            i = 0
            while i < runs:
                survivors = survivors + [[]]
                positions = positions + [[]]

                if len(target_sequence) > distance_between_primers*(i + 1) + 100:
                    targets = maketargets(target_sequence[distance_between_primers*(i + 1):distance_between_primers*(i + 1) + 100],primer_length)

                    j = 0
                    while j < len(targets):
                        flag = True

                        if threeprimetest(targets[j]) == False:
                            flag = False
                        if checkconsecutive(targets[j]):
                            flag = False
                        GC = countGC(targets[j])
                        if GC < 0.4 or GC > 0.6:
                            flag = False

                        if flag:
                            survivors[-1] = survivors[-1] + [targets[j]]
                            positions[-1] = positions[-1] + [distance_between_primers*(i + 1) + j + 1]

                        j = j + 1

                i = i + 1

        print "Your FWD sequencing primers are the following:"
        print ""

        i = 0
        while i < len(survivors):
            print "----------------------------------------------"
            print "Primer " + str(i + 1) + ":"
            if len(survivors[i]) == 0:
                print "No suitable primer could be found between bases " + str(distance_between_primers*i) + " and " + str(distance_between_primers*i + 100) + "."
                print ""
            else:
                print survivors[i][0]
                print "Positions " + str(positions[i][0]) + " to " + str(positions[i][0] + primer_length - 1) + "."
                print ""
            i = i + 1

        print ""
        print "For additional sequences, you may use the seqsurvivors and seqpositions arrays."

        seqsurvivors = survivors
        seqpositions = positions

# This subroutine calculates the percent GC content of an input sequence.

def GCperc(sequence):
    frac = countGC(sequence)
    percent = round(100.*frac, 1)
    print str(percent) + "% GC."

# This method calculates the melting temperature of a sequence with its
# reverse complement. (Only one sequence is needed as input.)

def findTm(sequence):
    if len(sequence) < 14:
        Tm = smalloligoTm(sequence)
    else:
        Tm = largeoligoTm(sequence)

    Tm = round(Tm, 1)

    print "Tm = " + str(Tm) + " degrees Celsius."

# This subroutine detects CpG islands in a submitted sequence.

def CpGislands(sequence):
    islands = CpGislanddetector(sequence)
    print "Sequence length is " + str(len(sequence)) + " bp."
    print ""
    if len(islands) == 0:
        print "No CpG islands were detected in your sequence."
    else:
        print "The following CpG islands were detected:"
        print ""
        i = 0
        while i < len(islands):
            print "Island " + str(i + 1) + ":"
            print "Positions " + str(islands[i][0] + 1) + " to " + str(islands[i][1] + 1) + "."
            print "Sequence:"
            print sequence[islands[i][0]:islands[i][1] + 1]
            print ""

            i = i + 1

        i = 0
        while i < len(islands):
            plot(islands[i],[0,0],'b-')
            i = i + 1

        xlabel('Nucleotide position')
        show()

# This subroutine detects CpG sites in a submitted sequence.

def CpGsites(sequence):
    sites = CpGsitedetector(sequence)
    print "Sequence length is " + str(len(sequence)) + " bp."
    if len(sites) == 0:
        print "No CpG sites were detected in your sequence."
        print ""
    else:
        seqstring = ''
        i = 0
        while i < len(sites):
            if i < len(sites) - 1:
                seqstring = seqstring + str(sites[i] + 1) + ', '
            else:
                seqstring = seqstring + str(sites[i] + 1) + '.'
            i = i + 1

        print "CpG sites were detected at positions " + seqstring

        plot(sites,[0]*len(sites),'b|')
        xlabel('Nucleotide position')
        show()
        
# This subroutine generates a list of activating RNA targets compatible with
# the rules laid out in Place et al., then BLASTs them against the human genome
# and eliminates off-target hits.

def make_activators(number,batch_name,length = 19):
    print "Generating sequences.."
    candidates = generate_sequences(number,length,True)
    i = 0
    survivors = []
    print "Validating candidates.."
    while i < len(candidates):
        if validate(candidates[i],0,-1,False):
            survivors = survivors + [candidates[i]]
        i = i + 1

    ans = ''
    print ""
    while ans != 'n' and ans != 'y':
        ans = raw_input(str(len(survivors)) + ' potential targets have been found.  Cross-validate with BLASTER? (y/n)')

        if ans == 'y':
            print "Saving data to your library.."
            i = 0
            survivornames = []
            while i < len(survivors):
                output = open(homefolder + 'microRNA_Sequences/' + batch_name + '_' + str(i + 1) + '.fasta', 'w')
                output.write('>' + batch_name + ' survivor ' + str(i + 1) + '.\n')
                output.write(survivors[i])
                output.close()
                
                survivorname = batch_name + '_' + str(i + 1)
                survivornames = survivornames + [survivorname]

                i = i + 1

            if length < 20:
                shortflag = True
            else:
                shortflag = False

            print "BLASTing sequences.."
            BLASTER(survivornames,shortflag,True,homefolder)

            print ""
            print "Scanning BLAST results for high-relevance seed sequence matches.."
            
            i = 0
            while i < len(survivors):
                f = open(homefolder + 'BLAST_Results/microRNA_' + batch_name + '_' + str(i + 1) + '.txt')
                blastresult = f.read()
                f.close()

                badmatch = seedfinder(blastresult)

                if badmatch:
                    survivornames[i] = ''
                    survivors[i] = ''

                i = i + 1

            i = 0
            while i < len(survivors):
                if survivors[i] == '':
                    gone = survivornames.pop(i)
                    gone = survivors.pop(i)

                i = i + 1

        output = open(homefolder + 'User Archives/Activators_' + batch_name + '.txt', 'w')
        output.write('saRNA candidate\t\tSequence\n')
        i = 0
        while i < len(survivors):
            output.write(str(i + 1) + '\t\t' + survivors[i] + '\n')
            i = i + 1

        output.close()

        print ""
        print "Candidate sequences may be found in the file 'Activators_" + batch_name + ".txt' in the User Archives folder."

# Finds off-target sequences by doing a BLAST search.  The sequences in
# the file miRNA_IDs should be those of the GUIDE strands of the relevant
# duplexes.

def off_targets(batch_name,under_20_bp = True):
    
    print "Saving data to your library.."
    
    sRNAids = genfromtxt(homefolder + 'User/miRNA IDs.txt', dtype = str, skip_header = 1, usecols = 0)
    sRNAseqs = genfromtxt(homefolder + 'User/miRNA IDs.txt', dtype = str, skip_header = 1, usecols = 1)
    
    i = 0
    while i < len(sRNAids):
        output = open(homefolder + 'microRNA_Sequences/' + sRNAids[i] + '.fasta', 'w')
        output.write('> Sequence ' + sRNAids[i] + '.\n')
        output.write(sRNAseqs[i])
        output.close()
        
        i = i + 1
                
    print "BLASTing sequences.."
    BLASTER(sRNAids,under_20_bp,True,homefolder)

    print ""
    print "Scanning BLAST results for high-relevance seed sequence matches.."
            
    i = 0
    badmatches = []
    while i < len(sRNAids):
        badmatches = badmatches + [[]]
        f = open(homefolder + 'BLAST_Results/microRNA_' + sRNAids[i] + '.txt')
        blastresult = f.read()
        f.close()

        badmatches[i] = seedidentifier(blastresult)

        i = i + 1

    output = open(homefolder + 'User Archives/Off-Targets_' + batch_name + '.txt', 'w')
    output.write('sRNA candidate\t\tSequence\t\tOff-Target Hits\n')
    i = 0
    while i < len(sRNAids):
        output.write(sRNAids[i] + '\t\t' + sRNAseqs[i] + '\t\t')
        
        if len(badmatches[i]) == 0:
            output.write('OK\n')
        else:
            
            output.write(badmatches[i][0])
            j = 1
            while j < len(badmatches[i]):
                output.write(',' + badmatches[i][j])
                j = j + 1
        
            output.write('\n')
                     
        i = i + 1
        
    output.close()

    print ""
    print "Off-target hits may be accessed in the file 'Off-Targets_" + batch_name + ".txt' in the User Archives folder."

# This method returns the off-target matches of a given sequence,
# provided that the sequence is named and in your library.

def off_target(sequence_name):
    i = 0
    badmatches = []
    f = open(homefolder + 'BLAST_Results/microRNA_' + sequence_name + '.txt')
    blastresult = f.read()
    f.close()

    badmatches = seedidentifier(blastresult)

    print "Off-target matches for " + sequence_name + ":"
    if len(badmatches) == 0:
        print "None."
    else:
        matches = badmatches[0]
        i = 1
        while i < len(badmatches):
            matches = matches + ',' + badmatches[i]
            i = i + 1
        matches = matches + '.'
        print matches

# This subroutine saves a user-submitted function to your library under the
# folder Saved_Functions/.  The inputs are (1) a string containing the function
# itself, as it would be written in Python, including x as the independent
# variable and parameters of the user's choice (for example, 'x**2/(a + b)')
# (2) a list of strings of parameter values (for example, ['a','b']),
# (3) a string including the name of the function (under which it is to be saved
# and addressed, and (4) a string containing whatever miscellaneous information
# the user chooses to include regarding the function (default = the empty
# string).
# This function checks to make sure that all parameters in the function have
# been defined by actually running the function silently in the shell with all
# parameters set to zero, and checking to see whether a NameError is thrown.
# (Whether we divide by zero is irrelevant; all that matters is whether
# Python returns a NameError.)  This subroutine also checks to make sure that
# no function of the same name already exists in your library; if it does,
# the user is prompted to authorize an overwrite.  Lastly, this subroutine
# calls load_function() to immediately bring the saved function into active
# memory so it can be used immediately. The function name should contain no
# spaces or special characters; the subroutine checks for this as well.

def save_function(function,parameter_list,variable_list,function_name,function_info = ''):
    i = 0
    while i < len(parameter_list):
        exec(parameter_list[i] + ' = 0.')
        i = i + 1

    i = 0
    while i < len(variable_list):
        exec(variable_list[i] + ' = 0.')
        i = i + 1
        
    try:
        exec('def ' + function_name + '():pass')
    except SyntaxError:
        print "Warning: '" + function_name + "' is not a valid function name.  Please avoid spaces and special characters other than underscores."
    else:
        try:
            exec('y = ' + function)
        except ZeroDivisionError:
            flag = True
        except SyntaxError:
            print "Warning: " + function_name + " is not a valid function and could not be saved.  Please correct the syntax and try again."
            flag = False
        except NameError:
            print "Warning: one or more parameters or variables in your function have not been included in your parameter or variable lists.  The function " + function_name + " is thus not fully defined and could not be saved."
            flag = False
        else:
            flag = True

        if flag:
            if len(variable_list) == 0:
                print "Warning: a function must have at least one variable."
            else:
                variable_list = list(set(variable_list))
                funcdescript = 'NAME$ ' + function_name + '\n'
                if function_info == '':
                    funcdescript = funcdescript + 'INFO$ <NONE>\n'
                else:
                    funcdescript = funcdescript + 'INFO$ ' + function_info + '\n'

                i = 0
                while i < len(variable_list):
                    funcdescript = funcdescript + 'VAR' + str(i) + '$ ' + variable_list[i] + '\n'
                    i = i + 1

                parameter_list = list(set(parameter_list))
                i = 0
                while i < len(parameter_list):
                    funcdescript = funcdescript + 'PAR' + str(i) + '$ ' + parameter_list[i] + '\n'
                    i = i + 1

                funcdescript = funcdescript + 'FUNC$ ' + function + '\n'
                    
                ans = ''
                if os.path.isfile(homefolder + 'Saved_Functions/' + function_name + '.txt'):
                    while ans != 'y' and ans != 'n':
                        ans = raw_input("A function by the name " + function_name + " already exists in your library.  Overwrite? (y/n)")

                    if ans == 'y':
                        output = open(homefolder + 'Saved_Functions/' + function_name + '.txt', 'w')
                        output.write(funcdescript)
                        output.close()
                        print "The function " + function_name + " has been saved to your library."
                        load_function(function_name)

                else:
                    output = open(homefolder + 'Saved_Functions/' + function_name + '.txt', 'w')
                    output.write(funcdescript)
                    output.close()
                    print "The function " + function_name + " has been saved to your library."
                    load_function(function_name)

    i = 0
    while i < len(parameter_list):
        exec('del ' + parameter_list[i])
        i = i + 1

    i = 0
    while i < len(variable_list):
        exec('del ' + variable_list[i])
        i = i + 1

# This subroutine loads a function from your library into active memory, where it can
# be used by the program.  The subroutine raises a flag if the function name you have
# called does not exist in your library, or if you have already assigned a variable
# to active memory with the same name as the function you would like to load.  The
# subroutine parses the function file for the name, function information, parameter
# list, and the function itself.  It also generates an new subroutine using the exec
# command, that has the same name as the function itself (and no arguments) that the
# user can call to display the function information in the shell in a meaningful
# format.  The function is also added to the list of active functions, which is
# then re-sorted in alphabetical order for listing purposes.

def load_function(function_name):
    global activefunctions
    flag = False
    
    if os.path.isfile(homefolder + 'Saved_Functions/' + function_name + '.txt') == False:
        print "There is no function by the name " + function_name + " in your library."
    else:
        try:
            exec(function_name)
        except NameError:
            flag = True
        else:
            ans = ''
            while ans != 'n' and ans != 'y':
                ans = raw_input("Warning: there is already a variable by the name " + function_name + " in your current shell session. Overwrite? (y/n)")

            if ans == 'y':
                flag = True

        if flag:
            f = open(homefolder + 'Saved_Functions/' + function_name + '.txt')
            rawfunc = f.read()
            namestart = rawfunc.find('$') + 2
            nameend = rawfunc.find('\n')
            name = rawfunc[namestart:nameend]

            rawfunc = rawfunc[nameend + 1:len(rawfunc)]
            infostart = rawfunc.find('$') + 2
            infoend = rawfunc.find('\n')
            info = rawfunc[infostart:infoend]

            rawfunc = rawfunc[infoend + 1:len(rawfunc)]

            variables = []
            i = 0
            while rawfunc.find('VAR' + str(i) + '$') != -1:
                varstart = rawfunc.find('$') + 2
                varend = rawfunc.find('\n')
                variable = rawfunc[varstart:varend]
                variables = variables + [variable]
                rawfunc = rawfunc[varend + 1:len(rawfunc)]
                i = i + 1
            
            parameters = []
            i = 0
            while rawfunc.find('PAR' + str(i) + '$') != -1:
                parstart = rawfunc.find('$') + 2
                parend = rawfunc.find('\n')
                parameter = rawfunc[parstart:parend]
                parameters = parameters + [parameter]
                rawfunc = rawfunc[parend + 1:len(rawfunc)]
                i = i + 1

            funcstart = rawfunc.find('$') + 2
            funcend = rawfunc.find('\n')
            function = rawfunc[funcstart:funcend]

            funcdefinition = 'def ' + function_name + '(array_of_inputs'
            for parameter in parameters:
                funcdefinition = funcdefinition + ',' + parameter

            funcdefinition = funcdefinition + '):\n'
            if len(variables) == 1:
                funcdefinition = funcdefinition + '\tif len(shape(array_of_inputs)) != 1:\n'
            else:
                funcdefinition = funcdefinition + '\tif len(array_of_inputs) != len(_' + function_name + '[2]):\n'
                
            funcdefinition = funcdefinition + '\t\tprint "Warning: the first dimension of the input array must match the number of input variables to the function."\n'
            funcdefinition = funcdefinition + '\telse:\n'

            if len(variables) == 1:
                funcdefinition = funcdefinition + '\t\t' + variables[0] + ' = array_of_inputs\n'
            else:
                i = 0
                while i < len(variables):
                    funcdefinition = funcdefinition + '\t\t' + variables[i] + ' = array_of_inputs[' + str(i) + ']\n'
                    i = i + 1
            funcdefinition = funcdefinition + '\t\ty = ' + function + '\n'
            funcdefinition = funcdefinition + '\t\treturn y'
            
            outfile = open(homefolder + 'Developer Toolbox/output.txt','w')#
            outfile.write(funcdefinition)#
            outfile.close()#
            
            exec(funcdefinition) in globals()
            
            exec('_' + function_name + ' = ["' + name + '","' + info + '",' + str(variables) + ',' + str(parameters) + ',"' + function + '"]') in globals()

            exec('activefunctions = activefunctions + ["' + function_name + '"]') in globals()
            
            activefunctions = list(sort(list(set(activefunctions))))[::-1]

            print "Function " + function_name + " has been successfully loaded into active memory."

# This subroutine allows the user to ping the shell for information
# on a particular function or dataset, which is then displayed conveniently.
# The function name, information, parameters, and form are all listed.

def ping(function_or_dataset):
    funcflag = False
    dataflag = False
    
    if function_or_dataset in activefunctions:
        funcflag = True
    elif function_or_dataset in activedatasets:
        dataflag = True
    else:
        print "That function or dataset is not in active memory.  Please use the active_functions() or active_datasets() methods to determine which functions and datasets are currently active, and the load_function() or load_dataset() methods to activate a library function or library dataset."

    if funcflag and dataflag:
        print "Warning: RNAeye has detected both a function and a dataset by the name " + function_or_dataset + ".  Information on both will be displayed below:"
        print ""
    
    if funcflag:
        exec('obj = _' + function_or_dataset)
        
        if obj[0] in activefunctions:
            print "Function name: " + obj[0]
            print "Function information: " + obj[1]

            i = 0
            while i < len(obj[2]):
                print "Variable " + str(i) + ": " + obj[2][i]
                i = i + 1
                
            i = 0
            while i < len(obj[3]):
                print "Parameter " + str(i) + ": " + obj[3][i]
                i = i + 1

            print "Function: " + obj[4]
        else:
            print "Warning: " + function_or_dataset + " is not currently an active function."
        print ""

    if dataflag:
        exec('obj = _' + function_or_dataset)
        
        if obj[0] in activedatasets:
            print "Dataset name: " + obj[0]
            print "Dataset information: " + obj[3]
            print ""
            
            varlist = ""
            i = 0
            while i < len(obj[1]):
                varlist = varlist + obj[1][i] + "\t"
                i = i + 1

            print varlist

            i = 0
            while i < len(obj[2][0]):
                numline = ""
                j = 0
                while j < len(obj[1]):
                    numline = numline + str(obj[2][j][i]) + "\t"
                    j = j + 1
                print numline
                i = i + 1

        else:
            print "Warning: " + function_or_dataset + " is not currently an active dataset."

# This subroutine lists all of the functions that are currently active,
# meaning that they can be accessed by the user for standard operations.

def active_functions():
    if len(activefunctions) == 0:
        print "There are currently no functions in active memory."
    else:
        print "Functions currently in active memory:"
        i = 0
        act = activefunctions[::-1]
        while i < len(activefunctions):
            print act[i]
            i = i + 1

# This subroutine lists all of the functions that are currently in
# your library, which means that they are available to be loaded
# into the shell memory using the load_function() subroutine.

def library_functions():
    libfunctions = os.listdir(homefolder + 'Saved_Functions/')
    if len(libfunctions) == 0 or libfunctions == ['.DS_Store']:
        print "There are currently no functions in your library."
    else:
        libfunctions = list(sort(libfunctions))
        libfunctions.remove('.DS_Store')
        print "Functions currently in your library:"
        i = 0
        while i < len(libfunctions):
            funcname = libfunctions[i].replace('.txt','')
            print funcname
            i = i + 1

# This subroutine loads all functions in your library to the active
# shell.  In the absence of any library functions, it returns an
# error message.

def load_functions():
    libfunctions = os.listdir(homefolder + 'Saved_Functions/')
    if len(libfunctions) == 0 or libfunctions == ['.DS_Store']:
        print "There are no functions saved in your library."
    else:
        libfunctions.remove('.DS_Store')
        
        for function in libfunctions:
            funcname = function.replace('.txt','')
            load_function(funcname)

# This subroutine saves the dataset input into the file 'External
# Data.txt' to your library.  The input file should be formatted
# as follows: the first line should include the variable names,
# separated by tabs; each tab-delimited column below should contain
# the measurements of the corresponding variable.  The file can
# be copy-pasted from an Excel spreadsheet.  The dataset comments
# should be a string containing any information the user would
# like to include regarding the dataset.  The dataset is then
# automatically loaded into active memory using the load_dataset()
# subroutine.

def save_dataset(dataset_name,dataset_comments = ''):
    try:
        exec(dataset_name + ' = 0.')
    except SyntaxError:
        print "Warning: The name of your dataset should include no spaces or special characters other than underscores."
    else:
        dataset = loadtxt(homefolder + 'User/External Data.txt',dtype = str)
        if len(dataset) == 0:
            print "Warning: The file 'External Data.txt' is empty."
        elif len(dataset) == 1:
            print "Warning: only the title row could be found; there is no other data in the file 'External Data.txt'."
        else:
            varnames = list(dataset[0])
            
            dataset = dataset[1:len(dataset)]
            dataset = dataset.T

            if dataset_comments == '':
                outtext = '<NONE>\n'
            else:
                outtext = dataset_comments + '\n'
                
            i = 0
            while i < len(varnames):
                outtext = outtext + varnames[i] + '\t'
                i = i + 1

            outtext = outtext + '\n'
            i = 0
            while i < len(dataset[0]):
                j = 0
                while j < len(varnames):
                    outtext = outtext + dataset[j][i] + '\t'
                    j = j + 1
                outtext = outtext + '\n'
                i = i + 1

            output = open(homefolder  + 'Saved_Datasets/' + dataset_name + '.txt','w')
            output.write(outtext)
            output.close()

            print "The dataset " + dataset_name + " has been saved to your library."

            load_dataset(dataset_name)

# This subroutine loads a given dataset into active memory so that
# it can be used in the shell by other methods.  It also updates the
# acvitedatasets list to include the loaded dataset.

def load_dataset(dataset_name):
    global activedatasets
    flag = False

    if os.path.isfile(homefolder + 'Saved_Datasets/' + dataset_name + '.txt') == False:
        print "There is no dataset by the name " + dataset_name + " in your library."
    else:
        try:
            exec(dataset_name)
        except NameError:
            flag = True
        else:
            ans = ''
            while ans != 'n' and ans != 'y':
                ans = raw_input("Warning: there is already a variable by the name " + dataset_name + " in your current shell session.  Overwrite? (y/n)")

            if ans == 'y':
                flag = True

        if flag:
            rawdata = loadtxt(homefolder + 'Saved_Datasets/' + dataset_name + '.txt',skiprows = 1,dtype = str)
            varnames = list(rawdata[0])
            numbers = rawdata[1:len(rawdata)]
            numbers = numbers.T
            numbers = numbers.astype(float)

            f = open(homefolder + 'Saved_Datasets/' + dataset_name + '.txt')
            rawdata = f.read()
            infoend = rawdata.find('\n')
            info = rawdata[0:infoend]

            dataset = [dataset_name,varnames,numbers,info]
            dataset = str(dataset)

            numbers = str([numbers])

            exec('_' + dataset_name + ' = ' + dataset) in globals()
            exec(dataset_name + ' = ' + numbers + '[0]') in globals()
            exec('activedatasets = activedatasets + ["' + dataset_name + '"]') in globals()

# The set() command removes duplicates in case some datasets are
# overwritten.
            activedatasets = list(sort(list(set(activedatasets))))

            print "Dataset " + dataset_name + " has been successfully loaded into active memory."

# This subroutine lists all of the datasets that are currently active,
# meaning that they can be accessed by the user for standard operations.

def active_datasets():
    if len(activedatasets) == 0:
        print "There are currently no datasets available in active memory."
    else:
        print "Datasets currently in active memory:"
        i = 0
        while i < len(activedatasets):
            print activedatasets[i]
            i = i + 1

# This subroutine lists all of the datasets that are currently in
# your library, which means that they are available to be loaded
# into the shell memory using the load_dataset() subroutine.

def library_datasets():
    libdatasets = os.listdir(homefolder + 'Saved_Datasets/')
    if len(libdatasets) == 0:
        print "There are currently no datasets in your library."
    else:
        libdatasets = list(sort(libdatasets))
        print "Datasets currently in your library:"
        i = 0
        while i < len(libdatasets):
            datname = libdatasets[i].replace('.txt','')
            print datname
            i = i + 1

# This subroutine loads all datasets in your library to the active
# shell.  In the absence of any library datasets, it returns an
# error message.

def load_datasets():
    libdatasets = os.listdir(homefolder + 'Saved_Datasets/')
    if len(libdatasets) == 0:
        print "There are no datasets saved in your library."
    else:
        for dataset in libdatasets:
            datname = dataset.replace('.txt','')
            load_dataset(datname)

# This subroutine takes a function and fits its parameters to a
# given dataset using scipy.optimize's curve_fit.  The success of
# this subroutine is largely contingent upon the user's correct
# choice of function and parameters, as well as on a sufficient
# number of data points being present in the dataset itself.
# Since RNAeye graphs the dataset and its fitted function together,
# the user can apply an intuitive judgement as to the correctness
# of the fit.  A least-squares method is used here, so the fit is
# qualitative and not necessarily consistent with the fully correct
# Bayesian procedure.  best_fit is also designed to automatically
# detect data columns that have the same names as the input variables
# to the function, and match them accordingly.  If it finds no
# dataset columns with the same name as a variable, it will
# prompt the user to choose from among the available columns, so
# there is no need for the user to input his choices a priori.
# best_fit also allows the user the option to input a column of
# standard deviations that can be used to assign weights to the
# dependent variable data according to the uncertainty inherent
# in each data point.

def best_fit(function,dataset,logplot = False):
    model_resolution = 1000.
    flag = False
    
    try:
        exec('func = _' + function)
        exec('fullfunc = ' + function)
    except NameError:
        print "Warning: The function " + function + " is not currently in active memory."
    else:
        flag = True

    try:
        exec('data = _' + dataset)
    except NameError:
        print "Warning : The dataset " + dataset + " is not currently in active memory."
    else:
        flag = True
    
    if len(func[2]) >= len(data[1]):
        print "Warning: there are too few data series in dataset " + data[0] + " to fit the function " + func[0] + "."
    elif flag:
        funcflags = []
        dataflags = [False]*len(data[1])
        stackflag = True
        
        for variable in func[2]:
            j = 0
            funcflag = False
            while j < len(data[1]):
                if variable == data[1][j]:
                    if stackflag:
                        vararray = data[2][j]
                        stackflag = False
                    else:
                        vararray = vstack((vararray,data[2][j]))
                        
                    funcflag = True
                    dataflags[j] = True
                j = j + 1

            funcflags = funcflags + [funcflag]
            
            if funcflag == False:
                if stackflag:
                    vararray = zeros(len(data[2][0]))
                    stackflag = False
                else:
                    vararray = vstack((vararray,zeros(len(data[2][0]))))

        while False in funcflags:
            missingfuncs = where(array(funcflags) == False)[0]
            okanswers = list(array(data[1])[where(array(dataflags) == False)[0]])

            ans = ''
            while ans not in okanswers:
                print "The variable " + func[2][missingfuncs[0]] + " does not have a name that corresponds to a column in the dataset " + data[0] + ".  Select one of the following columns to assign to the variable " + func[2][missingfuncs[0]] + "."
                for answer in okanswers:
                    print answer
                ans = raw_input("Please select one of the above:")

            selectedindex = where(array(data[1]) == ans)[0][0]
            
            if len(shape(vararray)) == 1:
                vararray = data[2][selectedindex]
            else:
                vararray[missingfuncs[0]] = data[2][selectedindex]
                
            funcflags[missingfuncs[0]] = True
            dataflags[selectedindex] = True

        okanswers = list(array(data[1])[where(array(dataflags) == False)[0]])

        if len(okanswers) == 1:
            selectedindex = where(array(data[1]) == okanswers[0])[0][0]
        else:
            ans = ''
            while ans not in okanswers:
                print "Please assign, from among the following choices, the column of datset " + data[0] + " that is to represent the dependent variable of function " + function[0] + ":"
                for answer in okanswers:
                    print answer
                ans = raw_input("Please select one of the above:")

            selectedindex = where(array(data[1]) == ans)[0][0]

        ydata = data[2][selectedindex]
        dataflags[selectedindex] = True
        yname = data[1][selectedindex]

        if len(data[1]) >= len(function[2]) + 2:
            ans = ''
            while ans != 'y' and ans != 'n':
                ans = raw_input("Would you like to weight the dependent variable " + data[1][selectedindex] + " according to its standard deviations? (y/n)")

            if ans == 'y':
                okanswers = list(array(data[1])[where(array(dataflags) == False)[0]])

                ans = ''
                while ans not in okanswers and ans != 'none':
                    print "Please assign, from among the following choices, the column of dataset " + data[0] + " that is to represent the standard deviation of dependent variable " + data[1][selectedindex] + ":"
                    for answer in okanswers:
                        print answer
                    ans = raw_input("Please select one of the above, or type 'none' to cancel standard deviation weighting:")

                if ans == 'none':
                    sigma = None
                else:
                    selectedindex = where(array(data[1]) == ans)[0][0]
                    sigma = data[2][selectedindex]
                    dataflags[selectedindex] = True
            else:
                sigma = None
        else:
            sigma = None

        try:
            warnings.simplefilter('ignore',RuntimeWarning)
            estimates = curve_fit(fullfunc,vararray,ydata,sigma = sigma)[0]
        except RuntimeError:
            print "Warning: the data in dataset " + dataset + " may be insufficient to reliably determine the optimal fitting parameters for the function " + function + "."
        else:
            if len(shape(vararray)) == 1:
                interval = (max(vararray) - min(vararray))/model_resolution
                start = min(vararray)
                end = max(vararray) + interval
                varlines = arange(start,end,interval)
            else:
                interval = (max(vararray[0]) - min(vararray[0]))/model_resolution
                start = min(vararray[0])
                end = max(vararray[0]) + interval
                varlines = arange(start,end,interval)

                i = 1
                while i < len(vararray):
                    interval = (max(vararray[i]) - min(vararray[i]))/model_resolution
                    start = min(vararray[i])
                    end = max(vararray[i]) + interval
                    varlines = vstack((varlines,arange(start,end,interval)))
                    i = i + 1

            callfunc = 'ylines = fullfunc(varlines'
            i = 0
            while i < len(func[3]):
                callfunc = callfunc + ',estimates[' + str(i) + ']'
                i = i + 1

            callfunc = callfunc + ')'
            exec(callfunc)

            callfuncrsq = 'yforrsq = fullfunc(vararray'
            i = 0
            while i < len(func[3]):
                callfuncrsq = callfuncrsq + ',estimates[' + str(i) + ']'
                i = i + 1

            callfuncrsq = callfuncrsq + ')'
            exec(callfuncrsq)
                
            print "RNAeye's best estimated parameter fits:"
            i = 0
            while i < len(func[3]):
                print func[3][i] + " = " + str(estimates[i])
                i = i + 1

            print ""

            datamean = sum(ydata)/len(ydata)
            SStot = sum((ydata - datamean)**2)
            modelmean = sum(yforrsq)/len(yforrsq)
            SSres = sum((ydata - yforrsq)**2)

            R2 = 1. - (SSres/SStot)

            print "R squared value of the fit: " + str(R2) + "."

            if logplot:
                    varlines = log(varlines)/log(10)
                    vararray = log(vararray)/log(10)
            
            if len(shape(vararray)) == 1:
                figure(1)
                xlabel(func[2][0])
                ylabel(yname)
                plot(varlines,ylines,'b-')
                errorbar(vararray,ydata,yerr = sigma,fmt = 'g.',markersize = 15)
                ylim(0,max(ylines)*1.1)
                show()
            else:
                i = 0
                while i < len(func[2]):
                    figure(i)
                    xlabel(func[2][i])
                    ylabel(yname)
                    plot(varlines[i],ylines,'b-',vararray[i],ydata,'g.',markersize = 15)

                    i = i + 1

                ylim(0,max(ylines)*1.1)
                show()

def checkprimer(candidate,auto = True):
    flag = True
                
    if threeprimetest(candidate) == False:
        if not auto:
            print "3' test failed."
        flag = False
    if checkconsecutive(candidate):
        if not auto:
            print "Consecutive test failed."
        flag = False
    GC = countGC(candidate)
    if GC < 0.4 or GC > 0.6:
        if not auto:
            print "GC test failed."
        flag = False

    return flag

def findprimers(sequence,primer_length):
    targets = maketargets(sequence,primer_length)
    
    survivors = []
    startlocations = []
    endlocations = []
    i = 0
    while i < len(targets):
        if checkprimer(targets[i]):
            survivors = survivors + [targets[i]]
            startlocations = startlocations + [i + 1]
            endlocations = endlocations + [i + primer_length + 1]
        i = i + 1

    print "Results:"
    
    i = 0
    while i < len(survivors):
        print "----------- Candidate " + str(i + 1) + " -----------"
        print "Bases " + str(startlocations[i]) + " to " + str(endlocations[i]) + "."
        print survivors[i]
        print ""

        i = i + 1
    

