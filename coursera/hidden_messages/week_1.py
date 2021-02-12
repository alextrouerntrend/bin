# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 2021

@author: Alexander Trouern-Trend

Coursera Bioinformatics Specialization: Finding Hidden Messages - Week 1
"""

def patternmatch(text, pattern):
    '''
    How many occurences of a pattern are there in a given string
    and where are they located?
    
    Parameters
    ----------
    text : str
        Nucleic acid sequence search space.
    pattern : str
        nucleotide pattern to search for in text.

    Returns
    -------
    count : int
        number of occurences for pattern.
    positions : list
        list of postions where pattern begins (0-based).

    '''
    positions = []
    count = 0 
    space = len(pattern)
    end = len(text) - 1
    for i in range(end-space):
        if text[i:i+space] == pattern:
            count += 1
            positions.append(i)
            
    return count, positions

def formatter(positions):
    '''
    Take a list of positions and convert it to a space-separated string
    for input as answer in Stepnik learning software.
    
    Parameters
    ----------
    positions : list
        List of integer positions to be con.

    Returns
    -------
    out : str
        string of pattern indices.

    '''
    out = ""
    for i in positions:
        out += (str(i) + " ")
    return out  
    

def freqTable_ds(text, k):
    '''
    Create a frequency 'table' (really just two lists with values related 
    by index) for all k-mers (overlapping) and their frequencies in a sequence.
    this version includes a search for reverse complements, so it accounts for
    both strands (ds = double stranded).

    Parameters
    ----------
    text : str
        Nucleic acid sequence search space..
    k : int
        length of k-mers.

    Returns
    -------
    kmer : list
        list of k-mers related to freq by index.
    freq : list
        list of freuencies related to kmer by index.

    '''
    kmer = []
    freq = []
    n = len(text)
    for i in range(n-(k)):
        ptt = text[i:i+(k)]
        revcom = reverse_complement(ptt)
        if ptt in kmer:
            indx = kmer.index(ptt)
            freq[indx] = freq[indx] + 1
        elif revcom in kmer:
            indx = kmer.index(revcom)
            freq[indx] = freq[indx] + 1
        else:
            kmer.append(ptt)
            freq.append(1)
    return kmer, freq


def freqTable_ss(text, k):
    '''
    Create a frequency 'table' (really just two lists with values related 
    by index) for all k-mers (overlapping) and their frequencies in a sequence.
    this version does not include a search for reverse complements, so it only 
    accounts for one strand (ss = single stranded).

    Parameters
    ----------
    text : str
        Nucleotide sequence search space..
    k : int
        length of k-mers.

    Returns
    -------
    kmer : list
        list of k-mers related to freq by index.
    freq : list
        list of freuencies related to kmer by index.

    '''
    kmer = []
    freq = []
    n = len(text)
    for i in range(n-(k)):
        ptt = text[i:i+(k)]
        if ptt in kmer:
            indx = kmer.index(ptt)
            freq[indx] = freq[indx] + 1
        else:
            kmer.append(ptt)
            freq.append(1)
    return kmer, freq


def BetterFrequentWords(text, k):
    '''
    Finds the most frequent kmers of a given length within a nucleic acid
    sequence search space. Calls freqTable_ss().
    

    Parameters
    ----------
    text : str
        Nucleic acid sequence search space.
    k : int
        lenth of kmers.

    Returns
    -------
    maxfreq: int
        highest kmer frequency.
    freqwords: list
        kmers with frequency == maxfreq.
    '''
    freqwords = []
    kmer, freq = freqTable_ss(text, k)
    maxfreq = max(freq)
    locs = []
    idx = 0
    for i in freq:
        if i == maxfreq:
            locs.append(idx)
        idx += 1
    for i in locs:
        freqwords.append(kmer[i])
    return maxfreq, list(set(freqwords))


def reverse_complement(string):
    '''
    Generate the reverse complement of a sequence.

    Parameters
    ----------
    string : str
        input sequence.

    Returns
    -------
    new : str
        reverse complement of input sequence.

    '''
    rev = string[::-1]
    new = ""
    for i in rev:
        if i == "A":
            new += "T"
        if i == "T":
            new += "A"
        if i == "G":
            new += "C"
        if i == "C":
            new += "G"
    return new


def clump_finder(window, threshold, kmer_length, genome):
    """
    Completely inefficient algorithm to locate 'clumps' of kmers in a genome
    which might indicate the Origin of Replication region.
    
    Parameters
    ----------
    window : int
        size of sliding window to use.
    threshold : int
        number of kmers that must occur within window for clump.
    kmer_length : int
        size of kmers to search for.
    genome : str
        input sequence.
    Returns
    -------
    captures : list
        kmers that occur in clumps.

    """
    captures = []
    size = len(genome)
    for i in range(size-window):
        end = i + window
        part = genome[i:end]
        kmer, freq= freqTable_norev(part, kmer_length)
        locs = []
        idx = 0
        for x in freq:
            if x >= threshold:
                locs.append(idx)
            idx += 1
        for z in locs:
            if kmer[z] not in captures:
                captures.append(kmer[z])
    return captures


# The following was written based on suggestions from user Roddie_Reventar


def kmer_index(kmer_length, genome):
    '''
    Index the (overlapping) kmers in a genome.

    Parameters
    ----------
    kmer_length : int
        length of kmer.
    genome : str
        imput sequence.

    Returns
    -------
    kmers : dict
        dictionary of kmers as keys and list of positions in genome as values.

    '''
    kmers = {}
    size = len(genome)
    for i in range(size-(kmer_length-1)):
        kmer = genome[i:i+kmer_length]
        if kmer in kmers:
            kmers[kmer].append(i)
        else:
            kmers[kmer] = [i]
    return kmers


def index_checker(window, threshold, kmers, kmer_length):
    '''
    

    Parameters
    ----------
    window : int
        size of window used to define 'clump' of repeated kmers.
    threshold : int
        number of repeated kmers in window to define 'clump'.
    kmers : dict
        Output of kmer_index() dictionary of kmers as keys and list of
        positions in genome as values.
    kmer_length : int
        length of kmer.

    Returns
    -------
    hitkmers : dict
        kmers representing clumps. Sequence as key, positions as values.
        includes all positions across genome.
    '''
    hitkmers = {}
    for kmer, loc_set in kmers.items():
        if len(loc_set) >= threshold:
            set_len = len(loc_set)
            for i in range(set_len-(threshold-1)):
                if loc_set[i+(threshold-1)] - loc_set[i] <= (window-kmer_length):
                    hitkmers[kmer] = loc_set
                    break
    return hitkmers           


def roddie_clumpfinder(window, threshold, kmer_length, genome):
    '''
    See documentation for kmer_index() and index_checker() functions.
    '''
    kmers = kmer_index(kmer_length, genome)
    hitkmers = index_checker(window, threshold, kmers, kmer_length)
    return hitkmers

    