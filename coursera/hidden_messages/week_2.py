# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 2021

@author: Alexander Trouern-Trend

Coursera Bioinformatics Specialization: Finding Hidden Messages - Week 2

"""
from itertools import product

def formatter(positions):
    out = ""
    for i in positions:
        out += (str(i) + " ")
    out = out[:-1]
    return out  


def gc_skew(string):
    '''
    Tracks the G-C skew while reading down a sequence.

    Parameters
    ----------
    string : str
        DNA sequence.

    Returns
    -------
    skew : list
        Array of skew values, one per nucleotide base.

    '''
    skew = [0]
    count = 0
    for i in string:
        if i == "G":
            count += 1
        if i == "C":
            count += -1
        skew.append(count)
    return skew


def min_skew(seq):
    '''
    Finds global minimum of skew array for given sequence. Calls gc_skew().

    Parameters
    ----------
    seq : str
        Input nucleotide sequence.

    Returns
    -------
    low : int
        skew minimum value.
    minima : list
        indices where skew = low.

    '''
    skew = gc_skew(seq)
    low = min(skew)
    minima = []
    for e, i in enumerate(skew):
        if i == low:
            minima.append(e)
    return low, minima


def hamming_distance(seq1, seq2):
    '''
    Calculate Hamming distance between two strings of equivalent length.

    Parameters
    ----------
    seq1 : str
        first sequence.
    seq2 : str
        second sequence.

    Returns
    -------
    dist : int
        Hamming distance.

    '''
    dist = 0
    if len(seq1) != len(seq2):
        return('strings have unequal length!')
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist


# Approximate pattern match to account for wobble between high-identity K-mers

def approximte_match(pattern, text, d):
    '''
    Finds matches within text which are 100% identical to pattern, or below the
    maximum Hamming distance set by d.

    Parameters
    ----------
    pattern : str
        Query sequence.
    text : str
        Target sequence.
    d : int
        Maximum acceptable hamming distance between pattern and match.

    Returns
    -------
    positions : list
        Indices where matches begin.

    '''
    positions = []
    pattern = pattern.strip()
    text = text.strip()
    space = len(pattern)
    end = len(text)
    for i in range(end-space+1):
        if hamming_distance(text[i:i+space], pattern) <= int(d.strip()):
            positions.append(i)
    return positions


# Modify approximte_match() to include a counter.

def count_d(pattern, text, d):
    '''
    See description for approximte_match(). Same inputs with one extra return,
    a count of matches.

    '''
    positions = []
    pattern = pattern.strip()
    text = text.strip()
    space = len(pattern)
    end = len(text)
    for i in range(end-space+1):
        if hamming_distance(text[i:i+space], pattern) <= int(d):
            positions.append(i)
    return positions, [len(positions)]


# Frequent Words with Mismatches problem

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


def ImmediateNeighbors(pattern):
    '''
    Returns set of input sequence of length k and all possible k-mers of Hamming
    distance = 1 (neighbors).

    Parameters
    ----------
    pattern : str
        Input nucleotide sequence.

    Returns
    -------
    neighborhood : set
        Input sequence and all neighbors.

    '''
    neighborhood = set()
    neighborhood.add(pattern)
    for i in range(len(pattern)):
        symbol = pattern[i]
        for x in "ATCG".replace(symbol, ""):
            neighbor = pattern[:i] + x + pattern[i+1:]
            neighborhood.add(neighbor)
    return neighborhood

# All possible nt combinations of length k


def arrangements(k):
    '''
    Returns all possible nucleotide sequences of length k.

    '''
    result = product('ACGT', repeat = k)
    result = map(list, result)
    words = []
    for item in result:
        word = ''
        for letter in item:
            word += letter
        words.append(word)
    return(words)

def neighbors(pattern, d):
    '''
    Find all neighbors of given pattern less than or equal to given hamming
    distancem, d.

    Parameters
    ----------
    pattern : str
        Input sequence.
    d : int
        Hamming distance threshold.

    Returns
    -------
    neighborhood: list
        Array of strings including input sequence and all neighbors within
        specified hamming distance.

    '''
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return ["A", "T", "C", "G"]
    neighborhood = set()
    # Obtain suffix
    suffix = pattern[1:]
    prefix = pattern[0]
        
    # recursive algorithm (repeat function for each suffix)
    suffixNeighbors = neighbors(suffix, d)
    for i in suffixNeighbors:
        if hamming_distance(suffix, i) < d:
            for x in "ACGT":
                neighborhood.add(x+i)
        else:
            neighborhood.add(prefix+i)
    return neighborhood


def FrequentWordsWithMismatches(text, k, d):
    '''
    A frequent k-mer finder which accounts for mismatches less than or equal to
    Hamming distance, d. A single-strand search function.

    Parameters
    ----------
    text : Sequence search space
        Input genome or genome fragment.
    k : k-mer length
        Length of k-mers to analyze.
    d : maximum hamming distance
        For each window, count will increase for all sequences within this
        dissimilarity threshold.

    Returns
    -------
    patterns : list
        Most frequent patterns of length k.

    '''
    patterns = []
    freqMap = {}
    n = len(text)
    for i in range(n-(k-1)):
        pattern = text[i:i+k]
        neighborhood = neighbors(pattern, d)
        for j in neighborhood:
            if j not in freqMap.keys():
                freqMap[j] = 1
            else:
                freqMap[j] = freqMap[j] + 1
    m = max(freqMap.values())
    for i in freqMap.keys():
        if freqMap[i] == m:
            patterns.append(i)
    return patterns


# Include reverse Complements

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

def fwwm_revcom(text, k, d):
    '''
    A frequent k-mer finder which accounts for mismatches less than or equal to
    Hamming distance, d and includes search for reverse complement sequences
    (double stranded search).

    Parameters
    ----------
    text : Sequence search space
        Input genome or genome fragment.
    k : k-mer length
        Length of k-mers to analyze.
    d : maximum hamming distance
        For each window, count will increase for all sequences within this
        dissimilarity threshold.

    Returns
    -------
    patterns : list
        Most frequent patterns of length k.

    '''
    patterns = []
    freqMap = {}
    n = len(text)
    for i in range(n-(k-1)):
        pattern = text[i:i+k]
        neighborhood = list(neighbors(pattern, d))
        revcoms = []
        for t in neighborhood:
            rev = reverse_complement(t)
            revcoms.append(rev)
        neighborhood = neighborhood + revcoms
        for j in neighborhood:
            if j not in freqMap.keys():
                freqMap[j] = 1
            else:
                freqMap[j] = freqMap[j] + 1
    m = max(freqMap.values())
    for i in freqMap.keys():
        if freqMap[i] == m:
            patterns.append(i)
    return patterns
