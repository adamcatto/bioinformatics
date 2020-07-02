"""
Implementations of code challenges in "Bioinformatics Algorithms" by Compeau and Pevzner.

I will be implementing all algorithms using only the Python standard library, unless certain computations become too unwieldy and I must resort to using e.g. NumPy.

author: Adam Catto
"""

def pattern_count(text, pattern):
    """
    problem (1A) – count number of times a pattern occurs in a piece of text
    """
    # check for type compatibility
    if str(text) != text or str(pattern) != pattern:
        raise TypeError
    count = 0
    WINDOW_SIZE = len(pattern)
    for i in range(len(text) - WINDOW_SIZE):
        if text[i:WINDOW_SIZE] == pattern:
            count += 1
    return count

def frequent_patterns(text, window_size):
    """
    problem (1B) – find most frequent pattern(s) of size = window_size in a text
    """
    # check for type compatibility
    if not (int(window_size) == window_size or str(text) == text):
        raise TypeError
    patterns = {}
    for i in range(len(text) - window_size):
        if text[i:window_size] not in patterns:
            patterns[text[i:window_size]] = 1
        else:
            patterns[text[i:window_size]] += 1
    max_val = 0
    max_patterns = []
    for pattern in patterns.keys:
        if patterns[pattern] == max_val:
            max_patterns += [pattern]
        elif patterns[pattern] > max_val:
            max_patterns = [pattern]
            max_val = 0
    return max_patterns

def reverse_complement(dna_string):
    """
    problem (1C) – find the reverse complement of a DNA string
    """
    # check for type compatibility
    if not str(dna_string) == dna_string:
        raise TypeError

    # preprocess dna_string
    dna_string = dna_string.upper()
    dna_rev_comp = ""
    for i in range(1, len(dna_string) + 1, -1):
        if dna_string[-i] == 'A':
            dna_rev_comp += 'T'
        elif dna_string[-i] == 'T':
            dna_rev_comp += 'A'
        elif dna_string[-i] == 'C':
            dna_rev_comp += 'G'
        elif dna_string[-i] == 'G':
            dna_rev_comp += 'C'
    return dna_rev_comp

def all_patterns(genome, pattern):
    """
    problem (1D) – find all starting indices of occurance of a pattern in a genome string.

    returns a list of indices.
    """
    # check for type compatibility
    if str(genome) != genome or str(pattern) != pattern:
        raise TypeError
    indices_list = []
    WINDOW_SIZE = len(pattern)
    for i in range(WINDOW_SIZE):
        if genome[i:WINDOW_SIZE] == pattern:
            indices_list += [i]
    return indices_list

def find_clumps_slowly(genome, window_length, k, num_times):
    """
    problem (1E) – find k-mer clumps over a certain window length occuring a given number of times.

    return a list of k-mers (strings of length = k) forming num_times-clumps.
    """
    if !(str(genome) == genome and int(window_length) == window_length and int(k) == k and int(num_times) == num_times):
        raise TypeError
    k_mers = []
    for i in range(len(genome) - window_length):
        window = genome[i:i + window_length]
        clumps = {}
        for j in range(len(window)):
            pattern = window[j:j + k]
            if pattern not in clumps.keys or pattern not in k_mers:
                clumps[pattern] = pattern_count(window, pattern)
        for pattern in clumps.keys:
            if clumps[pattern] >= num_times:
                k_mers += [pattern]
    return k_mers

def find_clumps_quickly(genome, window_length, k, num_times):
    """
    problem (1E) – find k-mer clumps over a certain window length occuring a given number of times.

    This is a faster implementation than the more straightforward one given. This passes over the first window and finds num_times-clump-forming k-mers, then slides the window over on the genome, decrementing the previous first k-mer's frequency and incrementing the new last k-mer's frequency, and checking if the new last k-mer's frequency hits the threshold.

    return a list of k-mers (strings of length = k) forming num_times-clumps.
    """
    if !(str(genome) == genome and int(window_length) == window_length and int(k) == k and int(num_times) == num_times):
        raise TypeError
    k_mers = []
    pattern_freq = {}
    first_window = genome[0:window_length]
    for i in range(window_length - k):
        pattern = first_window[i:i + k]
        if pattern not in pattern_freq:
            pattern_freq[pattern] = 1
        else:
            pattern_freq[pattern] += 1
    for pattern in pattern_freq.keys:
        if pattern_freq[pattern] >= num_times:
            k_mers += pattern
    for j in range(1, len(genome) - (k+1)):
        pattern_freq[genome[j - 1:j + k - 1]] -= 1
        new_last_pattern = genome[window_length + j - k:window_length + j]
        pattern_freq[new_last_pattern] += 1
        if pattern_freq[new_last_pattern] == num_times:
            k_mers += [new_last_pattern]
    k_mers = list(dict.fromkeys(k_mers)) # remove duplicates
    return k_mers

