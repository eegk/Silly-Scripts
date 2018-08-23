#!/usr/bin/python -O
import numpy as np
from numpy import array

A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19

int_to_char = {0:'A',1:'R',2:'N',3:'D',4:'C',5:'Q',6:'E',7:'G',8:'H',9:'I',
10:'L',11:'K',12:'M',13:'F',14:'P',15:'S',16:'T',17:'W',18:'Y',19:'V'}

indel = -5
#blosum62
scoring = array([[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
[-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
[-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
[-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
[0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
[-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
[-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
[0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
[-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
[-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
[-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
[-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],
[-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],
[-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],
[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],
[1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
[0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],
[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],
[-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],
[0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4]])

class AlignmentFinder(object):
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.D = None

    def find_gobal_alignment(self):
        self.D = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int16)
        self._compute_array()
        print self.D
        return self._traceback()

    def _compute_array(self):
        for i in xrange(self.seq1.size+1):
            self.D[i,0] = i*indel
        for j in xrange(self.seq2.size+1):
            self.D[0,j] = j*indel
        for i in xrange(1, self.seq1.size+1):
            for j in xrange(1, self.seq2.size+1):
                self.D[i,j] = max(  self.D[i-1, j-1] + self._get_score(i, j),
                                    self.D[i-1, j] + indel,
                                    self.D[i, j-1] + indel)
    def _get_score(self, i, j):
        ''' The indexing is quite tricky because the matrix as one more row & column.
        That causes a shift between the matrix index and the sequence indices.
        Therefore, to obtain the correct nucleotide in the sequence, we must
        substract 1 to the matrix index. '''
        return scoring[self.seq1[i-1], self.seq2[j-1]]
    
    def _get_aligned_pair(self, i, j):
        n1 = int_to_char[self.seq1[i-1]] if i>0 else '_'
        n2 = int_to_char[self.seq2[j-1]] if j>0 else '_'
        return (n1, n2)

    def _traceback(self):
        alignment= []
        i = self.seq1.size
        j = self.seq2.size
        while i >0 and j>0:
            if self.D[i-1, j-1] + self._get_score(i, j) == self.D[i,j]:
                alignment.append(self._get_aligned_pair(i, j))
                i -= 1
                j -= 1
            elif self.D[i-1, j] + indel == self.D[i,j]:
                alignment.append(self._get_aligned_pair(i, 0))
                i -= 1
            else:
                alignment.append(self._get_aligned_pair(0, j))
                j -= 1
        while i > 0:
            alignment.append(self._get_aligned_pair(i, 0))
            i -= 1
        while j > 0:
            alignment.append(self._get_aligned_pair(0, j))
            j -= 1
        alignment.reverse()
        return alignment

def print_sequences(pairs):
    top_seq = []
    bottom_seq = []
    for (b, t) in pairs:
        bottom_seq.append(b)
        top_seq.append(t)
    for n in top_seq:
        print n,
    print ' '
    for n in bottom_seq:
        print n,

if __name__ == "__main__":
    s1 = array([V,V,V,R,T,Y,C,S,S,A,A], dtype=np.int16)
    s2 = array([V,V,R,T,S,F], dtype=np.int16)
    print s1
    print s2
    aligner = AlignmentFinder(s1, s2)
    pairs = aligner.find_gobal_alignment()
    score = aligner._get_score
    print_sequences(pairs)
    print score

