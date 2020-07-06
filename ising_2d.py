#!/usr/bin/env python
  
""" 
Python program: 2d Ising model movie at various temperatures

No guarantees for this.  I looked at some Ising codes for python on
the web and adapted them to this very quickly, not having written
anything in python before.  This is in python3 which, in particular,
has a different print statement from python2.

This makes a sort of quick and dirty flip-card movie by writing
configurations to the screen.

For example, in Unix, pipe the output to a file
    python ising_2d.py > foo
and then type "less foo" and press the space bar to advance the movie
frame by frame.  Probably there is a similar way to do it on a Windows
machine if you open up an MSDOS window.  The screen showing the file
should have a height of
		HEIGHT = LENGTH + 3.
"""

from scipy import * 
from pylab import *

MCS = 300 
L = 37
#T = 2.2691
T = 1.0 
Tc = 2.2691853142130216092

def Energy(spins):
    E = 0 
    for i in range(len(spins)):
        for j in range(len(spins)):
            S = spins[i,j] 
            Neighbors = spins[(i+1)%L, j] + spins[i,(j+1)%L] + spins[(i-1)%L,j] + spins[i,(j-1)%L]
            E += -Neighbors*S
    return E/2.

def Infinite_temperature(L):
    spins = zeros((L,L), dtype=int)
    for i in range(L):
        for j in range(L):
            spins [i,j] = sign(2*rand()-1)
    return spins

def Zero_temperature(L):
    spins = zeros((L,L), dtype=int)
    for i in range(L):
        for j in range(L):
            spins [i,j] = 1
    return spins

def McMove (MCS, spins, Five_cases):
    E = Energy(spins) 
    N = L*L
    for k in range(N):
        itest = int(rand()*L)
        jtest = int(rand()*L)
        (i,j) = (itest,jtest)
        S = spins[i,j]
        Neighbors = spins[(i+1)%L, j] + spins[i,(j+1)%L] + spins[(i-1)%L,j] + spins[i,(j-1)%L]
        P = Five_cases[4+S*Neighbors]
        if P>rand(): 
            spins[i,j] = -S
        E += 2*S*Neighbors
    return spins

def frame(spins):
    for i in range(L):
        for j in range(L):
            S = spins[i,j]
            if S == 1:
                print("*", end = " ") 
            else: print(" ", end = " ")
        print("|")
    return i 
if __name__ == '__main__':

    spins = Infinite_temperature(L)
#   spins = Zero_temperature(L)

    Five_cases= zeros(9, dtype=float) 
    Five_cases[4+4] = exp(-4.*2/T)
    Five_cases[4+2] = exp(-2.*2/T)
    Five_cases[4+0] = exp(0.*2/T)
    Five_cases[4-2] = exp( 2.*2/T)
    Five_cases[4-4] = exp( 4.*2/T)

    print("Temperature in units of Tc is"), T/Tc
    print("Time is"), 0
    frame(spins)
    for timestep in range(MCS):
        print("Temperature in units of Tc is", T/Tc)
        print("Time is", timestep + 1)
        McMove(MCS, spins, Five_cases)
        frame (spins) 

