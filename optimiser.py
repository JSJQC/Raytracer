# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 11:19:11 2020

@author: jos18

A module that contains code for optimising a singlet lens for a given 
focal length. If run as __main__, contains code to plot various plots 
relating to this optimisation

"""

import System as sy

import scipy as sp
import matplotlib.pyplot as plt


n = 1.5168

def adjustor3(lens1, focallength, lensn, d):
    c1 = lens1.getcurve()
    
    # This is a rearrangement of the lensmaker's eqn
    top = (1 / (focallength * (lensn - 1))) - c1
    bottom = ((((lensn - 1) * d) / lensn) * c1 ) - 1
    
    c2 = top / bottom
    
    lens2 = sy.createLens(d, c2, lensn, 1, 5) # creates right face of lens
    
    return lens2


def lenslist(focallength, lensn, d):
    # Range of curvatures to be tested
    curves = sp.arange(0.01, 0.151, 0.01)
    
    lens1list = []
    lens2list = []
    
    for curve in curves:
        lens1 = sy.createLens(0, curve, 1, lensn, 5)
        # Creates list of left faces of varying curvature
        lens1list.append(lens1)
        lens2list.append(adjustor3(lens1, focallength, lensn, d))
        
    return lens1list, lens2list


def spherAb(lenses):
    # lenses is a list of lenses CONTAINING THE OUTPUT PLANE
    beam = sy.createCylBeam(20, 0.2) # Beam passed through test lenses
    sy.propagateSystem(lenses, beam)
    aberration = lenses[-1].rms(beam)
    # Gives RMS at output plane for each test lens
    return aberration


def selectBest(aberrations):
    # Finds index of best lens faces
    lowestAb = min(aberrations)
    
    bestIndex = aberrations.index(lowestAb)
    
    return bestIndex


def optimum(focallength, lensn, d):
    # Returns the best lens configuration
    # focallength is from the edge of the lens
    lenslist1, lenslist2 = lenslist(focallength, lensn, d)
    output = sy.createOutput(focallength + (d/2))
    
    aberrations = []
    
    for x in range(len(lenslist1)):
        lenses = [lenslist1[x], lenslist2[x], output]
        aberrations.append(spherAb(lenses))
        
    bestIndex = selectBest(aberrations)
    print (bestIndex)
    
    return lenslist1[bestIndex], lenslist2[bestIndex]


if __name__ == "__main__":
    focal = 35.487
    d = 10
    
    optimum(focal, n, d)
    
    lens1, lens2 = lenslist(focal, n, d)
    output = sy.createOutput(40.487)
    
    lens1curves = []
    lens2curves = []
    aberrations = []
    
    for x in range(len(lens1)):
        lenses = [lens1[x], lens2[x], output]
        lens1curves.append(lens1[x].getcurve())
        lens2curves.append(lens2[x].getcurve())
        
        aberrations.append(spherAb(lenses))
    
    # Graphs lens 1 curvatures only -----------------
    
    plt.figure('Aberrations')
    plt.xlabel('Lens 1 Curvature (mm)')
    plt.ylabel('RMS at focal point (mm)')
    plt.minorticks_on
    plt.grid()
    plt.plot(lens1curves, aberrations)
     
    # ------------------------------------------------
    
    # 3D scatter plot used in final report -----------
    
    sy.threeDscatter(lens1curves, lens2curves, aberrations)
    
    # ----------------------------------------------

    # 2D plot that uses colour instead of z-axis ------------------------------
        
    fig, ax = plt.subplots()
    scat = ax.scatter(lens1curves, lens2curves, c = aberrations, \
                      s=120, marker='.')
    fig.colorbar(scat)
    plt.xlabel('Lens 1 Curvature (mm-1)')
    plt.ylabel('Lens 2 Curvature (mm-1)')
    plt.minorticks_on
    plt.grid()

    plt.show()
    
    # Wasn't clear, so wasn't used --------------------------------------------
