# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 16:43:12 2020

@author: jos18

A module that carries out simulation pertinent to the lab report
If run as __main__ will also produce figures used in said report

"""
import System as sy
import optimiser as op

import scipy as sp
import matplotlib.pyplot as plt

    
def task15(lenses, beamarray):
    # takes a list of lenses, and a list of beams
    # output must be at index -1
    results = []
    
    for beam in beamarray:
        sy.propagateSystem(lenses, beam)
        results.append(lenses[-1].rms(beam))
        # Finds RMS at output for each beam

    return results


def rmsresults(maxraynumber, resolution, lenses):
    # carries out task 15 for a range of beam diameters
    raynumbers = sp.arange(1, maxraynumber + 1, 1)
    beams = []

    for number in raynumbers:
        beams.append(sy.createCylBeam(number, resolution, 0, 0, -20))
        
    result = task15(lenses, beams)
    
    return raynumbers, result

def runexperiment(maxraynumber, lenses):
    # Gives plottable results
    res = 0.2 # controls beam density - steps in radius, essentially
    
    results = rmsresults(maxraynumber, res, lenses)
    
    diameters = results[0] * res
    
    return diameters, results[1]
    

def spherabFig(lenses, rays):
    # Used to produces a figure showing spherical aberrations
    sy.propagateSystem(lenses, rays)
    sy.plotyz(rays, lenses, 'Spherical Aberrations of Planoconvex Lens')
    

if __name__ == '__main__':
    
    ## Focal distance values ---------------------------------------------
    # All in mm, measured from z = 0
    
    #97.8305 : convexleft
    #67.895 : planoconvex clpr : (| shape
    #74.497 : planoconvex plcr : |) shape
    #40.487 : biconvex
    
    # ---------------------------------------------------------
    
    # Lens parts used to produces RMS-beam diameter graph ----------
    
    convexleft = sy.createLens(0, 0.03, 1, 1.5168, 33.3)
    planeright = sy.createLens(10, 0, 1.5168, 1, 33.3)
    
    planeleft = sy.createLens(0, 0, 1, 1.5168, 33.3)
    convexright = sy.createLens(10, -0.03, 1.5168, 1, 33.3)
    
    lensface1, lensface2 = op.optimum(35.487, 1.5168, 10)
    
    lenses2 = [lensface1, lensface2, sy.createOutput(40.487)]
    
    # ---------------------------------------------------
    
    # Produces RMS-beam diameter data for four lenses
    
    cl = runexperiment(20, [convexleft, sy.createOutput(97.8305)])
    clpr = runexperiment(20, [convexleft, planeright, sy.createOutput(67.895)])
    bi = runexperiment(20, [convexleft, convexright, sy.createOutput(40.487)])
    bestbi = runexperiment(20, lenses2)
    
    # Plots RMS-beam diameter data ------------
    # Figures 6 and 8 in the report
    
    plt.figure('RMS vs diameter')
    plt.title('RMS at focal point vs beam diameter', fontsize = 25)
    plt.xticks(fontsize = 25)
    plt.yticks(fontsize = 25)
    plt.xlabel('Beam Diameter (mm)', fontsize = 25)
    plt.ylabel('RMS at Focal Point (mm)', fontsize = 25)
    plt.minorticks_on
    plt.grid(which = 'both')
    plt.plot(cl[0], cl[1], label = 'Single Surface')
    plt.plot(clpr[0], clpr[1], label = 'Planoconvex Lens')
    plt.plot(bi[0], bi[1], label = 'Biconvex Lens')
    plt.plot(bestbi[0], bestbi[1], label = 'Optimised', linewidth = 3.5)
    plt.legend(fontsize = 30)
    plt.show()
    
    # -------------------------------------------
    
    # Creates figure 5 of the report -----------------
    
    beam = sy.createCylBeam(10, 0.1)
    lens = sy.createLens(0, 0.03, 1, 1.5168, 33)
    lenses = [lens, sy.createOutput(97.835)]
    sy.propagateSystem(lenses, beam)
    sy.plotyz(beam, lenses, 'plot')
    sy.plotspotxy(beam, 'Spot Plot in the XY Plane at z=-20mm', -20)
    sy.plotspotxy(beam, 'Spot Plot in the XY Plane at the Output Plane', \
                  97.835)
    
    # -------------------------------------------------------------------------
    

