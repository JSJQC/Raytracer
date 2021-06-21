# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 16:44:14 2020

@author: jos18

A module that handles the behaviour of a system of lens and rays
This includes beam creation, passing of rays through several lens,
and plotting of ray paths and positions

"""

import raytracer as ryt
import opticalequipment as oe

import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def createLens(z0, curvature, n1, n2, apertureRadius):
    # Shields optical equipment from use by later modules
    lens = oe.SphericalRefraction(z0, curvature, n1, n2, apertureRadius)
    return lens

def createOutput(z0):
    out = oe.OutputPlane(z0)
    return out

def createBundle(number, spacing, point = True, xpos = 0, \
                 yangle = 0, zstart = -20):
    
    if yangle >= sp.pi / 2:
        print ()
        print ('Angle too great')
        print ()
        raise Exception
    
    rayarray = []
    
    for y in range(number):
        
        if point: # A point source, radiating out
            ydirection = spacing * y
        
            ray = ryt.Ray([xpos,0,zstart], [0,ydirection,100])
            
        elif not point: # collimated flat 'beam' - y-axis only
            ydirection = 100 * sp.tan(yangle)
            ypos = spacing * y
            ystart = -1 * abs(zstart) * sp.tan(yangle)
            
            ray = ryt.Ray([xpos, ypos + ystart, zstart], [0, ydirection, 100])
            
        rayarray.append(ray)
        
        
        if y != 0: # PRoduces a ray in the negative position/direction
            if point:
                ray = ryt.Ray([xpos,0,zstart], [0,-1 * ydirection,100])
                
            elif not point:
                ystart = -1 * abs(zstart) * sp.tan(yangle)
                ray = ryt.Ray([xpos, (-1 * ypos) + ystart, zstart], \
                               [0, ydirection, 100])           
            
            rayarray.append(ray) # Holds all rays
            
    return rayarray



def createRecBeam(number, spacing, yangle):
    # Beam that is square in cross-section
    totray = []
    for num in range(number):
        if num != 0:
            # Puts multiple flat bundles side-by-side
            line1 = createBundle(number, spacing, False, num * spacing, yangle)
            line2 = createBundle(number, spacing, False, \
                                 -1 * num * spacing, yangle)
            totray += line1 + line2
        elif num == 0:
            line = createBundle(number, spacing, False, num, yangle)
            totray += line
        
    return totray
    
        
def createCylBeam(number, resolution = 0.1, xangle = 0, \
                  yangle = 0, zstart = -20): 
    # Familiar cylindrical beam
    
    rayarray = []
    arclength = (sp.pi / 4) * resolution
    # Sets distance between points on the circle
    
    zdirection = 100
    xdirection = sp.tan(xangle) * zdirection
    ydirection = sp.tan(yangle) * zdirection # Handles off-axis angles
    
    xstart = -1 * abs(zstart) * sp.tan(xangle) # Adjusts start point to
    ystart = -1 * abs(zstart) * sp.tan(yangle) # account for angle change

    for num in range(number):
        # Need to do r = 0 special case
        if num == 0:
            
            z = zstart # Single, central ray
            x = xstart
            y = ystart 
        
            ray = ryt.Ray([x, y, z], [xdirection, ydirection, zdirection])       
            
            rayarray.append(ray)
        
        else:
            
            r = num * resolution
            
            points = int((2 * sp.pi) / (arclength / r)) 
            # number of points changes with r to give constant density
            
            for point in range(points):
                # Goes round 2 pi creating rays
                theta = (arclength / r) * point
                
                z = zstart
                x = r * sp.cos(theta) + xstart
                y = r * sp.sin(theta) + ystart
                
                ray = ryt.Ray([x, y, z], [xdirection, ydirection, zdirection])
            
                rayarray.append(ray)
    
    return rayarray

def propagateSystem(lenses, rays):
    # LENSES MUST BE LIST OF LENS OBJECTS -- followed
    
    for lens in lenses:
        lens.encounterlens(rays)
        
    return rays


def plotyz(rayarray, lenses, title):
    # plots rays' paths in yz plane
    plottingarray = []
    
    for ray in rayarray: # rayarray is a beam or bundle of rays
        smallarray = []
        for array in ray.vertices(): # Goes through list of past positions
            smallarray.append([array[2], array[1]])
        plottingarray.append(smallarray)
        
    plt.figure()
    plt.title(title, fontsize = 30)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.xlabel('z (mm)', fontsize = 20)
    plt.ylabel('y (mm)', fontsize = 20)
    plt.minorticks_on
    plt.grid()
    
    for group in plottingarray:
        
        z = []
        y = []
        for points in group:
            z.append(points[0])
            y.append(points[1])
        plt.plot(z, y, 'b')
        
    for lens in lenses:
        toplot = lens.plotlensyz(rayarray) # Plots lens outlines
        plt.plot(toplot[0], toplot[1], 'orange')
    
    plt.show()
        
    
def plotspotxy(rayarray, title, z0 = 0):
    # Plots ray locations in the xy plane
    plottingarray = []
    
    for ray in rayarray:
        for array in ray.vertices():
            if array[2] == z0:
                plottingarray.append([array[0], array[1]])
    
    if z0 == -20: # Corresponds to plane of instantiation
        plt.figure('2')
    else: # Some other location
        plt.figure('3')
    x = []
    y = []
    for group in plottingarray:
        
        x.append(group[0])
        y.append(group[1])
        
        
    plt.title(title, fontsize = 20)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.xlabel('x (mm)', fontsize = 20)
    plt.ylabel('y (mm)', fontsize = 20)
    plt.minorticks_on
    plt.grid()
    plt.plot(x, y, '.', color = 'b', markersize = 20)
    
    plt.show()
    

def threeDscatter(x, y, z): 
    # Only used for one task in optimiser.py
    # Fits best here with other plotting functions
    
    
    fig = plt.figure()
    
    ax = Axes3D(fig)
    
    ax.set_xlabel('Lens 1 Curvatures', fontsize = 20)
    ax.set_ylabel('Lens 2 Curvatures', fontsize = 20)
    ax.set_zlabel('RMS at Focal Point', fontsize = 20)
    ax.xaxis.labelpad = 35
    ax.yaxis.labelpad = 35
    ax.zaxis.labelpad = 35
    
    ax.tick_params(labelsize = 20, pad = 5)
    
    ax.scatter(x, y, z, c='k', marker = 'o')
    
    plt.show()
  
    
 