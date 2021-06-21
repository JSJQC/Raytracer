# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:24:43 2020

@author: jos18

A module that handles all ray properties and behaviours
"""

import scipy as sp

class Ray:
    
    """
    A class for a single ray of light
    Contains all methods that directly act on and/or modify the ray
    Warning: Attributes are hidden, and accessing them outside this class
    should not be attempted
    """
    
    
    def __init__(self, position = [0,0,-20], direction = [0,0,0]):
        
        self.__p = sp.array(position) # These correspond to current parameters
        self.__k = sp.array(direction)
        
        self.__allp = [sp.array(position)] 
        # These create a running record of past positions and directions
        self.__allk = [sp.array(direction)] 
        
    def __repr__(self):
        return "Ray([{0}], [{1}])".format(self.__p, self.__k)
    
    def __str__(self):
        return "A ray at position {0} with direction {1}".format( \
                                  self.__p, self.__k)
        
    def getp(self):
        return self.__p
    
    def getk(self):
        return self.__k
    
    def append(self, newp = [0,0,-20], newk = None):
        
        if newk is None:
            newk = self.getk()
            newk.tolist() # The default direction is the current direction
        
        self.__allp.append(newp) # Adds new direction and position to record
        self.__allk.append(newk) # immediately
        
        self.__p = sp.array(newp)
        self.__k = sp.array(newk)
        
        return self
    
    def vertices(self):
        return self.__allp
    
    def checkallk(self):
        return self.__allk
    
    def terminate(self, newp):
        # Acts to terminate a ray by making it 'stationary'
        self.append(newp, [0,0,0])  
        
        return self
