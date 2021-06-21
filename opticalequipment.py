# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:33:19 2020

@author: jos18

A module that handles the properties and behaviours of the lenses
"""

import scipy as sp

class OpticalElement:
    
    """
    Parent class for all optical elements
    Warning: Optical object should not be instantiated if it is 
    not a child class of this
    """
    
    def __init__(self, z0):
        self.__z0 = z0
    
    def propagate_ray(self, ray):
        ##propagate a ray through the optical element
        
        raise NotImplementedError()
        
    def encounterlens(self, rayarray):
        # General handling of multiple rays (ie. a beam) 
        # - inherited by child classes
        for ray in rayarray:
            self.propagate_ray(ray)
            
        return rayarray
    
    def rms(self, rays):
        # Inherited by child classes
        total = 0
        for ray in rays:
            
            points = ray.vertices()
            
            for point in points:
                # Interested in positions only at the element's location
                if point[2] == self.__z0: 
                    
                    r2 = (point[0] ** 2) + (point[1] ** 2)
                    
                    total += r2
                    
        meansquare = total / len(rays)
        
        rms = sp.sqrt(meansquare)
        
        return rms
    
    
        
class SphericalRefraction(OpticalElement):
    
    """
    A class for a spherical refracting surface
    Also handles behaviour of a planar surface (curvature = 0)
    Warning: Attributes should not be accessed directly outside of this class
    """
    
    def __init__(self, z0 = 0, curvature = 1, n1 = 1, \
                 n2 = 1, apertureRadius = 1):
        super().__init__(z0)        
        self.__z0 = z0 # z-intercept
        
        self.__curve = curvature # 1/radius of curvature
        # -ve curvature is concave in +ve z, +ve curvature is convex for +ve z
        
        self.__n1 = n1 # refractive index on left
        
        self.__n2 = n2 # refractive index on right
        
        self.__aprad = apertureRadius 
        # maximum extent of surface from optical axis
        
        if self.__curve != 0:
            if abs((1 / self.__curve)) < self.__aprad:
                print ()
                print ('Aperture radius cannot exceed physical extent of lens')
                print ()
                raise Exception
        
    def __repr__(self):
        return "SphericalRefraction({0}, {1}, {2}, {3}, {4})".format(
                self.__z0, self.__curve, self.__n1, self.__n2, self.__aprad)
        
    def getz(self):
        return self.__z0
    
    def getcurve(self):
        return self.__curve
        
    def intercept(self, ray):
        
        intersection = True
        
        # retrieves the variables we'll need
        P = ray.getp()
        K = ray.getk()
        
        # Checks for planar surface
        try:
            curRad = 1/(self.__curve)
            flat = False
        except ZeroDivisionError:
            flat = True
            
        
        if not flat:
            O = sp.array([0, 0, self.__z0 + curRad]) # Centre of the lens
            # The lens is centred on the z-axis
            
            # calculates the quadratic variables
            
            a = sp.dot(K, K)
            #print (a)
            b = 2 * sp.dot(K, P - O)
            #print (b)
            c = sp.dot(P - O, P - O) - (curRad) ** 2
            #print (c)
            
            # Performs the quadratic calculation
            lambdaPlus = ((-1 * b) + sp.sqrt(b*b - (4 * a * c))) / (2 * a)
            lambdaNeg = ((-1 * b) - sp.sqrt(b*b - (4 * a * c))) / (2 * a)

            # Correct intercept depends on whether convex/concave
            if self.__curve > 0:
                lam = min([lambdaPlus, lambdaNeg])
            elif self.__curve < 0:
                lam = max([lambdaPlus, lambdaNeg])
            
            # Puts lambda back into vector equation to find intersection
            final = P + (K * lam)
            
            #CHECKED MANUALLY FOR 1 AND 2D
            
            if final.dtype == 'complex128': 
                # complex solutions mean no intercept
                intersection = False
            
        
        elif flat:
            try: # Handles the case where the ray and surface are parallel
                mu = (self.__z0 - P[2]) / K[2]
            except ZeroDivisionError:
                intersection = False
            
            if intersection == True:
                final = P + (K * mu)
                
            else:
                print ('Ray parallel to surface')
                
        # Handles constraints on the size of lens
        if intersection == True:
            if abs(final[0]) > self.__aprad or abs(final[1]) > self.__aprad:
                intersection = False
                    
        if intersection == False:
            final = None
        
        return final # final == None corresponds to no intersection
            
    
    def getnormal(self, point):
        
        # returns a vector
        
        if self.__curve == 0:
            normal = sp.array([0, 0, -1])
            
        else:
            centre = sp.array([0, 0, self.__z0 + (1/self.__curve)])
            normal = point - centre
            
        if self.__curve < 0:
            normal = -1 * normal
            
        return normal
    
    
    def refract(self, khat, nhat):
        
        check = sp.sin( sp.arccos(sp.dot(khat, nhat)))
        
        
        ## checking for total internal reflection
        
        if check > (self.__n2 / self.__n1):
            k2hat = None
            
        else:
            # Vector form of Snell's law
            alongsurface = (self.__n1 / self.__n2) * (sp.cross(nhat,
                           sp.cross(-1 * nhat, khat)))
            
            alongnormal = sp.sqrt( 1 - (((self.__n1 / self.__n2)**2) *
                                        sp.dot(sp.cross(nhat, khat), \
                                               sp.cross(nhat, khat))) )
            
            k2hat = alongsurface - (nhat * alongnormal)
            
        
        # Has been manually checked in 2D
            
        return k2hat

    def propagate_ray(self, ray):
        '''
        To do:
            Handle None inputs -- Done, tested
        '''
        
        terminate = False
        
        K = ray.getk()
        
        if K[2] != 0:
            termPoint = ray.getp() + K * ((self.__z0 - ray.getp()[2])/ K[2])
        else:
            termPoint = [0,0,-20] # Hardcoded start point of rays
        
        intersection = self.intercept(ray) # May return none
        
        if intersection is None:
            terminate = True # If ray does not meet optical element
             
        else:
        
            surfaceNormal = self.getnormal(intersection)
            
        
        
            # NORMALISE BOTH DIRECTION AND NORMAL HERE
        
            k1norm = K / sp.sqrt( abs(sp.dot(K, K)) )
            
            nnorm = surfaceNormal / sp.sqrt( abs( sp.dot(surfaceNormal, \
                                                         surfaceNormal)) )
            
            k2norm = self.refract(k1norm, nnorm) # May return none
            
            if k2norm is None:
                terminate = True # If total internal reflection is observed
            
            else:
                newRay = ray.append(intersection, k2norm)
            
        if terminate == True:
            # Sets new position (@ intersection) with direction [0,0,0]
            newRay = ray.terminate(termPoint)
            
        return newRay
    
    
    def plotlensyz(self, rayarray): # Plots lens in figures
        yfine = sp.linspace(-1 * self.__aprad, self.__aprad, 101)
        
        if self.__curve != 0:
            radius = 1 / self.__curve
            if radius > 0:
                z =  self.__z0 + radius - sp.sqrt(radius**2 - yfine**2)
            elif radius < 0:
                z =  self.__z0 + radius + sp.sqrt(radius**2 - yfine**2)
        else:
            z = sp.full(len(yfine), self.__z0)
        
        return z, yfine


class OutputPlane(OpticalElement):
    
    """
    Class for the output 'screen'
    Rays terminate at this surface - do not use for further propagation
    """
    
    def __init__(self, z0 = 0):
        super().__init__(z0)        
        self.__z0 = z0 # z-intercept
        
    def __repr__(self):
        return "OutputPlane({0})".format(self.__z0)
        
    def intercept(self, ray):
        # Same as curvature == 0 case above
        P = ray.getp()
        K = ray.getk()
        
        intersection = True
        try:
            mu = (self.__z0 - P[2]) / K[2]
        except ZeroDivisionError:
            intersection = False
            
        if intersection == True:
            final = P + (K * mu)
            
        else:
            print ('Ray parallel to surface')
                    
            
        if intersection == False:
            final = None
            
        return final
    
        
    def propagate_ray(self, ray):
        intersection = self.intercept(ray)
        
        endRay = ray.terminate(intersection)
        # Ray ALWAYS terminates
        
        return endRay
    
    def plotlensyz(self, rayarray):
        # Plots output plane in figures
        # Extent given by most displaced ray
        rayposfinal = []
        
        for ray in rayarray:
            rayposfinal.append(ray.getp())
            
        aprad = max(row[1] for row in rayposfinal)
        yfine = sp.linspace(-1 * aprad, aprad, 101)
        z = sp.full(len(yfine), self.__z0)
        
        return z, yfine
