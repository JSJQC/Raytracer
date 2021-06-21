# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 13:00:49 2020

@author: jakes

Module containing test code that demonstrates accuracy and functionality
Should be run cell-by-cell as some tests are designed to raise exceptions
"""

import raytracer as ryt
import opticalequipment as oe
import System as sy
import optimiser as op


import scipy as sp

#%%

# Testing docstrings
help(ryt.Ray)
help(sy)

#%% Tests of Ray() attributes and methods

testray = ryt.Ray()
print (testray)

print ()

try:
    print (testray.__p)
except AttributeError:
    print ("Hidden attribute check successful")
    
print ()

testray.append([0,0,0], [1,1,1])
print (testray)

print ("After termination: ")
testray.terminate([1,1,1])
print (testray)

print ()
print ("All directions:")
print (testray.checkallk())
print ("All positions:")
print (testray.vertices())
print ()

#%% Tests of incorrect SphericalRefraction instantiation

# Testing aperture radius bigger than possible physical size of lens
testsphere = oe.SphericalRefraction(0, 0.03, 1, 1.5, 50)

#%% Tests of interception and normal calculations

testsphere = oe.SphericalRefraction()
print (testsphere)

testsphere = oe.SphericalRefraction(0, 0.03, 1, 1.5, 33)
print (testsphere)

print ()

ray1 = ryt.Ray([0,0,-20], [1,1,1])
ray2 = ryt.Ray([0,0,-20], [0,0.1,1])
ray3 = ryt.Ray([0,0,-20], [0,1,0])


# Should return None for no intersection
print (testsphere.intercept(ray1))

# Should return [0, 2.006, 0.06] as done by hand
intercept = testsphere.intercept(ray2)
print (intercept)

print ()
testplane = oe.SphericalRefraction(0, 0, 1, 1.5, 33)
print (testplane)

# Should return [20, 20, 0] - flat lens is a square
# If combined with a spherical lens, this extra extent does not matter, as
# the ray will terminate if it misses a preceding curved lens
print (testplane.intercept(ray1))

print ()

# Should both return None for a ray with no z-direction
print (testsphere.intercept(ray3))
print (testplane.intercept(ray3))

# Check that a flat surface always give a purely -ve z facing normal
print (testplane.getnormal([1, 1, 0]))
# Expect [0, 2.006, -33.273] as done by hand
print (testsphere.getnormal(intercept))


#%% Tests for concave lens

testconcave = oe.SphericalRefraction(0, -0.03, 1, 1.5, 33)
ray2 = ryt.Ray([0,0,-20], [0, 0.1, 1])

# Intercept should be different to the example above
intercept = testconcave.intercept(ray2)
print (intercept)

# Normal has a -ve y component, as lens is curving back round towards lower z
print (testconcave.getnormal(intercept))

#%% Test to check refraction and propagation

print ()
testsphere = oe.SphericalRefraction(0, 0.03, 1, 1.5, 33)

raygood = ryt.Ray([0,0,-20], [0,0.1,1]) # meets lens
raybad = ryt.Ray([0,0,-20], [1,1,1]) # misses lens
rayawful = ryt.Ray([0,0,-20], [0,1,0]) # no z-component

# Should carry on with direction [0, 0.046, 0.10] as checked by hand
# Could not easily be tested manually for 3D
print (testsphere.propagate_ray(raygood))
print (testsphere.propagate_ray(raybad)) # Should terminate at plane of lens
print (testsphere.propagate_ray(rayawful)) # Should terminate at [0,0,-20]


#%% Testing of output plane

# intercept is the same case as curvature = 0 previously

testray = ryt.Ray([0,0,-20], [1,1,1])

outputplane = oe.OutputPlane(0)

# Ray should have [0,0,0] direction @ position [20, 20, 0]
print (outputplane.propagate_ray(testray))



#%% Checking simple RMS calculation

ray1 = ryt.Ray([1,1,-20], [0,0,1])
ray2 = ryt.Ray([-1,1,-20], [0,0,1])
ray3 = ryt.Ray([0, -1 *sp.sqrt(3) + 1, -20], [0,0,1])

rays = [ray1, ray2, ray3]

outputplane = oe.OutputPlane(0)
outputplane.encounterlens(rays)


# Expect sqrt[ (1/3) * (2 + 2 + (sqrt3 - 1)^2) ] = 1.2296 as the answer
print (outputplane.rms(rays))


#%% Extreme example of terminating many rays

widebeam = sy.createBundle(20, 3, False) # Creates wide, flat beam

smalllens = sy.createLens(0, 0.1, 1, 1.5, 10)
output = sy.createOutput(30)

sy.propagateSystem([smalllens, output], widebeam)

sy.plotyz(widebeam, [smalllens, output], 'Wide beam, small lens')



#%% A demonstration of key functions in System.py
# All values are theoretically obtained given the lens parameters
# Agreement of graph with expectations signifies correct application of theory

    
# n = 1.5168 is the refractive index of glass

convexleft = sy.createLens(42, 0.03, 1, 1.5168, 17.794)
planeright = sy.createLens(47, 0, 1.5168, 1, 17.794)
    
planeleft = sy.createLens(77, 0, 1, 1.5168, 17.794)
convexright = sy.createLens(82, -0.03, 1.5168, 1, 17.794)
    
outputsurface = sy.createOutput(144)
    
lenses = [convexleft, planeright, planeleft, convexright, outputsurface]

bundle = sy.createBundle(20, 0.5)

finalrays = sy.propagateSystem(lenses, bundle)
    
sy.plotyz(finalrays, lenses, 'Planoconvex demonstration')
# Should plot all rays, along with all optical elements
    
    
#%% Tests of various lenses

# Converging lenses -----------------

output = sy.createOutput(150)
convexconverge = sy.createLens(0, 0.03, 1, 1.5168, 33) 
concaveconverge = sy.createLens(0, -0.03, 1.5168, 1, 33) # Note switched n
lenses1 = [convexconverge, output]
lenses2 = [concaveconverge, output]

# Both these cases should converge

rays = sy.createBundle(10, 1, False) # Creates flat 'beam'
sy.propagateSystem(lenses1, rays)
sy.plotyz(rays, lenses1, 'Convex converging')

rays = sy.createBundle(10, 1, False) # Creates flat 'beam'
sy.propagateSystem(lenses2, rays)
sy.plotyz(rays, lenses2, 'Concave converging')

# --------------------------------------------

# Diverging lenses ----------------------------

convexdiverge = sy.createLens(0, 0.03, 1.5168, 1, 33) # Note change in n
concavediverge = sy.createLens(0, -0.03, 1, 1.5168, 33)

lenses3 = [convexdiverge, output]
lenses4 = [concavediverge, output]

# Both these cases should diverge

rays = sy.createBundle(10, 1, False) # Creates flat 'beam'
sy.propagateSystem(lenses3, rays)
sy.plotyz(rays, lenses3, 'Convex diverging')

rays = sy.createBundle(10, 1, False) # Creates flat 'beam'
sy.propagateSystem(lenses4, rays)
sy.plotyz(rays, lenses4, 'Concave diverging')

# ---------------------------------------------


#%% Test of an angled bundle

angledbundle = sy.createBundle(10, 0.5, False, 0, 0.1) 
# final parameter is the angle in radians

convexleft = sy.createLens(0, 0.03, 1, 1.5, 17.55)
planeright = sy.createLens(5, 0, 1.5, 1, 17.55)

output = sy.createOutput(150)

sy.propagateSystem([convexleft, planeright, output], angledbundle)
sy.plotyz(angledbundle, [convexleft, planeright, output], 'Angled bundle')

#%% Testing a bundle at too great an angle

angledbundle = sy.createBundle(10, 0.5, False, 0, 1.6)
# final parameter is the angle in radians

# Code below here should not execute -------

convexleft = sy.createLens(0, 0.03, 1, 1.5, 17.55)
planeright = sy.createLens(5, 0, 1.5, 1, 17.55)

output = sy.createOutput(150)

sy.propagateSystem([convexleft, planeright, output], angledbundle)
sy.plotyz(angledbundle, [convexleft, planeright, output], 'Angled bundle')


#%% Testing the code that calculates lens 2 curvature given lens 1 curvature

# In optimiser.py module


# Change second parameter if you want to test other curvatures
testlens1 = sy.createLens(0, 0.03, 1, 1.5, 6.5)
testlens2 = sy.createLens(0, 0.06, 1, 1.5, 6.5)

outputlens1 = op.adjustor3(testlens1, 100, 1.5, 2)
outputlens2 = op.adjustor3(testlens2, 100, 1.5, 2)

output = sy.createOutput(110)

# Both configurations should converge to around 100mm

rays = sy.createBundle(10, 0.2, False)
sy.propagateSystem([testlens1, outputlens1, output], rays)
sy.plotyz(rays, [testlens1, outputlens1, output], 'C1 = 0.03')

rays = sy.createBundle(10, 0.2, False)
sy.propagateSystem([testlens2, outputlens2, output], rays)
sy.plotyz(rays, [testlens2, outputlens2, output], 'C1 = 0.06')




