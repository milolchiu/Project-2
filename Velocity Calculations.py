#Forces acting on the particles
import math

#Drag Force
#constants
pF = 1.2 #kg/m^3
velocityApp = 10 #m/s
diameter = 0.01 #m
Re = (pF * abs(velocityApp) * diameter) / (1.8 * 10**-5) #Reynolds number
Cd = 24/Re #drag coefficient


drag = .5 * pF * Cd * ((math.pi / 4) * (diameter)**2) * abs(velocityApp) * velocityApp