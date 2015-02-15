#!/usr/local/bin/python3 -i

# calc.py
# 
# Loads a list of user-defined functions for using the python interpreter as a
# calculator. Analogue to .bashrc, .zshrc, and the like.

import math

# Evaluate the quadratic formula for ax^2+bx+c=0
def quad(a,b,c):
    determinant = math.sqrt(b**2-4*a*c)
    denomimator = 2*a
    return [(-b+determinant)/denominator,(-b-determinant)/denominator]

# Convert a number to scientific notation, 5 significant figures
def sci(x): return "{:.4e}".format(x)

### Constants ###
# Gas constant
r_atm = 0.0820574614 # L*atm*mol^-1*K^-1
r_mmhg = 62.3636711 # L*mmHg*mol^-1*K^-1
r_joule = 8.314462175 # J*K^-1*mol^-1
kw = 1.01e-14 # Equilibrium constant for auto-ionization of water, kw
pi = math.pi
e = math.e

### Convenience Functions ###
def ln(x): return math.log(x)
def log(x): return math.log10(x)
def rad(x): return math.radians(x)
def deg(x): return math.degrees(x)
def e(x): return math.exp(x)
def logbase(x,y): return math.log(x,y)
def sqrt(x): return math.sqrt(x)
def cumsum(*x): return math.fsum(*x)
def abs(x): return math.fabs(x)
def fact(x): return math.factorial(x)
def gamma(x): return math.gamma(x)
def sin(x): return math.sin(x)
def cos(x): return math.cos(x)
def sec(x): return 1/math.cos(x)
def csc(x): return 1/math.sin(x)
def tan(x): return math.tan(x)
def cot(x): return 1/math.tan(x)
def asin(x): return math.asin(x)
def acos(x): return math.acos(x)
def atan(x): return math.atan(x)
def sinh(x): return math.sinh(x)
def cosh(x): return math.cosh(x)
def asinh(x): return math.asinh(x)
def acosh(x): return math.acosh(x)

# Compute approximate golden ratios using fibonacci
def gold(n):
    i = 0
    fib2 = 1
    fib1 = 1
    thisfib = 1
    while i < n:
        i += 1
        fib2 = fib1
        fib1 = thisfib
        thisfib = fib1 + fib2
    ratio = thisfib / fib1
    print(str(thisfib) + "/" + str(fib1) + " = " + str(ratio))

### MiFi Data Usage Statistics ###
from datetime import datetime

def daysInMonth(month):
    shortmonths = [4,6,9,11]
    if month in shortmonths:
        return 30
    elif month == 2: # Don't bother with leap year
        return 28
    else: # For simplicity, assume 31 if input is invalid.
        return 31

def data(gb,total):
    now = datetime.now()
    if now.day >= 13:
        totalDays = daysInMonth(now.month)
        cycleDay = now.day - 12
    else:
        totalDays = daysInMonth(now.month - 1)
        cycleDay = now.day + totalDays - 12
    
    cycleUsage = gb * 1024
    idealRate = total * 1024 / totalDays
    idealUsage = idealRate * cycleDay
    cycleRate = cycleUsage / cycleDay
    netUsage = cycleUsage - idealUsage
    netRate = cycleRate - idealRate
    coefficient = cycleRate / idealRate
    daysUsed = cycleUsage / idealRate
    
    print("Cycle Usage: %f MB" % cycleUsage)
    print("Ideal Usage: %f MB" % idealUsage)
    print("Net Usage: %f MB" % netUsage)
    print("Cycle Rate: %f MB/day" % cycleRate)
    print("Ideal Rate: %f MB/day" % idealRate)
    print("Net Rate: %f MB/day" % netRate)
    print("Use Coefficient: %f" % coefficient)
    print("Cycle Day: %d / %d" % (cycleDay, totalDays))
    print("Ideal Day: %d" % daysUsed)
    
    if netUsage > 0:
        daysBehind = netUsage / idealRate + 1
        print("Days to catch up: %d" % daysBehind)


print("Loaded calc.py")
### End of script ###
