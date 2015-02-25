#!/usr/local/bin/python3 -i

script_version = "15-225-a"
# calc.py
# 
# Loads a list set of functions and variables for everyday calculator
# functionality. Written for use with Python 3, but *should* work fine
# with Python 2.


# This is free and unencumbered software released into the public domain.
# 
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
# 
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


import math
from datetime import datetime
from fractions import Fraction

# Evaluate the quadratic formula for ax^2+bx+c=0
def quad_det(a,b,c):
    return b**2-4*a*c
def quad(a,b,c):
    return [(-b+quad_det(a,b,c))/(2*a),(-b-quad_det(a,b,c))/(2*a)]

# Convert a number to scientific notation, 5 significant figures
def sci(x):
    return "{:.4e}".format(x)

### Constants ###
# Gas constant
r_atm = 0.0820574614 # L*atm*mol^-1*K^-1
r_mmhg = 62.3636711 # L*mmHg*mol^-1*K^-1
r_joule = 8.314462175 # J*K^-1*mol^-1
kw = 1.01e-14 # Equilibrium constant for auto-ionization of water, kw
pi = math.pi
e = math.e

### Temperature Conversions ###
f_zero = -459.67
c_zero = -273.15
k_zero = 0
def valid_temp(temp,zero):
    if temp < zero:
        print("Impossibru!")
        return float('NaN')
    return round(temp,8)
def fc_temp(f):
    c = (f-32) * (5/9)
    return valid_temp(c,c_zero)
def cf_temp(c):
    f = c * (9/5) + 32
    return valid_temp(f,f_zero)
def ck_temp(c):
    k = c + 273.15
    return valid_temp(k,k_zero)
def kc_temp(k):
    c = k - 273.15
    return valid_temp(c,c_zero)
def fk_temp(f):
    k = (f+459.67) * (5/9)
    return valid_temp(k,k_zero)
def kf_temp(k):
    f = k * (9/5) - 459.67
    return valid_temp(f,f_zero)


### Fractions ###
def getfrac(x):
    return Fraction(x).limit_denominator()
def frac(x):
    print( getfrac(x) )
def mix(x):
    fraction = getfrac(x)
    numerator = fraction.numerator
    denominator = fraction.denominator
    if numerator > denominator:
        whole = math.floor(x)
        mixed_numerator = numerator - whole * denominator
        print("%d %d/%d" % (whole, mixed_numerator, denominator))
    else:
        print("%d/%d" % (numerator, denominator))

### Convenience Functions ###
def ln(x): return math.log(x)
def log(x): return math.log10(x)
def rad(x): return math.radians(x)
def deg(x): return math.degrees(x)
def e(x): return math.exp(x)
def logbase(x,y): return math.log(x,y)
def sqrt(x): return math.sqrt(x)
def sum(*x): return math.fsum(*x)
def abs(x): return math.fabs(x)
def fact(x): return math.factorial(x)
def gamma(x): return math.gamma(x)
def hypot(x): return math.hypot(x)
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
def tanh(x): return math.tanh(x)
def asinh(x): return math.asinh(x)
def acosh(x): return math.acosh(x)
def atanh(x): return math.atanh(x)
def floor(x): return math.floor(x)
def ceil(x): return math.ceil(x)


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
def daysInMonth(month):
    shortmonths = [4,6,9,11]
    if month in shortmonths:
        return 30
    elif month == 2: # Don't bother with leap year
        return 28
    else: # For simplicity, assume 31 if input is invalid.
        return 31

def data(gb,total):
    
    reset_day = 13 # Day of month on which billing month rolls over
    
    now = datetime.now()
    if now.day >= reset_day:
        totalDays = daysInMonth(now.month)
        cycleDay = now.day - (reset_day - 1)
    else:
        totalDays = daysInMonth(now.month - 1)
        cycleDay = now.day + totalDays - (reset_day - 1)
    
    cycleUsage = gb * 1024
    idealRate = total * 1024 / totalDays
    idealUsage = idealRate * cycleDay
    cycleRate = cycleUsage / cycleDay
    netUsage = cycleUsage - idealUsage
    netRate = cycleRate - idealRate
    coefficient = cycleRate / idealRate
    daysUsed = cycleUsage / idealRate
    
    print("     Cycle Usage: %f MB" % cycleUsage)
    print("     Ideal Usage: %f MB" % idealUsage)
    print("       Net Usage: %f MB" % netUsage)
    print("      Cycle Rate: %f MB/day" % cycleRate)
    print("      Ideal Rate: %f MB/day" % idealRate)
    print("        Net Rate: %f MB/day" % netRate)
    print(" Use Coefficient: %f" % coefficient)
    print("       Cycle Day: %d / %d" % (cycleDay, totalDays))
    print("       Ideal Day: %d" % daysUsed)
    
    if netUsage > 0:
        daysBehind = netUsage / idealRate + 1
        print("        Catch up: %d" % daysBehind)


print("Loaded calc.py")
### End of script ###
