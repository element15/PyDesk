#!/usr/local/bin/python3 -i

script_version = "15-914a"
# calc.py
#
# Loads a list set of functions and variables for everyday calculator
# functionality. Written for use with Python 3, but *should* work fine with
# Python 2.
#
#
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or distribute
# this software, either in source code form or as a compiled binary, for any
# purpose, commercial or non-commercial, and by any means.
#
# In jurisdictions that recognize copyright laws, the author or authors of this
# software dedicate any and all copyright interest in the software to the public
# domain. We make this dedication for the benefit of the public at large and to
# the detriment of our heirs and successors. We intend this dedication to be an
# overt act of relinquishment in perpetuity of all present and future rights to
# this software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import math
from datetime import datetime
from fractions import Fraction
from decimal import Decimal

### General Math ###
# Evaluate the quadratic formula for ax^2+bx+c=0
def quad_det(a,b,c):
    return b**2-4*a*c
def quad(a,b,c):
    return [(-b+sqrt(quad_det(a,b,c)))/(2*a),(-b-sqrt(quad_det(a,b,c)))/(2*a)]

# Convert x to scientific notation
sigfig = 4
def sci(x):
    string = "{:." + str(sigfig - 1) + "e}"
    return string.format(x)

# Convert x to engineering notation
def eng(x):
    return Decimal(str(x)).normalize().to_eng_string()

# Get the midpoint
def mid(a, b):
    return (a+b)/2
def mid2(x1, y1, x2, y2):
    return [mid(x1, x2), mid(y1, y2)]

# Distance
def dist(a, b):
    return b-a
def dist2(x1, y1, x2, y2):
    return (dist(x1, x2)**2 + dist(y1, y2)**2)**0.5



### Constants ###
# Gas constant
r_atm = 0.0820574614 # L*atm*mol^-1*K^-1
r_mmhg = 62.3636711 # L*mmHg*mol^-1*K^-1
r_joule = 8.314462175 # J*K^-1*mol^-1
kw = 1.01e-14 # Equilibrium constant for auto-ionization of water, kw
pi = math.pi
e = math.e
g = 9.807 # Acceleration due to gravity in m*s^-1

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
def logbase(x, y): return math.log(x, y)
def exp(x): return math.exp(x)
def pow(x, y): return math.pow(x, y)
def sqrt(x): return math.sqrt(x)
def nrt(x, y): return math.pow(x, 1/y)
def abs(x): return math.fabs(x)
def fact(x): return math.factorial(x)
def gamma(x): return math.gamma(x)
def hypot(x): return math.hypot(x)
def floor(x): return math.floor(x)
def ceil(x): return math.ceil(x)
def sum(*x): return math.fsum(x)
# Trigonometry (degrees functions begin with 'd')
def rad(x): return math.radians(x)
def deg(x): return math.degrees(x)
def sin(x): return math.sin(x)
def dsin(x): return sin(rad(x))
def cos(x): return math.cos(x)
def dcos(x): return cos(rad(x))
def sec(x): return 1/math.cos(x)
def dsec(x): return sec(rad(x))
def csc(x): return 1/math.sin(x)
def dcsc(x): return csc(rad(x))
def tan(x): return math.tan(x)
def dtan(x): return tan(rad(x))
def cot(x): return 1/math.tan(x)
def dcot(x): return cot(rad(x))
def asin(x): return math.asin(x)
def dasin(x): return deg(asin(x))
def acos(x): return math.acos(x)
def dacos(x): return deg(acos(x))
def atan(x): return math.atan(x)
def datan(x): return deg(atan(x))
def sinh(x): return math.sinh(x)
def dsinh(x): return sinh(rad(x))
def cosh(x): return math.cosh(x)
def dcosh(x): return cosh(rad(x))
def tanh(x): return math.tanh(x)
def dtanh(x): return tanh(rad(x))
def asinh(x): return math.asinh(x)
def dasinh(x): return deg(asinh(x))
def acosh(x): return math.acosh(x)
def dacosh(x): return deg(acosh(x))
def atanh(x): return math.atanh(x)
def datanh(x): return deg(atanh(x))

### Specialty Math ###
# Economics
def elastic(q1, q2, p1, p2):
    q = dist(q1, q2) / mid(q2, q1)
    p = dist(p1, p2) / mid(p2, p1)
    return q/p

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

# Pythagorean Theorem. Nuff said.
def pyth(a, b):
    return sqrt(a**2 + b**2)
def pythleg(c, a):
    return sqrt(c**2 - a**2)

### Cellular Data Statistics ###
def days_in_month(month):
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
        totalDays = days_in_month(now.month)
        cycleDay = now.day - (reset_day - 1)
    else:
        totalDays = days_in_month(now.month - 1)
        cycleDay = now.day + totalDays - (reset_day - 1)
    
    cycleUsage = gb * 1024
    idealRate = total * 1024 / totalDays
    idealUsage = idealRate * cycleDay
    cycleRate = cycleUsage / cycleDay
    netUsage = idealUsage - cycleUsage
    netRate = idealRate - cycleRate
    coefficient = cycleRate / idealRate
    daysUsed = cycleUsage / idealRate
    
    print("     Cycle Usage: %d MiB" % cycleUsage)
    print("     Ideal Usage: %d MiB" % idealUsage)
    print("       Net Usage: %d MiB" % netUsage)
    print("      Cycle Rate: %d MiB/day" % cycleRate)
    print("      Ideal Rate: %d MiB/day" % idealRate)
    print("        Net Rate: %d MiB/day" % netRate)
    print(" Use Coefficient: %f" % coefficient)
    print("       Cycle Day: %d / %d" % (cycleDay, totalDays))
    print("       Ideal Day: %d" % daysUsed)
    
    if netUsage < 0:
        daysBehind = -netUsage / idealRate + 1
        print("        Catch up: %d" % daysBehind)


print("Loaded calc.py")
### End of script ###
