#!/usr/local/bin/python3 -i

script_version = "17-228a"
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

from collections import abc
from datetime import datetime
from decimal import Decimal
from fractions import Fraction
from random import randint
import math
import string

### General Math ###
# Evaluate the quadratic formula for ax^2+bx+c=0
def quad_det(a, b, c):
    return b**2-4*a*c
def quad(a, b, c):
    return [(-b+sqrt(quad_det(a, b, c)))/(2*a),
            (-b-sqrt(quad_det(a, b, c)))/(2*a)]

# Convert x to scientific notation
sigfig = 6 # if this value is < 1, it is assumed to be 1
def sci(x):
    global sigfig
    if sigfig < 1:
        sigfig = 1
    string = "{:." + str(sigfig - 1) + "e}"
    return string.format(x)
# Convert x to engineering notation
def eng(x):
    sci_string = sci(x)
    components = sci_string.split('e')
    exponent = int(components[1])
    offset = exponent % 3
    head = components[0].replace('.', '')
    is_negative = False
    if head[0] == '-':
        is_negative = True
        head = head[1:]
    head = head[:(offset + 1)] + ('.' if sigfig > 1 else '') + \
            head[(offset + 1):]
    new_exponent = exponent - offset
    out = ('-' if is_negative else '') + head + 'e' + \
            ('+' if new_exponent >= 0 else '-') + \
            ('0' if abs(new_exponent) < 10 else '') + str(abs(new_exponent))
    return out

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
# General
pi = math.pi
e = math.e
# Chemistry and Thermodynamics
r_atm = 0.0820574614 # Gas constant (L*atm/mol/K)
r_mmhg = 62.3636711 # Gas constant (L*mmHg/mol/K)
r_joule = 8.314462175 # Gas constant (J/mol/K)
kw = 1.01e-14 # Equilibrium constant for auto-ionization of water
avo = 6.022e23 # Avogadro constant (/mol)
# Mechanics
g = 9.80665 # Acceleration due to gravity (m/s)
g_ft = 32.174049 # Acceleration due to gravity (ft/s)
G = 6.673e-11 # Gravitational constant (N*m^2/kg^2)
# Electromagnetism
ele = 1.602e-19 # Elementary charge (C)
mu0 = math.pi / 2.5e+6 # Permiability of a vacuum (N/A^2)
e0 = 8.854e-12 # Permittivity of a vacuum (F/m)
ke = 8.988e9 # Coulomb constant (N*m^2/C^2)
c = 2.998e8 # Speed of light in a vacuum (m/s)
# Quantum Physics
h = 6.626e-34 # Planck constant (J*s)
hbar = 1.054e-34 # h-bar constant (J*s)
me = 9.109e-31 # Electron mass (kg)
mp = 1.673e-27 # Proton mass (kg)
mn = 1.675e-27 # Neutron mass (kg)
ev = 1.602e-19 # Electron-volt (J)


### Temperature Conversions ###
abs_zero_f = -459.67
abs_zero_c = -273.15
abs_zero_k = 0
def temp_check_zero(temp, zero):
    if temp < zero:
        print("Invalid. Result is below absolute zero.")
        return float('NaN')
    return round(temp, 8)
def temp_fc(f):
    c = (f-32) * (5/9)
    return temp_check_zero(c, abs_zero_c)
def temp_cf(c):
    f = c * (9/5) + 32
    return temp_check_zero(f, abs_zero_f)
def temp_ck(c):
    k = c + 273.15
    return temp_check_zero(k, abs_zero_k)
def temp_kc(k):
    c = k - 273.15
    return temp_check_zero(c, abs_zero_c)
def temp_fk(f):
    k = (f+459.67) * (5/9)
    return temp_check_zero(k, abs_zero_k)
def temp_kf(k):
    f = k * (9/5) - 459.67
    return temp_check_zero(f, abs_zero_f)


### Fractions ###
def getfrac(x):
    return Fraction(x).limit_denominator()
def frac(x):
    print(getfrac(x))
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
def log(x): return math.log(x)
def log10(x): return math.log10(x)
def logbase(x, y): return math.log(x, y)
def exp(x): return math.exp(x)
def exp10(x, n): return x * 10**n
def sqrt(x): return math.sqrt(x)
def nrt(x, y): return math.pow(x, 1/y)
def rec(x): return 1/x
def fact(x): return math.factorial(x)
def gamma(x): return math.gamma(x)
def hypot(x): return math.hypot(x)
def floor(x): return math.floor(x)
def ceil(x): return math.ceil(x)
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
# Convert base unit to SI prefix
def to_yotta(x): return exp10(x, -24)
def to_zetta(x): return exp10(x, -21)
def to_exa(x): return exp10(x, -18)
def to_peta(x): return exp10(x, -15)
def to_tera(x): return exp10(x, -12)
def to_giga(x): return exp10(x, -9)
def to_mega(x): return exp10(x, -6)
def to_kilo(x): return exp10(x, -3)
def to_centi(x): return exp10(x, 2)
def to_milli(x): return exp10(x, 3)
def to_micro(x): return exp10(x, 6)
def to_nano(x): return exp10(x, 9)
def to_angstrom(x): return exp10(x, 10)
def to_pico(x): return exp10(x, 12)
def to_femto(x): return exp10(x, 15)
def to_atto(x): return exp10(x, 18)
def to_zepto(x): return exp10(x, 21)
def to_yocto(x): return exp10(x, 24)
# Convert SI prefix to base unit
def from_yotta(x): return exp10(x, 24)
def from_zetta(x): return exp10(x, 21)
def from_exa(x): return exp10(x, 18)
def from_peta(x): return exp10(x, 15)
def from_tera(x): return exp10(x, 12)
def from_giga(x): return exp10(x, 9)
def from_mega(x): return exp10(x, 6)
def from_kilo(x): return exp10(x, 3)
def from_centi(x): return exp10(x, -2)
def from_milli(x): return exp10(x, -3)
def from_micro(x): return exp10(x, -6)
def from_nano(x): return exp10(x, -9)
def from_angstrom(x): return exp10(x, -10)
def from_pico(x): return exp10(x, -12)
def from_femto(x): return exp10(x, -15)
def from_atto(x): return exp10(x, -18)
def from_zepto(x): return exp10(x, -21)
def from_yocto(x): return exp10(x, -24)


# Generate diceware values
def diceware(n):
    for i in range(0, n):
        print(str(i + 1) + ": ", end="")
        for i in range(0, 5):
            print(randint(1, 6), end="")
        print()

## Vectors ##
# Note that all of these functions are intended for 3-dimensional vectors
def vcross(a, b): # cross (vector) product of vector a and vector b
    return [a[1] * b[2] - a[2] * b[1], \
            a[2] * b[0] - a[0] * b[2], \
            a[0] * b[1] - a[1] * b[0]]
def vadd(a, b): # add vector a to vector b
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
def vneg(a): # negate vector a
    return [-a[0], -a[1], -a[2]]
def vsub(a, b): # subtract vector b from vector a
    return vadd(a, vneg(b))
def vdot(a, b): # dot (scalar) product of vector a and vector b
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def vscale(a, alpha): # scale vector a by scalar alpha
    return [a[0] * alpha, a[1] * alpha, a[2] * alpha]
def vlen(a): # get absolute value
    return sqrt(a[0]**2 + a[1]**2 + a[2]**2)
def vproj(a, b): # projection of a onto b
    return vscale(b, vdot(a, b) / vdot(b, b))
def vunit(a): # makes a unit vector
    return vscale(a, 1 / vlen(a))
def vtheta(a, b): # find the angle between two vectors in radians
    return acos(vdot(a, b) / (vlen(a) * vlen(b)))
def dvtheta(a, b): # find the angle between two vectors in degrees
    return deg(vtheta(a, b))

## Lists ##
def flatten_list(*x):
    n = x
    m = []
    is_flat = False
    while not is_flat:
        is_flat = True
        for i in n:
            if isinstance(i, abc.Iterable):
                # for every list, add its elements to n, even if those elements
                # are themselves lists. The while loop will repeat until the
                # heirarchy is flat.
                is_flat = False
                for j in i:
                    m.append(j)
            else:
                m.append(i)
        n = m
        m = []
    return n
def to_float_list(*x):
    n = flatten_list(x)
    m = []
    for i in n:
        m.append(float(i))
    return m
def to_int_list(*x):
    n = flatten_list(x)
    m = []
    for i in n:
        m.append(int(i))
    return m
def fsum(*x): # Floating-point sum of a list
    return math.fsum(to_float_list(x))
def isum(*x): # Integer sum of a list
    n = to_int_list(x)
    sum = 0
    for i in n:
        sum += i
    return int(sum)
def mean(*x): # Arithmetic mean
    n = flatten_list(x)
    return fsum(n) / len(n)
def stdDev(*x): # Population Standard Deviation
    n = to_float_list(x)
    avg = mean(n)
    total_deviation = 0
    for i in n:
        total_deviation += (i - avg) ** 2
    return sqrt(1 / len(n) * total_deviation)
def pctRSD(*x): # %RSD
    try:
        return stdDev(x) / mean(x) * 100
    except ZeroDivisionError:
        return float('NaN')

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

# Pythagorean theorem
def pyth(a, b):
    return sqrt(a**2 + b**2)
def pythleg(c, a):
    return sqrt(c**2 - a**2)

# Convert an int n to an arbitrary base b (as string)
def to_base(n, b):
    num_dict = string.digits + string.ascii_lowercase
    if n < 0:
        sign = -1
    elif n == 0:
        return num_dict[0]
    else:
        sign = 1
    n *= sign
    digits = []
    while n:
        digits.append(num_dict[n % b])
        n //= b # Double-slash forces integer division
    if sign < 0:
        digits.append('-')
    digits.reverse()
    return ''.join(digits)

### Cellular Data Statistics ###
def days_in_month(month):
    shortmonths = [4, 6, 9, 11]
    if month in shortmonths:
        return 30
    elif month == 2: # Don't bother with leap year
        return 28
    else: # For simplicity, assume 31 if input is invalid.
        return 31

def data(gb, total):
    reset_day = 11 # Day of month on which billing month rolls over
    
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


### Exit functions ###
def exit():
    import sys
    sys.exit()
def quit(): exit()
def bye(): exit()

print("Loaded calc.py")
