#!/usr/bin/env python3 -i

script_version = "20-706a"
# calc.py
#
# Loads a list set of functions and variables for everyday calculator
# functionality.
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
from math import *
from sys import exit
import random
import re
import string
import subprocess

import numpy as np

#################
### Constants ###
#################

# Chemistry and Thermodynamics
r_atm = 0.0820574614 # Gas constant (L*atm/mol/K)
r_mmhg = 62.3636711 # Gas constant (L*mmHg/mol/K)
r_joule = 8.314462175 # Gas constant (J/mol/K)
r_btu = 1.98588 # Gas constant (BTU/lbmol/R)
r_psia = 10.7316 # Gas constant (psia*ft^3/lbmol/R)
kw = 1.01e-14 # Equilibrium constant for auto-ionization of water
avo = 6.022e23 # Avogadro constant (mol^-1)

# Mechanics
g = 9.80665 # Acceleration due to gravity (m/s)
g_ft = 32.174049 # Acceleration due to gravity (ft/s)
G = 6.673e-11 # Gravitational constant (N*m^2/kg^2)

# Electromagnetism
ele = 1.602e-19 # Elementary charge (C)
mu_0 = pi / 2.5e+6 # Permiability of a vacuum (N/A^2)
epsilon_0 = 8.854e-12 # Permittivity of a vacuum (F/m)
k_e = 8.988e9 # Coulomb constant (N*m^2/C^2)
c = 2.998e8 # Speed of light in a vacuum (m/s)

# Quantum Physics
h = 6.626e-34 # Planck constant (J*s)
hbar = 1.054e-34 # h-bar constant (J*s)
m_e = 9.109e-31 # Electron mass (kg)
m_p = 1.673e-27 # Proton mass (kg)
m_n = 1.675e-27 # Neutron mass (kg)
eV = 1.602e-19 # Electron-volt (J)

#############################
### Convenience Functions ###
#############################

# Trigonometry in degrees
sind = lambda x : sin(rad(x))
cosd = lambda x : cos(rad(x))
tand = lambda x : tan(rad(x))
asind = lambda x : deg(asin(x))
acosd = lambda x : deg(acos(x))
atand = lambda x : deg(atan(x))
sinhd = lambda x : sinh(rad(x))
coshd = lambda x : cosh(rad(x))
asinhd = lambda x : deg(asinh(x))
acoshd = lambda x : deg(acosh(x))
atanhd = lambda x : deg(atanh(x))

# Convert base unit to SI prefix
exp10 = lambda x, y : x * 10**y
to_yotta = lambda x : exp10(x, -24)
to_zetta = lambda x : exp10(x, -21)
to_exa = lambda x : exp10(x, -18)
to_peta = lambda x : exp10(x, -15)
to_tera = lambda x : exp10(x, -12)
to_giga = lambda x : exp10(x, -9)
to_mega = lambda x : exp10(x, -6)
to_kilo = lambda x : exp10(x, -3)
to_centi = lambda x : exp10(x, 2)
to_milli = lambda x : exp10(x, 3)
to_micro = lambda x : exp10(x, 6)
to_nano = lambda x : exp10(x, 9)
to_angstrom = lambda x : exp10(x, 10)
to_pico = lambda x : exp10(x, 12)
to_femto = lambda x : exp10(x, 15)
to_atto = lambda x : exp10(x, 18)
to_zepto = lambda x : exp10(x, 21)
to_yocto = lambda x : exp10(x, 24)

# Convert SI prefix to base unit
from_yotta = lambda x : exp10(x, 24)
from_zetta = lambda x : exp10(x, 21)
from_exa = lambda x : exp10(x, 18)
from_peta = lambda x : exp10(x, 15)
from_tera = lambda x : exp10(x, 12)
from_giga = lambda x : exp10(x, 9)
from_mega = lambda x : exp10(x, 6)
from_kilo = lambda x : exp10(x, 3)
from_centi = lambda x : exp10(x, -2)
from_milli = lambda x : exp10(x, -3)
from_micro = lambda x : exp10(x, -6)
from_nano = lambda x : exp10(x, -9)
from_angstrom = lambda x : exp10(x, -10)
from_pico = lambda x : exp10(x, -12)
from_femto = lambda x : exp10(x, -15)
from_atto = lambda x : exp10(x, -18)
from_zepto = lambda x : exp10(x, -21)
from_yocto = lambda x : exp10(x, -24)

########################
### Unit Conversions ###
########################

gnu_units_output = re.compile(
    r'(?:\t?(?P<reci_note>reciprocal conversion)?\n?'
    r'\t\* (?P<normal>[\d\.\-\+e]+)\n'
    r'\t/ (?P<reciprocal>[\d\.\-\+e]+))|'
    r'(?:(?P<conform_note>conformability error)\n'
    r'\t[\d\.\-\+e]+ (?P<in_unit>[^\n]+)\n'
    r'\t[\d\.\-\+e]+ (?P<out_unit>[^\n]+))')
gnu_units_executable = 'gunits'

# Evaluate a query using [GNU Units](en.wikipedia.org/wiki/GNU_Units),
# returning a tuple containing the direct conversion, and the reciprocal
# conversion, respectively.
#   v: Number to convert
#   a: Unit to convert from
#   b: Unit to convert to
# If any errors occur, or if the output does not match the expected format, the
# output of GNU Units will be returned directly as a string.
def units(v, a, b):
    result = subprocess.run([gnu_units_executable, f'{v}{a}', b],
        stdout=subprocess.PIPE).stdout.decode('utf-8')
    m = gnu_units_output.match(result)
    if m:
        if m.group('conform_note'):
            return (m.group('conform_note'), m.group('in_unit'),
                m.group('out_unit'))
        if m.group('reci_note'):
            return (float(m.group('normal')), float(m.group('reciprocal')),
                m.group('reci_note'))
        return (float(m.group('normal')), float(m.group('reciprocal')))
    return result
u = units # shorthand

# Perform a unit conversion, but only return the normal conversion instead of a
# tuple. If the conversion results in a conformability or unknown error, `None`
# is returned. In addition, the full output of `units()` is printed unles
# `silent=True` is specified.
def uu(v, a, b, silent=False):
    conv = units(v, a, b)
    if not silent:
        print(conv)
    if not isinstance(conv, tuple) or conv[0] == 'conformability error':
        return None
    return conv[0]

####################
### General Math ###
####################

# Evaluate the quadratic formula for ax^2+bx+c=0
quad_det = lambda a, b, c : b**2 - 4*a*c
def quad(a, b, c):
    soln_mean = -b/(2*a)
    soln_radius = sqrt(quad_det(a, b, c))/(2*a)
    return soln_mean+soln_radius, soln_mean-soln_radius

mid = lambda a, b : (a+b)/2
mid2d = lambda x1, y1, x2, y2 : (mid(x1, x2), mid(y1, y2))
dist2d = lambda x1, y1, x2, y2 : ((x2-x1)**2 + (y2-y1)**2)**0.5

# Linear interpolation
lint = lambda x1, xn, x2, y1, y2 : (y2 - y1) / (x2 - x1) * (xn - x1) + y1

# Pythagorean leg
leg = lambda a, c : sqrt(abs(c**2 - a**2))

######################
### Thermodynamics ###
######################

# Calculate enthalpy specific heat on a mole basis
def heat_cp_mol(T, *coeff):
    coeff = flatten_list(coeff)
    c_p = 0;
    for i in range(0, 3):
        c_p += coeff[i] * T**i
    return c_p

# Calculate internal energy specific heat on a mole basis
heat_cv_mol = lambda R, T, *coeff : heat_cp_mol(T, coeff) - R

# Calculate enthalpy on a mole basis
def heat_h_mol(T, *coeff):
    coeff = flatten_list(coeff)
    h = 0
    for i in range(0, 3):
        h += coeff[i] * T**(i + 1) / (i + 1)
    return h

# Calculate internal energy on a mole basis
def heat_u_mol(R, T, *coeff):
    coeff = flatten_list(coeff)
    h =  heat_h_mol(T, coeff)
    return h - R * T

#########################
### Number Formatting ###
#########################

# Scientific notation
def sci(x, sigfig=6):
    sigfig
    if sigfig < 1:
        sigfig = 1
    string = '{:.' + str(sigfig - 1) + 'e}'
    return string.format(x)

# Engineering notation
def eng(x, sigfig=6):
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

###############################
### Temperature Conversions ###
###############################

# Absolute zero checks
abs_zero_f = -459.67
abs_zero_c = -273.15
abs_zero_k = 0
def temp_check_zero(temp, zero):
    if temp < zero:
        print('Result is below absolute zero')
        return float('NaN')
    return round(temp, 8)

# Conversions
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

#################
### Fractions ###
#################

getfrac = lambda x : Fraction(x).limit_denominator()

frac = lambda x : print(getfrac(x))

def mix(x):
    fraction = getfrac(x)
    numerator = fraction.numerator
    denominator = fraction.denominator
    if numerator > denominator:
        whole = floor(x)
        mixed_numerator = numerator - whole * denominator
        print("%d %d/%d" % (whole, mixed_numerator, denominator))
    else:
        print("%d/%d" % (numerator, denominator))

###############
### Vectors ###
###############

# Convert an n-dimensional vector to a 3-dimensional one
def vector_to_3d(a):
    n = len(a)
    if n == 1:
        return (a[0], 0, 0)
    if n == 2:
        return (a[0], a[1], 0)
    if n >= 3:
        return (a[0], a[1], a[2])
    return (0, 0, 0) # This should only happen if n == 0

def vcross(a, b): # cross (vector) product of vector a and vector b
    a_ = vector_to_3d(a)
    b_ = vector_to_3d(b)
    return (a_[1] * b_[2] - a_[2] * b_[1], \
            a_[2] * b_[0] - a_[0] * b_[2], \
            a_[0] * b_[1] - a_[1] * b_[0])
def vadd(a, b): # add vector a to vector b
    a_ = vector_to_3d(a)
    b_ = vector_to_3d(b)
    return (a_[0] + b_[0], a_[1] + b_[1], a_[2] + b_[2])
def vneg(a): # negate vector a
    a_ = vector_to_3d(a)
    return (-a_[0], -a_[1], -a_[2])
def vsub(a, b): # subtract vector b from vector a
    return vadd(a, vneg(b))
def vdot(a, b): # dot (scalar) product of vector a and vector b
    a_ = vector_to_3d(a)
    b_ = vector_to_3d(b)
    return a_[0] * b_[0] + a_[1] * b_[1] + a_[2] * b_[2]
def vscale(a, alpha): # scale vector a by scalar alpha
    a_ = vector_to_3d(a)
    return (a_[0] * alpha, a_[1] * alpha, a_[2] * alpha)
def vlen(a): # get absolute value
    a_ = vector_to_3d(a)
    return sqrt(a_[0]**2 + a_[1]**2 + a_[2]**2)
def vproj(a, b): # projection of a onto b
    return vscale(b, vdot(a, b) / vdot(b, b))
def vunit(a): # makes a unit vector
    return vscale(a, 1 / vlen(a))
def vtheta(a, b): # find the angle between two vectors in radians
    return acos(vdot(a, b) / (vlen(a) * vlen(b)))
def dvtheta(a, b): # find the angle between two vectors in degrees
    return deg(vtheta(a, b))

#############
### Lists ###
#############

# Given a list containing some combination of (possibly deeply nested) Iterables
# and non-iterables, produce a single list of non-iterables. No guarentees are
# made regarding the order of the output values with respect to the input
# structure.
def flatten_list(*x):
    n = x
    m = []
    is_flat = False
    while not is_flat:
        is_flat = True
        for i in n:
            if isinstance(i, abc.Iterable):
                is_flat = False
                for j in i:
                    m.append(j)
            else:
                m.append(i)
        n = m
        m = []
    return n

# Flatten using `flatten_list()` and cast all values to `dtype`.
to_type_list = lambda dtype, *x : [dtype(i) for i in flatten_list(x)]

# Integer sum of a list
isum = lambda *x : sum(to_type_list(int, x))

# Arithmetic mean of a list
def mean(*x):
    n = flatten_list(x)
    return fsum(n) / len(n)

# Population Standard Deviation of a list
def stdDev(*x):
    n = to_type_list(float, x)
    avg = mean(n)
    total_deviation = 0
    for i in n:
        total_deviation += (i - avg) ** 2
    return sqrt(1 / len(n) * total_deviation)

# %RSD of a list
def pctRSD(*x):
    try:
        return stdDev(x) / mean(x) * 100
    except ZeroDivisionError:
        return float('NaN')

################################
### Cellular Data Statistics ###
################################

# Given an integer from 1 to 12 (inclusive) representing a month, or the name
# of a month return the number of days in that month, ignoring leap years.
def days_in_month(month):
    month_names = { 'jan': 1, 'feb': 2, 'mar': 3, 'apr': 4, 'may': 5, 'jun': 6,
            'jul': 7, 'aug': 8, 'sep': 9, 'oct': 10, 'nov': 11, 'dec': 12}
    shortmonths = [4, 6, 9, 11]
    try:
        month_number = int(month)
    except ValueError:
        try:
            month_number = month_names[str(month).lower()[:3]]
        except KeyError:
            month_number = 0
    if month_number in shortmonths:
        return 30
    elif month_number == 2:
        return 28
    else: # For simplicity, assume 31 days if input is invalid.
        return 31


# Given the current amount of data used (in Gigabytes), and the total data
# allowance for each month (also in Gigabytes), calculate statistics for how
# much data should be used to yield a uniform usage pattern throughout the
# month.
def data(gb, total, reset_day=11):
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

    print(
        f'     Cycle Usage: {cycleUsage:d} MiB\n'
        f'     Ideal Usage: {idealUsage:d} MiB\n'
        f'       Net Usage: {netUsage:d} MiB\n'
        f'      Cycle Rate: {cycleRate:d} MiB/day\n'
        f'      Ideal Rate: {idealRate:d} MiB/day\n'
        f'        Net Rate: {netRate:d} MiB/day\n'
        f' Use Coefficient: {coefficient:f}\n'
        f'       Cycle Day: {cycleDay:d} / {totalDays:d}\n'
        f'       Ideal Day: {daysUsed:d}')

    if netUsage < 0:
        daysBehind = -netUsage / idealRate + 1
        print(f'        Catch up: {daysBehind:d}')

######################
### Exit functions ###
######################

# Alternatives for sys.exit.
quit = exit
bye = exit

print(f'PyDesk, version {script_version}')
