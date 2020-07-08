#!/usr/bin/env python3 -i

script_version = "20-708a"
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
rad = radians
deg = degrees
to_rad = lambda f : lambda x : f(rad(x))
to_deg = lambda f : lambda x : deg(f(x))

sind = to_rad(sin)
cosd = to_rad(cos)
tand = to_rad(tan)
asind = to_deg(asin)
acosd = to_deg(acos)
atand = to_deg(atan)
sinhd = to_rad(sinh)
coshd = to_rad(cosh)
asinhd = to_deg(asinh)
acoshd = to_deg(acosh)
atanhd = to_deg(atanh)

# Convert base unit to SI prefix
exp10 = lambda x, y : x * 10**y

to_exp = lambda n : lambda x : exp10(x, n)

to_yotta = to_exp(-24)
to_zetta = to_exp(-21)
to_exa = to_exp(-18)
to_peta = to_exp(-15)
to_tera = to_exp(-12)
to_giga = to_exp(-9)
to_mega = to_exp(-6)
to_kilo = to_exp(-3)
to_centi = to_exp(2)
to_milli = to_exp(3)
to_micro = to_exp(6)
to_nano = to_exp(9)
to_angstrom = to_exp(10)
to_pico = to_exp(12)
to_femto = to_exp(15)
to_atto = to_exp(18)
to_zepto = to_exp(21)
to_yocto = to_exp(24)

# Convert SI prefix to base unit

from_yotta = to_exp(24)
from_zetta = to_exp(21)
from_exa = to_exp(18)
from_peta = to_exp(15)
from_tera = to_exp(12)
from_giga = to_exp(9)
from_mega = to_exp(6)
from_kilo = to_exp(3)
from_centi = to_exp(-2)
from_milli = to_exp(-3)
from_micro = to_exp(-6)
from_nano = to_exp(-9)
from_angstrom = to_exp(-10)
from_pico = to_exp(-12)
from_femto = to_exp(-15)
from_atto = to_exp(-18)
from_zepto = to_exp(-21)
from_yocto = to_exp(-24)

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

# Convert a decimal value to (degrees, minutes, seconds)
def to_dms(n):
    hemisphere = copysign(1, n)
    n_pos = abs(n)
    degrees = int(n_pos)
    minutes_full = 60 * (n_pos-degrees)
    minutes = int(min_full)
    seconds = 60 * (minutes_full-minutes)
    return copysign(degrees, hemisphere), minutes, seconds
# Convert a (degrees, minutes, seconds) value to decimal
def from_dms(n):
    deg = abs(n[0])
    dec = deg + n[1]/60 + n[2]/3600
    return copysign(dec, n[0])
# Convert a pair of decimal degree values to a pretty DMS string
def pretty_dms(lat, lon):
    latd, latm, lats = to_dms(lat)
    lond, lonm, lons = to_dms(lon)
    ns_hemisphere = 'N' if lat >= 0 else 'S'
    ew_hemisphere = 'E' if lon >= 0 else 'W'
    return (
    	f'{abs(latd):.0f}˚ {latm:.0f}\' {lats:.3f}" {ns_hemisphere}, '
        f'{abs(lond):.0f}˚ {lonm:.0f}\' {lons:.3f}" {ew_hemisphere}')

###############################
### Temperature Conversions ###
###############################

# Absolute zero checks
abs_zero = {
    'f': -459.67,
    'c': -273.15,
    'k': 0,
    'r': 0,
}
def temp_check_zero(temp, scale):
    if temp < abs_zero[scale]:
        print('Result is below absolute zero')
        return float('NaN')
    return round(temp, 8)

temp_fc = lambda f : temp_check_zero((f-32)/1.8, 'c')
temp_cf = lambda c : temp_check_zero(c*1.8 + 32, 'f')
temp_ck = lambda c : temp_check_zero(c+273.15, 'k')
temp_kc = lambda k : temp_check_zero(k-273.15, 'c')
temp_fk = lambda f : temp_check_zero((f+459.67)/1.8, 'k')
temp_kf = lambda k : temp_check_zero(k*1.8 - 459.67, 'f')
temp_fr = lambda f : temp_check_zero(f+459.67, 'r')
temp_rf = lambda r : temp_check_zero(r-459.67, 'f')
temp_kr = lambda k : temp_check_zero(k*1.8, 'r')
temp_rk = lambda r : temp_check_zero(r/1.8, 'k')
temp_cr = lambda c : temp_check_zero((c+273.15)*1.8, 'r')
temp_rc = lambda r : temp_check_zero(r/1.8 - 273.15, 'c')

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
    return tuple((list(a) + [0] * 3)[:3])

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
