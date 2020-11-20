#!/usr/bin/env python3 -i

script_version = "20-925c"
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
import cmath
import random
import re
import string
import subprocess

import numpy as np
import pyproj

import config

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
c = 299_792_458 # Speed of light in a vacuum (m/s)

# Quantum Physics
h = 6.626e-34 # Planck constant (J*s)
hbar = 1.054e-34 # h-bar constant (J*s)
m_e = 9.109e-31 # Electron mass (kg)
m_p = 1.673e-27 # Proton mass (kg)
m_n = 1.675e-27 # Neutron mass (kg)
eV = 1.602e-19 # Electron-volt (J)

# Coordinate transformations
EPSG_ECEF = 4978 # ID for Earth Centered, Earth Fixed coordinates
EPSG_LLA = 4326 # ID for Latitude, Longitude, Altitude coordinates
earth_dia = 12_742_018 # Mean diameter of the earth (m)

#############################
### Convenience Functions ###
#############################

ln = log

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

exp10 = lambda x, y : x * 10**y

si_prefixes = {
    'Y': 24, 'Z': 21, 'E': 18, 'P': 15, 'T': 12,
    'G': 9, 'M': 6, 'k': 3, 'h': 2, 'da': 1,
    'd': -1, 'c': -2, 'm': -3, 'µ': -6, 'n': -9,
    'Å': -10, 'p': -12, 'f': -15, 'a': -18, 'z': -21, 'y': -24,
    }
si_prefix_names = {
    'yotta': 'Y', 'zetta': 'Z', 'exa': 'E', 'peta': 'P', 'tera': 'T',
    'giga': 'G', 'mega': 'M', 'kilo': 'k', 'hecto': 'h', 'deca': 'da',
    'deci': 'd', 'centi': 'c', 'milli': 'm', 'micro': 'µ', 'nano': 'n',
    'angstrom': 'Å', 'pico': 'p', 'femto': 'f', 'atto': 'a', 'zepto': 'z',
    'yocto': 'y',
    }
def to_prefix(val, prefix):
    """Convert a base unit to an SI prefix."""
    if prefix.lower() in si_prefix_names:
        prefix = si_prefix_names[prefix.lower()]
    return exp10(val, -si_prefixes[prefix])
def from_prefix(val, prefix):
    """Convert an SI prefix to base units."""
    if prefix.lower() in si_prefix_names:
        prefix = si_prefix_names[prefix.lower()]
    return exp10(val, si_prefixes[prefix])

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

def units(v, a, b):
    """Convert units using GNU Units.

    Parameters:
        - v: Value to convert
        - a: Current unit
        - b: Desired unit

    Returns:
        Tuple containing the direct conversion ('*') and the recriprocal
        conversion ('/'), respectively. If any errors occur, the output of
        GNU Units will be returned directly as a string.
    """
    result = subprocess.run([config.gnu_units_executable, f'{v}{a}', b],
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

def uu(v, a, b, silent=False):
    """Convert using `units()` and unpack the resulting tuple.

    Parameters:
        - v: Value to convert
        - a: Current unit
        - b: Desired unit
        - silent: Do not print the full conversion result to stdout
            (default: False)

    Returns:
        Direct conversion ('*'). If any errors occur, None is returned.
    """
    conv = units(v, a, b)
    if not silent:
        print(conv)
    if not isinstance(conv, tuple) or conv[0] == 'conformability error':
        return None
    return conv[0]

def temp(t, scales):
    """Convert temperature `t` between scales specified by `scales`."""
    conv = { # (absolute_zero, kelvin_per_degree)
        'f': (-459.67, 5/9),
        'c': (-273.15, 1),
        'r': (0, 5/9),
        'k': (0, 1),
    }
    s1, s2 = scales.lower()
    z1, f1 = conv[s1]
    z2, f2 = conv[s2]
    t_kelvin = (t-z1)*f1
    if t_kelvin < 0:
        print('Result is below absolute zero')
        return float('NaN')
    t_out = t_kelvin/f2 + z2
    return round(t_out, 8)

####################
### General Math ###
####################

quad_det = lambda a, b, c : b**2 - 4*a*c
def quad(a, b, c):
    """Quadratic formula for ax**2 + bx + c = 0."""
    soln_mean = -b/(2*a)
    soln_radius = cmath.sqrt(quad_det(a, b, c))/(2*a)
    return soln_mean+soln_radius, soln_mean-soln_radius

mid = lambda a, b : (a+b)/2
mid2d = lambda x1, y1, x2, y2 : (mid(x1, x2), mid(y1, y2))
dist2d = lambda x1, y1, x2, y2 : ((x2-x1)**2 + (y2-y1)**2)**0.5

# Linear interpolation
def lint(x1, xn, x2, y1, y2):
    """Linear interpolate points.

    Calculate some point `yn` between `y1` and `y2` which is proportional to
    the point `xn` between `x1` and `x2`.
    """
    return (y2-y1)/(x2-x1)*(xn-x1) + y1

# Pythagorean leg
leg = lambda a, c : sqrt(abs(c**2 - a**2))

######################
### Thermodynamics ###
######################

def heat_cp_mol(T, *coeff):
    """Enthalpy specific heat on a mole basis."""
    coeff = flatten_list(coeff)
    c_p = 0;
    for i in range(0, 3):
        c_p += coeff[i] * T**i
    return c_p

def heat_cv_mol(R, T, *coeff):
    """Internal energy specific heat on a mole basis."""
    return heat_cp_mol(T, coeff) - R

def heat_h_mol(T, *coeff):
    """Enthalpy on a mole basis."""
    coeff = flatten_list(coeff)
    h = 0
    for i in range(0, 3):
        h += coeff[i] * T**(i + 1) / (i + 1)
    return h

def heat_u_mol(R, T, *coeff):
    """Internal energy on a mole basis."""
    coeff = flatten_list(coeff)
    h =  heat_h_mol(T, coeff)
    return h - R * T

#########################
### Number Formatting ###
#########################

getfrac = lambda x, denom=config.max_frac_denom : (
    Fraction(x).limit_denominator(denom))
frac = lambda x, denom=config.max_frac_denom : print(getfrac(x, denom))

def mix(x, denom=config.max_frac_denom):
    fraction = getfrac(x, denom)
    numerator = fraction.numerator
    denominator = fraction.denominator
    if numerator > denominator:
        whole = floor(x)
        mixed_numerator = numerator - whole * denominator
        print("%d %d/%d" % (whole, mixed_numerator, denominator))
    else:
        print("%d/%d" % (numerator, denominator))

def sci(x, sigfig=6):
    """Format a number in scientific notation."""
    sigfig
    if sigfig < 1:
        sigfig = 1
    string = '{:.' + str(sigfig - 1) + 'e}'
    return string.format(x)

def eng(x, sigfig=6):
    """Format a number in engineering notation."""
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
    """Convert a number to an arbitrary base."""
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

def to_dms(n):
    """Convert a decimal degree value to (degrees, minutes, seconds)."""
    hemisphere = copysign(1, n)
    n_pos = abs(n)
    degrees = int(n_pos)
    minutes_full = 60 * (n_pos-degrees)
    minutes = int(minutes_full)
    seconds = 60 * (minutes_full-minutes)
    return copysign(degrees, hemisphere), minutes, seconds
def from_dms(n):
    """Convert a (degrees, minutes, seconds) value to decimal degrees."""
    deg = abs(n[0])
    dec = deg + n[1]/60 + n[2]/3600
    return copysign(dec, n[0])
def pretty_dms(lat, lon):
    """Convert a pair of decimal degree values to a pretty DMS string."""
    latd, latm, lats = to_dms(lat)
    lond, lonm, lons = to_dms(lon)
    ns_hemisphere = 'N' if lat >= 0 else 'S'
    ew_hemisphere = 'E' if lon >= 0 else 'W'
    return (
    	f'{abs(latd):.0f}˚ {latm:.0f}\' {lats:.3f}" {ns_hemisphere}, '
        f'{abs(lond):.0f}˚ {lonm:.0f}\' {lons:.3f}" {ew_hemisphere}')

ecef_lla = pyproj.Transformer.from_crs(EPSG_ECEF, EPSG_LLA).transform
lla_ecef = pyproj.Transformer.from_crs(EPSG_LLA, EPSG_ECEF).transform

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

def flatten_list(*x):
    """Flatten nested lists.

    No guarentees are made regarding the order of the output values with
    respect to the input structure.
    """
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

def mean(*x):
    """Arithmetic mean."""
    n = flatten_list(x)
    return fsum(n) / len(n)

def stdDev(*x):
    """Population standard deviation."""
    n = to_type_list(float, x)
    avg = mean(n)
    total_deviation = 0
    for i in n:
        total_deviation += (i - avg) ** 2
    return sqrt(1 / len(n) * total_deviation)

def pctRSD(*x):
    """Percent relative standard deviation (using population stdDev)."""
    try:
        return stdDev(x) / mean(x) * 100
    except ZeroDivisionError:
        return float('NaN')

################################
### Cellular Data Statistics ###
################################

def days_in_month(month):
    """Number of days in the month (leap years notwithstanding)."""
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

def data(gb, total, reset_day=config.cellular_reset_day):
    """Ration limited cellular data.

    Given the current amount of data used (in GiB) and the total data allowance
    for each month (in GiB), calculate statistics for how much data should be
    used to yield a uniform usage pattern throughout the month.
    """
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
