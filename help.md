# PyDesk

PyDesk is a library of functions that can be loaded into a Python interpreter for use as a genera-purpose scientific and engineering calculator.

To invoke PyDesk, run `python3 -i /path/to/calc.py`

## Command Categories
Functions in PyDesk fall into the following categories:
   * Constants
   * Fundamentals
   * Trigonometry
   * SI prefixes
   * General-purpose
   * Number formatting
   * Temperature
   * Fractions
   * Vectors
   * Lists
   * Cellular data statistics
   * Exit
   * Help

## Constants
PyDesk includes several mathematical and physical constants which may be helpful for calculations. These constants and their units are outlined below.

| Variable | Value | Units | Description |
| --- | --- | --- | --- |
| pi | 3.141592653589793 | [N/A] | The ratio between a circle's circumference and its diameter |
| e | 2.718281828459045 | [N/A] | Euler's constant |
| r_atm | 0.0820574614 | L atm/(mol K) | The universal gas constant |
| r_mmhg | 62.3636711 | L mmHg/(mol K) | The universal gas constant |
| r_joule | 8.314462175 | J/(mol K) | The universal gas constant |
| kw | 1.01e-14 | [N/A] | The equilibrium constant for the auto-ionization of water |
| avo | 6.022e23 | 1/mol | Avogardo's constant |
| g | 9.80665 | m/s | Acceleration due to earth's gravity |
| g_ft | 32.174049 | ft/s | Acceleration due to earth's gravity |
| G | 6.673e-11 | N m<sup>2</sup>/kg<sup>2</sup> | Gravitational constant |
| ele | 1.602e-19 | C | Elementary charge |
| mu_0 | pi / 2.5e+6 | N/A<sup>2</sup> | Permeability of a vacuum |
| epsilon_0 | 8.854e-12 | F/m | Permittivity of a vacuum |
| k_e | 8.988e9 | N m<sup>2</sup>/C<sup>2</sup> | Coulomb's constant |
| c | 2.998e8 | m/s | Speed of light in a vacuum |
| h | 6.626e-34 | J s | Planck's constant |
| hbar | 1.054e-34 | J s | H-bar constant |
| m_e | 9.109e-31 | kg | Electron mass |
| m_p | 1.673e-27 | kg | Proton mass |
| m_n | 1.675e-27 | kg | Neutron mass |
| eV | 1.602e-19 | J | Electron-Volt |

## Fundamentals
The following fundamental mathematical functions delegate to built-in Python functions.

| Function | Arguments | Returns | Description |
| --- | --- | --- | --- |
| ln(x) | float x | float | Natural logarithm of x |
| log(x) | float x | float | Natural logarithm of x |
| log10(x) | float x | float | Base-10 logarithm of x |
| logbase(x, y) | float x, float y | float | Base-y logarithm of x |
| exp(x) | float x | float | e<sup>x</sup> |
| exp10(x) | float x | float | 10<sup>x</sup> |
| sqrt(x) | float x | float | Square root of x |
| nrt(x, y) | float x, float y | float | x<sup>1/y</sup> |
| rec(x) | float x | float | Recripocal of x |
| fact(x) | int x | int | Factorial of x |
| gamma(x) | float x | float | Gamma function of x |
| hypot(x, y) | float x, float y | float | Euclidian norm of x and y, or (x<sup>2</sup> + y<sup>2</sup>)<sup>1/2</sup> |
| floor(x) | float x | int | Round x down |
| ceil(x) | float x | int | Round x up |

## Trigonometry
The following functions delegate to built-in Python functions.

### Angle Conversion
| Function | Arguments | Returns | Description |
| --- | --- | --- | --- |
| rad(x) | float x | float | Convert angle x (degrees) to radians |
| deg(x) | float x | float | Convert angle x (radians) to degrees |

### Trigonometric Functions
All of the following functions take arguments in radians. For any of them, a degree version can be obtained by adding a `d` to the front of the function name. For example, `sin()` becomes `dsin()` and `atan()` becomes `datan()`.

| Function | Arguments | Returns | Description |
| --- | --- | --- | --- |
| sin(x) | float x | float | Sine of x |
| cos(x) | float x | float | Cosine of x |
| sec(x) | float x | float | Secant of x |
| csc(x) | float x | float | Cosecant of x |
| tan(x) | float x | float | Tangent of x |
| cot(x) | float x | float | Cotangent of x |
| asin(x) | float x | float | Arcsine of x |
| acos(x) | float x | float | Arccosine of x |
| atan(x) | float x | float | Arctangent of x |
| sinh(x) | float x | float | Hyperbolic sine of x |
| cosh(x) | float x | float | Hyperbolic cosine of x |
| tanh(x) | float x | float | Hyperbolic tangent of x |
| asinh(x) | float x | float | Hyperbolic arcsine of x |
| acosh(x) | float x | float | Hyperbolic arccosine of x |
| atanh(x) | float x | float | Hyperbolic arctangent of x |

## SI Prefixes
These functions convert to a base unit from the specified SI prefix. To convert to prefixed units from base units, substitute `to` for `from`. For example, `from_kilo()` becomes `to_kilo()`.

| Function | Arguments | Returns | Description |
| --- | --- | --- | --- |
| from_yotta(x) | float x | float | x \* 10<sup>24</sup> |
| from_zetta(x) | float x | float | x \* 10<sup>21</sup> |
| from_exa(x) | float x | float | x \* 10<sup>18</sup> |
| from_peta(x) | float x | float | x \* 10<sup>15</sup> |
| from_tera(x) | float x | float | x \* 10<sup>12</sup> |
| from_giga(x) | float x | float | x \* 10<sup>9</sup> |
| from_mega(x) | float x | float | x \* 10<sup>6</sup> |
| from_kilo(x) | float x | float | x \* 10<sup>3</sup> |
| from_centi(x) | float x | float | x \* 10<sup>-2</sup> |
| from_milli(x) | float x | float | x \* 10<sup>-3</sup> |
| from_micro(x) | float x | float | x \* 10<sup>-6</sup> |
| from_nano(x) | float x | float | x \* 10<sup>-9</sup> |
| from_angstrom(x) | float x | float | x \* 10<sup>-10</sup> |
| from_pico(x) | float x | float | x \* 10<sup>-12</sup> |
| from_femto(x) | float x | float | x \* 10<sup>-15</sup> |
| from_atto(x) | float x | float | x \* 10<sup>-18</sup> |
| from_zepto(x) | float x | float | x \* 10<sup>-21</sup> |
| from_yocto(x) | float x | float | x \* 10<sup>-24</sup> |

## General-purpose
General math functions not included in Python's standard library

| Function | Arguments | Returns | Description |
| --- | --- | --- | --- |
| quad_det(a, b, c) | float a, float b, float c | float | Calculate the quadratic determinant (b<sup>2</sup> - 4ac) of <br /> a parabola in the form ax<sup>2</sup> + bx + c = 0 |
| quad(a, b, c) | float a, float b, float c | (float,<br />float) | Find the real roots of a parabola in the form <br /> ax<sup>2</sup> + bx + c = 0. Parabolae with complex roots will <br /> throw a domain exception. |
| mid(a, b) | float a, float b | float | Find the midpoint between two 1-dimensional points |
| mid2(x1, y1, x2, y2) | float x1, float y1, <br /> float x2, float y2 | (float,<br />float) | Find the midpoint between two 2-dimensional points |
| dist(a, b) | float a, float b | float | Find the distance from point a to point b |
| dist2(x1, y1, x2, y2) | float x1, float y1, <br /> float x2, float y2 | (float,<br />flaot) | Find the distance from point (x1, y1) to point (x2, y2) |
| lint<br />(x1, xn, x2, y1, y2) | float x1, float xn, <br /> float x2, float y1, <br /> float y2 | float | Use linear interpolation to find a point between <br /> y1 and y2, given point xn between x1 and x2 |
| pythleg(c, a) | float c, float a | float | Calculate the length of the remaining leg of a <br /> right triangle with leg a and hypotenuse c |
| diceware(n = 5) | float n | void | Generate a specified number (defaults to 5) of <br /> values for lookup in the Diceware table of words <br /> and print them to stdout |

## Number Formatting

| Function | Arguments | Returns | Description |
| --- | --- | --- | --- |
| sci(x, sigfig = 6) | float x, float sigfig | string | Format a given number in scientific notation. The <br /> value is rounded to 6 significant figures unless <br /> sigfig is specified. |
| eng(x, sigfig = 6) | float x, float sigfig | string | Format a given number in engineering notation. The <br /> value is rounded to 6 significant figures unless <br /> sigfig is specified. |
| to_base(n, b) | int n, int b | string | Convert an integer n to arbitrary base b |

## Temperature Conversions

| Function | Arguments | Returns | Description |
| --- | --- | --- | --- |
| temp_fc(f) | float f | float | Convert from degrees Fahrenheit to degrees Celsius |
| temp_cf(c) | float f | float | Convert from degrees Celsius to degrees Fahrenheit |
| temp_ck(c) | float f | float | Convert from degrees Celsius to Kelvins |
| temp_kc(k) | float f | float | Convert from Kelvins to degrees Celsius |
| temp_fk(f) | float f | float | Convert from degrees Fahrenheit to Kelvins |
| temp_kf(k) | float f | float | Convert from Kelvins to degrees Fahrenheit |
