# PyDesk

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
