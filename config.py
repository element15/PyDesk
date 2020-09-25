"""Configuration for PyDesk"""

# Maximum allowable for fraction denominators. This prevents unreasonably
# large values from resulting when converting a float to a fraction.
max_frac_denom = 1024

# Executable path for GNU Units. If the executable is in your $PATH, then
# a full filepath need not be specified.
gnu_units_executable = 'gunits'

# The day on which rationing statistics reset.
cellular_reset_day = 11