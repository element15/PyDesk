#!/usr/bin/osascript

# launch.sh
# 
# This script allows access to a Python interpreter window without the need to
# create and maintain multiple Terminal profiles. Modification of the 'target'
# identifier may be necessary. Once the script suits your needs, you may use
# Automator to create a thin wrapper application to launch PyDesk form the Dock.
# This script may be run as a native AppleScript or as a UNIX shell script.


# Equivilant of '~'
set home to (the POSIX path of (path to home folder))

# Modify here if necessary. Points to ~/git/PyDesk/calc.py by default.
set target to (home & "/git/PyDesk/calc.py")

try
    python_script
on error
    set python_script to POSIX file target as alias
end try

tell application "Terminal" to open python_script
