#!/usr/bin/osascript

# launch.sh
# 
# This script allows access to a Python interpreter window without the need to
# create and maintain multiple Terminal profiles. Modification of the 'target'
# identifier may be necessary. Once the script suits your needs, you may use
# Automator to create a thin wrapper application to launch PyDesk form the Dock.
# This script may be run as a native AppleScript or as a UNIX shell script.


# Modify here if necessary
set target to "~/git/PyDesk/calc.py"

tell application "Terminal"
    do script "clear; " & target & "; exit"
    activate
end tell
