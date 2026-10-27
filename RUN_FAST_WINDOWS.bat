@echo off
cd /d "%~dp0"
py -3 run.py --quality quick --figures all --workers 4
if errorlevel 1 echo Program ended with errors. See output above.
pause
