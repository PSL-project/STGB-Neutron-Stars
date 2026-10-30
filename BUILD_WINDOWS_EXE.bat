@echo off
setlocal
cd /d "%~dp0"

echo [1/4] Installing build tools...
py -m pip install --upgrade pyinstaller pillow
if errorlevel 1 goto :fail

echo [2/4] Generating neutron-star icon...
py tools\generate_neutron_star_icon.py
if errorlevel 1 goto :fail

echo [3/4] Building STGB neutron star.exe...
pyinstaller --noconfirm --clean --onefile --windowed ^
  --name "STGB neutron star" ^
  --icon "assets\neutron_star.ico" ^
  --add-data "data;data" ^
  --paths "src" ^
  --collect-data matplotlib ^
  "gui_run.py"
if errorlevel 1 goto :fail

echo [4/4] Done.
echo Executable: dist\STGB neutron star.exe
exit /b 0

:fail
echo Build failed.
exit /b 1
