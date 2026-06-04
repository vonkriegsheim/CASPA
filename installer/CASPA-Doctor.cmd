@echo off
rem Verify the CASPA install (Python pkgs, R pkgs, tool resolution).
setlocal
set "ROOT=%~dp0"
set "PATH=%ROOT%miniforge3;%ROOT%miniforge3\Library\bin;%ROOT%miniforge3\Scripts;%PATH%"
if exist "%ROOT%R\bin\Rscript.exe" set "PATH=%ROOT%R\bin;%PATH%"
"%ROOT%miniforge3\python.exe" "%ROOT%caspa\doctor.py"
echo.
pause
