@echo off
rem Launch the CASPA setup GUI using the bundled Python/R (no system PATH needed).
setlocal
set "ROOT=%~dp0"
set "PATH=%ROOT%miniforge3;%ROOT%miniforge3\Library\bin;%ROOT%miniforge3\Scripts;%PATH%"
if exist "%ROOT%R\bin\Rscript.exe" set "PATH=%ROOT%R\bin;%PATH%"
"%ROOT%miniforge3\python.exe" "%ROOT%caspa\gui.py"
