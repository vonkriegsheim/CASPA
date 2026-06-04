@echo off
rem Open a console with the bundled CASPA Python/R on PATH, for running the pipeline.
set "ROOT=%~dp0"
set "PATH=%ROOT%miniforge3;%ROOT%miniforge3\Library\bin;%ROOT%miniforge3\Scripts;%PATH%"
if exist "%ROOT%R\bin\Rscript.exe" set "PATH=%ROOT%R\bin;%PATH%"
cd /d "%ROOT%"
echo CASPA environment ready. Examples:
echo.
echo   python caspa\init.py --workdir %%USERPROFILE%%\MyExperiment --pg-matrix C:\data\report.pg_matrix.tsv --species human
echo   python caspa\run.py  --workdir %%USERPROFILE%%\MyExperiment --cores 8
echo.
echo (or use the "CASPA Setup (GUI)" shortcut to build the config by form.)
echo.
cmd /k
