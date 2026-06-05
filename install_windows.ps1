<#
.SYNOPSIS
  Install CASPA's Python and R dependencies on Windows, then verify.

.DESCRIPTION
  bioconda has no win-64 Bioconductor builds, so on Windows CASPA uses miniforge
  Python (pip) + native CRAN R (BiocManager). This script installs both and runs
  caspa\doctor.py. It does NOT install miniforge or R themselves — install those
  first (links printed if missing).

.PARAMETER PythonExe
  Path to miniforge/conda python.exe. Auto-detected if omitted.
.PARAMETER RscriptExe
  Path to Rscript.exe. Auto-detected if omitted.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File .\install_windows.ps1
#>
param(
  [string]$PythonExe,
  [string]$RscriptExe
)

$ErrorActionPreference = "Stop"
$repo = $PSScriptRoot

function Find-Python {
  if ($PythonExe -and (Test-Path $PythonExe)) { return $PythonExe }
  $cands = @(
    "$env:USERPROFILE\miniforge3\python.exe",
    "$env:USERPROFILE\mambaforge\python.exe",
    "$env:USERPROFILE\miniconda3\python.exe",
    "$env:USERPROFILE\anaconda3\python.exe",
    "C:\ProgramData\miniforge3\python.exe"
  )
  foreach ($c in $cands) { if (Test-Path $c) { return $c } }
  $cmd = Get-Command python -ErrorAction SilentlyContinue
  if ($cmd) { return $cmd.Source }
  return $null
}

function Find-Rscript {
  if ($RscriptExe -and (Test-Path $RscriptExe)) { return $RscriptExe }
  $cmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue
  if ($cmd) { return $cmd.Source }
  $rdirs = Get-ChildItem "C:\Program Files\R" -Directory -ErrorAction SilentlyContinue |
           Where-Object { Test-Path "$($_.FullName)\bin\Rscript.exe" } |
           Sort-Object Name -Descending
  if ($rdirs) { return "$($rdirs[0].FullName)\bin\Rscript.exe" }
  return $null
}

Write-Host "=== CASPA Windows installer ===" -ForegroundColor Cyan

$py = Find-Python
if (-not $py) {
  Write-Error ("No miniforge/conda Python found. Install miniforge from " +
    "https://github.com/conda-forge/miniforge/releases then re-run this script.")
  exit 1
}
Write-Host "Python : $py"

$rs = Find-Rscript
if (-not $rs) {
  Write-Error ("No R found. Install the pinned R 4.5.2 (Bioconductor 3.21 - the tested combo) from " +
    "https://cran.r-project.org/bin/windows/base/old/4.5.2/ then re-run this script.")
  exit 1
}
Write-Host "Rscript: $rs"

Write-Host "`n[1/3] Installing Python packages (pip)..." -ForegroundColor Cyan
& $py -m pip install -r "$repo\requirements-windows.txt"
if ($LASTEXITCODE -ne 0) { Write-Error "pip install failed."; exit 1 }

Write-Host "`n[2/3] Installing R / Bioconductor packages (BiocManager)..." -ForegroundColor Cyan
& $rs "$repo\pipeline\scripts\R\install_r_packages.R"
if ($LASTEXITCODE -ne 0) { Write-Error "R package install reported missing packages."; exit 1 }

# Pin the R that got the packages, so run.py uses exactly this one even if the
# machine has other CRAN R versions (run.py honours CASPA_RSCRIPT first).
Write-Host "`nPinning Rscript via CASPA_RSCRIPT (user env): $rs"
[Environment]::SetEnvironmentVariable("CASPA_RSCRIPT", $rs, "User")
$env:CASPA_RSCRIPT = $rs   # also for this session's doctor check

Write-Host "`n[3/3] Verifying (caspa doctor)..." -ForegroundColor Cyan
& $py "$repo\caspa\doctor.py" --rscript $rs
$doctor = $LASTEXITCODE

$pyDir = Split-Path $py
$rsDir = Split-Path $rs
Write-Host "`n=== Done ===" -ForegroundColor Cyan
Write-Host "Add these to PATH so python/snakemake/Rscript resolve when you run the pipeline:"
Write-Host "  $pyDir"
Write-Host "  $pyDir\Scripts"
Write-Host "  $rsDir"
Write-Host "`nThen run:" -ForegroundColor Cyan
Write-Host "  python caspa\run.py --workdir <workdir> --cores N"
exit $doctor
