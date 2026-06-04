<#
  CASPA installer bootstrap.

  Installs a private miniforge (Python) and native CRAN R *inside the install
  folder*, then the CASPA Python + R/Bioconductor packages. Self-contained: no
  admin rights, nothing added to the system PATH, nothing else on the machine is
  touched. Run by the Inno Setup installer after the files are extracted.
#>
param([string]$Root = (Split-Path -Parent $PSScriptRoot))

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"   # much faster Invoke-WebRequest downloads

function Say($m) { Write-Host "[caspa-setup] $m" -ForegroundColor Cyan }
function Die($m) { Write-Host "[caspa-setup] ERROR: $m" -ForegroundColor Red; exit 1 }

$mf      = Join-Path $Root "miniforge3"
$py      = Join-Path $mf   "python.exe"
$rdir    = Join-Path $Root "R"
$tmp     = Join-Path $env:TEMP "caspa-setup"
New-Item -ItemType Directory -Force -Path $tmp | Out-Null

Say "Installing CASPA into: $Root"

# --- 1. Miniforge (Python) -------------------------------------------------
if (-not (Test-Path $py)) {
    try {
        Say "Downloading Miniforge (Python)..."
        $mfexe = Join-Path $tmp "Miniforge3.exe"
        Invoke-WebRequest "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe" -OutFile $mfexe
        Say "Installing Miniforge to $mf ..."
        # NSIS silent install; /D must be last and unquoted.
        Start-Process -Wait -FilePath $mfexe -ArgumentList "/InstallationType=JustMe", "/RegisterPython=0", "/AddToPath=0", "/S", "/D=$mf"
    } catch { Die "Miniforge download/install failed: $($_.Exception.Message)" }
    if (-not (Test-Path $py)) { Die "Miniforge install did not produce $py" }
} else { Say "Miniforge already present." }

# --- 2. Native CRAN R (bioconda has no win-64 Bioconductor) -----------------
$rscript = (Get-ChildItem -Path $rdir -Recurse -Filter Rscript.exe -ErrorAction SilentlyContinue | Select-Object -First 1)
if (-not $rscript) {
    try {
        Say "Finding the latest R release on CRAN..."
        $base = "https://cran.r-project.org/bin/windows/base/"
        $page = Invoke-WebRequest $base -UseBasicParsing
        $rel  = ($page.Links | Where-Object { $_.href -match 'R-\d+\.\d+\.\d+-win\.exe' } | Select-Object -First 1).href
        if (-not $rel) { Die "Could not find the R installer link on CRAN." }
        Say "Downloading $rel ..."
        $rexe = Join-Path $tmp "R-win.exe"
        Invoke-WebRequest ($base + $rel) -OutFile $rexe
        Say "Installing R to $rdir ..."
        Start-Process -Wait -FilePath $rexe -ArgumentList "/VERYSILENT", "/SUPPRESSMSGBOXES", "/DIR=$rdir"
    } catch { Die "R download/install failed: $($_.Exception.Message)" }
    $rscript = (Get-ChildItem -Path $rdir -Recurse -Filter Rscript.exe -ErrorAction SilentlyContinue | Select-Object -First 1)
    if (-not $rscript) { Die "R install did not produce Rscript.exe under $rdir" }
} else { Say "R already present." }
$rscript = $rscript.FullName

# --- 3. Python packages ----------------------------------------------------
Say "Installing Python packages (pip) - a few minutes..."
& $py -m pip install --no-input --no-warn-script-location -r (Join-Path $Root "requirements-windows.txt")
if ($LASTEXITCODE -ne 0) { Die "pip install failed." }

# --- 4. R / Bioconductor packages ------------------------------------------
Say "Installing R/Bioconductor packages (BiocManager) - this can take 20-40 minutes..."
& $rscript (Join-Path $Root "pipeline\scripts\R\install_r_packages.R")
if ($LASTEXITCODE -ne 0) { Die "R/Bioconductor package install reported missing packages." }

# --- 5. Verify -------------------------------------------------------------
Say "Verifying the installation..."
& $py (Join-Path $Root "caspa\doctor.py") --rscript $rscript

Say "CASPA install complete. Launch it from the Start Menu: 'CASPA Setup (GUI)'."
