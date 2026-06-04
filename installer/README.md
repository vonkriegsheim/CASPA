# CASPA Windows installer

A no-admin bootstrapper `.exe` (built with [Inno Setup](https://jrsoftware.org/isinfo.php)).
On double-click it extracts CASPA and then `caspa-bootstrap.ps1` installs a **private
miniforge (Python) and native CRAN R inside the install folder** plus all CASPA
Python + R/Bioconductor packages — nothing is added to the system or PATH, and no
admin rights are needed. It then adds Start-Menu shortcuts.

> The R/Bioconductor stack can't be bundled (no win-64 bioconda), so the installer
> **downloads** miniforge + R + packages at install time: needs internet, pulls
> several GB, and takes ~20–40 minutes (the Bioconductor build is the slow part).

## What the user gets

Start Menu → **CASPA**:
- **CASPA Setup (GUI)** — the form-based config builder (`caspa gui`).
- **CASPA Console** — a console with the bundled Python/R on PATH, to run
  `python caspa\run.py --workdir … --cores N`.
- **CASPA Doctor** — verifies the install.

Installs to `%LOCALAPPDATA%\CASPA` (per-user, no admin). Uninstall removes it,
including the bundled `miniforge3\` and `R\`.

## Building it

CI does this automatically — the **build-installer** workflow compiles
`caspa-setup.iss` on a Windows runner and uploads **`caspa-setup.exe`** as an
artifact (and attaches it to GitHub Releases). Trigger it by pushing the
`installer` branch or via the Actions tab.

Locally (needs Inno Setup 6):

```
ISCC.exe installer\caspa-setup.iss      ->  installer\Output\caspa-setup.exe
```

## Files

| File | Role |
|------|------|
| `caspa-setup.iss` | Inno Setup script (layout, shortcuts, runs the bootstrap) |
| `caspa-bootstrap.ps1` | Installs miniforge + R + packages into the install folder |
| `CASPA-GUI.cmd` / `CASPA-Console.cmd` / `CASPA-Doctor.cmd` | Launchers (bundled env on PATH) |
