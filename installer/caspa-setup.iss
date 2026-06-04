; CASPA Windows installer (Inno Setup).
; A no-admin bootstrapper: extracts CASPA, then caspa-bootstrap.ps1 installs a
; private miniforge + native R *inside the install folder* and all packages.
; Compile with: ISCC.exe caspa-setup.iss  (CI does this on a Windows runner).

#define MyAppName "CASPA"
#define MyAppVersion "0.1.0"
#define RepoRoot ".."

[Setup]
AppId={{A1B2C3D4-E5F6-47A8-9B0C-1D2E3F405162}}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher=von Kriegsheim lab, University of Edinburgh
DefaultDirName={localappdata}\CASPA
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=Output
OutputBaseFilename=caspa-setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible

[Files]
; CASPA source tree -> {app}
Source: "{#RepoRoot}\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion; \
  Excludes: "\.git\*,\.github\*,\installer\Output\*,\installer\caspa-setup.iss,*.pyc,*.pyo,__pycache__\*,\.snakemake\*,_*.log,*.exe"
; launchers at the install root (so %~dp0 resolves to {app})
Source: "CASPA-GUI.cmd";     DestDir: "{app}"; Flags: ignoreversion
Source: "CASPA-Console.cmd"; DestDir: "{app}"; Flags: ignoreversion
Source: "CASPA-Doctor.cmd";  DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\CASPA\CASPA Setup (GUI)"; Filename: "{app}\CASPA-GUI.cmd";     WorkingDir: "{app}"
Name: "{autoprograms}\CASPA\CASPA Console";     Filename: "{app}\CASPA-Console.cmd"; WorkingDir: "{app}"
Name: "{autoprograms}\CASPA\CASPA Doctor";      Filename: "{app}\CASPA-Doctor.cmd";  WorkingDir: "{app}"

[Run]
; Install Python + R + all packages into {app}. Console is visible so the user
; sees progress; Inno waits for it to finish.
Filename: "powershell.exe"; \
  Parameters: "-NoProfile -ExecutionPolicy Bypass -File ""{app}\installer\caspa-bootstrap.ps1"" -Root ""{app}"""; \
  StatusMsg: "Installing Python and R dependencies (downloads several GB; can take 20-40 minutes)..."; \
  Flags: waituntilterminated

[UninstallDelete]
Type: filesandordirs; Name: "{app}\miniforge3"
Type: filesandordirs; Name: "{app}\R"
