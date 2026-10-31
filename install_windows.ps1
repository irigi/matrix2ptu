# Installs matrix2ptu on Windows into a private Python environment (.venv)
# inside this folder. Safe to run again: it updates the existing installation.
#
#   Double-click install_windows.bat, or in PowerShell:
#   powershell -ExecutionPolicy Bypass -File install_windows.ps1 [-Dev] [-NoValidate] [-Recreate]
#
#   -Dev         also install pytest (to run the tests)
#   -NoValidate  install only the converter (numpy + pyyaml), not the validation readers
#   -Recreate    delete .venv and build it from scratch

param(
    [switch]$Dev,
    [switch]$NoValidate,
    [switch]$Recreate
)

$ErrorActionPreference = "Stop"
$MinMajor = 3
$MinMinor = 12
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $Root

function Fail($msg) {
    Write-Host ""
    Write-Host "ERROR: $msg" -ForegroundColor Red
    exit 1
}

# Returns the command (as an array) of the first Python >= 3.12 found, or $null.
# "py" is the official Python launcher; "python" may be the Microsoft Store stub,
# which prints nothing and fails, so every candidate is actually run.
function Find-Python {
    $candidates = @(@("py", "-3"), @("python"), @("python3"))
    foreach ($c in $candidates) {
        $exe = $c[0]
        $pre = @($c | Select-Object -Skip 1)
        if (-not (Get-Command $exe -ErrorAction SilentlyContinue)) { continue }
        try {
            $v = & $exe @pre -c "import sys; print('%d.%d' % sys.version_info[:2])" 2>$null
        } catch { continue }
        if ($LASTEXITCODE -ne 0 -or -not $v) { continue }
        $parts = "$v".Trim().Split(".")
        $maj = [int]$parts[0]; $min = [int]$parts[1]
        if ($maj -gt $MinMajor -or ($maj -eq $MinMajor -and $min -ge $MinMinor)) {
            Write-Host "Found Python $v ($($c -join ' '))"
            return ,$c
        }
        Write-Host "Skipping Python $v ($($c -join ' ')): version $MinMajor.$MinMinor or newer is needed"
    }
    return $null
}

Write-Host "=== matrix2ptu installer ===" -ForegroundColor Cyan
Write-Host "Folder: $Root"

if (-not (Test-Path (Join-Path $Root "pyproject.toml"))) {
    Fail "pyproject.toml not found. Run this script from the matrix2ptu folder."
}

$venv = Join-Path $Root ".venv"
$venvPy = Join-Path $venv "Scripts\python.exe"

if ($Recreate -and (Test-Path $venv)) {
    Write-Host "Removing old .venv ..."
    Remove-Item -Recurse -Force $venv
}

if (-not (Test-Path $venvPy)) {
    $py = Find-Python
    if (-not $py) {
        Fail ("Python $MinMajor.$MinMinor or newer was not found.`n" +
              "Install it from https://www.python.org/downloads/windows/ " +
              "(tick 'Add python.exe to PATH' in the installer), then run this script again.")
    }
    $pyExe = $py[0]
    $pyArgs = @($py | Select-Object -Skip 1)
    Write-Host "Creating the environment in .venv ..."
    & $pyExe @pyArgs -m venv $venv
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path $venvPy)) { Fail "Could not create .venv" }
} else {
    Write-Host "Using the existing environment in .venv (use -Recreate to rebuild it)"
}

Write-Host "Updating pip ..."
& $venvPy -m pip install --upgrade pip --disable-pip-version-check -q
if ($LASTEXITCODE -ne 0) { Fail "pip update failed (is the internet connection working?)" }

$extras = @()
if (-not $NoValidate) { $extras += "validate" }
if ($Dev) { $extras += "dev" }
$target = "."
if ($extras.Count -gt 0) { $target = ".[" + ($extras -join ",") + "]" }

Write-Host "Installing matrix2ptu $target ..."
& $venvPy -m pip install --disable-pip-version-check -e $target
if ($LASTEXITCODE -ne 0) { Fail "Installation failed, see the messages above." }

$exe = Join-Path $venv "Scripts\matrix2ptu.exe"
& $exe --help | Out-Null
if ($LASTEXITCODE -ne 0) { Fail "matrix2ptu was installed but does not start." }

Write-Host ""
Write-Host "matrix2ptu is installed." -ForegroundColor Green
Write-Host "Start it by double-clicking matrix2ptu_terminal.bat, which opens a terminal"
Write-Host "in this folder where the 'matrix2ptu' command works, for example:"
Write-Host "    matrix2ptu convert configs\my_data.yaml --dry-run"
Write-Host "If you move this folder, run the installer again with -Recreate."
