@echo off
rem Double-click to open a terminal in this folder with matrix2ptu ready to use.
cd /d "%~dp0"
if not exist ".venv\Scripts\activate.bat" (
    echo matrix2ptu is not installed yet. Double-click install_windows.bat first.
    pause
    exit /b 1
)
cmd /k "call .venv\Scripts\activate.bat && echo matrix2ptu is ready. Type: matrix2ptu --help"
