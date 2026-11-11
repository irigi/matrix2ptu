@echo off
rem Double-click to open a terminal in this folder with the matrix2ptu command ready.
cd /d "%~dp0"
set "PATH=%~dp0;%PATH%"
cmd /k "echo matrix2ptu is ready. Type: matrix2ptu --help"
