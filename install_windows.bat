@echo off
rem Double-click to install or update matrix2ptu (runs install_windows.ps1).
rem Options are passed on, e.g.:  install_windows.bat -Dev
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install_windows.ps1" %*
set RC=%ERRORLEVEL%
echo.
pause
exit /b %RC%
