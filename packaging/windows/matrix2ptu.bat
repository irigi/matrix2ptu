@echo off
rem Runs matrix2ptu with the Python bundled in this folder, e.g.:
rem   matrix2ptu convert configs\my_data.yaml --dry-run
"%~dp0python\python.exe" -m matrix2ptu %*
