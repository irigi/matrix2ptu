matrix2ptu - convert Abberior MATRIX .raw / streamed .obf photon streams to PicoQuant PTU
==========================================================================================

Nothing to install: this folder contains its own Python. It can be anywhere
(also a USB stick); it is only read, never changed.

START
  Double-click matrix2ptu_terminal.bat. A terminal opens with the matrix2ptu
  command ready. If Windows says "Windows protected your PC", click
  "More info" and "Run anyway" (the files are not signed).

CONVERT (example for data in D:\Data\2026-11-11)
  matrix2ptu inspect D:\Data\2026-11-11\Matrix.raw --msr D:\Data\2026-11-11\Matrix.msr
  copy configs\matrix.yaml D:\Data\2026-11-11\matrix.yaml
  notepad D:\Data\2026-11-11\matrix.yaml
      set input.raw (or input.obf), input.msr and output.directory as
      absolute paths in single quotes, e.g.  raw: 'D:\Data\2026-11-11\Matrix.raw'
  matrix2ptu convert D:\Data\2026-11-11\matrix.yaml --dry-run
  matrix2ptu convert D:\Data\2026-11-11\matrix.yaml --validate

  Example configs for other kinds of data are in configs\. All options,
  output files and how to open the PTUs: README.md in this folder.

FROM OTHER SCRIPTS
  Call matrix2ptu.bat with its full path, e.g.
  "C:\matrix2ptu\matrix2ptu-<version>-windows-x64\matrix2ptu.bat" convert my.yaml
  Exit code: 0 = ok, 1 = error, 2 = a validation check failed.

PIXEL REASSIGNMENT
  Experimental and disabled by default (correctness not established yet).
  convert never applies it.

Source, updates: https://github.com/irigi/matrix2ptu/releases
