@echo off
setlocal
set "ROOT=D:\visual\GSE121212_RASGRP1_repro"
set "RSCRIPT=D:\Program Files\R\R-4.5.1\bin\Rscript.exe"
if not exist "D:\" exit /b 1
if not exist "%ROOT%\temp" mkdir "%ROOT%\temp"
set "TMP=%ROOT%\temp"
set "TEMP=%ROOT%\temp"
set "TMPDIR=%ROOT%\temp"
set "R_LIBS_USER=%ROOT%\R_library"
set "R_USER_CACHE_DIR=%ROOT%\package_cache"
if not exist "%RSCRIPT%" (
  echo Set RSCRIPT in this file to your Rscript.exe path.
  exit /b 1
)
"%RSCRIPT%" --vanilla "%~dp0install_dependencies.R"
if errorlevel 1 exit /b 1
"%RSCRIPT%" --vanilla "%~dp0RASGRP1_three_groups.R"
if errorlevel 1 exit /b 1
endlocal

