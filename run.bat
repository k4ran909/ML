@echo off
title Tubulin Drug Discovery Pipeline - Master Runner
echo ===================================================================
echo           TUBULIN DRUG DISCOVERY PIPELINE: MASTER RUNNER
echo ===================================================================
echo.
echo Starting the pipeline execution...
echo.

:: Navigate to the script directory and run the master python runner
cd /d "%~dp0src"
py all.py

echo.
echo ===================================================================
echo Execution finished.
echo ===================================================================
echo.
pause
