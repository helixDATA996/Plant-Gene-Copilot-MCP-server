@echo off
setlocal
cd /d "%~dp0"

where py >nul 2>nul
if errorlevel 1 (
    echo Python was not found. Install Python 3.10 or later and enable Add Python to PATH.
    pause
    exit /b 1
)

echo Checking and installing Python dependencies...
py -m pip install -r requirements.txt
if errorlevel 1 (
    echo Dependency installation failed. Check your network or run: py -m pip install -r requirements.txt
    pause
    exit /b 1
)

echo.
echo Evaluation webpage: http://localhost:5000
echo Close this window to stop the server.
echo.
py server.py
pause
