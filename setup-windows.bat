@echo off
setlocal EnableExtensions EnableDelayedExpansion

echo ===============================================
echo Top6Meta Modeling Environment Setup (Windows)
echo ===============================================
echo.

:: ------------------------------------------------
:: Detect Conda
:: ------------------------------------------------

set CONDA_FOUND=0

:: 1) Conda on PATH
where conda >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    echo Conda found on PATH.
    set CONDA_FOUND=1
    goto create_env
)

:: 2) Miniconda (default)
if exist "%UserProfile%\miniconda3\Scripts\conda.bat" (
    echo Found Miniconda installation.
    call "%UserProfile%\miniconda3\Scripts\activate.bat"
    set CONDA_FOUND=1
    goto create_env
)

:: 3) Anaconda (default)
if exist "%UserProfile%\anaconda3\Scripts\conda.bat" (
    echo Found Anaconda installation.
    call "%UserProfile%\anaconda3\Scripts\activate.bat"
    set CONDA_FOUND=1
    goto create_env
)

:: ------------------------------------------------
:: Conda not found → install Miniconda
:: ------------------------------------------------
:install_conda

echo Conda not found. Installing Miniconda...
echo.

powershell -Command "Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser -Force"
powershell -Command ^
 "Invoke-WebRequest -Uri 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe' ^
  -OutFile 'Miniconda3-installer.exe'"

if not exist "Miniconda3-installer.exe" (
    echo ERROR: Failed to download Miniconda installer.
    pause
    exit /b 1
)

echo Installing Miniconda silently...
Miniconda3-installer.exe ^
 /InstallationType=JustMe ^
 /RegisterPython=1 ^
 /AddToPath=0 ^
 /S ^
 /D=%UserProfile%\miniconda3

del Miniconda3-installer.exe

echo Initializing Conda...
"%UserProfile%\miniconda3\Scripts\conda.exe" init

echo.
echo Conda installed successfully.
echo PLEASE CLOSE this window, open a NEW Command Prompt,
echo then re-run this script.
echo.
pause
exit /b 0

:: ------------------------------------------------
:: Create / activate environment
:: ------------------------------------------------
:create_env

echo.
echo Checking for requirements.yml...

if not exist "requirements.yml" (
    echo ERROR: requirements.yml not found.
    echo Run this script from the project root directory.
    pause
    exit /b 1
)

echo.
echo Checking if environment "modeling_studio" exists...

conda env list | findstr /C:"modeling_studio" >nul
if %ERRORLEVEL% EQU 0 (
    echo Environment already exists.
) else (
    echo Creating environment from requirements.yml...
    conda env create -f requirements.yml
    if %ERRORLEVEL% NEQ 0 (
        echo ERROR: Failed to create environment.
        pause
        exit /b 1
    )
)

echo.
echo Activating environment...
call conda activate modeling_studio

if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Failed to activate environment.
    pause
    exit /b 1
)

echo.
echo ===============================================
echo Setup Complete!
echo ===============================================
echo.
echo Environment "modeling_studio" is ready.
echo.
echo To run the application:
echo   cd gui
echo   python top6meta.py
echo.
echo For future use:
echo   conda activate modeling_studio
echo ===============================================
echo.

pause