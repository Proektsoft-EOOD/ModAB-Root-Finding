@echo off
setlocal
cd /d "%~dp0"
dotnet build tests\ModAB.Switching.Verification -c Release -m:1 -nr:false
if errorlevel 1 exit /b 1
dotnet tests\ModAB.Switching.Verification\bin\Release\net10.0\ModAB.Switching.Verification.dll --output results\verification.json
