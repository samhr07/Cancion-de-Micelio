@echo off
rem uso: compilar.bat ruta_fuente.cpp salida.exe [include_dir]
set AQUI=%~dp0
set BASE=%AQUI%inc
if not "%~3"=="" set BASE=%~3
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cd /d "%AQUI%"
if not exist build mkdir build
cl /nologo /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%BASE%" %1 /Fo"build\\" /Fe"build\%2" /link /STACK:1073741824
