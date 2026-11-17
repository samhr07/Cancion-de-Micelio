@echo off
rem uso: compilar.bat fuente.cpp salida.exe
set AQUI=%~dp0
pushd "%AQUI%..\..\..\..\..\ark_cxx\base_auditada\include"
set BASE=%CD%
popd
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cd /d "%AQUI%"
if not exist build mkdir build
cl /nologo /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%BASE%\ark\metrologia" /I"%BASE%" %1 /Fo"build\\" /Fe"build\%2"
