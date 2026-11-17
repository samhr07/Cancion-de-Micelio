@echo off
rem compilar.bat FUENTE [INCLUDE]  -> build\<fuente>.exe  (por omision la base nueva; pila de 1 GB, la de las baterias de Lilyth)
rem rutas normalizadas con %%~f: sin los ..\ el include pasaba de 260 caracteres y no abria
set AQUI=%~dp0
for %%I in ("%AQUI%..\..\..\..\ark_cxx\base_auditada\include") do set BASE=%%~fI
if not "%~2"=="" for %%I in ("%~2") do set BASE=%%~fI
for %%I in ("%AQUI%..\..\red") do set RED=%%~fI
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
cd /d "%AQUI%"
if not exist build mkdir build
cl /nologo %EXTRA% /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%BASE%" /I"%RED%\inc" /I"%RED%" "%~1" /Fo"build\\" /Fe"build\%~n1.exe" /link /STACK:1073741824
