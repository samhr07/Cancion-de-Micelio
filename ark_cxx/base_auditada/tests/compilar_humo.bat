@echo off
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1
if not exist build mkdir build
cl /nologo /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"..\include" humo_base_auditada.cpp /Fo"build\\" /Fe"build\humo_base_auditada.exe" && build\humo_base_auditada.exe
