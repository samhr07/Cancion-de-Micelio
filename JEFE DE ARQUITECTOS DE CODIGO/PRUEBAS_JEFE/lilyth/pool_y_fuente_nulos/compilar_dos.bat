@echo off
set EXTRA=/DBASE_NUEVA
call "%~dp0compilar.bat" volcado_previa.cpp "%~dp0base_previa\include"
call "%~dp0compilar.bat" volcado_nueva.cpp
