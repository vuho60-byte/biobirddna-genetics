@echo off
set ROOT=%~dp0
python -B -m unittest discover -s "%ROOT%tests" -v
if errorlevel 1 exit /b %errorlevel%
pushd "%ROOT%npm"
call npm test
set TEST_EXIT=%errorlevel%
popd
exit /b %TEST_EXIT%
