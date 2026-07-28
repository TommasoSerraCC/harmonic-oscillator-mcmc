@echo off
REM Script per compilare il documento LaTeX due volte
REM in modo da risolvere tutti i riferimenti incrociati

echo Compilazione 1/2...
pdflatex -interaction=nonstopmode relazione.tex

echo.
echo Compilazione 2/2...
pdflatex -interaction=nonstopmode relazione.tex

echo.
echo Compilazione completata! I riferimenti dovrebbero ora essere corretti.
pause
