# Script PowerShell per compilare il documento LaTeX due volte
# in modo da risolvere tutti i riferimenti incrociati

Write-Host "Compilazione 1/2..." -ForegroundColor Cyan
pdflatex -interaction=nonstopmode relazione.tex

Write-Host "`nCompilazione 2/2..." -ForegroundColor Cyan
pdflatex -interaction=nonstopmode relazione.tex

Write-Host "`nCompilazione completata! I riferimenti dovrebbero ora essere corretti." -ForegroundColor Green
