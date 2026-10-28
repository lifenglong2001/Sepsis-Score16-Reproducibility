[CmdletBinding()]
param([string]$Rscript = "Rscript")
$ErrorActionPreference = "Continue"
$root = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$failed = 0
foreach ($t in @("test_static_integrity.R", "test_checks_integrity.R")) {
  Write-Output "===== $t"
  & $Rscript (Join-Path $root "tests/$t")
  if ($LASTEXITCODE -ne 0) { $failed++ }
}
if ($failed -gt 0) { Write-Output "TESTS FAILED: $failed"; exit 1 }
Write-Output "ALL TESTS PASSED"
exit 0
