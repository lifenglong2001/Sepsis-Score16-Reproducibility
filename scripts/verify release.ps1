[CmdletBinding()]
param(
 [Parameter(Mandatory=$true)][string]$Config,
 [Parameter(Mandatory=$true)][string]$Rscript,
 [Parameter(Mandatory=$true)][string]$OutputDirectory,
 [ValidateSet('core','benchmarks','pbmc_score','pbmc','scoring_portable','cellular','cross_task','clinical_task','cts_reconstruction')]
 [string[]]$Stages=@('core','benchmarks','pbmc_score','pbmc','scoring_portable','cellular')
)
# Verification performs the analysis again from configured patient/donor inputs.
# It cannot certify a release by comparing already exported result CSV files.
& (Join-Path $PSScriptRoot 'run release.ps1') -Config $Config -Rscript $Rscript -OutputDirectory $OutputDirectory -Stages $Stages
exit $LASTEXITCODE
