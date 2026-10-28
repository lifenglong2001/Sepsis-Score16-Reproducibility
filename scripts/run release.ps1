[CmdletBinding()]
param(
 [Parameter(Mandatory=$true)][string]$Config,
 [Parameter(Mandatory=$true)][string]$Rscript,
 [Parameter(Mandatory=$true)][string]$OutputDirectory,
 [ValidateSet('core','benchmarks','pbmc_score','pbmc','scoring_portable','cellular','cross_task','clinical_task','cts_reconstruction','pbmc_replication','neutrophil_paired','granulopoiesis')]
 [string[]]$Stages=@('core','benchmarks','pbmc_score','pbmc','scoring_portable','cellular')
)
$ErrorActionPreference='Stop'
try {
 $codeRoot=(Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
 $configPath=(Resolve-Path -LiteralPath $Config).Path
 if (-not (Test-Path -LiteralPath $Rscript -PathType Leaf)) { throw "Rscript executable not found: $Rscript" }
 # Read as UTF-8 explicitly: Windows PowerShell 5.1 otherwise assumes the ANSI
 # code page for BOM-less files and mangles non-ASCII (e.g. Chinese) input paths.
 $settings=Get-Content -LiteralPath $configPath -Raw -Encoding UTF8 | ConvertFrom-Json
 $adult=@($settings.discovery)+@($settings.external)
 $needed=@()
 foreach ($stage in $Stages) {
  switch ($stage) {
   'cts_reconstruction' { $needed += @('cts_core_expression','cts_core_labels','cts_target_expression','cts_public117_labels') }
   'pbmc_replication' { $needed += @('gse290679_pseudobulk_counts','gse290679_pseudobulk_metadata','gse290679_direction_registry','gmt') }
   'neutrophil_paired' { $needed += 'gse216009_neutrophil_maturity_scores' }
   'granulopoiesis' {
    $needed += $adult
    $needed += @($adult | ForEach-Object { 'expression_'+$_ })
    $needed += @('gse216009_pseudobulk_expression','mars_joined')
   }
   'cross_task' { $needed += @('cross_task_GSE236713','cross_task_GSE63042_diagnosis','cross_task_GSE63042_mortality') }
   'clinical_task' { $needed += @('clinical_GSE185263','clinical_GSE154918') }
   'core' { $needed += $adult }
   'benchmarks' { $needed += @('VANISH','GSE65682','vanish_labels','mars_labels') }
   'pbmc' { $needed += 'pbmc_scores' }
   'pbmc_score' { $needed += @('pbmc_object','pbmc_scores','gmt') }
   'cellular' { $needed += @('mixture_object','deconv_estimates','GSE185263') }
   'scoring_portable' {
    $needed += $adult; $needed += 'gmt'
    $needed += @($adult | ForEach-Object { 'expression_'+$_ })
    $needed += @($settings.discovery | ForEach-Object { 'ssgsea_expression_'+$_ })
   }

  }
 }
 foreach ($key in @($needed | Select-Object -Unique)) {
  $value=$settings.inputs.$key
  if (-not $value -or -not (Test-Path -LiteralPath $value -PathType Leaf)) { throw "Required configured input missing: $key = $value" }
 }
 $out=[IO.Path]::GetFullPath($OutputDirectory)
 if ((Test-Path -LiteralPath $out) -and @(Get-ChildItem -LiteralPath $out -Force).Count -gt 0) { throw "Output directory must be new or empty: $out" }
 New-Item -ItemType Directory -Path $out -Force | Out-Null
 New-Item -ItemType Directory -Path (Join-Path $out 'logs') | Out-Null
 $settings.output=Join-Path $out 'results'
 $effective=Join-Path $out 'effective_config.json'
 # Write BOM-less UTF-8 so R's jsonlite can parse it directly.
 $effectiveJson=$settings | ConvertTo-Json -Depth 30
 [IO.File]::WriteAllText($effective,$effectiveJson,(New-Object System.Text.UTF8Encoding($false)))
 $env:SEPSIS_CODE_ROOT=$codeRoot
 # Do not force a locale that may be absent on this host: R then keeps its own
 # default and emits no spurious "Setting LC_* failed" lines on stderr.
 [Environment]::SetEnvironmentVariable('LC_ALL',$null,'Process')
 [Environment]::SetEnvironmentVariable('LANG',$null,'Process')
 if (@($settings.r_libraries).Count -gt 0) { $env:R_LIBS=@($settings.r_libraries) -join [IO.Path]::PathSeparator }
 $records=@()
 foreach ($stage in $Stages) {
  $started=(Get-Date).ToUniversalTime().ToString('o')
  $clock=[Diagnostics.Stopwatch]::StartNew()
  $log=Join-Path $out "logs/$stage.log"
  # Windows PowerShell 5.1 turns any native stderr line into a terminating
  # NativeCommandError under ErrorActionPreference=Stop. Capture both streams
  # with the preference relaxed, then restore it and use the exit code.
  $savedEAP=$ErrorActionPreference
  $ErrorActionPreference='Continue'
  & $Rscript (Join-Path $codeRoot 'scripts/run release.R') $effective $stage *> $log
  $nativeExit=$LASTEXITCODE
  $ErrorActionPreference=$savedEAP
  $clock.Stop()
  $records += [pscustomobject]@{stage=$stage;start_utc=$started;elapsed_seconds=$clock.Elapsed.TotalSeconds;exit_code=$nativeExit;log=$log}
  $records | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $out 'execution.json') -Encoding utf8
  if ($nativeExit -ne 0) { throw "Reproduction stage failed: $stage (exit $nativeExit); see $log" }
 }
 Write-Output "PASS: fresh analysis executed for stages $($Stages -join ','); results: $out"
 exit 0
} catch {
 Write-Error -Message $_ -ErrorAction Continue
 exit 1
}
