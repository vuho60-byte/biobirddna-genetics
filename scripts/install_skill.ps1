param(
  [Parameter(Mandatory = $true)]
  [string]$DestinationRoot,
  [string]$PythonExecutable = "python"
)

$ErrorActionPreference = "Stop"
$SkillRoot = [IO.Path]::GetFullPath(
  (Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path))
)
$DestinationRoot = [IO.Path]::GetFullPath($DestinationRoot)
$Target = Join-Path $DestinationRoot "biobirddna-genetics"
$Stage = Join-Path $DestinationRoot (
  ".biobirddna-genetics.install-" + [Guid]::NewGuid().ToString("N")
)
$LockPath = "$Target.install-lock"
$LockStream = $null

New-Item -ItemType Directory -Path $DestinationRoot -Force | Out-Null

try {
  try {
    $LockStream = [IO.File]::Open(
      $LockPath,
      [IO.FileMode]::CreateNew,
      [IO.FileAccess]::Write,
      [IO.FileShare]::None
    )
  }
  catch {
    throw "Another BioBirdDNA Genetics installation is active: $LockPath"
  }
  if (Test-Path -LiteralPath $Target) {
    throw "Target already exists: $Target. Move or version the existing skill before installing."
  }

  New-Item -ItemType Directory -Path $Stage | Out-Null
  $Prefix = $SkillRoot.TrimEnd("\", "/") + [IO.Path]::DirectorySeparatorChar
  Get-ChildItem -LiteralPath $SkillRoot -Recurse -File -Force |
    ForEach-Object {
      $Relative = $_.FullName.Substring($Prefix.Length)
      $Parts = $Relative -split "[\\/]"
      if (
        $Parts -contains "__pycache__" -or
        $_.Extension -in @(".pyc", ".pyo")
      ) {
        return
      }
      $Destination = Join-Path $Stage $Relative
      $Parent = Split-Path -Parent $Destination
      New-Item -ItemType Directory -Path $Parent -Force | Out-Null
      Copy-Item -LiteralPath $_.FullName -Destination $Destination
    }

  & $PythonExecutable -X utf8 -B (Join-Path $Stage "scripts\verify_package.py") --root $Stage
  if ($LASTEXITCODE -ne 0) {
    throw "Package integrity validation failed in staging: $Stage"
  }

  [IO.Directory]::Move($Stage, $Target)
  & $PythonExecutable -X utf8 -B (Join-Path $Target "scripts\verify_package.py") --root $Target
  if ($LASTEXITCODE -ne 0) {
    throw "Package integrity validation failed after final rename: $Target"
  }
  Write-Host "Installed verified BioBirdDNA Genetics skill to $Target"
}
finally {
  if (Test-Path -LiteralPath $Stage) {
    Remove-Item -LiteralPath $Stage -Recurse -Force
  }
  if ($null -ne $LockStream) {
    $LockStream.Dispose()
    [IO.File]::Delete($LockPath)
  }
}
