$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$blockedExtensions = @('.fastq','.fq','.bam','.sam','.sif','.rds','.RData','.docx','.pdf','.key','.pem')
$files = Get-ChildItem -LiteralPath $repo -Recurse -File -Force | Where-Object { $_.FullName -notmatch '\\.git\\' }
$bad = @()
foreach ($file in $files) {
    $lower = $file.Name.ToLowerInvariant()
    if ($file.Length -gt 50MB) { $bad += "File larger than 50 MB: $($file.FullName)" }
    if ($blockedExtensions -contains $file.Extension) { $bad += "Blocked file type: $($file.FullName)" }
    if ($lower.EndsWith('.fastq.gz') -or $lower.EndsWith('.fq.gz')) { $bad += "Raw sequence file: $($file.FullName)" }
}
$patterns = 'hy18164053837|ZhuanZ.DESKTOP-PH97BKO|AKIA[0-9A-Z]{16}|ghp_[A-Za-z0-9]{20,}|github_pat_|BEGIN (RSA |OPENSSH )?PRIVATE KEY|password\s*=' 
$textFiles = $files | Where-Object { $_.Extension -in @('.R','.Rmd','.py','.sh','.slurm','.sbatch','.ps1','.md','.txt','.tsv','.yml','.yaml','.json') -and $_.Name -ne 'validate_public_repo.ps1' }
$hits = $textFiles | Select-String -Pattern $patterns -CaseSensitive:$false -ErrorAction SilentlyContinue
if ($hits) { $bad += ($hits | ForEach-Object { "Sensitive pattern: $($_.Path):$($_.LineNumber)" }) }
if ($bad.Count) { $bad | ForEach-Object { Write-Host $_ -ForegroundColor Red }; throw 'Public-release safety check failed.' }
Write-Host "PASS: $($files.Count) files checked; no raw data, oversized files, private keys, tokens, or known personal paths detected." -ForegroundColor Green
