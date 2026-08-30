$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$log = Join-Path $repo 'github_upload_log.txt'

function Invoke-Git {
    param([Parameter(ValueFromRemainingArguments = $true)][string[]]$Arguments)
    Write-Host ("git " + ($Arguments -join ' ')) -ForegroundColor Cyan
    & $script:git @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Git command failed (exit code $LASTEXITCODE): git $($Arguments -join ' ')"
    }
}

try {
    Start-Transcript -LiteralPath $log -Force | Out-Null
    & (Join-Path $repo 'tools/validate_public_repo.ps1')
    $gitCommand = Get-Command git -ErrorAction SilentlyContinue
    if ($gitCommand) {
        $script:git = $gitCommand.Source
    } else {
        $gitCandidates = @(
            (Join-Path $env:USERPROFILE '.cache\codex-runtimes\codex-primary-runtime\dependencies\native\git\cmd\git.exe'),
            (Join-Path $env:LOCALAPPDATA 'GitHubDesktop\app-*\resources\app\git\cmd\git.exe'),
            'C:\Program Files\Git\cmd\git.exe',
            'C:\Program Files\Git\bin\git.exe'
        )
        $script:git = $gitCandidates | ForEach-Object { Get-Item $_ -ErrorAction SilentlyContinue } |
            Sort-Object FullName -Descending | Select-Object -First 1 -ExpandProperty FullName
    }
    if (-not $script:git -or -not (Test-Path -LiteralPath $script:git)) {
        throw 'Git was not found. Install Git for Windows and run this script again.'
    }
    Write-Host ('Using Git: ' + $script:git) -ForegroundColor DarkGray

    Add-Type -AssemblyName Microsoft.VisualBasic
    $remote = [Microsoft.VisualBasic.Interaction]::InputBox(
        'Paste the URL of the new GitHub repository, then click OK.',
        'GitHub repository URL',
        'https://github.com/username/morchella-soil-microbiome-analysis.git'
    ).Trim().TrimEnd('/')
    if (-not $remote) { throw 'No GitHub repository URL was provided.' }
    if ($remote -notmatch '^https://github\.com/[^/]+/[^/]+(?:\.git)?$' -and $remote -notmatch '^git@github\.com:[^/]+/[^/]+(?:\.git)?$') {
        throw 'Invalid GitHub URL. Example: https://github.com/username/repository.git'
    }

    Push-Location $repo
    try {
        if (-not (Test-Path '.git')) { Invoke-Git init }
        Invoke-Git branch -M main
        Invoke-Git add --all

        $staged = & $script:git diff --cached --name-only
        if ($LASTEXITCODE -ne 0) { throw 'Unable to inspect staged files.' }
        if ($staged) {
            Invoke-Git commit -m 'Release reproducible analysis code'
        } elseif (-not (& $script:git rev-parse --verify HEAD 2>$null)) {
            throw 'No files are available to commit.'
        } else {
            Write-Host 'The local code is already committed.' -ForegroundColor Yellow
        }

        $remotes = @(& $script:git remote)
        if ($remotes -contains 'origin') { Invoke-Git remote set-url origin $remote }
        else { Invoke-Git remote add origin $remote }

        & $script:git ls-remote --exit-code --heads origin main *> $null
        $remoteMainExists = ($LASTEXITCODE -eq 0)
        if ($remoteMainExists) {
            Write-Host 'A remote main branch was detected. Checking whether the repository was pre-populated.' -ForegroundColor Yellow
            Invoke-Git fetch origin main
            & $script:git merge-base HEAD origin/main *> $null
            if ($LASTEXITCODE -ne 0) {
                throw 'The remote repository is not empty (usually because README, .gitignore, or License was selected during creation). To avoid overwriting remote content, create a completely empty repository and run this script again.'
            }
            Invoke-Git pull --rebase origin main
        }

        Invoke-Git push -u origin main
        Write-Host 'GitHub upload completed.' -ForegroundColor Green
    }
    finally { Pop-Location }
}
catch {
    Write-Host ''
    Write-Host ('Upload did not complete: ' + $_.Exception.Message) -ForegroundColor Red
    Write-Host ('Full log: ' + $log) -ForegroundColor Yellow
    exit 1
}
finally {
    try { Stop-Transcript | Out-Null } catch {}
}
