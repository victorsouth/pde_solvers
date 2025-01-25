# Получаем абсолютный путь к директории, в которой находится этот скрипт
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

# Создаем директорию для файлов журнала
$LogDir = Join-Path $ScriptDir "..\pipeline_result\build_mingw_result"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

# ----------------------------------------------------------------------------------------------------------------------------------------------
# Сборка и установка pde_solvers
# ----------------------------------------------------------------------------------------------------------------------------------------------
Set-Location -Path "../../"
New-Item -ItemType Directory -Force -Path "build" | Out-Null
Set-Location -Path "build"
# Собираем pde_solvers
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-w -Wno-#pragma-messages" -DCMAKE_INSTALL_PREFIX="../Libs/MinGW" > "$LogDir/pde_solvers_build.log" 2>&1
cmake --build . -j 8 >> "$LogDir/pde_solvers_build.log" 2>&1
$MakeResult = $LASTEXITCODE

# Выводим ошибки из файла pde_solvers_build.log, содержащие 'error:' или 'fatal error:'
Select-String -Pattern "error:|fatal error:|Error" -Path "$LogDir/pde_solvers_build.log" | ForEach-Object { $_.Line }
if ($MakeResult -ne 0) {
    Write-Host "--------------- ERROR ---------------"
    Write-Host "Error: pde_solvers [$env:PDE_SOLVERS_BRANCH] build failed"
    Write-Host "-------------------------------------"
    exit 1
}

Write-Host "pde_solvers [$env:PDE_SOLVERS_BRANCH] successfully assembled"

# Устанавливаем pde_solvers
cmake --install . >> "$LogDir/pde_solvers_build.log" 2>&1

if ($LASTEXITCODE -ne 0) {
    Write-Host "--------------- ERROR ---------------"
    Write-Host "Error: pde_solvers [$env:PDE_SOLVERS_BRANCH] installation failed"
    Write-Host "-------------------------------------"
    exit 1
}
Write-Host "pde_solvers [$env:PDE_SOLVERS_BRANCH] installed successfully"