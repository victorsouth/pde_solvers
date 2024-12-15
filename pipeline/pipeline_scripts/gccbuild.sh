#!/bin/sh

#Получаем абсолютный путь к директории, в которой находится этот скрипт
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#Создаем директорию для файлов журнала
LOG_DIR="${SCRIPT_DIR}/../pipeline_result/build_gcc_result"
mkdir -p ${LOG_DIR}

#Сборка и установка pde_solvers
mkdir -p build
cd build 
#Собираем pde_solvers
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-w -Wno-#pragma-messages" -DCMAKE_INSTALL_PREFIX=../Libs/GCC > ${LOG_DIR}/pde_solvers_build.log 2>&1
make -j$(nproc) >> ${LOG_DIR}/pde_solvers_build.log 2>&1
MAKE_RESULT=$?
# Выводим ошибки из файла pde_solvers_build.log, содержащие 'error:' или 'fatal error:'
grep -E "error:|fatal error:|Error" ${LOG_DIR}/pde_solvers_build.log
if [ $MAKE_RESULT -ne 0 ]; then
    echo "--------------- Result ---------------"
	echo "Ошибка: сборка pde_solvers [$PDE_SOLVERS_BRANCH] не удалась"
    echo "--------------------------------------"
	exit 1
fi
echo "pde_solvers [$PDE_SOLVERS_BRANCH] успешно собран"
#Устанавливаем pde_solvers
make install >> ${LOG_DIR}/pde_solvers_build.log 2>&1
if [ $? -ne 0 ]; then
    echo "--------------- Result ---------------"
	echo "Ошибка: установка pde_solvers [$PDE_SOLVERS_BRANCH] не удалась"
	echo "--------------------------------------"
    exit 1
fi
echo "pde_solvers [$PDE_SOLVERS_BRANCH] успешно установлен"