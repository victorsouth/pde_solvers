#!/bin/sh

# Получаем абсолютный путь к директории, в которой находится этот скрипт
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RES_DIR="${SCRIPT_DIR}/../pipeline_result/tests_out"
mkdir -p ${RES_DIR}
cd "${SCRIPT_DIR}/../../Libs/${1}/bin/"
TEST_NAME="$2"

if [ ! -f "./$TEST_NAME" ]; then
    echo "--------------- FAIL ---------------"
    echo "File $TEST_NAME NOT FOUND"
    echo "--------------------------------------"
    exit 1
fi

./${TEST_NAME} --gtest_output=xml:${RES_DIR}/${TEST_NAME}_out.xml > ${RES_DIR}/${TEST_NAME}_out.txt 2>&1
TEST_RETURN=$?

# Проверка на наличие сообщений об аварийном завершении
CRASH_ERRORS=$(grep -E 'terminate called after throwing|std::logic_error|Segmentation fault|Aborted' ${RES_DIR}/${TEST_NAME}_out.txt)
GLOBAL_TEAR_DOWN_PRESENT=$(grep -q '\[----------\] Global test environment tear-down' ${RES_DIR}/${TEST_NAME}_out.txt && echo "yes" || echo "no")

# Проверка на неудачные тесты или сбои
failed_tests=$(grep -E 'FAILED|ERROR|Failure|exception' ${RES_DIR}/${TEST_NAME}_out.txt)

if [[ -n "$failed_tests" || -n "$CRASH_ERRORS" || $TEST_RETURN -ne 0 || "$GLOBAL_TEAR_DOWN_PRESENT" = "no" ]]; then
    echo "---------- Failed tests or crashes found ----------"
    if [[ "$GLOBAL_TEAR_DOWN_PRESENT" = "yes" ]]; then
		tac ${RES_DIR}/${TEST_NAME}_out.txt | awk '/\[----------\] Global test environment tear-down/{exit} 1' | tac
	fi
    if [[ -n "$CRASH_ERRORS" ]]; then
        echo "Test failure detected:"
        echo "$CRASH_ERRORS"
		echo "More: ${TEST_NAME}_out.txt"
    fi
	echo "------------------------------------------------------"
    echo "--------------- List of tests ---------------"
    ./${TEST_NAME} --gtest_list_tests
    echo "-------------------------------------------"
    exit 1
else
    echo "---------- All tests were successful ---------"
    tac ${RES_DIR}/${TEST_NAME}_out.txt | awk '/\[----------\] Global test environment tear-down/{exit} 1' | tac
    echo "-------------------------------------------"
    echo "--------------- List of tests ---------------"
    ./${TEST_NAME} --gtest_list_tests
    echo "-------------------------------------------"
fi