#!/bin/bash
# Скрипт для статического анализа всей кодовой базы

#Получаем абсолютный путь к директории, в которой находится этот скрипт
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#Определяем куда идет ветка
if [ "${CI_PIPELINE_SOURCE}" = "merge_request_event" ]; then
  echo -e "\033[32m---------------MR---------------"
  echo -e "\033[32morigin/${CI_MERGE_REQUEST_SOURCE_BRANCH_NAME} to origin/${CI_MERGE_REQUEST_TARGET_BRANCH_NAME}"
else
  echo -e "\033[32m-------------COMMIT-------------"
fi

declare repoPath="${CI_PROJECT_DIR}"
declare path=$repoPath"/**/*"

shopt -s globstar
for i in $path
do
    if [ -f "$i" ];
    then
        echo "${i%/*}""/""${i##*/}"| grep -owP "[^\s]+\.cpp|[^\s]+\.h" >> files
    fi
done

for line in $(cat files);
do
	echo $line
	clang-tidy $line -checks='modernize-use-default-member-init,modernize-use-default,cppcoreguidelines-virtual-class-destructor,modernize-use-using,readability-convert-member-functions-to-static' >> clang_full_log.txt 2>&1
done

cat "clang_full_log.txt" | grep "modernize-use-default-member-init\|modernize-use-default\|cppcoreguidelines-virtual-class-destructor\|modernize-use-using\|readability-convert-member-functions-to-static"| egrep -v "clang-diagnostic" >> clang_log.txt

cp "clang_log.txt" "${SCRIPT_DIR}/../pipeline_result/clang_log.txt"

if [ -s clang_log.txt ]; then
	echo -e "\033[33m----------Warnings----------"
	tac clang_log.txt
    exit 1;
else
	echo -e "\033[32m--------No-warnings---------"
	exit 0
fi
