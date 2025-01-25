#!/bin/bash

# Скрипт предназначен для поиска и выбора одноименных веток в используемых репозиториях.  

shopt -s nocasematch  # Включить режим без учета регистра

#Получение всех веток репозиториев
FIXED_SOLVERS_BRANCHES="$(git ls-remote --heads https://gitlab-ci-token:${CI_FIXED_SOLVERS_TOKEN}@git.ktk.gubkin.ru/solvers/fixed_solvers.git | awk '{print $2}' | cut -d'/' -f3-)"

#Задаем стандартные ветки
DEFAULT_FIXED_SOLVERS_BRANCH="main"

#Определяем PDE_SOLVERS_BRANCH (ветку в репозитории PDE_SOLVERS)
if [ "$CI_PIPELINE_SOURCE" == "merge_request_event" ]; then
	export PDE_SOLVERS_BRANCH="$CI_MERGE_REQUEST_SOURCE_BRANCH_NAME"
else
	export PDE_SOLVERS_BRANCH="$CI_COMMIT_BRANCH"
fi

#Удаляем файл если он есть. Данный файл используется для передачи переменных между стадиями в CI.
filename='branches.txt'
if [[ "$filename" != "" && -f "$filename" ]]; then
	rm -v "$filename"
fi

#Функция для определения веток
defining_branch() {
    local branches="$1" # Первый аргумент: перечень всех веток репозитория
    local target_branch="$2" # Второй аргумент: название ветки, которую мы ищем
    local variable_name="$3" # Третий аргумент: название переменной в которой будем хранить название ветки
    local default_branch="$4" # Четвертый аргумент: название стандартной ветки
	local mathfold_branch="math/$target_branch" # Доп. переменная для поиска ветки в папке pdelib
	local found=0  # Переменная для отслеживания, было ли соответствие найдено

    for branch in $branches; do
        if [[ "$branch" == "$mathfold_branch" ]]; then
            echo "$variable_name=$mathfold_branch" >> "$filename"
			found=1
			break
		elif [[ "$branch" == "$target_branch" ]]; then
			echo "$variable_name=$target_branch" >> "$filename"
			found=1
			break
		fi
    done
    if [[ $found -eq 0 ]]; then
        echo "$variable_name=$default_branch" >> "$filename"
    fi
}

# Определяем название ветки в репозитории fixed_solvers и записываем его в файл $filename
defining_branch "$FIXED_SOLVERS_BRANCHES" "$PDE_SOLVERS_BRANCH" "FIXED_SOLVERS_BRANCH" "$DEFAULT_FIXED_SOLVERS_BRANCH"
# Записываем название ветки в репозитории PDE_SOLVERS в файл $filename
echo "PDE_SOLVERS_BRANCH=$PDE_SOLVERS_BRANCH" >> "$filename"

shopt -u nocasematch # Выключить режим без учета регистра

cat "$filename" # Выводим содержимое файла $filename
