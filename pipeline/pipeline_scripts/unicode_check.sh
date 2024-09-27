#!/bin/bash

error_list="pipeline/pipeline_result/encoding_error_list.txt"

# проверяем все файлы в репозитории с расширением .m .cpp .h
for file in $(find . -type f \( -name "*.m" -o -name "*.cpp" -o -name "*.h" \)); do
    # проверяем кодировку файла
    file_encoding=$(file -b --mime-encoding "$file")
    if [ "$file_encoding" != "utf-8" ] && [ "$file_encoding" != "us-ascii" ]; then
        # записываем путь и название файла в error_list
        echo "$file" >> "$error_list"
    fi
done

echo "--------------- Result ---------------"

# Если файлы с неправильной кодировкой найдены, выводим ошибку и список файлов в отдельный файл
if [ -s "$error_list" ]; then
    echo "Error: some files are not in UTF-8 encoding. See the list of affected files in pipeline/pipeline_result/encoding_error_list.txt"
    cat "$error_list"
    echo "--------------------------------------"
    exit 1
else
    echo "All files are in Unicode encoding."
    echo "--------------------------------------"
fi
