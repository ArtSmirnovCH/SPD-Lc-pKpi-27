#!/bin/bash

# Создаём массив для хранения сумм
declare -a sums=(0 0 0 0 0 0 0 0 0 0 0 0 0 0) # 13 + 1

# Папка для итогового файла
input_directory="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data/signal/txt"
output_dir="/afs/jinr.ru/user/a/asmirnov/CutFlow"
output_file="${output_dir}/CutFlow_signal.txt"

mkdir -p "$output_dir"

if ls *.txt 1> /dev/null 2>&1; then
    rm *.txt
fi

for file in "$input_directory"/CutFlow_*.txt; do
    if [[ -f "$file" ]]; then
        i=0
        while IFS= read -r line; do
            sums[$i]=$((sums[$i] + line))
            ((i++))
        done < "$file"
    fi
done

for sum in "${sums[@]}"; do
    echo "$sum" >> "$output_file"
done

echo "Суммирование завершено. Результаты записаны в $output_file"

input_file="$output_file"
output_file="${output_dir}/CutFlow_ratios_signal.txt"

previous_value=""

while IFS= read -r line; do
    if [[ -n "$previous_value" ]]; then
        if [[ "$line" =~ ^-?[0-9]+(\.[0-9]+)?$ && "$previous_value" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
            ratio=$(echo "scale=6; $line / $previous_value" | bc)
            echo "$ratio" >> "$output_file"
        else
            echo "Ошибка: \"$line\" или \"$previous_value\" не является числом." >> "$output_file"
        fi
    fi
    previous_value="$line"
done < "$input_file"

echo "Результаты записаны в '$output_file'."


declare -a sums=(0 0 0 0 0 0 0 0 0 0 0 0 0 0) # 13 + 1

# II
# Папка для итогового файла
input_directory="/eos/user/a/asmirnov/reco_data_lambda_c/analysed_data/background/txt/"
output_dir="/afs/jinr.ru/user/a/asmirnov/CutFlow"
output_file="${output_dir}/CutFlow_background.txt"

mkdir -p "$output_dir"

for file in "$input_directory"/CutFlow_*.txt; do
    if [[ -f "$file" ]]; then
        i=0
        while IFS= read -r line; do
            sums[$i]=$((sums[$i] + line))
            ((i++))
        done < "$file"
    fi
done

for sum in "${sums[@]}"; do
    echo "$sum" >> "$output_file"
done

echo "Суммирование завершено. Результаты записаны в $output_file"

input_file="$output_file"
output_file="${output_dir}/CutFlow_ratios_background.txt"

# Инициализируем предыдущую переменную
previous_value=""

# Читаем файл построчно
while IFS= read -r line; do
    # Если предыдущая строка не пустая, вычисляем отношение
    if [[ -n "$previous_value" ]]; then
        # Проверяем, что оба значения являются числами
        if [[ "$line" =~ ^-?[0-9]+(\.[0-9]+)?$ && "$previous_value" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
            ratio=$(echo "scale=6; $line / $previous_value" | bc)
            echo "$ratio" >> "$output_file"
        else
            echo "Ошибка: \"$line\" или \"$previous_value\" не является числом." >> "$output_file"
        fi
    fi
    # Запоминаем текущее значение для следующей итерации
    previous_value="$line"
done < "$input_file"

echo "Результаты записаны в '$output_file'."

exit
