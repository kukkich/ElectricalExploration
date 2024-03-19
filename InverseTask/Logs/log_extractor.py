import re
import argparse

def extract_numbers_and_save_to_file(log_file_path, pattern, output_file_path):
    """
    Извлекает числовые значения из логов, соответствующие заданному паттерну и сохраняет их в файл.

    Parameters:
    log_file_path (str): Путь к файлу логов.
    pattern (str): Регулярное выражение для поиска в логах.
    output_file_path (str): Путь к файлу для сохранения результатов.
    """

    compiled_pattern = re.compile(pattern)

    with open(log_file_path, 'r', encoding='utf-8') as log_file, \
         open(output_file_path, 'w', encoding='utf-8') as output_file:
        for line in log_file:
            matches = compiled_pattern.findall(line)
            for match in matches:
                output_file.write(match + '\n')

def main():
    parser = argparse.ArgumentParser(description="Извлекает числа из логов и сохраняет их в файл.")
    parser.add_argument("log_file_path", help="Путь к файлу логов.")
    parser.add_argument("pattern", help="Регулярное выражение для поиска числовых значений.")
    parser.add_argument("output_file_path", help="Путь к файлу для сохранения результатов.")

    args = parser.parse_args()

    extract_numbers_and_save_to_file(args.log_file_path, args.pattern, args.output_file_path)

if __name__ == "__main__":
    main()
