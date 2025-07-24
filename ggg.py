# convert_egm2008.py
import os
print("Текущая рабочая директория:", os.getcwd())

print("Файлы в директории:", os.listdir())


input_file = "egm2008.dat"         # Входной файл (исходник EGM2008)
output_file = "111111egm2008_clean.dat"  # Выходной файл (очищенный)

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if not line.startswith("RECOEF"):
            continue
        parts = line.split()
        if len(parts) >= 4:
            n = parts[1]
            m = parts[2]
            c = parts[3]
            s = parts[4] if len(parts) > 4 else "0.0"
            outfile.write(f"{n} {m} {c} {s}\n")

print(f"Готово! Файл сохранён как {output_file}")
