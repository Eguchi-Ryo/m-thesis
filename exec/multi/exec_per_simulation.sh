abs_path=$(cd $(dirname exec_list.csv); pwd)/$(basename exec_list.csv)
python3 exec_per_simulation.py $abs_path