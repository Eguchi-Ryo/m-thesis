python3 make_answer_data.py $1
python3 make_3D_answer_data.py $1 $2
python3 make_elastix_data.py $1 
python3 matching.py $1 $3
python3 make_folders.py $1
python3 calc_3D_elastix.py $1
python3 exec_elastix_simulation.py $1

