import subprocess
import os
import glob
import sys
import numpy as np
import pandas as pd


def exec_simulation(csv_path: str):
    base_path = os.path.abspath("make_answer_data.py")
    base_path = base_path[:base_path.rfind("/") - 11]

    trans_file_base = "making_data/transform3d.txt"

    # ここは使用するデータでパスを変える必要がある
    ###########################
    ii_path = "../../image/"
    ii_name = "IMG%04d.tiff"
    il_path = "../../label/"
    il_name = "IMG%04d.tiff"
    ###########################

    out_img = "tiff_img/"
    out_scl = "max_tiff_img/"
    out_csv_path = "../data/"
    bin = base_path + "/mpm/Bin/membrane"
    ai_path = "../../making_data/max_tiff_img/Answer/"
    ai_name = "%04d.tiff"
    dice_file = "../data/dice_score.csv"
    info_path = "../../info.txt"
    unconverged_path = "../data/unconvrged_list.csv"
    divergence_path = "../data/divergence_list.csv"

    # シミュレーションのリストを取得．

    while (1):

        # シミュレーションのファイルを探す
        data_file = pd.read_csv(csv_path)
        print(data_file.head())
        try:
            if (len(data_file) == 0):
                break
            path = data_file.head(1)["path"].values[0]
            cs_x = str(data_file.head(1)["x"].values[0])
            cs_y = str(data_file.head(1)["y"].values[0])
            cs_z = str(data_file.head(1)["z"].values[0])
            data_file[1:].to_csv(csv_path, index=False)
            print(path)
            print(cs_x, type(cs_x))
            print(cs_y, type(cs_y))
            print(cs_z, type(cs_z))
        except:
            data_file[1:].to_csv(csv_path, index=False)
            continue

        os.chdir(path + "simulation")
        dir_name = "%s_%s_%s" % (cs_x, cs_y, cs_z)
        print(dir_name)
        os.chdir(dir_name)

        out_csv = out_csv_path + "/data%s.csv" % dir_name
        out_log = os.getcwd() + "/log.log"
        trans_file = "../../" + trans_file_base
        tar_csv = "../../making_data/answer_data.csv"


        cmd = bin + " -outlog " + out_log + " -iipath " + ii_path + " -iiname " + ii_name + " -ilpath " + il_path + " -ilname " + il_name + " -outimg " + out_img + " -outscl " + out_scl + " -outcsv " + out_csv + " -trans " + trans_file + " -cs_x " + cs_x + " -cs_y " + \
            cs_y + " -cs_z " + cs_z + " -tar 1234 -mode 4 -output 0" + " -aipath " + ai_path + " -ainame " + ai_name + " -dfile " + \
            dice_file + " -info " + info_path + " -unconverged " + unconverged_path + \
            " -divergence " + divergence_path + " -tarcsv " + tar_csv

        print("exec_cmd : " + cmd)
        subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    args = sys.argv
    if len(args) == 2:
        exec_simulation(args[1])
    else:
        print("invalid input")
