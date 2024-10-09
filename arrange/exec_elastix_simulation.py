import subprocess
import os
import glob
import sys
import numpy as np

import pandas as pd

# 0行めはmm単位、1行めはpixel単位


def exec_simulation(cwd):
    base_path = os.path.abspath("calc_3D_elastix.py")
    base_path = base_path[:base_path.rfind("/") - 8]

    os.chdir(cwd)

    # ここは使用するデータでパスを変える必要がある
    ###########################
    ii_path = "../../image/"
    ii_name = "IMG%04d.tiff"
    il_path = "../../label/"
    il_name = "IMG%04d.tiff"
    ###########################

    out_csv_path = "../data/"

    os.chdir("simulation/elastix")

    out_csv = out_csv_path + "elastix_data.csv"
    out_log = os.getcwd() + "/log.log"
    print(out_csv)
    trans_file = "../../3D_elastix/dir/TransformParameters.0.R2.txt"

    bin = base_path + "/mpm_mac/Bin/membrane"
    ai_path = "../../making_data/max_tiff_img/Answer/"
    ai_name = "IMG%04d.tiff"  ##<--hara's way : %04d.tiff, myway : IMG%04d.tiff
    dice_file = "../data/dice_score.csv"
    elx_img = "../elastix/"
    info_path = "../../info.txt"

    cmd = bin + " -outlog " + out_log + " -iipath " + ii_path + " -iiname " + ii_name + " -ilpath " + il_path + " -ilname " + il_name + " -mode 3 -tar 1234 -output 1 -3Dtrans " + \
        trans_file + " -outElximg " + elx_img + " -out3Dcsv " + out_csv + " -aipath " + \
        ai_path + " -ainame " + ai_name + " -dfile " + dice_file + " -info " + info_path
    subprocess.run(cmd, shell=True)
    os.chdir("../../")


def calc_elastix_error():

    # spacingについて
    spacing = [0, 0, 0]
    with open('info.txt') as f:
        for line in f:
            if line.find("pixel dimension") != -1:
                print(line.split())
                spacing[0] = float(line.split()[3])
                spacing[1] = float(line.split()[4])
                spacing[2] = float(line.split()[5])

    print(spacing[0], spacing[1], spacing[2])
    os.chdir("simulation/data")
    column_list = ["index", "x", "y", "z", "cs", "dcm", "type"]

    # データのインポート
    elx_df = pd.read_csv("elastix_data.csv", header=None)
    ans_df = pd.read_csv("../../making_data/answer_data.csv", header=None)

    ans_df.columns = column_list
    elx_df.columns = column_list

    print(elx_df.head())
    print(ans_df.head())

    ans_df_list = []
    elx_df_list = []
    ans_df_list.append(ans_df)
    ans_df_list.append(ans_df.query("type == 1"))
    ans_df_list.append(ans_df.query("type == 2"))
    ans_df_list.append(ans_df.query("type == 3"))
    ans_df_list.append(ans_df.query("type == 4"))

    elx_df_list.append(elx_df)
    elx_df_list.append(elx_df.query("type == 1"))
    elx_df_list.append(elx_df.query("type == 2"))
    elx_df_list.append(elx_df.query("type == 3"))
    elx_df_list.append(elx_df.query("type == 4"))

    # 出力用のdfを定義
    out_df = pd.DataFrame()
    out_index = ["all", "panc", "stom", "duo", "lkid"]

    idx = -1
    for (ans, elx) in zip(ans_df_list, elx_df_list):
        idx = idx + 1
        sub = pd.DataFrame()  # データ整形用の一時的な配列を定義
        # mm単位にする場合
        sub["each_error_mm"] = spacing[0]*spacing[0]*pow(ans["x"] - elx["x"], 2) + spacing[1]*spacing[1]*pow(
            ans["y"] - elx["y"], 2) + spacing[2]*spacing[2]*pow(ans["z"] - elx["z"], 2)
        # pixel単位にする場合
        sub["each_error_px"] = pow(ans_df["x"] - elx["x"], 2) + pow(
            ans_df["y"] - elx["y"], 2) + pow(ans_df["z"] - elx["z"], 2)
        each_error_val_mm = sub["each_error_mm"].values  # 誤差を計算できるようにnumpyに伝達
        each_error_val_mm = np.sqrt(each_error_val_mm)  # 平方根を計算
        error_mm = np.sum(each_error_val_mm) / \
            float(np.size(each_error_val_mm))  # 誤差平均
        print(error_mm, out_index[idx])

        each_error_val_px = sub["each_error_px"].values  # 誤差を計算できるようにnumpyに伝達
        each_error_val_px = np.sqrt(each_error_val_px)  # 平方根を計算
        error_px = np.sum(each_error_val_px) / \
            float(np.size(each_error_val_px))  # 誤差平均
        print(error_px, out_index[idx])

        # dfに出力
        out_df[out_index[idx]] = [error_mm, error_px]

    print(out_df)
    out_df.to_csv('elastix_error.csv')


def calc_init_error():
    os.chdir("../../")
    print(os.getcwd())

    # spacingについて
    spacing = [0, 0, 0]
    with open('info.txt') as f:
        for line in f:
            if line.find("pixel dimension") != -1:
                print(line.split())
                spacing[0] = float(line.split()[3])
                spacing[1] = float(line.split()[4])
                spacing[2] = float(line.split()[5])

    print(spacing[0], spacing[1], spacing[2])
    os.chdir("simulation/data")
    column_list = ["index", "x", "y", "z", "cs", "dcm", "type"]

    # データのインポート
    init_df = pd.read_csv("../../making_data/init_position.csv", header=None)
    ans_df = pd.read_csv("../../making_data/answer_data.csv", header=None)

    ans_df.columns = column_list
    init_df.columns = column_list

    print(init_df.head())
    print(ans_df.head())

    ans_df_list = []
    init_df_list = []
    ans_df_list.append(ans_df)
    ans_df_list.append(ans_df.query("type == 1"))
    ans_df_list.append(ans_df.query("type == 2"))
    ans_df_list.append(ans_df.query("type == 3"))
    ans_df_list.append(ans_df.query("type == 4"))

    init_df_list.append(init_df)
    init_df_list.append(init_df.query("type == 1"))
    init_df_list.append(init_df.query("type == 2"))
    init_df_list.append(init_df.query("type == 3"))
    init_df_list.append(init_df.query("type == 4"))

    # 出力用のdfを定義
    out_df = pd.DataFrame()
    out_index = ["all", "panc", "stom", "duo", "lkid"]

    idx = -1
    for (ans, init) in zip(ans_df_list, init_df_list):
        idx = idx + 1
        sub = pd.DataFrame()  # データ整形用の一時的な配列を定義
        # mm単位にする場合
        sub["each_error_mm"] = spacing[0]*spacing[0]*pow(ans["x"] - init["x"], 2) + spacing[1]*spacing[1]*pow(
            ans["y"] - init["y"], 2) + spacing[2]*spacing[2]*pow(ans["z"] - init["z"], 2)
        # pixel単位にする場合
        sub["each_error_px"] = pow(ans_df["x"] - init["x"], 2) + pow(
            ans_df["y"] - init["y"], 2) + pow(ans_df["z"] - init["z"], 2)
        each_error_val_mm = sub["each_error_mm"].values  # 誤差を計算できるようにnumpyに伝達
        each_error_val_mm = np.sqrt(each_error_val_mm)  # 平方根を計算
        error_mm = np.sum(each_error_val_mm) / \
            float(np.size(each_error_val_mm))  # 誤差平均
        print(error_mm, out_index[idx])

        each_error_val_px = sub["each_error_px"].values  # 誤差を計算できるようにnumpyに伝達
        each_error_val_px = np.sqrt(each_error_val_px)  # 平方根を計算
        error_px = np.sum(each_error_val_px) / \
            float(np.size(each_error_val_px))  # 誤差平均
        print(error_px, out_index[idx])

        # dfに出力
        out_df[out_index[idx]] = [error_mm, error_px]

    print(out_df)
    out_df.to_csv('init_error.csv')


if __name__ == '__main__':
    args = sys.argv
    if len(args) == 2:
        exec_simulation(args[1])
        calc_elastix_error()
        calc_init_error()
    else:
        print("invalid input")
