import re
import subprocess
import os
import sys
from tracemalloc import start

import pandas as pd
from pyparsing import col


def matching(dir_path, cs_list, exec_flg):

    base_path = os.path.abspath("matching.py")
    print(base_path)
    base_path = base_path[:base_path.rfind("/") - 8]
    print(base_path)
    # elastixの位置合わせによるデータを保管するフォルダに移動
    os.chdir(dir_path + "/elastix_data")

    # フォルダを作り、各断面でelastixを実行する
    for i in cs_list:
        dir_name = "%04d" % i
        # もしすでにディレクトリが存在していたら削除し、新規作成
        if os.path.exists(dir_name):
            subprocess.run("rm -r "+dir_name, shell=True)
        subprocess.run("mkdir "+dir_name, shell=True)

        # そのディレクトリに移動し、elastix実行に必要なディレクトリを作成する
        os.chdir(dir_name)
        subprocess.run("mkdir dir", shell=True)
        subprocess.run("mkdir dirres", shell=True)

        # 全ての断層において、elastixを実行したい場合は以下のコメントアウトを解除

        # elastixをコマンド以外で実行するとき、明示的にパスを指定する必要がある
        dyld_path = "LD_LIBRARY_PATH=" + base_path + "/elastix/lib"
        elastix_path = base_path + "/elastix/bin/elastix"
        transformix_path = base_path + "/elastix/bin/transformix"
        param_path = base_path + "/param/2D_elastix.txt"

        if exec_flg == 1:
            # シミュレーションによる答えの位置
            mv_path = "../../ans_data/IMG%04d.tiff" % i

            # base位置
            fx_path = "../../base_data/IMG%04d.tiff" % i
            el_cmd = dyld_path + " " + elastix_path + " -f " + fx_path + \
                " -m " + mv_path + " -p " + param_path + " -out dir"

            subprocess.run(el_cmd, shell=True)

            with open('dir/TransformParameters.0.R2.txt') as f:
                s = f.read()
            with open('dir/TransformParameters.0.R2.txt', "w") as f:
                f.write(s.replace("mhd", "tiff"))
            tr_cmd = dyld_path + " " + transformix_path + " -in "+mv_path + \
                " -out dirres -tp dir/TransformParameters.0.R2.txt -def all"
            subprocess.run(tr_cmd, shell=True)

        os.chdir("../")


def find_calc_all_cs(cwd):
    start_num = -10
    end_num = -10
    with open(cwd+"making_data/condition.txt") as f:
        lines = f.readlines()
        for line in lines:
            # print(line)
            if "csnum" in line:
                start_num = int(re.findall(r"\d+", line)[0])
                end_num = int(re.findall(r"\d+", line)[1])
    print(start_num, end_num)
    cs_list = [i for i in range(start_num, end_num+1)]

    return cs_list


if __name__ == '__main__':
    args = sys.argv
    if len(args) == 3:
        cs_list = list()
        cs_list = find_calc_all_cs(args[1])
        matching(args[1], cs_list, int(args[2]))
    else:
        print("invalid output")
