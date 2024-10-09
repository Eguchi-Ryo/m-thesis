import subprocess
import os
import glob
import sys


rm_flag = False


def make_folders(cwd):

    os.chdir(cwd)

    cs_list = []
    dir_number = sorted(glob.glob("elastix_data/*/"))

    tar_dirname = ["all", "panc", "panc-stom", "panc-duo",
                   "panc-lkid", "panc-stom-duo", "panc-stom-lkid", "panc-duo-lkid"]

    # シミュレーションのリストを取得．
    for num_buf in dir_number:
        splitted = num_buf.split("/")
        cs_list.append(splitted[-2])
    print(cs_list)

    os.chdir("simulation/")

    # elastix用のディレクトリを作成
    subprocess.run("mkdir elastix", shell=True)
    subprocess.run("mkdir elastix/tiff_img", shell=True)
    subprocess.run("mkdir elastix/max_tiff_img", shell=True)

    if os.path.exists("data") and rm_flag == True:
        subprocess.run("rm -r data", shell=True)

    subprocess.run("mkdir data", shell=True)
    for i in tar_dirname:
        subprocess.run("mkdir data/" + i, shell=True)

    index = 0
    for i in cs_list:

        if os.path.exists(i) and rm_flag == True:
            subprocess.run("rm -r "+i, shell=True)

        subprocess.run("mkdir "+i, shell=True)
        os.chdir(i)
        for k in tar_dirname:
            subprocess.run("mkdir " + k, shell=True)
            os.chdir(k)
            subprocess.run("mkdir tiff_img", shell=True)
            subprocess.run("mkdir max_tiff_img", shell=True)

            os.chdir("../")

        index = index + 1
        os.chdir("../")


if __name__ == '__main__':
    args = sys.argv
    if len(args) == 2:
        make_folders(args[1])
    else:
        print("invalid input")
