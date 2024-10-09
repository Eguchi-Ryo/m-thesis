import os
import subprocess
import sys
import random
import shutil
import glob
import pandas as pd
from tifffile import TiffFile, imwrite
import re


def make_ans_data(cwd):
    base_path = os.path.abspath("make_3D_answer_data.py")
    base_path = base_path[:base_path.rfind("/") - 8]
    os.chdir(cwd)
    copy_folders_path = os.path.join(cwd, "answer")
    organ_list = ['All', 'Pancreas', 'Duodenum', 'Stomach', 'Left_Kidney']
    copy_folders = []
    for organ in organ_list:
        copy_path  = os.path.join(copy_folders_path, organ)
        copy_folders.append(copy_path)
    print("copy folders:"+str(copy_folders))
    os.chdir("./making_data/tiff_img/")
    #os.mkdir("Answer")
    answer_dir = os.path.join(cwd, "making_data", "tiff_img", "Answer")
    for folder in copy_folders:
        duplicate_path = os.path.join(answer_dir, folder.split('/')[-1])
        #shutil.copytree(folder, duplicate_path)
    os.chdir(cwd)
    os.chdir("making_data")
    os.chdir("max_tiff_img")
    os.mkdir("Answer")
    duplicate_path = None
    scl_dir = os.path.join(cwd, "making_data", "max_tiff_img", "Answer")
    for folder in copy_folders:
        duplicate_path = os.path.join(scl_dir, folder.split('/')[-1])
        shutil.copytree(folder, duplicate_path)

    original_csv_path = os.path.join(cwd, "making_data", "answer_data_geometric.csv")
    copy_csv_path = os.path.join(cwd, "making_data", "answer_data.csv")
    shutil.copy(original_csv_path, copy_csv_path)

if __name__ == "__main__":
    args = sys.argv
    if len(args) == 2:
        make_ans_data(args[1])
    else:
        print("invalid input")