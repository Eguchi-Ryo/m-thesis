import math
import re
import sys
import os
import subprocess
import cv2
import glob
import numpy as np
import random
from PIL import Image
import pandas as pd
from tifffile import TiffFile, imwrite
import cv2


def pre_exec_elastix(cwd):
    base_path = os.path.abspath("make_answer_data.py")
    base_path = base_path[:base_path.rfind("/") - 8]
    os.chdir(cwd)
    subprocess.run("mkdir buf", shell=True)
    os.chdir("buf")

    # elastix?????????????????
    subprocess.run("mkdir dir", shell=True)
    subprocess.run("mkdir dirres", shell=True)

    os.chdir(cwd)

    os.chdir("making_data")
    subprocess.run("mkdir tiff_img", shell=True)
    subprocess.run("mkdir tiff_img/Geometric", shell=True)

    base_tar = "./tiff_img/Geometric/"
    base_img_path = "../image/"
    answer_label_path = "../answer/"

    img_size = [0, 0, 0]

    # ???????

    base_img_files = sorted(glob.glob(base_img_path + "*.tiff"))
    label_files = []
    label_files.append(sorted(glob.glob(answer_label_path + "Pancreas/*.tiff")))
    label_files.append(sorted(glob.glob(answer_label_path + "Stomach/*.tiff")))
    label_files.append(sorted(glob.glob(answer_label_path + "Duodenum/*.tiff")))
    label_files.append(sorted(glob.glob(answer_label_path + "Left_Kidney/*.tiff")))

    out_df = pd.DataFrame()
    init_df= pd.DataFrame()

    out_data_x = list()
    out_data_y = list()
    out_data_z = list()
    out_data_cs = list()
    out_data_val = list()
    out_data_type = list()

    init_data_x = list()
    init_data_y = list()
    init_data_z = list()
    init_data_cs = list()
    init_data_val = list()
    init_data_type = list()

    vol_img = list()
    out_img = list()

    # ??????????????
    for file in base_img_files:
        print(file)
        splitted = file.split("/")
        cs = int(re.sub(r"\D", "", splitted[-1]))
        print(cs)
        with TiffFile(file) as tif:
            img = tif.asarray()
            if img_size[0] == 0:
                img_size[0] = img.shape[0]
                img_size[1] = img.shape[1]
                img_size[2] = len(base_img_files)
                print(img_size)

            vol_img.append(img)

    # ??
    print(len(vol_img), vol_img[0].shape)

    # ??????????
    for z in range(len(vol_img)):
        buf = np.zeros(vol_img[z].shape)
        out_img.append(buf)

    name_list = list()
    org_type = 0

    # elastix????????????????????????????????
    exist_flag = False
    cross_sectional = [0, 0]

    for labels in label_files:
        org_type = org_type + 1
        for label in labels:
            print(label)
            splitted = label.split("/")
            name_list.append(splitted[-1])
            cs = int(re.sub(r"\D", "", splitted[-1]))
            # cs = 100
            # TIFF?????????
            with TiffFile(label) as tif:
                # np???????
                img = tif.asarray()

                print(img_size)

                # ??????????????????????
                if np.amax(img) != 0 and org_type == 1 and exist_flag == False:
                    cross_sectional[0] = cs
                    exist_flag = True
                elif exist_flag and cross_sectional[1] == 0 and np.amax(img) == 0:
                    cross_sectional[1] = cs - 1
                    exist_flag = False

                for col in range(len(img[0])):
                    for row in range(len(img)):
                        if img[row][col] != 0:
                            # print(row,col)
                            # ?????????
                            # ?????????????

                            init_vec = np.array(
                                [col, (img_size[1]-row), 1], dtype='float')
                            vec = np.array(
                                [col, (img_size[1]-row), 1], dtype='float')
 
                            # print(vec)
                            if vec[0] >= img_size[0] or vec[0] < 0 or vec[1] >= img_size[1] or vec[1] < 0:
                                print(cs, vec)

                            # ?????

                            out_img[cs-1][img_size[1] -
                                          int(vec[1])][int(vec[0])] = vol_img[cs-1][row][col]

                            # ??????
                            out_data_x.append(vec[0])
                            out_data_y.append(vec[1])
                            out_data_z.append(cs-1)
                            out_data_cs.append(cs)
                            out_data_val.append(vol_img[cs-1][row][col])
                            out_data_type.append(org_type)

                            init_data_x.append(init_vec[0])
                            init_data_y.append(init_vec[1])
                            init_data_z.append(cs-1)
                            init_data_cs.append(cs)
                            init_data_val.append(vol_img[cs-1][row][col])
                            init_data_type.append(org_type)

    # ????df???
    out_df["x"] = out_data_x
    out_df["y"] = out_data_y
    out_df["z"] = out_data_z
    out_df["cs"] = out_data_cs
    out_df["val"] = out_data_val
    out_df["type"] = out_data_type

    init_df["x"] = init_data_x
    init_df["y"] = init_data_y
    init_df["z"] = init_data_z
    init_df["cs"] = init_data_cs
    init_df["val"] = init_data_val
    init_df["type"] = init_data_type

    print(out_df.head())
    out_df.to_csv("answer_data_geometric.csv", header=False)
    init_df.to_csv("init_position.csv", header=False)
    
    # TIFF???????????
    tiff_info = list()
    tiff_mode = list()
    for labels in label_files:
        for label in labels:
            with Image.open(label) as img:
                tiff_info.append(img.info)
                tiff_mode.append(img.mode)

    # ?????
    for i, img in enumerate(out_img):
        img = img.astype(np.uint16)

        newimg = Image.fromarray(img)
        newimg.save(base_tar+name_list[i], **tiff_info[i])

if __name__ == "__main__":
    args = sys.argv
    if len(args) == 2:
        pre_exec_elastix(args[1])
    else:
        print("invalid input")