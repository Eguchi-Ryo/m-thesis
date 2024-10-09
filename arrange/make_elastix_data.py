from asyncio import subprocess
import glob
import imp
import os
import sys
from tifffile import TiffFile, imwrite
from PIL import Image
import numpy as np


def make_elastix_data(cwd):

    # 作業ディレクトリに移動
    os.chdir(cwd)

    # simulationデータがあるディレクトリに移動
    os.chdir("making_data")

    base_tar = "../base_data/"

    ##########
    # ここは今後変える必要がある
    base_img_path = "../image/"
    #########
    ans_tar = "../ans_data/"
    # 画像を読み込む
    ans_files = sorted(glob.glob("tiff_img/Answer/All/*.tiff"))
    base_files = sorted(glob.glob(base_img_path + "*.tiff"))
    # data_files = ["tiff_img/200/0100.tiff"]
    print(base_img_path + "tiff_img/*.tiff")
    print(len(base_files))
    print(len(ans_files))
    # 最大値を保存
    max_val = 0
    for file in ans_files:
        n_list = file.split("/")
        print(n_list)
        with TiffFile(file) as tif:
            img = tif.asarray()
            max_val = max(max_val, np.amax(img))
            # imwrite(tar+n_list[-1],img)

    # 画像の情報を保存
    with Image.open(base_files[0]) as img:
        tiff_info = img.info

    # 全ての画像に倍率をかけて保存する
    scale = int(65535/max_val)

    for i, file in enumerate(ans_files):
        n_list = file.split("/")
        print(n_list)
        with TiffFile(file) as tif:
            img = tif.asarray()
            img = img * scale
            newimg = Image.fromarray(img)
            fname = "IMG%04d.tiff" % (i+1)
            newimg.save(ans_tar + fname, **tiff_info)
            # imwrite(ans_tar+fname,img)

    for i, file in enumerate(base_files):
        n_list = file.split("/")
        print(n_list)
        with TiffFile(file) as tif:
            img = tif.asarray()
            img = img * scale
            newimg = Image.fromarray(img)
            fname = "IMG%04d.tiff" % (i+1)
            newimg.save(base_tar + fname, **tiff_info)
            # imwrite(base_tar+fname,img)


if __name__ == "__main__":
    args = sys.argv
    if len(args) == 2:
        make_elastix_data(args[1])
    else:
        print("ivnvalid input")
