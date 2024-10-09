import glob
import os
import subprocess
import sys
from PIL import Image
import numpy as np
from tifffile import TiffFile

base_path = os.path.abspath("calc_3D_elastix.py")
base_path = base_path[:base_path.rfind("/")- 8]

def make_3D_tiff(cwd):
    os.chdir(cwd)

    # 対象のフォルダに移動
    os.chdir("3D_elastix")

    # まとめる画像データ
    base_img_path = "../base_data/"
    ans_img_path = "../ans_data/"

    # 出力先
    ans_name = "ans_stack.tiff"
    base_name = "base_stack.tiff"


    # 画像を読み込む
    ans_files = sorted(glob.glob(ans_img_path + "*.tiff"))
    base_files = sorted(glob.glob(base_img_path + "*.tiff"))

    # 画像の情報を保存
    with Image.open(base_files[0]) as img:
        tiff_info = img.info

    ans_stack = []
    # 3次元画像を作成
    for i,file in enumerate(ans_files):
        n_list = file.split("/")
        print(n_list)
        print(i)
        with TiffFile(file) as tif:
            img = tif.asarray()
            newimg = Image.fromarray(img)
            ans_stack.append(newimg)

    ans_stack[0].save(ans_name , save_all = True, append_images=ans_stack[1:], **tiff_info)

    base_stack = []
    # 3次元画像を作成
    for i,file in enumerate(base_files):
        n_list = file.split("/")
        print(n_list)
        print(i)
        with TiffFile(file) as tif:
            img = tif.asarray()
            newimg = Image.fromarray(img)
            base_stack.append(newimg)

    base_stack[0].save(base_name , save_all = True, append_images=base_stack[1:], **tiff_info)

def calc_3D_elastix():
    # elastixの位置合わせによるデータを保管するフォルダに移動
    # 上の関数により該当フォルダにいることが前提

    subprocess.run("mkdir dir",shell = True)
    subprocess.run("mkdir dirres",shell = True)

    # elastixをコマンド以外で実行するとき、明示的にパスを指定する必要がある
    dyld_path = "LD_LIBRARY_PATH=" + base_path + "/elastix/lib"
    elastix_path = base_path + "/elastix/bin/elastix"
    param_path = base_path + "/param/3D_elastix.txt"
    transformix_path = base_path + "/elastix/bin/transformix"

    # シミュレーションによる答えの位置
    mv_path = "ans_stack.tiff"


    # base位置
    fx_path = "base_stack.tiff"
    el_cmd = dyld_path + " " + elastix_path +  " -f " + fx_path + " -m " + mv_path + " -p " + param_path + " -out dir"
    
    
    #print(el_cmd)
    subprocess.run(el_cmd,shell=True)

    
    with open('dir/TransformParameters.0.R2.txt') as f:
        s = f.read()
    with open('dir/TransformParameters.0.R2.txt',"w") as f:
        f.write(s.replace("mhd","tiff"))
    
    tr_cmd = dyld_path + " " + transformix_path + " -in "+mv_path + " -out dirres -tp dir/TransformParameters.0.R2.txt -def all"
    subprocess.run(tr_cmd,shell = True)



    

    

if __name__ == "__main__":
    args = sys.argv
    if len(args) == 2:
        make_3D_tiff(args[1])
        calc_3D_elastix()
    else :
        print("invalid output")