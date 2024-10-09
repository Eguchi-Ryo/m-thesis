import subprocess
import os
import glob
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def exec_simulation(csv_path):
    base_path = os.path.abspath("make_answer_data.py")
    base_path = base_path[:base_path.rfind("/") - 11]
    trans_file_base = "making_data/transform3d.txt"

    # ??????????????????????
    ###########################
    ii_path = "../../../image/"
    ii_name = "IMG%04d.tiff"
    il_path = "../../../label/"
    il_name = "IMG%04d.tiff"
    ###########################

    out_img = "tiff_img/"
    out_scl = "max_tiff_img/"
    out_csv_path = "../../data/"
    bin = base_path + "/mpm_mac/Bin/membrane"
    ai_path = "../../../making_data/max_tiff_img/Answer/"
    ##modified  (prev "%04d.tiff"##)
    ai_name = "IMG%04d.tiff"
    #################
    dice_file = "../../data/dice_score.csv"
    info_path = "../../../info.txt"
    unconverged_path = "../../data/unconvrged_list.csv"
    divergence_path = "../../data/divergence_list.csv"
    df_csv = "../../../registration_opt/deformation_field_csv/analysis_mhd_file_affine.csv"
    cs_list = "../../../organ_cs_list.txt"
    affine3D_path = "../../../registration/TransformParameters.0.txt" ##後で決める,コマンドライン引数の"-3Daffine"使えるか要チェック
    df_ratio_list = "../../../deformation_ratio_list.txt"

    deformationField_path = "../../../registration/deformationField_3D.raw"

    tar_dict = {"all":"1234", "panc":"1","panc-stom":"12","panc-duo":"13","panc-lkid":"14","panc-stom-duo":"123","panc-stom-lkid":"124","panc-duo-lkid":"134"}

    # ????????????????

    while(1):

        #????????????????
        data_file = pd.read_csv(csv_path)
        print(data_file.head())
        try:
            if(len(data_file) == 0) :
                break
            path = data_file.head(1)["path"].values[0]
            tar_dirname = data_file.head(1)["tar"].values[0]
            num = "%04d" % (data_file.head(1)["num"].values[0])
            data_file[1:].to_csv(csv_path,index=False)
            print(path)
            print(tar_dirname)
            print(num)
        except:
            data_file[1:].to_csv(csv_path,index=False)
            continue

        

        print(path,num,tar_dirname)
        os.chdir(path + "simulation")
        os.chdir(num)
        os.chdir(tar_dirname)

        out_csv = out_csv_path + tar_dirname + "/data%04d.csv" % int(num)
        out_log = os.getcwd() + "/log.log"
        # print(out_csv)
        trans_file = "../../../" + trans_file_base
        tar = tar_dict[tar_dirname]

        cmd = bin + " -outlog " + out_log + " -iipath " + ii_path + " -iiname " + ii_name + " -ilpath " + il_path +  " -ilname " + il_name + " -outimg " + out_img + " -outscl " + out_scl + " -outcsv " + out_csv + " -trans " + trans_file + " -cs " +  num  + " -tar " + tar + " -mode 2 -output 1" + " -aipath " + ai_path + " -ainame " + ai_name + " -dfile " + dice_file + " -info " + info_path + " -unconverged " + unconverged_path  + " -divergence " + divergence_path + " -VolumeAffine " + affine3D_path + " -deformationField_path " + deformationField_path
        print("exec : " + num + "," + tar_dirname)
        # print("exec_cmd : " + cmd)
        subprocess.run(cmd,shell=True)
        #plot_errors(out_log)

def plot_errors(log_file_path):
    plt.ion()
    fig, ax = plt.subplots()
    ax.set_xlabel("Step")
    ax.set_ylabel("Error")
    ax.set_title("Real-time Error Plot")

    steps = []
    errors = []

    with open(log_file_path, 'r') as log_file:
        lines = log_file.readlines()
        for line in lines:
            if 'error =' in line:
                step = int(line.split(',')[0])
                error = float(line.split('=')[1])
                steps.append(step)
                errors.append(error)
                ax.plot(steps, errors, 'r-')
                plt.draw()
                plt.pause(0.01)

    plt.ioff()
    plt.show()
         

if __name__ == '__main__':
    args = sys.argv
    if len(args) == 2:
        exec_simulation(args[1])
    else:
        print("invalid input")

