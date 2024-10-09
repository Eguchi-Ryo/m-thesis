import os
import subprocess
import sys

def initialize(cwd):
    os.chdir(cwd)
    os.getcwd()
    subprocess.run("mkdir ans_data",shell=True)
    subprocess.run("mkdir base_data",shell=True)
    subprocess.run("mkdir elastix_data",shell=True)
    subprocess.run("mkdir simulation",shell=True)
    subprocess.run("mkdir making_data",shell=True)
    subprocess.run("mkdir image",shell = True)
    subprocess.run("mkdir label",shell = True)
    subprocess.run("mkdir 3D_elastix",shell = True)
    subprocess.run("mkdir answer", shell = True)

    # making_dataのしたにはサブディレクトリが存在
    os.chdir("making_data")
    subprocess.run("mkdir tiff_img",shell=True)
    subprocess.run("mkdir max_tiff_img",shell=True)

    os.chdir("../simulation")
    subprocess.run("mkdir data",shell=True)
    subprocess.run("mkdir elastix",shell=True)    

if __name__ == "__main__":
    args = sys.argv
    if len(args) == 2:
        initialize(args[1])
    else:
        print("invalid input")


