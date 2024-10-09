import SimpleITK as sitk
import numpy as np
import csv 

# MHDファイルのパス
path = "/Users/lelab/Downloads/exam/patient2/registration/deformationField_3D.mhd"

# 画像データの読み込み
image = sitk.ReadImage(path)
image_array = sitk.GetArrayFromImage(image)

# 配列の次元を取得
print("image array shape : "+str(image_array.shape))
depth, width, height, dim = image_array.shape

"""
# CSVファイルのパス
csv_path = "/Users/lelab/Downloads/exam/patient1/debug/deformationField3d.csv"
# CSVファイルに書き込み
with open(csv_path, 'w', newline='') as file:
    writer = csv.writer(file)
    for z in range(depth):
        for x in range(width):
            for y in range(height):
                writer.writerow([z, x, y, image_array[z, x, y]])
"""

print("Array shape: " + str(image_array.shape))
print("Displacement at [0,0,0]: " + str(image_array[0, 0, 0]))
print("Displacement at [232,575,575]: " + str(image_array[232, 575, 575]))
print("Displacement at [232,0,0]: " + str(image_array[232, 0, 0]))
print("Displacement at [116,287,287]: " + str(image_array[116, 287, 287]))
