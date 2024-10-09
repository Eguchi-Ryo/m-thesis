#nrrdを可視化するコード

import SimpleITK as sitk
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# NRRDファイルの読み込み
itk_image = sitk.ReadImage("/Users/lelab/Downloads/exam/patient2/image/new_3d_object_planes.nrrd")
#itk_image = sitk.ReadImage("/Users/lelab/Downloads/exam/patient1/image/combined_object.nrrd")

# SimpleITKの画像をnumpy配列に変換
image_array = sitk.GetArrayFromImage(itk_image)

# 非ゼロのボクセルの座標を取得
nonzero_coords = np.nonzero(image_array)

# 3Dプロットの準備
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# 非ゼロのボクセルをプロット
ax.scatter(nonzero_coords[2], nonzero_coords[1], nonzero_coords[0], marker='.', color='b')

# 軸ラベルの設定
ax.set_xlabel('X軸')
ax.set_ylabel('Y軸')
ax.set_zlabel('Z軸')

# プロットの表示
plt.title('3D Object from NRRD File')
plt.show()
