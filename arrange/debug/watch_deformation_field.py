import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt

# MHDファイルのパス
mhd_file_path = "/Users/lelab/Downloads/exam/patient2/registration/axial/deformationField.mhd"
# MHDファイルの読み込み
image = sitk.ReadImage(mhd_file_path)
size = image.GetSize()
aspect_ratio = size[1] / size[0]

# MHDイメージからベクトルの取得
velocity_field = sitk.GetArrayFromImage(image)

# ベクトルのxとy成分
x_velocity = velocity_field[..., 0]  # x方向の成分
y_velocity = velocity_field[..., 1]  # y方向の成分

# ベクトルのサイズと位置
sample_interval = 15  # サンプリング間隔
y = np.arange(0, y_velocity.shape[0], sample_interval)
x = np.arange(0, x_velocity.shape[1], sample_interval)
X, Y = np.meshgrid(x, y)

# サンプリングされたベクトル
sampled_x_velocity = x_velocity[::sample_interval, ::sample_interval]
sampled_y_velocity = y_velocity[::sample_interval, ::sample_interval] 

# ベクトルの大きさを計算
magnitude = np.sqrt(sampled_x_velocity**2 + sampled_y_velocity**2)

# 矢印の描画と色の設定
plt.quiver(X, Y, sampled_x_velocity, sampled_y_velocity, magnitude, cmap='coolwarm', scale=700)

# カラーバーの表示
plt.colorbar(label='magnitude(pixel)', orientation='vertical')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('deformation field')

plt.gca().invert_yaxis()  # Y軸を反転
plt.gca().set_aspect(aspect_ratio)
plt.grid(True)
plt.show()


