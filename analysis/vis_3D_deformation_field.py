import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
import matplotlib.cm as cm

# MHDファイルのパス
mhd_file_path = "/Users/lelab/Downloads/exam/patient1/registration/deformationField_3D.mhd"

# MHDファイルの読み込み
image = sitk.ReadImage(mhd_file_path)

# MHDイメージからベクトルの取得
velocity_field = sitk.GetArrayFromImage(image)

# ベクトルのx、y、z成分
x_velocity = velocity_field[..., 0]  # x方向の成分
y_velocity = velocity_field[..., 1]  # y方向の成分
z_velocity = velocity_field[..., 2]  # z方向の成分

# ベクトルのサイズと位置
sample_interval = 50  # サンプリング間隔
y = np.arange(0, y_velocity.shape[0], sample_interval)
x = np.arange(0, x_velocity.shape[1], sample_interval)
z = np.arange(0, z_velocity.shape[2], sample_interval)
Y, X, Z = np.meshgrid(y, x, z, indexing='ij')  # メッシュグリッドの生成

# サンプリングされたベクトル
sampled_x_velocity = x_velocity[::sample_interval, ::sample_interval, ::sample_interval]
sampled_y_velocity = y_velocity[::sample_interval, ::sample_interval, ::sample_interval]
sampled_z_velocity = z_velocity[::sample_interval, ::sample_interval, ::sample_interval]

# ベクトルの大きさを計算
magnitude = np.sqrt(sampled_x_velocity**2 + sampled_y_velocity**2 + sampled_z_velocity**2)

# 正規化して色マップを適用
norm = mcolors.Normalize(vmin=magnitude.min(), vmax=magnitude.max())
cmap = cm.coolwarm

# 色の設定 (RGBA)
colors = cmap(norm(magnitude))

# 3Dプロットのセットアップ
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# ベクトル場の描画（矢印）
quiver = ax.quiver(X, Y, Z, sampled_x_velocity, sampled_y_velocity, sampled_z_velocity, color=colors.reshape(-1, 4), length=15, normalize=True)

# カラーバーの表示
mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array(magnitude)
plt.colorbar(mappable, ax=ax, label='Magnitude')
"""
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
"""
ax.set_title('3D deformation field with magnitude color')

plt.show()

