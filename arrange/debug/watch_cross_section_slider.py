import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import nrrd
from skimage import measure

# NRRDファイルを読み込む
stack1, header1 = nrrd.read("/Users/lelab/Downloads/exam/patient1/debug/shifted_stack.nrrd")  # 赤（正解）"D:/3_skelton/patient1/day2/3d_object_stack.nrrd"
stack2, header2 = nrrd.read("/Users/lelab/Downloads/exam/patient1/image/3d_object_stack.nrrd")  # 青"D:/3_skelton/patient1/day1/transformed_3d_object_stack.nrrd"

# スタックのサイズを取得する
size_z, size_y, size_x = stack1.shape

# 重なったエリアを抽出する
intersection = np.logical_and(stack1, stack2)

# カラーマップを設定する
cmap_stack1 = np.array([[1, 0, 0, 0.5]])  # 赤
cmap_stack2 = np.array([[0, 0, 1, 0.5]])  # 青
cmap_intersection = np.array([[1, 0, 1, 0.5]])  # マゼンタ
cmap_contour = np.array([[0, 0, 0, 1]])  # 黒

# 初期断面を設定する
z_idx = 120
y_idx = size_y // 2
x_idx = size_x // 2

# 図を作成する
fig, ax = plt.subplots(1, 3, figsize=(15, 5))

# Z方向の断面を表示する
rgb_image_z = np.zeros((size_y, size_x, 4))
im_z = ax[0].imshow(rgb_image_z)
ax[0].set_title('Coronal')
ax[0].axis('off')

# Y方向の断面を表示する
rgb_image_y = np.zeros((size_z, size_x, 4))
im_y = ax[1].imshow(rgb_image_y)
ax[1].set_title('Sagittal')
ax[1].axis('off')

# X方向の断面を表示する
rgb_image_x = np.zeros((size_z, size_y, 4))
im_x = ax[2].imshow(rgb_image_x)
ax[2].set_title('Axial')
ax[2].axis('off')

# スライダーを作成する
ax_slider_z = plt.axes([0.2, 0.02, 0.6, 0.03])
slider_z = Slider(ax_slider_z, 'Coronal', 0, size_z - 1, valinit=z_idx, valstep=1)
ax_slider_y = plt.axes([0.2, 0.06, 0.6, 0.03])
slider_y = Slider(ax_slider_y, 'Sagittal', 0, size_y - 1, valinit=y_idx, valstep=1)
ax_slider_x = plt.axes([0.2, 0.1, 0.6, 0.03])
slider_x = Slider(ax_slider_x, 'Axial', 0, size_x - 1, valinit=x_idx, valstep=1)

# スライダーの値が変更されたときの処理
def update_z(val):
    z_idx = int(slider_z.val)
    rgb_image_z = np.zeros((size_y, size_x, 4))
    rgb_image_z[stack1[z_idx] > 0] = cmap_stack1
    rgb_image_z[stack2[z_idx] > 0] = cmap_stack2
    rgb_image_z[intersection[z_idx] > 0] = cmap_intersection

    # stack1の輪郭を黒く塗る
    contours = measure.find_contours(stack1[z_idx], 0.5)
    for contour in contours:
        contour = contour.astype(int)
        rgb_image_z[contour[:, 0], contour[:, 1]] = cmap_contour

    im_z.set_data(rgb_image_z)
    fig.canvas.draw_idle()

def update_y(val):
    y_idx = int(slider_y.val)
    rgb_image_y = np.zeros((size_z, size_x, 4))
    rgb_image_y[stack1[:, y_idx, :] > 0] = cmap_stack1
    rgb_image_y[stack2[:, y_idx, :] > 0] = cmap_stack2
    rgb_image_y[intersection[:, y_idx, :] > 0] = cmap_intersection

    # stack1の輪郭を黒く塗る
    contours = measure.find_contours(stack1[:, y_idx, :], 0.5)
    for contour in contours:
        contour = contour.astype(int)
        rgb_image_y[contour[:, 0], contour[:, 1]] = cmap_contour

    im_y.set_data(rgb_image_y)
    fig.canvas.draw_idle()

def update_x(val):
    x_idx = int(slider_x.val)
    rgb_image_x = np.zeros((size_z, size_y, 4))
    rgb_image_x[stack1[:, :, x_idx] > 0] = cmap_stack1
    rgb_image_x[stack2[:, :, x_idx] > 0] = cmap_stack2
    rgb_image_x[intersection[:, :, x_idx] > 0] = cmap_intersection

    # stack1の輪郭を黒く塗る
    contours = measure.find_contours(stack1[:, :, x_idx], 0.5)
    for contour in contours:
        contour = contour.astype(int)
        rgb_image_x[contour[:, 0], contour[:, 1]] = cmap_contour

    im_x.set_data(rgb_image_x)
    fig.canvas.draw_idle()

slider_z.on_changed(update_z)
slider_y.on_changed(update_y)
slider_x.on_changed(update_x)

# ラベルを表示
fig.text(0.3, 0.95, 'shifted volume', color='red')
fig.text(0.6, 0.95, 'volume', color='blue')

# 初期断面を表示
update_z(z_idx)
update_y(y_idx)
update_x(x_idx)

# 図を表示する
plt.tight_layout()
plt.show()