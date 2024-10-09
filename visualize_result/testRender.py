from paraview.simple import *

# PLYファイルのパス
ply_file = "/Users/lelab/Downloads/exam/animation/ply_lkid/output0888.ply"

# PLYリーダーでファイルを読み込む
ply_reader = PLYReader(FileName=ply_file)

# レンダリング用ビューを取得または作成
render_view = GetActiveViewOrCreate('RenderView')

# データを表示する
display = Show(ply_reader, render_view)

# 表示タイプをポイントとして設定
display.SetRepresentationType('Points')

# RGB属性を使ってカラーマップを設定する
ColorBy(display, ('POINTS', 'RGB'))

# カラーマップのスケールを設定
rgb_lut = GetColorTransferFunction('RGB')
rgb_lut.RescaleTransferFunction(0.0, 255.0)

# レンダリング実行
Render()

# スクリーンショット保存（オプション）
SaveScreenshot('/Users/lelab/Downloads/exam/animation/ply_lkid/output0888_visualization.png', render_view)
"""
"""   
import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from scipy.spatial.transform import Rotation as R

def load_stack(folder):
    files = sorted([f for f in os.listdir(folder) if f.endswith('.tiff')])
    stack = []
    for file in files:
        img = Image.open(os.path.join(folder, file))
        stack.append(np.array(img))
    return np.stack(stack, axis=2)

def create_colored_stack(stack1, stack2):
    stack1_colored = np.zeros((*stack1.shape, 3), dtype=np.uint8)
    stack2_colored = np.zeros((*stack2.shape, 3), dtype=np.uint8)
    
    stack1_colored[..., 2] = stack1  # Blue channel
    stack2_colored[..., 0] = stack2  # Red channel
    
    combined_stack = stack1_colored + stack2_colored
    overlap = (stack1 > 0) & (stack2 > 0)
    combined_stack[overlap, :] = [255, 255, 0]  # Yellow for overlapping regions
    
    return combined_stack

def rotate_points(points, rotation):
    return rotation.apply(points)

def init():
    ax.scatter(rotated_points[:, 0], rotated_points[:, 1], rotated_points[:, 2], c=colors, marker='o', s=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    return fig,

def animate(i):
    ax.view_init(elev=elevation, azim=i)
    return fig,

folder1 = '/Users/lelab/Downloads/exam/patient1/ans_data'
folder2 = '/Users/lelab/Downloads/exam/patient1/image'

stack1 = load_stack(folder1)
stack2 = load_stack(folder2)

colored_stack = create_colored_stack(stack1, stack2)

# Define the spacing for each axis
spacing = np.array([1.2, 0.8, 0.8])  # Adjust this array based on actual spacing

# Get indices of non-zero voxels
non_zero_indices = np.argwhere(colored_stack[..., 0] + colored_stack[..., 1] + colored_stack[..., 2] > 0)

# Sample the indices
sample_size = 10000
if len(non_zero_indices) > sample_size:
    sampled_indices = non_zero_indices[np.random.choice(non_zero_indices.shape[0], sample_size, replace=False)]
else:
    sampled_indices = non_zero_indices

# Apply spacing to the points
points = sampled_indices[:, [2, 1, 0]] * spacing  # x, y, z with spacing applied
colors = colored_stack[sampled_indices[:, 0], sampled_indices[:, 1], sampled_indices[:, 2]] / 255.0  # Normalize RGB values

# Rotate the point
rotation = R.from_euler('xyz', [45, 125, 45], degrees=True)  # Rotate by 45 degrees around each axis
rotated_points = rotate_points(points, rotation)

# Create the plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Initial plot settings
elevation = 30  # Set elevation for a more isometric view
frames = 360  # Number of frames for full rotation

# Create animation
ani = FuncAnimation(fig, animate, init_func=init, frames=frames, interval=70, blit=False)

# Save or display the animation
# ani.save('rotation_animation.gif', writer='imagemagick')  # Save as gif
plt.show()  # Display the animation
