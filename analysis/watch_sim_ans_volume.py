import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from mpl_toolkits.mplot3d import Axes3D

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
    combined_stack[overlap, :] = [255, 255, 0] 
    
    return combined_stack

def plot_stack(stack, sample_size):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # 非ゼロボクセルのインデックスを取得
    non_zero_indices = np.argwhere(stack[..., 0] + stack[..., 1] + stack[..., 2] > 0)
    
    # サンプルサイズ以上のポイントがある場合はランダムサンプル
    if len(non_zero_indices) > sample_size:
        sampled_indices = non_zero_indices[np.random.choice(non_zero_indices.shape[0], sample_size, replace=False)]
    else:
        sampled_indices = non_zero_indices
    
    # x, y, zの座標を取得（xとyを入れ替え）
    y, x, z = sampled_indices[:, 0], sampled_indices[:, 1], sampled_indices[:, 2]  # xとyを入れ替え
    colors = stack[y, x, z] / 255.0  # RGB値を正規化
    
    # スキャッタープロット
    ax.scatter(x, y, z, c=colors, marker='o', s=2)
    
    # 原点のプロット（赤い大きな点）
    ax.scatter(0, 0, 0, c='red', marker='o', s=100, label="Origin")
    
    # X, Y, Z軸方向の矢印をプロット
    arrow_length = max(stack.shape) // 10  # 矢印の長さを適切に設定
    ax.quiver(0, 0, 0, arrow_length, 0, 0, color='r', linewidth=1.5, arrow_length_ratio=0.2)
    ax.quiver(0, 0, 0, 0, arrow_length, 0, color='g', linewidth=1.5, arrow_length_ratio=0.2)
    ax.quiver(0, 0, 0, 0, 0, arrow_length, color='b', linewidth=1.5, arrow_length_ratio=0.2)
    
    # ラベル付け
    ax.text(arrow_length, 0, 0, 'X', color='red', fontsize=12)
    ax.text(0, arrow_length, 0, 'Y', color='green', fontsize=12)
    ax.text(0, 0, arrow_length, 'Z', color='blue', fontsize=12)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    
    # 凡例を追加して原点を示す
    ax.legend()
    
    plt.show()



# フォルダパスを指定してください
folder1 = '/Users/lelab/Downloads/exam/patient2/ans_data'
#folder2 = '/Users/lelab/Downloads/exam/patient1/image'
folder2 = '/Users/lelab/Downloads/exam/patient2/image'

stack1 = load_stack(folder1)
stack2 = load_stack(folder2)

colored_stack = create_colored_stack(stack1, stack2)

plot_stack(colored_stack, 10000)
