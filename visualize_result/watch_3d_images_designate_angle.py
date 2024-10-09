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
    combined_stack[overlap, :] = [255, 255, 0]  # Yellow for overlapping regions
    
    return combined_stack

def plot_stack(stack, sample_size, elevation, azimuth):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Get indices of non-zero voxels
    non_zero_indices = np.argwhere(stack[..., 0] + stack[..., 1] + stack[..., 2] > 0)
    
    # Randomly sample the indices if there are more than sample_size points
    if len(non_zero_indices) > sample_size:
        sampled_indices = non_zero_indices[np.random.choice(non_zero_indices.shape[0], sample_size, replace=False)]
    else:
        sampled_indices = non_zero_indices
    
    x, y, z = sampled_indices[:, 2], sampled_indices[:, 1], sampled_indices[:, 0]
    colors = stack[z, y, x] / 255.0  # Normalize RGB values
    
    ax.scatter(x, y, z, c=colors, marker='o', s=1)
    """
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    """
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    ax.view_init(elev=elevation, azim=azimuth)  # Set the view angle

    plt.show()

# フォルダパスを指定してください
folder1 = '/Users/lelab/Downloads/exam/patient1/ans_data'
folder2 = '/Users/lelab/Downloads/exam/patient1/image'
#folder2 = '/Users/lelab/Downloads/exam/patient2/simulation/0097/all/tiff_img/Estimation/All'
#folder2 = '/Users/lelab/Downloads/exam/result_hara/patient2/cs_97_dice_0.732/All'

stack1 = load_stack(folder1)
stack2 = load_stack(folder2)

colored_stack = create_colored_stack(stack1, stack2)

plot_stack(colored_stack, 10000, -85, 0)  # 描画する点の数を指定

