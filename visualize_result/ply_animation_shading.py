import os
import open3d as o3d
import numpy as np
import imageio

# ディレクトリのパス
ply_directory_panc = "/Volumes/KIOXIA/eguchi/seminar_result/animation/ply_panc"
ply_directory_stomach = "/Volumes/KIOXIA/eguchi/seminar_result/animation/ply_stom"
ply_directory_duo = "/Volumes/KIOXIA/eguchi/seminar_result/animation/ply_duo"
ply_directory_lkid = "/Volumes/KIOXIA/eguchi/seminar_result/animation/ply_lkid/"

# 出力アニメーションファイル名とパス
#output_animation_path = "/Volumes/KIOXIA/eguchi/seminar_result/animation/output_animation.gif"

# アニメーション生成の設定
vis = o3d.visualization.Visualizer()
vis.create_window()

# 回転角度（ラジアン）
rotation_angle_x = np.pi / 2 - np.pi/5
rotation_angle_y = 0*np.pi / 2 #+ np.pi/5
rotation_angle_z = 0*np.pi / 2 

# 膵臓の連番PLYファイルをソートして読み込み
ply_files_panc = sorted([file for file in os.listdir(ply_directory_panc) if file.endswith(".ply")])
print(ply_files_panc)

# 胃の連番PLYファイルをソートして読み込み
ply_files_stomach = sorted([file for file in os.listdir(ply_directory_stomach) if file.endswith(".ply")])

ply_files_duodenum = sorted([file for file in os.listdir(ply_directory_duo) if file.endswith(".ply")])

ply_files_lkid = sorted([file for file in os.listdir(ply_directory_lkid) if file.endswith(".ply")])
# 両方のファイル数が同じであることを確認
#assert len(ply_files_panc) == len(ply_files_stomach), "The number of PLY files for pancreas and stomach must be the same."

images = []

# ウィンドウサイズを指定
width = 800
height = 600

# アニメーション生成の設定
vis = o3d.visualization.Visualizer()
vis.create_window(width=width, height=width)

opt = vis.get_render_option()
opt.light_on = True
opt.point_size = 5

for i in range(len(ply_files_panc)):
    print("number: "+str(i))
    ply_path_panc = os.path.join(ply_directory_panc, ply_files_panc[i])
    ply_path_stomach = os.path.join(ply_directory_stomach, ply_files_stomach[i])
    ply_path_duo = os.path.join(ply_directory_duo, ply_files_duodenum[i])
    ply_path_lkid = os.path.join(ply_directory_lkid, ply_files_lkid[i])
    
    # PLYファイルを読み込む
    pcd_panc = o3d.io.read_point_cloud(ply_path_panc)
    pcd_stomach = o3d.io.read_point_cloud(ply_path_stomach)
    pcd_duo = o3d.io.read_point_cloud(ply_path_duo)
    pcd_lkid = o3d.io.read_point_cloud(ply_path_lkid)
    #法線の設定（影を付けるため）
    pcd_panc.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=0.1, max_nn=30))
    pcd_stomach.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=0.1, max_nn=30))
    pcd_duo.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=0.1, max_nn=30))
    pcd_lkid.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=0.1, max_nn=30))

    # 回転行列を作成
    rotation_matrix_x = np.array([[1, 0, 0],
                                   [0, np.cos(rotation_angle_x), -np.sin(rotation_angle_x)],
                                   [0, np.sin(rotation_angle_x), np.cos(rotation_angle_x)]])
    
    rotation_matrix_y = np.array([[np.cos(rotation_angle_y), 0, np.sin(rotation_angle_y)],
                                   [0, 1, 0],
                                   [-np.sin(rotation_angle_y), 0, np.cos(rotation_angle_y)]])
    
    rotation_matrix_z = np.array([[np.cos(rotation_angle_z), -np.sin(rotation_angle_z), 0],
                                   [np.sin(rotation_angle_z), np.cos(rotation_angle_z), 0],
                                   [0, 0, 1]])

    rotation_matrix = rotation_matrix_x @ rotation_matrix_y @ rotation_matrix_z

    # 膵臓のポイントクラウドを回転させる
    pcd_panc.points = o3d.utility.Vector3dVector(np.dot(np.asarray(pcd_panc.points), rotation_matrix))
    
    # 胃のポイントクラウドを回転させる
    pcd_stomach.points = o3d.utility.Vector3dVector(np.dot(np.asarray(pcd_stomach.points), rotation_matrix))

    pcd_duo.points = o3d.utility.Vector3dVector(np.dot(np.asarray(pcd_duo.points), rotation_matrix))
    print("duo points: " + str(pcd_duo.points))
    pcd_lkid.points = o3d.utility.Vector3dVector(np.dot(np.asarray(pcd_lkid.points), rotation_matrix))
    #print("lkidpoints: " + str(pcd_lkid.points))
    # 膵臓の中心座標を計算
    center_panc = np.asarray(pcd_panc.points).mean(axis=0)

    # 胃の中心座標を計算
    center_stomach = np.asarray(pcd_stomach.points).mean(axis=0)

    center_duo = np.asarray(pcd_duo.points).mean(axis=0)

    center_lkid = np.asarray(pcd_lkid.points).mean(axis=0)

    view = vis.get_view_control()
    view.set_zoom(0.5)

    # 膵臓の中心座標を注視点に設定
    view.set_lookat(center_panc)

    # 膵臓のポイントクラウドを可視化
    vis.add_geometry(pcd_panc)

    # 胃のポイントクラウドを可視化
    vis.add_geometry(pcd_stomach)

    vis.add_geometry(pcd_duo)

    vis.add_geometry(pcd_lkid)

    # マウス操作を受け付ける
    vis.poll_events()
    vis.update_renderer()

    # ウィンドウの内容を画像にキャプチャして保存
    img = vis.capture_screen_float_buffer()
    
    images.append(((np.asarray(img) * 255).astype(np.uint8) * 255).astype(np.uint8))

# アニメーションをGIF形式で保存
#imageio.mimsave(output_animation_path, images, duration=0.1)

# ウィンドウを破棄
vis.destroy_window()