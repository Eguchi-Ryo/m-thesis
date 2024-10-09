import subprocess
import re
from PIL import Image
import os
import numpy as np
import csv
import itk
from tifffile import TiffFile
import pandas as pd
import sys
import glob
import cv2
import SimpleITK as sitk

#####elastix environment#####
dyld_path= 'LD_LIBRARY_PATH=/home/ec2-user/Desktop/research_handover/elastix/lib'
elastix_path = '/home/ec2-user/Desktop/research_handover/elastix/bin/elastix'
transformix_path = '/home/ec2-user/Desktop/research_handover/elastix/bin/transformix'  #shell = True <-- be remind
####path####
start_panc_cs_num = 1 #いつでも1にするように（全断層の倍率を出したいから）
formatted_start_num = str(start_panc_cs_num).zfill(4)
driven_cs_num = 94 #100 OK
formatted_driven_cs_num = str(driven_cs_num).zfill(4)

driven_img_fixed = "/home/ec2-user/Desktop/patient2_debug/image/IMG0094.tiff"#"/home/ec2-user/Desktop/tes2/image/IMG0038.tiff"  #folder_path + "IMG" + formatted_driven_cs_num + ".tiff"
driven_img_moved = "/home/ec2-user/Desktop/patient2_debug/ans_data/IMG0094.tiff"  #"/home/ec2-user/Desktop/tes2/ans_data/" + "IMG" + formatted_driven_cs_num + ".tiff"

mask = "/home/ec2-user/Downloads/mask/patient2/mask_plus_15pix.mhd"

Cor_fixed_img = "/home/ec2-user/Downloads/cor_patient2_ume/day2_IMG0338.tiff"
Cor_moved_img = "/home/ec2-user/Downloads/cor_patient2_ume/day1_IMG0338.tiff"


def make_folder(cwd):

    os.chdir(cwd)
    os.mkdir("registration_opt")
    os.chdir("registration_opt")
    os.mkdir("deformation_field")
    os.mkdir("transform_parameter")
    os.mkdir("registration_output_cor")
    os.mkdir("deformation_field_csv")
    os.mkdir("registration_Bspline")
    base = os.path.join(cwd, "registration_opt")
    return base


def registration_Axial_bspline(cwd):

    bspline_registration_dir = os.path.join(cwd, "registration_opt", "registration_Bspline")
    registration_param_path_Axi = "/home/ec2-user/Downloads/mask/Bspline_first.txt"
    registration_cmd = dyld_path + " " + elastix_path + " -f " + driven_img_fixed + " -fMask " + mask + " -m " + driven_img_moved + " -mMask " + mask + " -out " + bspline_registration_dir + " -p " + registration_param_path_Axi
    #registration_cmd = dyld_path + " " + elastix_path + " -f " + driven_img_moved + " -m " + driven_img_fixed + " -out " + bspline_registration_dir + " -p " + registration_param_path_Axi
    subprocess.run(registration_cmd, shell=True)   
    #transformix
    transform_parameter = os.path.join(bspline_registration_dir, "TransformParameters.0.txt")#bspline_registration_dir + "TransformParameters.0.R0.txt"
    transformix_Bspline_out_dir = bspline_registration_dir
    transformix_cmd = dyld_path + " " + transformix_path + " -def all" + " -in " + driven_img_fixed + " -out " + transformix_Bspline_out_dir + " -tp " + transform_parameter
    subprocess.run(transformix_cmd, shell=True)    

    displacement_file = os.path.join(transformix_Bspline_out_dir, "deformationField.mhd")#transformix_Bspline_out_dir + "deformationField.mhd"
    displacement_image = sitk.ReadImage(displacement_file)
    size = displacement_image.GetSize()
    print("Deformation Field size: ", size)
    displacement_image = sitk.Cast(displacement_image, sitk.sitkVectorFloat64)
    displacement_array = sitk.GetArrayViewFromImage(displacement_image)

    deformationField_csv = os.path.join(base, "deformation_field_csv", "analysis_mhd_file_affine.csv")
    with open(deformationField_csv, "w", newline='') as f:  # newline='' ???????????
        writer = csv.writer(f)

        height, width = displacement_array.shape[:2]

        # ????????
        writer.writerow(["x", "y", "dx", "dy"])
        """
        # ??????????????
        for y in range(height):
            for x in range(width):
                RightHand_y = height -1 -y
                dx, pre_dy = displacement_array[x, RightHand_y]  #?mhd????????????????
                dy = - pre_dy
                writer.writerow([x, y, dx, dy])
        """
        for y in range(height):
            for x in range(width):
                RightHand_y = height -1 -y
                dx = displacement_array[y, x, 0]
                dy = displacement_array[y, x, 1]
                writer.writerow([x, RightHand_y, dx, dy])

    # ????NumPy?????????????????????????????
    disp = np.loadtxt(deformationField_csv, delimiter=",", skiprows=1, dtype=float)

    # ?????OpenCV?????
    with TiffFile(driven_img_fixed) as tif:
        # np???????
        img = tif.asarray()
        img = img.astype(np.uint16)

    height, width = img.shape
    
    bspline_transformed_image = np.zeros_like(img)

    spacing = [1.0, 1.0]

    for col in range(width):
        for row in range(height):                       
            #??????
            x = col
            y = height - row - 1  # ??????? 
            dx = displacement_array[col, row, 0] * spacing[0] 
            dy = displacement_array[col, row, 1] * spacing[1] 

            init_vec = np.array(
                [x, y, 1], dtype='float')
            vec = np.array(
                [x, y, 1], dtype='float')
            # print(vec)
            # ??????????spacing????
            vec[0] = vec[0]*spacing[0]
            vec[1] = vec[1]*spacing[1]
            # print(vec)
            init_vec[0] = init_vec[0]*spacing[0]
            init_vec[1] = init_vec[1]*spacing[1]

            #??
            vec[0] += dx
            vec[1] -= dy

            # spacing?????
            init_vec[0] = init_vec[0]/spacing[0]
            init_vec[1] = init_vec[1]/spacing[1]

            vec[0] = vec[0]/spacing[0]
            vec[1] = vec[1]/spacing[1]

            int_vec0 = int(round(vec[0]))
            int_vec1 = int(round(vec[1]))

            int_init_vec0 = int(round(init_vec[0]))
            int_init_vec1 = int(round(init_vec[1]))

            if 0 <= int_vec0 < width and 0 <= height - 1 - int_vec1 < height:
                bspline_transformed_image[int_vec1, int_vec0] = img[row, col]
            else:
                # ???????
                int_vec0 = max(0, min(int_vec0, width - 1))
                int_vec1 = max(0, min(int_vec1, height - 1))
                bspline_transformed_image[int_vec1, int_vec0] = img[row, col]
            
        registration_image_dir = os.path.join(cwd, "registration_opt", "registration_Bspline", "deformed_image")
        os.makedirs(registration_image_dir, exist_ok=True)

        filename = os.path.join(registration_image_dir, "deformed_image.tiff")
        newimg = Image.fromarray(bspline_transformed_image.astype(np.uint8))
        newimg.save(filename)

        res_img = cv2.imread(filename, cv2.IMREAD_ANYDEPTH)
        bspline_transformed_image_path = filename

    return bspline_transformed_image_path, displacement_file, res_img

##全断層の倍率を求められるようにゼロピクセルもカウントするようにした
def count_non_black_images(folder_path):
    valid_file = 0
    for filename in os.listdir(folder_path):
        if filename.endswith(('.png', '.jpg', '.jpeg', '.tif', '.tiff')):
            image_path = os.path.join(folder_path, filename)
            image = Image.open(image_path)
            image_array = np.array(image)
            if np.count_nonzero(image_array) >= 0: 
                valid_file += 1
    return valid_file


def registration_Axial(grid_spacing_x, grid_spacing_y, registration_param_path, base, fix_image):

    driven_param_path = os.path.join(base, "transform_parameter")

    with open(registration_param_path, 'r') as regi_param_file:
        read_regi_param_lines = regi_param_file.readlines()

    grid_spacing = f"{grid_spacing_x} {grid_spacing_y}"
    new_line = f"(GridSpacingSchedule {grid_spacing} )\n"
    print("Axi:" + str(grid_spacing))

    for i, line in enumerate(read_regi_param_lines):
        if line.startswith("(GridSpacingSchedule "):
            read_regi_param_lines[i] = f"(GridSpacingSchedule {grid_spacing})\n" 
            break
    
    with open(registration_param_path, 'w') as file:
        file.writelines(read_regi_param_lines)

    registration_cmd = dyld_path + " " + elastix_path + " -f " + driven_img_moved + " -m " + fix_image + " -out " + driven_param_path + " -p " + registration_param_path
    subprocess.run(registration_cmd, shell=True)


def registration_Cornal(grid_spacing_x, grid_spacing_z, Cor_registration_path, Cor_registration_out):
    
    with open(Cor_registration_path, 'r') as cor_param_file:
        read_cor_param_lines = cor_param_file.readlines()

    grid_spacing = f"{grid_spacing_x} {grid_spacing_z}"
    print("Cor:" + str(grid_spacing))

    for i, line in enumerate(read_cor_param_lines):
        if line.startswith("(GridSpacingSchedule "):
            read_cor_param_lines[i] = f"(GridSpacingSchedule {grid_spacing} )\n"
            break
    
    with open(Cor_registration_path, 'w') as file:
        file.writelines(read_cor_param_lines)

    registration_cmd = dyld_path + " " + elastix_path + " -f " + Cor_moved_img + " -m " + Cor_fixed_img + " -out " + Cor_registration_out + " -p " + Cor_registration_path
    subprocess.run(registration_cmd, shell=True)


def extract_number_of_parameters(transform_param_path): #transform_param_path?Cor

    with open(transform_param_path, 'r') as trans_file:
        content = trans_file.read()
    pattern = r'\(NumberOfParameters (\d+)\)'
    match = re.search(pattern, content)
    if match:
        return int(match.group(1))


def correspond_cs_to_control_points_low(number_of_lows_of_control_points, ratio_of_deformation_list):
    
    each_cross_section_ratio = np.arange(0, spacing_z * number_of_CrossSection, spacing_z)
    each_control_points_low_ratio = np.arange(0, spacing_z * number_of_CrossSection, spacing_z * number_of_CrossSection / number_of_lows_of_control_points)
    corresponding_indices = []
    corresponding_ratio = []

    for cs_value in each_cross_section_ratio:
        min_index = np.argmin(np.abs(each_control_points_low_ratio - cs_value))
        corresponding_indices.append(min_index)

    corresponding_ratio = [ratio_of_deformation_list[index] for index in corresponding_indices]

    return corresponding_ratio


def determine_weight(number_of_control_point_per_low, transform_param_path, cwd):
    
    with open(transform_param_path, 'r') as file:  #Cornal????????
        content = file.read()

    pattern = r'TransformParameters\s*((?:-?\d+\.\d+\s*)+)'
    match = re.search(pattern, content)

    if match:
        TransformParameters = match.group(1)
        Pre_TransformParameters_list = [float(param) for param in TransformParameters.split()]
        TransformParameters_list = Pre_TransformParameters_list[:int(len(Pre_TransformParameters_list)/2)]

        TransformParameters_per_low = [TransformParameters_list[i:i + number_of_control_point_per_low] for i in range(0, len(TransformParameters_list), number_of_control_point_per_low)]

    #????????????
    degree_of_transform_parameters = [np.sum(np.abs(transform_parameter)) for transform_parameter in TransformParameters_per_low]
    #?????????????
    modified_ratio_list = correspond_cs_to_control_points_low(len(degree_of_transform_parameters), degree_of_transform_parameters)
    print("how largely move each low:" + str(modified_ratio_list))
    #????????
    standard_cs_num = driven_cs_num - start_panc_cs_num
    standard = modified_ratio_list[standard_cs_num]
    print("standard : " + str(standard))
    #????????
    ratio_of_deformation_list = np.array(modified_ratio_list) / standard
    print("ratio list shape: " + str(ratio_of_deformation_list.shape))
    output_file_path = os.path.join(cwd, "deformation_ratio_list.txt")
    np.savetxt(output_file_path, ratio_of_deformation_list, delimiter=",")
    return ratio_of_deformation_list

def sort_csv():
    csv_file_path = "/home/ec2-user/Desktop/target_data_geometric_proposal_way.csv"
    data = pd.read_csv(csv_file_path, header=None)  # ??????????header=None?????
    data.columns = ["index", "x", "y", "z", "cs", "val", "type", "process order", "order"]
    sorted_data = data.sort_values(by="process order")
    sorted_csv = "/home/ec2-user/Desktop/target_data_geometric.csv"
    sorted_data.to_csv(sorted_csv, index=False, header=False)

def dice_GT_init(image1, image2):
    # ???????
    image1 = cv2.imread(image1, cv2.IMREAD_GRAYSCALE)
    image2 = cv2.imread(image2, cv2.IMREAD_GRAYSCALE)
    _, image1_bin = cv2.threshold(image1, 0, 255, cv2.THRESH_BINARY)
    _, image2_bin = cv2.threshold(image2, 0, 255, cv2.THRESH_BINARY)

    # ????????
    intersection = np.logical_and(image1_bin, image2_bin)
    intersection_count = np.count_nonzero(intersection)
    total_count = np.count_nonzero(image1_bin) + np.count_nonzero(image2_bin)
    dice = 2 * intersection_count / total_count

    print("GT_init Dice: " + str(dice))

def dice_GT_result(image1, image2):
    # ???????
    image1 = cv2.imread(image1, cv2.IMREAD_GRAYSCALE)
    image2 = cv2.imread(image2, cv2.IMREAD_GRAYSCALE)
    _, image1_bin = cv2.threshold(image1, 0, 255, cv2.THRESH_BINARY)
    _, image2_bin = cv2.threshold(image2, 0, 255, cv2.THRESH_BINARY)

    # ????????
    intersection = np.logical_and(image1_bin, image2_bin)
    intersection_count = np.count_nonzero(intersection)
    total_count = np.count_nonzero(image1_bin) + np.count_nonzero(image2_bin)
    dice = 2 * intersection_count / total_count

    print("GT_result Dice: " + str(dice))

def is_black(image):
    # ?????????????
    return np.all(image == 0)

def find_nonzero_range(image_files, organ_path):
    start_index = None
    end_index = None

    for i, image_file in enumerate(image_files):
        image_path = os.path.join(organ_path, image_file)
        img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

        if not is_black(img):
            if start_index is None:
                start_index = i
            end_index = i

    return start_index, end_index

def search_cs_each_organ(cwd):
    output_file = os.path.join(cwd, "organ_cs_list.txt")
    with open(output_file, 'w') as f:
        f.write(f"start deformation list cs num {start_panc_cs_num}\n")  # Pancreas?????
    with open(output_file, 'a') as f:
        for organ_dir in os.listdir(os.path.join(cwd, "label")):
            organ_name = organ_dir.split('_')[0]  # ???????????????
            if organ_name == "Left":
                organ_name = "Left_Kidney"  # Left?Left_Kidney???
            organ_path = os.path.join(cwd, "label", organ_dir)
            image_files = sorted([f for f in os.listdir(organ_path) if f.endswith('.tiff')])

            start_index, end_index = find_nonzero_range(image_files, organ_path)

            if start_index is not None and end_index is not None:
                start_image = image_files[start_index]
                end_image = image_files[end_index]

                start_num = int(start_image.split('IMG')[1].split('.')[0])
                end_num = int(end_image.split('IMG')[1].split('.')[0])

                f.write(f"{organ_name} {start_num} {end_num}\n")
            else:
                print(f"No non-black images found for {organ_name}")

                
if __name__ == '__main__':
    args = sys.argv
    if len(args) == 2:
        cwd = args[1]

        base = make_folder(cwd)

        transformed_bsp, bsp_deformation_field, image_for_input_dice = registration_Axial_bspline(cwd)

        info_path = os.path.join(cwd, "info.txt")
        with open(info_path, 'r') as read_info:
            lines = read_info.read()
        pattern = r'\(Space (\d+\.\d+) (\d+\.\d+) (\d+\.\d+)\)'
        match = re.search(pattern, lines)
        if match:
            spacing_x = float(match.group(1))
            spacing_y = float(match.group(2))
            spacing_z = float(match.group(3))
        else:
            print("couldn't find spacing")

        registration_param_path = "/home/ec2-user/Desktop/research_handover_eguchi/param/2D_elastix_Bspline.txt"#base + "input/2D_elastix_Bspline.txt"  ##???????????
        with open(registration_param_path, 'r') as regi_param_file:
            read_regi_param_lines = regi_param_file.readlines()
        grid_spacing = f"{spacing_z} {spacing_z}"
        for i, line in enumerate(read_regi_param_lines):
            if line.startswith("(FinalGridSpacingInVoxels "):               #1????????
                read_regi_param_lines[i] = f"(FinalGridSpacingInVoxels {grid_spacing})\n" 
                break
        with open(registration_param_path, 'w') as file:
            file.writelines(read_regi_param_lines)

        folder_path = os.path.join(cwd, "base_data")
        number_of_CrossSection = count_non_black_images(folder_path)
        print("number of cross section is " + str(number_of_CrossSection))
        #GridSize??
        with open(registration_param_path, 'r') as regi_param_file:
            read_regi_param_lines_GridSize = regi_param_file.readlines()
        grid_size = f"{number_of_CrossSection} {number_of_CrossSection}"
        for i, line in enumerate(read_regi_param_lines_GridSize):
            if line.startswith("(GridSize "):
                read_regi_param_lines_GridSize[i] = f"(GridSice {grid_size})\n" 
                break
        with open(registration_param_path, 'w') as file:
            file.writelines(read_regi_param_lines_GridSize)

        #Axi??????
        with Image.open(driven_img_fixed) as Axi_img:
            print("Fixed image: ", Axi_img.size)
            Axi_img_W, Axi_img_H = Axi_img.size
        Axi_grid_spacing_W = Axi_img_W / number_of_CrossSection
        Axi_grid_spacing_H = Axi_img_H / number_of_CrossSection
        registration_Axial(Axi_grid_spacing_W, Axi_grid_spacing_H, registration_param_path, base, transformed_bsp)

        #Cor??????
        #?????????
        with Image.open(Cor_fixed_img) as Cor_img:
            Cor_img_H, Cor_img_W = Cor_img.size
        Cor_grid_spacing_W = Cor_img_W / number_of_CrossSection
        Cor_grid_spacing_H = Cor_img_H / number_of_CrossSection
        registration_Cor_output = os.path.join(base, "registration_output_cor")
        registration_Cornal(Cor_grid_spacing_W, Cor_grid_spacing_H, registration_param_path, registration_Cor_output)
        Cornal_transform_parameters = os.path.join(registration_Cor_output, "TransformParameters.0.R0.txt")
        number_of_Bsp_param = extract_number_of_parameters(Cornal_transform_parameters)
        number_of_control_point_per_low = int(np.sqrt(number_of_Bsp_param / 2))
        #????
        deformation_ratio_list = determine_weight(number_of_control_point_per_low, Cornal_transform_parameters, cwd)
        print("deformation ratio list is" + str(deformation_ratio_list)) 
        search_cs_each_organ(cwd)
        #sort_csv()
        dice_GT_init(driven_img_moved, driven_img_fixed)
        dice_GT_result(driven_img_moved, transformed_bsp)
        
    else:
        print("invalid input")

