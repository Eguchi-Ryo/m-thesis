import SimpleITK as sitk
import os
import numpy as np
from scipy.ndimage import binary_dilation, generate_binary_structure
import subprocess
import sys
import re
from PIL import Image
import shutil

#####select cross-section to make skelton############  
axial_plane = 106  #zの位置を指定（体の高さ方向の位置）
sagittal_plane = 170 #xの位置を指定（体の左右方向の位置）
coronal_plane = 210 #yの位置を指定（体の奥行き方向の位置）
######################

###elastix path ###
dyld_path = "DYLD_LIBRARY_PATH=/Users/lelab/Downloads/skelton_registration_mac-main/elastix/lib/libANNlib-5.1.1.dylib"
elastix_path = "/Users/lelab/Downloads/skelton_registration_mac-main/elastix/bin/elastix"#"/home/ec2-user/Desktop/research_handover/elastix/bin/elastix"
registration_param_path = "/Users/lelab/Downloads/skelton_registration_mac-main/param/2D_elastix.txt"
registration_param_path_mask = "/Users/lelab/Downloads/3d_registratoion_from_3cs/param/registration_3d_mask.txt"
###################

def make_dir(wd):
    os.chdir(wd)
    if os.path.exists("registration"):
        shutil.rmtree("registration")
    os.mkdir("registration")
    os.chdir("registration")
    os.mkdir("param")
    os.mkdir("axial")
    os.mkdir("coronal")
    os.mkdir("sagittal")

def read_image(image_path):
    """Reads a single image and returns it as a numpy array."""
    image = sitk.ReadImage(image_path)
    return sitk.GetArrayFromImage(image)

def stack_images(directory):
    #Stacks all images in a directory in order of slice number.
    # Get list of files in directory
    files = os.listdir(directory)
    
    # Filter files to keep only those matching the format IMGxxxx.tiff
    pattern = re.compile(r'IMG(\d{4})\.tiff$', re.IGNORECASE)
    files = [file for file in files if pattern.match(file)]
    
    # Sort files based on slice number in filename
    files.sort(key=lambda x: int(pattern.match(x).group(1)))
    
    # Read each image and store in a list
    image_stack = []
    for file in files:
        image_path = os.path.join(directory, file)
        image_array = read_image(image_path)
        image_stack.insert(0, image_array)
  
    # Convert list to numpy array for 3D stack
    image_stack = np.array(image_stack)
    print("stack image shape: "+str(image_stack.shape))
    
    return image_stack

def construct_3d_object(image_stack):
    #Constructs a 3D object from a stack of images using original pixel values.
    # Get shape of the image stack
    depth, height, width = image_stack.shape
    
    # Initialize 3D object array
    object_3d = np.zeros((depth, height, width), dtype=image_stack.dtype)
    image_type = image_stack.dtype
    
    # Iterate through each slice and copy pixel values to 3D object array
    for z in range(depth):
        object_3d[z, :, :] = image_stack[z, :, :]
    
    return object_3d, image_type

def extract_plane_from_3d_object(object_3d, plane, position=0):
    #Extracts a plane from a 3D object based on specified axis and position.
    #object3d[z,x,y]の順番になっている
    if plane == 'xy':
        plane_data = object_3d[position, :, :] #axial
    elif plane == 'xz':
        plane_data = object_3d[:, position, :] #coronal
    elif plane == 'yz':
        plane_data = object_3d[:, :, position] #sagittal
    else:
        raise ValueError(f"Unsupported plane: {plane}. Choose from 'xy', 'xz', or 'yz'.")
    
    #return np.flipud(plane_data)  # Y????  反転させてはいけない。要注意!!!!!
    return plane_data

def make_directory_AxiCorSag(base):
    os.makedirs(os.path.join(base, "axial"), exist_ok=True)
    os.makedirs(os.path.join(base, "coronal"), exist_ok=True)
    os.makedirs(os.path.join(base, "sagittal"), exist_ok=True)
    
def save_as_image(path, plane_data, type):
    plane_data = plane_data.astype(type)
    
    # Pillow で画像に変換して保存します
    image = Image.fromarray(plane_data)
    file_full_path = os.path.join(path, "cross_section_image.tiff")
    image.save(file_full_path)
    
def process_directory(base_dir, axial_plane, sagittal_plane, coronal_plane):
    #Processes images in the given directory to create and save 3D objects and planes.
    #before, afterの二つを指定
    axial_dir = os.path.join(base_dir) #base_dir でpatient/image, patient/ans_dataまで行く
    
    #axial, coronal, sagittalのディレクトリを作る
    make_directory_AxiCorSag(base_dir)
    
    # Stack images
    image_stack = stack_images(axial_dir) 
    
    # Construct 3D object
    object_3d, image_type = construct_3d_object(image_stack)
    
    # Save 3D object
    save_path = os.path.join(base_dir, '3d_object_stack.nrrd')
    itk_image = sitk.GetImageFromArray(object_3d)
    itk_image.SetDirection([1, 0, 0, 0, 1, 0, 0, 0, 1])  # ??????????
    itk_image.SetOrigin((0, 0, 0))
    itk_image.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(itk_image, save_path)
    print(f"3D object created and saved to {save_path}")
    
    # Extract planes
    axi_pos = axial_plane
    cor_pos = object_3d.shape[2] - coronal_plane  #y軸（体の奥行き方向） 
    sag_pos = object_3d.shape[1] - sagittal_plane   #x軸（体の左右方向）
    plane_xy = extract_plane_from_3d_object(object_3d, plane='xy', position=axi_pos)
    axi_image_path = os.path.join(base_dir, "axial")
    save_as_image(axi_image_path, plane_xy, image_type)
    plane_xz = extract_plane_from_3d_object(object_3d, plane='xz', position=cor_pos)
    cor_image_path = os.path.join(base_dir, "coronal")
    save_as_image(cor_image_path, plane_xz, image_type)
    plane_yz = extract_plane_from_3d_object(object_3d, plane='yz', position=sag_pos)
    sag_image_path = os.path.join(base_dir, "sagittal")
    save_as_image(sag_image_path, plane_yz, image_type)
    
    # Create new 3D object using the extracted planes
    new_3d_object = np.zeros(object_3d.shape, dtype=object_3d.dtype)
    new_3d_object[axi_pos, :, :] = plane_xy
    print("plane_xy shape: "+str(new_3d_object[axi_pos, :, :].shape))
    new_3d_object[:, cor_pos, :] = plane_xz
    print("plane_xz shape: "+str(new_3d_object[:, cor_pos, :].shape))
    new_3d_object[:, :, sag_pos] = plane_yz
    
def load_nrrd_for_mask(file_path):
    image = sitk.ReadImage(file_path)
    array = sitk.GetArrayFromImage(image)
    return image, array

def save_mask(mask, reference_image, output_file_path):
    mask_image = sitk.GetImageFromArray(mask.astype(np.uint8))
    mask_image.CopyInformation(reference_image)
    sitk.WriteImage(mask_image, output_file_path)

def make_mask(wd):
    nrrd1_path = os.path.join(wd, "image", "combined_object.nrrd")
    nrrd2_path = os.path.join(wd, "ans_data", "combined_object.nrrd")
    image1, array1 = load_nrrd_for_mask(nrrd1_path)
    image2, array2 = load_nrrd_for_mask(nrrd2_path)

    # ????????
    mask1 = array1 != 0
    mask2 = array2 != 0

    # OR??????????
    combined_mask = np.logical_or(mask1, mask2)

    # 3????????????
    structure = generate_binary_structure(rank=3, connectivity=1)

    # ?????????????????????????
    dilated_mask = binary_dilation(combined_mask, structure=structure, iterations=15)

    # ?????
    image_shape = combined_mask.shape
    dilated_mask[0, :, :] = 0
    dilated_mask[-1, :, :] = 0
    dilated_mask[:, 0, :] = 0
    dilated_mask[:, -1, :] = 0
    dilated_mask[:, :, 0] = 0
    dilated_mask[:, :, -1] = 0

    # ????????
    mask_dest_path = os.path.join(wd, "registration", "param", "mask.nrrd")
    save_mask(dilated_mask, image1, mask_dest_path)

def make_deformation_field(cs, wd):
    os.environ["DYLD_LIBRARY_PATH"] = "/usr/local/lib"
    fx = os.path.join(wd, "ans_data", cs, "cross_section_image.tiff")
    mv = os.path.join(wd, "image", cs, "cross_section_image.tiff")
    mask = os.path.join(wd, "registration", "param", "mask.nrrd")
    out = os.path.join(wd, "registration", cs)

    # Create command
    el_cmd = [
        "/Users/lelab/Downloads/skelton_registration_mac-main/elastix/bin/elastix",
        "-f", fx,
        "-m", mv,
        "-p", registration_param_path,
        "-out", out
    ]

    el_cmd_use_mask = [
        "/Users/lelab/Downloads/skelton_registration_mac-main/elastix/bin/elastix",
        "-f", fx,
        "-fMask", mask,
        "-m", mv,
        "-mMask", mask,
        "-p", registration_param_path_mask,
        "-out", out
    ]
    # Execute command
    subprocess.run(el_cmd, check=True)
    
    transformix_param = os.path.join(out, "TransformParameters.0.txt")
    tr_cmd = [
        "/Users/lelab/Downloads/skelton_registration_mac-main/elastix/bin/transformix",
        "-def", "all",
        "-in", mv,
        "-out", out,
        "-tp", transformix_param
    ]
    subprocess.run(tr_cmd, check=True)

def load_nrrd(file_path):
    return sitk.ReadImage(file_path)

def save_nrrd(image, file_path):
    sitk.WriteImage(image, file_path)

def calculate_dice_score(image1, image2):
    array1 = sitk.GetArrayFromImage(image1) > 0
    array2 = sitk.GetArrayFromImage(image2) > 0
    
    intersection = np.logical_and(array1, array2)
    volume_sum = array1.sum() + array2.sum()
    
    if volume_sum == 0:
        return 1.0
    
    return 2.0 * intersection.sum() / volume_sum

def make_df_voxels(initialize, cs, wd):
    width, height, depth =  initialize.GetSize()
    df_voxel = np.zeros((depth, width, height, 2), dtype=np.float64)
    df_path = os.path.join(wd, "registration", cs, "deformationField.mhd")
    df_image = sitk.ReadImage(df_path)
    df_image = sitk.Cast(df_image, sitk.sitkVectorFloat64)
    df_array = sitk.GetArrayFromImage(df_image)
    if cs == "axial":
        print("axial df shape: "+str(df_array.shape))
        for i in range(depth):
            df_voxel[i, :, :, :] = df_array
    elif cs == "coronal":
        print("coronal df shape: "+str(df_array.shape))
        for i in range(width):
            df_voxel[:, i, :, :] = df_array
    else:
        print("sagittal df shape: "+str(df_array.shape))
        for i in range(height):
            df_voxel[:, :, i, :] = df_array
            
    df_image_sitk = sitk.GetImageFromArray(df_voxel)
    df_image_sitk.SetSpacing(initialize.GetSpacing())
    df_image_sitk.SetOrigin(initialize.GetOrigin())
    df_image_sitk.SetDirection(initialize.GetDirection())
    output_path = os.path.join(wd, "registration", cs, "deformationField_stack.mhd")
    sitk.WriteImage(df_image_sitk, output_path)
            
def compose_voxels(initialize, wd):
    width, height, depth =  initialize.GetSize()
    df_voxel_3D = np.zeros((depth, width, height, 3), dtype=np.float64)
    df_axi_stack_path = os.path.join(wd, "registration", "axial", "deformationField_stack.mhd")
    df_cor_stack_path = os.path.join(wd, "registration", "coronal", "deformationField_stack.mhd")
    df_sag_stack_path = os.path.join(wd, "registration", "sagittal", "deformationField_stack.mhd")

    df_axi = sitk.GetArrayFromImage(sitk.ReadImage(df_axi_stack_path))
    df_voxel_3D[:, :, :, 0] += df_axi[:, :, :, 0]  # X 方向の変位
    df_voxel_3D[:, :, :, 1] += df_axi[:, :, :, 1]  # Y 方向の変位
    df_cor = sitk.GetArrayFromImage(sitk.ReadImage(df_cor_stack_path))
    df_voxel_3D[:, :, :, 1] += df_cor[:, :, :, 0]  # Y 方向の変位
    df_voxel_3D[:, :, :, 2] += df_cor[:, :, :, 1]  # Z 方向の変位
    df_sag = sitk.GetArrayFromImage(sitk.ReadImage(df_sag_stack_path))
    df_voxel_3D[:, :, :, 0] += df_sag[:, :, :, 0]  # X 方向の変位
    df_voxel_3D[:, :, :, 2] += df_sag[:, :, :, 1]  # Z 方向の変位
    
    df_voxel_3D /= 2

    df_image_sitk = sitk.GetImageFromArray(df_voxel_3D)
    df_image_sitk.SetSpacing(initialize.GetSpacing())
    df_image_sitk.SetOrigin(initialize.GetOrigin())
    df_image_sitk.SetDirection(initialize.GetDirection())
    output_path = os.path.join(wd, "registration", "deformationField_3D.mhd")
    sitk.WriteImage(df_image_sitk, output_path)
    print(f"Deformation field saved as {output_path}")
    
def adopt_3D_deformation_field(before_image_path, wd):
    df3d_path = os.path.join(wd, "registration", "deformationField_3D.mhd")
    image = sitk.ReadImage(before_image_path)
    df3d = sitk.ReadImage(df3d_path, sitk.sitkVectorFloat64)

    img_ratio = (0.00173611, 0.00173611, 0.00263889)
    spacing = (0.78, 0.78, 1.2)

    #image.SetSpacing(spacing)
    #image.SetSpacing(img_ratio)
    #df3d.SetSpacing(img_ratio)

    df_array = sitk.GetArrayFromImage(df3d)
    disp1 = df_array[0,0,0]
    disp2 = df_array[232, 0, 0]

    df3d_array = sitk.GetArrayFromImage(df3d)  # (z, y, x, vector)
    displacement = df3d_array[66, 244, 283]  # [dx, dy, dz] の順番で取得
    dx, dy, dz = displacement
    print(f"Displacement at (244, 283, 66): (dx, dy, dz) = ({dx}, {dy}, {dz})")

    print(f"Image Direction: {image.GetDirection()}")
    print(f"Image Spacing: {image.GetSpacing()}")
    print(f"Image Origin: {image.GetOrigin()}")

    print(f"Displacement Field Direction: {df3d.GetDirection()}")
    print(f"Displacement Field Spacing: {df3d.GetSpacing()}")
    print(f"Displacement Field Origin: {df3d.GetOrigin()}")
    print(f"Displacement Field shape: {df3d.GetSize()}")

    
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(image)
    #resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetTransform(sitk.DisplacementFieldTransform(df3d))
    
    deformed_image = resampler.Execute(image)

    
    output_path = os.path.join(wd, "registration", "transfromed.nrrd")
    sitk.WriteImage(deformed_image, output_path)
    
    return deformed_image
    
if __name__ == '__main__':
    args = sys.argv
    if len(args) == 2:
        working_dir = args[1]
        print("working dir: " + str(working_dir))
    else:
        print("invalid input")

    make_dir(working_dir)

    # Directories for day1 and day2
    base_dirs = [os.path.join(working_dir, "image"), os.path.join(working_dir, "ans_data")]
    
    # Process each directory
    for base_dir in base_dirs:
        process_directory(base_dir, axial_plane, sagittal_plane, coronal_plane)

    #make mask
    make_mask(working_dir)

    cs_list = ["axial", "coronal", "sagittal"]
    for cs in cs_list:
        make_deformation_field(cs, working_dir)
    
    #apply affine and calc dice
    transform_params_path = os.path.join(working_dir, "registration", "TransformParameters.0.txt")#"/home/ec2-user/Desktop/3d_skelton/test_regi/TransformParameters.0.txt"
    fixed_image_path = os.path.join(working_dir, "ans_data", "3d_object_stack.nrrd")#"/home/ec2-user/Desktop/3d_skelton/patient1/day2/3d_object_stack.nrrd"
    moving_image_path = os.path.join(working_dir, "image", "3d_object_stack.nrrd")#"/home/ec2-user/Desktop/3d_skelton/patient1/day1/3d_object_stack.nrrd"
    fixed_image = load_nrrd(fixed_image_path)
    print("fx image(initialize) origin: "+str(fixed_image.GetOrigin()))
    print("fx image(initialize) direction: "+str(fixed_image.GetDirection()))
    moving_image = load_nrrd(moving_image_path)
    #print("3d voxel size: " + str(fixed_image.GetSize()))  #<- (576, 576, 233)
    
    #変位場のstackと合成
    for cs in cs_list:
        make_df_voxels(fixed_image, cs, working_dir)
    compose_voxels(fixed_image, working_dir)
    transformed_image = adopt_3D_deformation_field(moving_image_path, working_dir)


    # Calculate the Dice score between the fixed image and the transformed moving image
    ini_dice = calculate_dice_score(fixed_image, moving_image)
    print(f"Dice Score init-ans: {ini_dice:.4f}")

    dice_score = calculate_dice_score(fixed_image, transformed_image)
    print(f"Dice Score trans-ans: {dice_score:.4f}")

