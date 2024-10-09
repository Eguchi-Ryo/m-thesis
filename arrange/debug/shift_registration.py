import SimpleITK as sitk
import os
import numpy as np
from scipy.ndimage import binary_dilation, generate_binary_structure, shift
import subprocess
import sys
import re
from PIL import Image
import nrrd


#####select cross-section to make skelton############  
axial_plane = 90  #zの位置を指定（体の高さ方向の位置）
sagittal_plane = 150 #xの位置を指定（体の左右方向の位置）
coronal_plane = 170 #yの位置を指定（体の奥行き方向の位置）
######################

###elastix path ###
dyld_path = "DYLD_LIBRARY_PATH=/Users/lelab/Downloads/skelton_registration_mac-main/elastix/lib/libANNlib-5.1.1.dylib"
elastix_path = "/Users/lelab/Downloads/skelton_registration_mac-main/elastix/bin/elastix"#"/home/ec2-user/Desktop/research_handover/elastix/bin/elastix"
registration_param_path = "/Users/lelab/Downloads/skelton_registration_mac-main/param/2D_elastix.txt"
registration_param_path_mask = "/Users/lelab/Downloads/3d_registratoion_from_3cs/param/registration_3d_mask.txt"
###################

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

def shift_vol(stack, wd):
    shift_vec = [20, 0, 0]  #(z,y,x) ここで平行移動の移動量を定義
    shifted_stack = shift(stack, shift_vec)
    save_path = os.path.join(wd, "debug", "shifted_stack.nrrd")
    itk_image = sitk.GetImageFromArray(shifted_stack)
    itk_image.SetDirection([1, 0, 0, 0, 1, 0, 0, 0, 1])  # ??????????
    itk_image.SetOrigin((0, 0, 0))
    itk_image.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(itk_image, save_path)
    
def load_nrrd(file_path):
    """Load a NRRD file and return the image."""
    return sitk.ReadImage(file_path)

def load_nrrd_for_mask(file_path):
    image = sitk.ReadImage(file_path)
    array = sitk.GetArrayFromImage(image)
    return image, array

def save_mask(mask, reference_image, output_file_path):
    mask_image = sitk.GetImageFromArray(mask.astype(np.uint8))
    mask_image.CopyInformation(reference_image)
    sitk.WriteImage(mask_image, output_file_path)
    
def make_mask(wd, fx, mv):
    image1 = Image.open(fx)
    array1 = np.array(Image.open(fx))
    array2 = np.array(Image.open(mv))

    # ????????
    mask1 = array1 != 0
    mask2 = array2 != 0

    # OR??????????
    combined_mask = np.logical_or(mask1, mask2)

    # 3????????????
    structure = generate_binary_structure(rank=2, connectivity=1)

    # ?????????????????????????
    dilated_mask = binary_dilation(combined_mask, structure=structure, iterations=15)

    # ?????
    image_shape = combined_mask.shape
    dilated_mask[0, :] = 0
    dilated_mask[-1, :] = 0
    dilated_mask[:, 0] = 0
    dilated_mask[:, -1] = 0
    
    reference_image = sitk.GetImageFromArray(array1)
    # ????????
    mask_dest_path = os.path.join(wd, "registration", "param", "mask.nrrd")
    save_mask(dilated_mask, reference_image, mask_dest_path)
    
def make_deformation_field(cs, wd):
    os.environ["DYLD_LIBRARY_PATH"] = "/usr/local/lib"
    fx = os.path.join(wd, "debug", cs, "cross_section_image.tiff")
    mv = os.path.join(wd, "image", cs, "cross_section_image.tiff")
    make_mask(wd, fx, mv)
    mask = os.path.join(wd, "registration", "param", "mask.nrrd")
    out = os.path.join(wd, "debug","registration", cs)

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
        "-p", registration_param_path,#registration_param_path_mask,
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
    
def save_as_image(path, plane_data, type):
    plane_data = plane_data.astype(type)
    
    # Pillow で画像に変換して保存します
    image = Image.fromarray(plane_data)
    file_full_path = os.path.join(path, "cross_section_image.tiff")
    image.save(file_full_path)
    
def make_df_voxels(initialize, cs, wd):
    width, height, depth =  initialize.GetSize()
    df_voxel = np.zeros((depth, width, height, 2), dtype=np.float64)
    df_path = os.path.join(wd, "debug","registration", cs, "deformationField.mhd")
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
    output_path = os.path.join(wd, "debug","registration", cs, "deformationField_stack.mhd")
    sitk.WriteImage(df_image_sitk, output_path)
    
def compose_voxels(initialize, wd):
    width, height, depth =  initialize.GetSize()
    df_voxel_3D = np.zeros((depth, width, height, 3), dtype=np.float64)
    df_axi_stack_path = os.path.join(wd, "debug", "registration", "axial", "deformationField_stack.mhd")
    df_cor_stack_path = os.path.join(wd, "debug","registration", "coronal", "deformationField_stack.mhd")
    df_sag_stack_path = os.path.join(wd, "debug","registration", "sagittal", "deformationField_stack.mhd")

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
    output_path = os.path.join(wd, "debug", "registration", "deformationField_3D.mhd")
    sitk.WriteImage(df_image_sitk, output_path)
    print(f"Deformation field saved as {output_path}")
    
def adopt_3D_deformation_field(before_image_path, wd):
    df3d_path = os.path.join(wd, "debug", "registration", "deformationField_3D.mhd")
    image = sitk.ReadImage(before_image_path)
    df3d = sitk.ReadImage(df3d_path, sitk.sitkVectorFloat64)
    
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(image)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetTransform(sitk.DisplacementFieldTransform(df3d))
    
    deformed_image = resampler.Execute(image)
    
    output_path = os.path.join(wd, "debug", "registration", "transfromed.nrrd")
    sitk.WriteImage(deformed_image, output_path)
    
    return deformed_image

def calculate_dice_score(image1, image2):
    #"""Calculate the Dice score between two images."""
    array1 = sitk.GetArrayFromImage(image1) > 0
    array2 = sitk.GetArrayFromImage(image2) > 0
    
    intersection = np.logical_and(array1, array2)
    volume_sum = array1.sum() + array2.sum()
    
    if volume_sum == 0:
        return 1.0
    
    return 2.0 * intersection.sum() / volume_sum
    
if __name__ == '__main__':
    args = sys.argv
    if len(args) == 2:
        working_dir = args[1]
        print("working dir: " + str(working_dir))
    else:
        print("invalid input")
    image_path = os.path.join(working_dir, "image")
    image_stack = stack_images(image_path)
    shift_vol(image_stack, working_dir)
    
    before_image_path = os.path.join(working_dir, "image", "3d_object_stack.nrrd")
    shifted_image_path = os.path.join(working_dir, "debug", "shifted_stack.nrrd")
    
    before_image = load_nrrd(before_image_path)
    shifted_image = load_nrrd(shifted_image_path)

    # Construct 3D object
    object_3d, image_type = construct_3d_object(image_stack)
    print("object3d(image) shape: "+str(object_3d.shape))
    itk_image = sitk.GetImageFromArray(object_3d)
    itk_image.SetDirection([1, 0, 0, 0, 1, 0, 0, 0, 1])  # ??????????
    itk_image.SetOrigin((0, 0, 0))
    itk_image.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(itk_image, before_image_path)

    # Extract planes
    axi_pos = axial_plane
    cor_pos = object_3d.shape[2] - coronal_plane  #y軸（体の奥行き方向） 
    sag_pos = object_3d.shape[1] - sagittal_plane   #x軸（体の左右方向）
    plane_xy = extract_plane_from_3d_object(object_3d, plane='xy', position=axi_pos)
    axi_image_path = os.path.join(working_dir, "image","axial")
    save_as_image(axi_image_path, plane_xy, image_type)
    plane_xz = extract_plane_from_3d_object(object_3d, plane='xz', position=cor_pos)
    cor_image_path = os.path.join(working_dir, "image","coronal")
    save_as_image(cor_image_path, plane_xz, image_type)
    plane_yz = extract_plane_from_3d_object(object_3d, plane='yz', position=sag_pos)
    sag_image_path = os.path.join(working_dir, "image", "sagittal")
    save_as_image(sag_image_path, plane_yz, image_type)
    
    #shift imageのtypeと配列を取得
    shift_data, header = nrrd.read(shifted_image_path)
    shift_data = np.transpose(shift_data, (2, 1, 0))
    print("shift data shape: "+str(shift_data.shape))
    shift_array = shift_data
    shift_type = shift_data.dtype
    #extract plane in shift data
    axi_pos = axial_plane
    cor_pos = object_3d.shape[2] - coronal_plane  #y軸（体の奥行き方向） 
    sag_pos = object_3d.shape[1] - sagittal_plane   #x軸（体の左右方向）
    plane_xy = extract_plane_from_3d_object(shift_data, plane='xy', position=axi_pos)
    axi_image_path = os.path.join(working_dir, "debug","axial")
    save_as_image(axi_image_path, plane_xy, shift_type)
    plane_xz = extract_plane_from_3d_object(shift_data, plane='xz', position=cor_pos)
    cor_image_path = os.path.join(working_dir, "debug","coronal")
    save_as_image(cor_image_path, plane_xz, shift_type)
    plane_yz = extract_plane_from_3d_object(shift_data, plane='yz', position=sag_pos)
    sag_image_path = os.path.join(working_dir, "debug", "sagittal")
    save_as_image(sag_image_path, plane_yz, shift_type)
    
    cs_list = ["axial", "coronal", "sagittal"]
    for cs in cs_list:
        make_deformation_field(cs, working_dir)
        make_df_voxels(before_image, cs, working_dir)
        
    compose_voxels(before_image, working_dir)
    transformed_image = adopt_3D_deformation_field(before_image_path, working_dir)

    # Calculate the Dice score between the before image and the transformed moving image
    ini_dice = calculate_dice_score(before_image, shifted_image)
    print(f"Dice Score init-ans: {ini_dice:.4f}")

    dice_score = calculate_dice_score(shifted_image, transformed_image)
    print(f"Dice Score trans-ans: {dice_score:.4f}")

    