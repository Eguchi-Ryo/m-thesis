import SimpleITK as sitk
import os
import numpy as np
from scipy.ndimage import binary_dilation, generate_binary_structure
import subprocess
import sys
import re


#####select cross-section to make skelton############  
axial_plane = 115  #zの位置を指定（体の高さ方向の位置）
sagittal_plane = 170 #xの位置を指定（体の左右方向の位置）
coronal_plane = 180 #yの位置を指定（体の奥行き方向の位置）
######################

###elastix path ###
dyld_path = "DYLD_LIBRARY_PATH=/Users/lelab/Downloads/skelton_registration_mac-main/elastix/lib/libANNlib-5.1.1.dylib"
elastix_path = "/Users/lelab/Downloads/skelton_registration_mac-main/elastix/bin/elastix"#"/home/ec2-user/Desktop/research_handover/elastix/bin/elastix"
registration_param_path = "/Users/lelab/Downloads/skelton_registration_mac-main/param/3D_elastix.txt"
registration_param_path_mask = "/Users/lelab/Downloads/3d_registratoion_from_3cs/param/registration_3d_mask.txt"
###################

def make_dir(wd):
    os.chdir(wd)
    if not os.path.exists("registration"):
        os.mkdir("registration")
        os.chdir("registration")
        os.mkdir("param")

def read_image(image_path):
    """Reads a single image and returns it as a numpy array."""
    image = sitk.ReadImage(image_path)
    return sitk.GetArrayFromImage(image)

def stack_images(directory):
    """Stacks all images in a directory in order of slice number."""
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
        image_stack.append(image_array)
    
    # Convert list to numpy array for 3D stack
    image_stack = np.array(image_stack)
    
    return image_stack

def construct_3d_object(image_stack):
    #Constructs a 3D object from a stack of images using original pixel values.
    # Get shape of the image stack
    depth, height, width = image_stack.shape
    
    # Initialize 3D object array
    object_3d = np.zeros((depth, height, width), dtype=image_stack.dtype)
    
    # Iterate through each slice and copy pixel values to 3D object array
    for z in range(depth):
        object_3d[z, :, :] = image_stack[z, :, :]
    
    return object_3d

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

def process_directory(base_dir, axial_plane, sagittal_plane, coronal_plane):
    #Processes images in the given directory to create and save 3D objects and planes.
    # Axial directory
    axial_dir = os.path.join(base_dir)
    
    # Stack images
    image_stack = stack_images(axial_dir)
    
    # Construct 3D object
    object_3d = construct_3d_object(image_stack)
    
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
    plane_xz = extract_plane_from_3d_object(object_3d, plane='xz', position=cor_pos)
    plane_yz = extract_plane_from_3d_object(object_3d, plane='yz', position=sag_pos)
    
    # Create new 3D object using the extracted planes
    new_3d_object = np.zeros(object_3d.shape, dtype=object_3d.dtype)
    new_3d_object[axi_pos, :, :] = plane_xy
    new_3d_object[:, cor_pos, :] = plane_xz
    new_3d_object[:, :, sag_pos] = plane_yz
    
    # Save new 3D object
    new_save_path = os.path.join(base_dir, 'new_3d_object_planes.nrrd')
    itk_new_image = sitk.GetImageFromArray(new_3d_object)
    itk_new_image.SetDirection([1, 0, 0, 0, 1, 0, 0, 0, 1]) 
    itk_new_image.SetOrigin((0, 0, 0))
    itk_new_image.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(itk_new_image, new_save_path)
    print(f"New 3D object created and saved to {new_save_path}")

    #Make and save combined 3D object(from 3 cross-section)
    combined_3d_object = extrude_and_combine_planes(plane_xy, plane_xz, plane_yz, object_3d) 
    combined_obj_path = os.path.join(base_dir, 'combined_object.nrrd')
    itk_new_image = sitk.GetImageFromArray(combined_3d_object)
    itk_new_image.SetDirection([1, 0, 0, 0, 1, 0, 0, 0, 1]) 
    itk_new_image.SetOrigin((0, 0, 0))
    itk_new_image.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(itk_new_image, combined_obj_path)
    print(f"New 3D object created and saved to {combined_obj_path}")

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

def registration(wd):
    os.environ["DYLD_LIBRARY_PATH"] = "/usr/local/lib"
    fx = os.path.join(wd, "ans_data", "combined_object.nrrd")
    mv = os.path.join(wd, "image", "combined_object.nrrd")
    mask = os.path.join(wd, "registration", "param", "mask.nrrd")
    out = os.path.join(wd, "registration")

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

def load_nrrd(file_path):
    """Load a NRRD file and return the image."""
    return sitk.ReadImage(file_path)

def save_nrrd(image, file_path):
    """Save an image as a NRRD file."""
    sitk.WriteImage(image, file_path)

def apply_affine_transform(image, transform_parameters):
    """Apply an affine transform to the input image."""
    # Create the affine transform
    transform = sitk.AffineTransform(3)
    transform.SetParameters(transform_parameters)
    
    # Resample the image using the transform
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(image)
    resampler.SetTransform(transform)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(0)
    
    return resampler.Execute(image)

def calculate_dice_score(image1, image2):
    """Calculate the Dice score between two images."""
    array1 = sitk.GetArrayFromImage(image1) > 0
    array2 = sitk.GetArrayFromImage(image2) > 0
    
    intersection = np.logical_and(array1, array2)
    volume_sum = array1.sum() + array2.sum()
    
    if volume_sum == 0:
        return 1.0
    
    return 2.0 * intersection.sum() / volume_sum

def extract_transform_parameters(file_path):
    """Extract transform parameters from the elastix parameter file."""
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("(TransformParameters"):
                parameters = line.strip().split(" ")[1:]  # Split and remove the initial label
                parameters[-1] = parameters[-1].rstrip(')')  # Remove trailing parenthesis from the last parameter
                return list(map(float, parameters))
    return None

#3平面のAND算でオブジェクトを構成
"""
def extrude_and_combine_planes(plane_xy, plane_xz, plane_yz, object3d):
    #Combines extruded planes into a single 3D object using logical AND.
    # Initialize extruded planes
    extruded_xy = np.repeat(plane_xy[np.newaxis, :, :], object3d.shape[0], axis=0)
    extruded_xz = np.repeat(plane_xz[:, np.newaxis, :], object3d.shape[1], axis=1)
    extruded_yz = np.repeat(plane_yz[:, :, np.newaxis], object3d.shape[2], axis=2)

    # Combine using logical AND
    combined_object_mask = np.logical_and(extruded_xy, extruded_xz)
    combined_object_mask = np.logical_and(combined_object_mask, extruded_yz)

    combined_object = np.where(combined_object_mask, object3d, 0)

    return combined_object.astype(object3d.dtype) 
"""
"""
#2平面のAND算でオブジェクトを構成
def extrude_and_combine_planes(plane_xy, plane_xz, plane_yz, object3d):
    # Initialize extruded planes with pixel values, filling with zeros where there are no values
    extruded_xy = np.repeat(plane_xy[np.newaxis, :, :], object3d.shape[0], axis=0)
    extruded_xz = np.repeat(plane_xz[:, np.newaxis, :], object3d.shape[1], axis=1)
    extruded_yz = np.repeat(plane_yz[:, :, np.newaxis], object3d.shape[2], axis=2)

    # Initialize a sum array and a count array to calculate the average
    sum_values = np.zeros(object3d.shape, dtype=float)
    count_values = np.zeros(object3d.shape, dtype=int)

    # Add values from each extruded plane and count overlaps
    non_zero_xy = extruded_xy != 0
    non_zero_xz = extruded_xz != 0
    non_zero_yz = extruded_yz != 0

    sum_values += np.where(non_zero_xy, extruded_xy, 0)
    count_values += non_zero_xy.astype(int)

    sum_values += np.where(non_zero_xz, extruded_xz, 0)
    count_values += non_zero_xz.astype(int)

    sum_values += np.where(non_zero_yz, extruded_yz, 0)
    count_values += non_zero_yz.astype(int)

    # Calculate the average where there is at least two overlaps
    valid_overlap = count_values >= 2
    average_values = np.zeros(object3d.shape, dtype=float)
    average_values[valid_overlap] = sum_values[valid_overlap] / count_values[valid_overlap]

    return average_values
"""

def extrude_and_combine_planes(plane_xy, plane_xz, plane_yz, object3d):
    # Initialize extruded planes
    extruded_xy = np.repeat(plane_xy[np.newaxis, :, :], object3d.shape[0], axis=0)
    extruded_xz = np.repeat(plane_xz[:, np.newaxis, :], object3d.shape[1], axis=1)
    extruded_yz = np.repeat(plane_yz[:, :, np.newaxis], object3d.shape[2], axis=2)

    # Calculate logical AND for each pair
    overlap_xy_xz = np.logical_and(extruded_xy, extruded_xz)
    overlap_xy_yz = np.logical_and(extruded_xy, extruded_yz)
    overlap_xz_yz = np.logical_and(extruded_xz, extruded_yz)

    # Combine using logical OR to include regions where at least two overlap
    combined_object_mask = np.logical_or(overlap_xy_xz, overlap_xy_yz)
    combined_object_mask = np.logical_or(combined_object_mask, overlap_xz_yz)
    
    object3d_before = os.path.join(working_dir, "image", "3d_object_stack.nrrd")
    print("object3d before path: "+object3d_before)
    object3d_before_image = sitk.ReadImage(object3d_before)
    object3d_before_array = sitk.GetArrayFromImage(object3d_before_image)
    
    combined_object = np.where(combined_object_mask, object3d_before_array, 0)

    return combined_object.astype(object3d.dtype)

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

    #registration
    registration(working_dir)
    
    #apply affine and calc dice
    transform_params_path = os.path.join(working_dir, "registration", "TransformParameters.0.txt")#"/home/ec2-user/Desktop/3d_skelton/test_regi/TransformParameters.0.txt"
    fixed_image_path = os.path.join(working_dir, "ans_data", "3d_object_stack.nrrd")#"/home/ec2-user/Desktop/3d_skelton/patient1/day2/3d_object_stack.nrrd"
    moving_image_path = os.path.join(working_dir, "image", "3d_object_stack.nrrd")#"/home/ec2-user/Desktop/3d_skelton/patient1/day1/3d_object_stack.nrrd"

    # Load the transform parameters
    transform_parameters = extract_transform_parameters(transform_params_path)
    if transform_parameters is None:
        raise ValueError("Transform parameters not found in the specified file.")

    # Load the images
    fixed_image = load_nrrd(fixed_image_path)
    moving_image = load_nrrd(moving_image_path)

    # Apply the affine transform to the moving image
    transformed_moving_image = apply_affine_transform(moving_image, transform_parameters)

    # Save the transformed image (optional)
    transformed_image_path = os.path.join(working_dir, "image", "transformed_3d_object_stack.nrrd")#"/home/ec2-user/Desktop/3d_skelton/patient1/day1/transformed_3d_object_stack.nrrd"
    save_nrrd(transformed_moving_image, transformed_image_path)

    # Calculate the Dice score between the fixed image and the transformed moving image
    ini_dice = calculate_dice_score(fixed_image, moving_image)
    print(f"Dice Score init-ans: {ini_dice:.4f}")

    dice_score = calculate_dice_score(fixed_image, transformed_moving_image)
    print(f"Dice Score trans-ans: {dice_score:.4f}")
