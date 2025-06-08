import glob
import subprocess
import numpy as np


# Delete AOA directories function
def remove_directories(pattern: str):
    """
    Removes directories matching the given pattern.
    """
    # Find all directories matching the pattern
    directories = glob.glob(pattern)

    # If there are directories to remove, run the subprocess
    if directories:
        subprocess.run(["rm", "-r"] + directories, check=True)


def _rotate_points(xy_array, alpha_degrees):
    alpha_rad = np.radians(alpha_degrees)
    rotation_matrix = np.array([
        [np.cos(alpha_rad), -np.sin(alpha_rad)],
        [np.sin(alpha_rad),  np.cos(alpha_rad)]
    ])
    return xy_array @ rotation_matrix.T  # Transpose for correct orientation

def rotate_points(df, alpha_degrees):
    xy_array = df[['x', 'y']].to_numpy()
    xy_array_rot = _rotate_points(xy_array, alpha_degrees)
    df.x = xy_array_rot[:, 0]
    df.y = xy_array_rot[:, 1]
    return df



