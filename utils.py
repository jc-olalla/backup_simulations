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

def rotate_points(df, alpha_degrees, x_ref=0.0, y_ref=0.0):
    # Step 1: Translate points to origin
    xy_array = df[['x', 'y']].to_numpy()
    xy_centered = xy_array - np.array([x_ref, y_ref])

    # Step 2: Rotate
    xy_rotated = _rotate_points(xy_centered, alpha_degrees)

    # Step 3: Translate back
    xy_final = xy_rotated + np.array([x_ref, y_ref])

    # Update dataframe
    df['x'] = xy_final[:, 0]
    df['y'] = xy_final[:, 1]

    return df


