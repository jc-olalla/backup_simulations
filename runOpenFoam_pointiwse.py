"""Run simulation for each AOA and create necessary post-processing"""
import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt
from rdp import rdp

def compute_alpha_pointwise(Re, alpha, previous_alpha, root_path, profile_name):
    parallel_comp = True
    cp_plot = True
    force_plot = True
    residuals_plot = True   # keep on otherwise cp_plot does not work without last iteration value
    steady = True
    unsteady = False
    initialise_sim = False
    VTK_paraview = True

    if steady:
        ### Create the working directory ###
        subprocess.run(["mkdir", "-p", "openfoam_runs_pointwise/AOA_" + str(alpha)])   # Create a directory AOA_{AOA}
        subprocess.run(["cp", "-r", "openfoam_runs_pointwise/Base_AOA/constant", "openfoam_runs_pointwise/Base_AOA/system", "openfoam_runs_pointwise/AOA_" + str(alpha)])    # Copying Base_AOA folders and content to the last mentioned directory AOA
        subprocess.run(["cp", "-r", "openfoam_runs_pointwise/Base_AOA/polyMesh", "openfoam_runs_pointwise/AOA_" + str(alpha) + "/constant"]) # Copies mesh and places it in the last mentioned directory
        subprocess.run(["cp", "-r", "openfoam_runs_pointwise/Base_AOA/0", "openfoam_runs_pointwise/AOA_" + str(alpha)])     # Copying Base_AOA/0 folderw and content to the last mentioned directory AOA


        ### Change AOA in includeDict  ###
        f = open("openfoam_runs_pointwise/AOA_" + str(alpha) + '/system/includeDict', 'r+')
        t = f.readlines()
        new_lines = []
        for line in t:
            if 'alpha' in line.split():
                new_lines = np.append(new_lines, 'alpha               ' + str(alpha) + ';\n')
            elif 'Re' in line.split():
                new_lines = np.append(new_lines, 'Re                  ' + str(Re) + ';\n')
            else:
                new_lines = np.append(new_lines, line)

        f.close()
        f = open("openfoam_runs_pointwise/AOA_" + str(alpha) + '/system/includeDict', 'w')
        f.writelines(new_lines)
        f.close()


        ### Renumber Pointwise mesh to Openfoam mesh ###
        case_dir = "openfoam_runs_pointwise/AOA_" + str(alpha)
        env = os.environ.copy()
        env["PWD"] = os.path.abspath(case_dir)
        subprocess.run(["renumberMesh"], cwd=case_dir, env=env, check=True)


        # Initialise next sim with previous converged simulation field
        if initialise_sim and alpha != previous_alpha:
            subprocess.run(f'mapFields ../AOA_{previous_alpha} -sourceTime latestTime', shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))


        ### Run simpleFoam ###
        if parallel_comp:
            subprocess.run(['decomposePar'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))
            subprocess.run(['mpirun -np 4 simpleFoam -parallel > log'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))
            subprocess.run(['reconstructPar'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))
        else:
            subprocess.run(['simpleFoam > log'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))


    ### Unsteady sim ###
    if unsteady:
        # Changing working directory
        subprocess.run(["cp", "-r", "openfoam_runs_pointwise/Base_AOA_unsteady/system", "openfoam_runs_pointwise/AOA_" + str(alpha)])    # copying the folders and content to the last mentioned directory

        # Delete decomposed subdomains
        subprocess.run(["rm -rf processor*"], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))

        # Delete postProcessing
        subprocess.run(["rm -r postProcessing"], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))

        # Check if float
        def is_float(string):
            try:
                float(string)
                return True
            except ValueError:
                return False

        # List all directories that are float values
        time_dirs_all = [d for d in os.listdir("openfoam_runs_pointwise/AOA_" + str(alpha)) if is_float(d)]

        # Convert to floats and find the maximum value
        time_dirs_numbers = [float(d) for d in time_dirs_all if d != "1"]
        max_value = str(max(time_dirs_numbers)) if time_dirs_numbers else None

        # Delete all directories excpet "1" and last iteration
        time_dirs = [d for d in time_dirs_all if d != "1" and d != max_value]
        for time_dir in time_dirs:
            subprocess.run([f'rm -r {time_dir}'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))
            print(f"Deleted directory: {time_dir}")


        # Run the shell command to get the latest time directory
        result = subprocess.run(
            ['ls -d [0-9]* | sort -n | tail -1'],
            shell=True,
            cwd="openfoam_runs_pointwise/AOA_" + str(alpha),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if result.returncode == 0:  # Check if the command was successful
            latest_time = result.stdout.strip()  # Get the latest time directory
            print(f"Latest time directory: {latest_time}")

            # Copy the latest time directory to 0
            subprocess.run([f'cp -r {latest_time} 0'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))
            print(f"Copied {latest_time} to 0 directory.")

            # List all numerical directories and delete them except "0"
            time_dirs = [d for d in os.listdir("openfoam_runs_pointwise/AOA_" + str(alpha)) if d.isdigit() and d != '0']
            for time_dir in time_dirs:
                subprocess.run([f'rm -r {time_dir}'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))
                print(f"Deleted directory: {time_dir}")
        else:
            print(f"Error finding latest time directory: {result.stderr}")



        ### Change AOA in includeDict  ###
        f = open("openfoam_runs_pointwise/AOA_" + str(alpha) + '/system/includeDict', 'r+')
        t = f.readlines()
        new_lines = []
        for line in t:
            if 'alpha' in line.split():
                new_lines = np.append(new_lines, 'alpha               ' + str(alpha) + ';\n')
            elif 'Re' in line.split():
                new_lines = np.append(new_lines, 'Re                  ' + str(Re) + ';\n')
            else:
                new_lines = np.append(new_lines, line)

        f.close()
        f = open("openfoam_runs_pointwise/AOA_" + str(alpha) + '/system/includeDict', 'w')
        f.writelines(new_lines)
        f.close()


        ### Run pimpleFoam ###
        if parallel_comp:
            subprocess.run(['decomposePar'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))
            subprocess.run(['mpirun -np 4 pimpleFoam -parallel > log'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))
            subprocess.run(['reconstructPar'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))
        else:
            subprocess.run(['pimpleFoam > log'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))



    ### Output surfaces postProcess ###
    subprocess.run(['postProcess -func surfaces'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))

    ### Output yPlus postProcess ###
    subprocess.run(['simpleFoam -postProcess -func yPlus'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))

    ### Output VTK paraview field ###
    if VTK_paraview:
        subprocess.run(['foamToVTK -time -latestTime'], shell=True, cwd="openfoam_runs_pointwise/AOA_" + str(alpha))


    ### Residual convergence plot ###
    if residuals_plot:
        residual_file_path = "openfoam_runs_pointwise/AOA_" + str(alpha) + "/postProcessing/solverInfo/1/solverInfo.dat"

        # Check if the file exists
        if not os.path.isfile(residual_file_path):
            raise FileNotFoundError(f"The file '{residual_file_path}' does not exist. Please check the file path.")

        # Initialize variables
        iteration = []
        Ux_final = []
        Uy_final = []
        k_final = []
        omega_final = []
        p_final = []

        with open(residual_file_path, "r") as file:
            lines = file.readlines()[2:]

            for line in lines:
                line = line.strip()
                # Split the line into parts
                parts = line.split()

                # Extract relevant columns
                iteration.append(float(parts[0]))
                Ux_final.append(float(parts[3]))
                Uy_final.append(float(parts[6]))
                k_final.append(float(parts[11]))
                p_final.append(float(parts[16]))
                omega_final.append(float(parts[21]))

        min_y = min(
            min(Ux_final), min(Uy_final), min(k_final),
            min(omega_final), min(p_final))


        # Plot residual values
        plt.figure(figsize=(12, 7))
        plt.plot(iteration, Ux_final, label='Ux_final Residuals', marker='o', markersize=1)
        plt.plot(iteration, Uy_final, label='Uy_final Residuals', marker='x', markersize=1)
        plt.plot(iteration, k_final, label='k_final Residuals', marker='s', markersize=1)
        plt.plot(iteration, omega_final, label='omega_final Residuals', marker='d', markersize=1)
        plt.plot(iteration, p_final, label='p_final Residuals', marker='^', markersize=1)
        plt.plot([0, iteration[-1]], [10**-6, 10**-6], "--", color='k')

        plt.xlabel('Iteration Number')
        plt.ylabel('Residual Value')
        plt.title(f'AOA: {alpha}, profile: {profile_name}', fontsize=8)
        plt.yscale('log')
        plt.ylim(bottom=min_y, top=0.1)
        plt.yticks([10**i for i in range(-1, -8, -1)])
        plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.savefig("openfoam_runs_pointwise/AOA_" + str(alpha) + "/postProcessing/solverInfo/1/residuals" + str(alpha) + ".png", format='png', dpi=300)  # Saving figure to results
        plt.close()


    ### Force coefficient residual plot ###
    if force_plot:
        force_file_path = "openfoam_runs_pointwise/AOA_" + str(alpha) + "/postProcessing/forceCoeffs/1/coefficient.dat"

        # Check if the file exists
        if not os.path.isfile(force_file_path):
            raise FileNotFoundError(f"The file '{force_file_path}' does not exist. Please check the file path.")

        # Initialize variables
        iteration = []
        residual_Cd = []
        residual_Cl = []
        residual_CmPitch = []

        with open(force_file_path, "r") as file:
            lines = file.readlines()[13:]

            # Initial values for residuals
            previous_Cd = 0.0
            previous_Cl = 0.0
            previous_CmPitch = 0.0

            for line in lines:
                # Split the line into parts
                parts = line.split()

                # Extract columns
                iteration.append(float(parts[0]))
                Cd_value = float(parts[1])
                Cl_value = float(parts[3])
                CmPitch_value = float(parts[5])


                # Calculate residuals
                if len(iteration) > 1:
                    residual_Cd.append(abs(Cd_value - previous_Cd))  # Absolute difference for Cd residual
                    residual_Cl.append(abs(Cl_value - previous_Cl))  # Absolute difference for Cl residual
                    residual_CmPitch.append(abs(CmPitch_value - previous_CmPitch))  # Absolute difference for CmPitch residual
                else:
                    residual_Cd.append(0.0)  # First iteration, no residual
                    residual_Cl.append(0.0)  # First iteration, no residual
                    residual_CmPitch.append(0.0)  # First iteration, no residual

                # Update previous values for the next iteration
                previous_Cd = Cd_value
                previous_Cl = Cl_value
                previous_CmPitch = CmPitch_value

        # Plot force residual plot
        plt.figure(figsize=(12, 7))
        plt.plot(iteration, residual_Cd, label='Residual Cd', marker='o', markersize=1)
        plt.plot(iteration, residual_Cl, label='Residual Cl', marker='x', markersize=1)
        plt.plot(iteration, residual_CmPitch, label='Residual Cm', marker='s', markersize=1)

        plt.xlabel('Iteration Number')
        plt.ylabel('Residual')
        plt.title(f'Res Cl, Cd and Cm at AOA: {alpha}deg, profile: {profile_name}', fontsize=8)
        plt.yscale('log')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("openfoam_runs_pointwise/AOA_" + str(alpha) + "/postProcessing/forceCoeffs/1/force_residuals" + str(alpha) + ".png", format='png', dpi=300)
        plt.close()


    ### Pressure coefficient plot ###
    if cp_plot:
        # File path
        pressure_file_path = "openfoam_runs_pointwise/AOA_" + str(alpha) + "/postProcessing/surfaces/" + str(int(iteration[-1])) + "/p_airfoilSurface.raw"
        cp_output_file_path = "openfoam_runs_pointwise/AOA_" + str(alpha) + "/postProcessing/surfaces/" + str(int(iteration[-1])) + "/cp_AOA_" + str(alpha) + ".dat"

        data = []
        with open(pressure_file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                # Split each line into its components (x, y, z, Cp)
                parts = line.split()

                # Check if the line has at least 4 columns (x, y, z, Cp)
                if len(parts) >= 4:
                    x = float(parts[0])
                    y = float(parts[1])
                    z = float(parts[2])
                    Cp = float(parts[3])

                    # Only include the data where z is not equal to 1
                    if z != 1:
                        data.append((x, y, Cp))

        data = np.array(data)


        # Apply RDP to x and Cp only
        reduced_points = rdp(data[:, [0, 2]], epsilon=0.006)  # Get reduced indices, modify eps to desired data resolution, smaller equals more resolution

        # Add Y value to reduced data again to plot the according airfoil shape with it
        reduced_data = []
        for point in reduced_points:
            # Find indices where x and Cp match the reduced point
            indices = np.where((np.isclose(data[:, 0], point[0])) & (np.isclose(data[:, 2], point[1])))[0]

            # Check if at least one match is found
            if len(indices) > 0:
                # Use the first match (in case there are duplicates)
                idx = indices[0]
                reduced_data.append(data[idx])

        reduced_data = np.array(reduced_data)

        # Write reduced data to dat file
        with open(cp_output_file_path, 'w') as f:
            f.write('# x y Cp\n')  # Header
            for x, y, Cp in reduced_data:
                f.write(f"{x} {y} {Cp}\n")

        x_vals = data[:, 0]
        y_vals = data[:, 1]
        Cp_vals = data[:, 2]
        x_vals_reduced = reduced_data[:, 0]
        Cp_vals_reduced = reduced_data[:, 2]

        # Plot Cp plot
        plt.figure(figsize=(10, 6))
        plt.plot(x_vals, Cp_vals, marker='o', markersize=0.5, linestyle='-', linewidth=0.5)
        plt.plot(x_vals_reduced, Cp_vals_reduced, marker='o', markersize=0.5, linestyle='-', linewidth=0.5)
        plt.plot(x_vals, -y_vals, marker='o', markersize=0.5, linestyle='-', color='k', linewidth=0.5)
        plt.xlabel('x/c [-]')
        plt.ylabel('Cp [-]')
        plt.yticks(np.arange(np.round(min(Cp_vals) * 1.2, decimals=1), np.round(max(Cp_vals) * 1.2, decimals=1), 0.2))
        plt.gca().invert_yaxis()
        plt.title(f'cp plot, AOA: {alpha}deg, profile: {profile_name}', fontsize=8)
        plt.grid(True)
        plt.savefig("openfoam_runs_pointwise/AOA_" + str(alpha) + "/postProcessing/surfaces/" + str(int(iteration[-1])) + "/cp_AOA_" + str(alpha) + ".png", format='png', dpi=500)  # Saving figure to results
        plt.close()


    ### Write last line of forceCoeffs.dat and residuals (solverInfo.dat) into polar_pointwise.dat ###
    force_coeffs_path = "openfoam_runs_pointwise/AOA_" + str(alpha) + "/postProcessing/forceCoeffs/1/coefficient.dat"
    residuals_path = "openfoam_runs_pointwise/AOA_" + str(alpha) + "/postProcessing/solverInfo/1/solverInfo.dat"
    yplus_path = "openfoam_runs_pointwise/AOA_" + str(alpha) + "/postProcessing/yPlus/" + str(int(iteration[-1])) + "/yPlus.dat"
    output_file = root_path + "/results/polar_pointwise.dat"

    # Define the header row
    header = "# Alpha Cd Cs Cl CmRoll CmPitch CmYaw Cd(f) Cd(r) Cs(f) Cs(r) Cl(f) Cl(r) Ux_final Uy_final k_final p_final omega_final yplus"

    # Check if the output file already exists and has content
    try:
        with open(output_file, 'r') as datfile:
            file_content = datfile.read()
            file_has_header = header in file_content
    except FileNotFoundError:
        file_has_header = False  # File does not exist, so no header yet

    # If the file is new or doesn't have a header, write the header
    if not file_has_header:
        with open(output_file, 'w') as datfile:
            datfile.write(header + "\n")

    # Read the last line of the force coefficients file
    with open(force_coeffs_path, 'r') as f:
        force_lines = f.readlines()
        last_force_line = force_lines[-1].strip().split()

        # Replace the first value with alpha
        last_force_line[0] = str(alpha)

    # Read the last line of the residual/solverInfo file
    try:
        with open(residuals_path, 'r') as f:
            residuals_lines = f.readlines()
            last_residuals_line = residuals_lines[-1].strip().split()

            # Extract residual values
            Ux_final = last_residuals_line[3]
            Uy_final = last_residuals_line[6]
            k_final = last_residuals_line[11]
            p_final = last_residuals_line[16]
            omega_final = last_residuals_line[21]


            # Combine the residuals into a string
            residuals_data = f"{Ux_final} {Uy_final} {k_final} {p_final} {omega_final}"

    except FileNotFoundError:
        residuals_data = "N/A N/A N/A N/A N/A"  # Placeholder if solverInfo.dat is missing

    # Read the last line of the yplus file and extract the max yplus value
    try:
        with open(yplus_path, 'r') as f:
            yplus_lines = f.readlines()
            last_yplus_line = yplus_lines[-1].strip().split()
            yplus_max = last_yplus_line[3]
            yplus_data = f"{yplus_max}"

    except FileNotFoundError:
        yplus_data = "N/A"  # Placeholder if yPlus.dat is missing

    # Combine forces, residuals and yplus max
    combined_line = " ".join(last_force_line) + " " + residuals_data + " " + yplus_data

    # Append the combined data to the polar output file
    with open(output_file, 'a') as datfile:
        datfile.write(combined_line + "\n")
    return ()