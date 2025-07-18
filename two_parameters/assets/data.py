from odbAccess import openOdb
from abaqusConstants import *
import os
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import csv

def getPeak(nodeData):
    """ Can define custom peak extraction here """
    return nodeData[-1]*((4*322.58)/410)*100
    # returns the peak from the csv

def extractDatathree(odbFile):
    """Extract displacement data from ODB file and return peak value"""
    step_name = 'Step-1'
    frame_number = -1  # Last frame

    try:
        # === Open ODB and access data ===
        odb = openOdb(path=odbFile)
        frame = odb.steps[step_name].frames[frame_number]

        # --- Displacement Field Output ---
        disp_field = frame.fieldOutputs['RF']
        disp_subset = disp_field.getSubset(position=NODAL)

        # === Collect node data ===
        node_data = []

        for v in disp_subset.values:
            #node_label = v.nodeLabel
            displacement = v.data[-1]  # Tuple: (RF1, RF2, RF3)
            #node = v.instance.getNodeFromLabel(node_label)
            #coordinates = node.coordinates  # Tuple: (X, Y, Z)
            
            #row = [node_label] + list(coordinates) + list(displacement)
            node_data.append(displacement)

        # Sort by node label
        #node_data.sort(key=lambda x: x[-1])
        node_data.sort()

        val = getPeak(node_data)
        odb.close()
        
        return val
        
    except Exception as e:
        print(f"Error extracting data from {odbFile}: {e}")
        if 'odb' in locals():
            odb.close()
        raise

def extractDatatwo(odbFile):
    """Extract RF3 vs Time for Node 84, scale RF3, and return max scaled RF3"""
    step_name = 'Step-1'
    node_label = 84
    instance_name = 'SPHEREINDENT-1'
    scale_factor = (4 * 322.58 / 410) * 0.001

    try:
        odb = openOdb(path=odbFile)
        step = odb.steps[step_name]
        instance = odb.rootAssembly.instances[instance_name]

        scaled_rf3_values = []

        for frame in step.frames:
            time = frame.frameValue
            rf_field = frame.fieldOutputs['RF']
            rf_values = rf_field.getSubset(position=NODAL).values
            
            for val in rf_values:
                if val.nodeLabel == node_label and val.instance.name == instance_name:
                    rf3 = val.data[2]
                    scaled_rf3 = rf3 * scale_factor
                    scaled_rf3_values.append(scaled_rf3)
                    break

        odb.close()
        return max(scaled_rf3_values) if scaled_rf3_values else None

    except Exception as e:
        print(f"Error extracting data from {odbFile}: {e}")
        if 'odb' in locals():
            odb.close()
        raise

def extractDataone(odbFile):
    """Extract RF3 vs Time for Node 84 and U3 peak from mirrored plot"""
    step_name = 'Step-1'
    node_label = 84
    instance_name = 'SPHEREINDENT-1'
    scale_factor = (4 * 322.58 / 410) * 0.001

    try:
        odb = openOdb(path=odbFile)
        step = odb.steps[step_name]
        instance = odb.rootAssembly.instances[instance_name]

        # === Part 1: Extract and scale RF3 values ===
        scaled_rf3_values = []
        for frame in step.frames:
            rf_field = frame.fieldOutputs['RF']
            rf_values = rf_field.getSubset(position=NODAL).values
            for val in rf_values:
                if val.nodeLabel == node_label and val.instance.name == instance_name:
                    rf3 = val.data[2]
                    scaled_rf3 = rf3 * scale_factor
                    scaled_rf3_values.append(scaled_rf3)
                    break  # Only one match expected per frame

        # === Part 2: Extract U3 field from last frame and compute mirrored plot peak ===
        frame = step.frames[-1]
        disp_field = frame.fieldOutputs['U']
        disp_subset = disp_field.getSubset(position=NODAL).values

        u3_data = []
        for v in disp_subset:
            if v.instance.name == instance_name:
                node = v.instance.getNodeFromLabel(v.nodeLabel)
                x, y, z = node.coordinates
                u1, u2, u3 = v.data
                if round(z, 4) == 20:  # Top surface only
                    u3_data.append([x, y, u3])

        if not u3_data:
            raise ValueError("No U3 data found on top surface (Z=20)")

        df = pd.DataFrame(u3_data, columns=['X', 'Y', 'U3'])

        # Interpolate to regular grid
        x = df['X'].values
        y = df['Y'].values
        u3 = df['U3'].values
        xi = np.linspace(min(x), 10, 100)
        yi = np.linspace(10, 20, 100)
        xi_grid, yi_grid = np.meshgrid(xi, yi)
        ui_grid = griddata((x, y), u3, (xi_grid, yi_grid), method='linear')

        # Mirror the grid in 4 quadrants
        original = ui_grid
        mirror_y = np.fliplr(original)
        mirror_x = np.flipud(original)
        mirror_both = np.flipud(np.fliplr(original))
        top_half = np.concatenate([mirror_y, original], axis=1)
        bottom_half = np.concatenate([mirror_both, mirror_x], axis=0)
        full_plot = np.concatenate([top_half, bottom_half], axis=0)

        # Get U3 values along center vertical line (X=0)
        x_index_center = full_plot.shape[1] // 2
        u3_along_y = -full_plot[:, x_index_center]

        # Peak value in U3 along Y at X=0
        peak_u3 = np.nanmax(u3_along_y)

        odb.close()
        return (max(scaled_rf3_values) if scaled_rf3_values else None, peak_u3)

    except Exception as e:
        print(f"Error extracting data from {odbFile}: {e}")
        if 'odb' in locals():
            odb.close()
        raise
        
def extractData(odbFile):
    """
    Extract RF3 vs Time for Node 84 and U3 peak from mirrored plot and 
    displacement data from ODB file and return peak and pile up
    """
    output_csv = "displacement_data.csv"
    step_name = 'Step-1'
    frame_number = -1  # Last frame

    try:
        # === Open ODB and access data ===
        odb = openOdb(path=odbFile)
        frame = odb.steps[step_name].frames[frame_number]

        # --- Displacement Field Output ---
        disp_field = frame.fieldOutputs['U']
        disp_subset = disp_field.getSubset(position=NODAL)

        # === Collect node data ===
        node_data = []

        for v in disp_subset.values:
            node_label = v.nodeLabel
            displacement = v.data  # Tuple: (U1, U2, U3)
            node = v.instance.getNodeFromLabel(node_label)
            coordinates = node.coordinates  # Tuple: (X, Y, Z)
            
            row = [node_label] + list(coordinates) + list(displacement)
            node_data.append(row)

        node_data.sort(key=lambda x: x[0])

        # === Write to CSV ===
        with open(output_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Node Label', 'X', 'Y', 'Z', 'U1', 'U2', 'U3'])
            writer.writerows(node_data)

        # Load CSV data
        df = pd.read_csv("displacement_data.csv")

        # Filter top surface (Z = 20)
        df_top = df[df['Z'] == 20]

        # Extract coordinates and U3 values
        x = df_top['X'].values
        y = df_top['Y'].values
        u3 = df_top['U3'].values

        # Create regular grid for interpolation
        xi = np.linspace(min(x), 10, 100)
        yi = np.linspace(10, 20, 100)
        xi_grid, yi_grid = np.meshgrid(xi, yi)
        ui_grid = griddata((x, y), u3, (xi_grid, yi_grid), method='nearest')

        # Create the 4 mirrored versions
        # Original (top-right quadrant)
        original = ui_grid

        # Mirror about y=10 (top-left quadrant)
        mirror_y = np.fliplr(original)

        # Mirror about x=0 (bottom-right quadrant)
        mirror_x = np.flipud(original)

        # Mirror about both x=0 and y=10 (bottom-left quadrant)
        mirror_both = np.flipud(np.fliplr(original))

        # Combine all 4 quadrants
        # Top half: mirror_y (left) + original (right)
        top_half = np.concatenate([mirror_y, original], axis=1)

        # Bottom half: mirror_both (left) + mirror_x (right)
        bottom_half = np.concatenate([mirror_both, mirror_x], axis=1)

        # Full combined plot: bottom_half (bottom) + top_half (top)
        full_plot = np.concatenate([top_half, bottom_half], axis=0)

        x_index_center = full_plot.shape[1] // 2
        u3_along_y = -full_plot[:, x_index_center]

        # === RF3 extraction ===
        node_label = 84
        instance_name = 'SPHEREINDENT-1'
        scale_factor = (4 * 322.58 / 410) * 0.001

        step = odb.steps[step_name]
        instance = odb.rootAssembly.instances[instance_name]

        scaled_rf3_values = []

        for frame in step.frames:
            rf_field = frame.fieldOutputs['RF']
            rf_values = rf_field.getSubset(position=NODAL).values

            for val in rf_values:
                if val.nodeLabel == node_label and val.instance.name == instance_name:
                    rf3 = val.data[2]
                    scaled_rf3 = rf3 * scale_factor
                    scaled_rf3_values.append(scaled_rf3)
                    break

        odb.close()

        # === Return both peak U3 and max scaled RF3 ===
        if len(u3_along_y) == 0:
            peak_u3 = 0.0
        else:
            # Remove NaN values before finding max
            u3_clean = u3_along_y[~np.isnan(u3_along_y)]
            peak_u3 = np.max(u3_clean) if len(u3_clean) > 0 else 0.0
        
        peak_rf3 = max(scaled_rf3_values) if scaled_rf3_values else 0.0

        # Clean up temporary file
        if os.path.exists(output_csv):
            os.remove(output_csv)

        return peak_rf3, peak_u3*1000

    except Exception as e:
        print(f"Error extracting data from {odbFile}: {e}")
        if 'odb' in locals():
            odb.close()
        # Clean up temporary file on error
        if os.path.exists(output_csv):
            os.remove(output_csv)
        raise