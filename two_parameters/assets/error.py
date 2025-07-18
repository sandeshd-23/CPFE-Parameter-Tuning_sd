import os

def calculateResiduals(data, trueData):
    """Calculate squared residual between simulated and experimental data"""
    w_rf = 1.0
    w_u3 = 1.0
    rf, u3 = data
    rf_true, u3_true = trueData
    rf_error = (rf - rf_true) / rf_true
    u3_error = (u3 - u3_true) / u3_true

    residual = w_rf * rf_error**2 + w_u3 * u3_error**2

    return residual, rf_error, u3_error

def modify(parameter):
    """Modify the kCRSS.f and kMaterialParam.f files to update in real time"""
    peak_parameter, pileUp_parameter = parameter
    
    # Modify kMaterialParam.f
    try:
        # Check if file exists before trying to read
        if not os.path.exists(r'umat_directory\kMaterialParam.f'):
            raise FileNotFoundError("kMaterialParam.f not found")
            
        # Read entire file as a single string
        with open(r'umat_directory\kMaterialParam.f', 'r') as f:
            content = f.read()

        # Split the file by newline character
        lines = content.split('\n')

        # Replace the 152nd line (index 151)
        if len(lines) < 152:
            raise ValueError(f"kMaterialParam.f has only {len(lines)} lines, need at least 152")
        
        lines[151] = f"		XTAUC1 = {peak_parameter}"

        # Join lines back into a single string
        modified_content = '\n'.join(lines)

        # Write back to file
        with open(r'umat_directory\kMaterialParam.f', 'w') as f:
            f.write(modified_content)

    except Exception as e:
        print(f"Error modifying kMaterialParam.f: {e}")
        raise
    
    # Modify kCRSS.f
    try:
        # Check if file exists before trying to read
        if not os.path.exists(r'umat_directory\kCRSS.f'):
            raise FileNotFoundError("kCRSS.f not found")
            
        # Read entire file as a single string
        with open(r'umat_directory\kCRSS.f', 'r') as f:
            content = f.read()

        # Split the file by newline character
        lines = content.split('\n')

        # Replace the 84th line (index 83)
        if len(lines) < 84:
            raise ValueError(f"kCRSS.f has only {len(lines)} lines, need at least 84")
        
        lines[83] = f"          tauc = tauc + {pileUp_parameter}*G12*(burgerv(1))*sqrt(gndtot)"

        # Join lines back into a single string
        modified_content = '\n'.join(lines)

        # Write back to file
        with open(r'umat_directory\kCRSS.f', 'w') as f:
            f.write(modified_content)

    except Exception as e:
        print(f"Error modifying kCRSS.f: {e}")
        raise