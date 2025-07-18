def calculateResiduals(data, trueData):
    """Calculate squared residual between simulated and experimental data"""
    return (data - trueData)**2

def modify(parameter):
    try:
        # Read entire file as a single string
        with open(r'umat_directory\kMaterialParam.f', 'r') as f:
            s = f.read()

        # Split the file by newline character
        lines = s.split('\n')

        # Replace the 152th line (index 151)
        if len(lines) < 152:
            raise ValueError("File has fewer than 99 lines.")
        lines[151] = f"		XTAUC1 = {parameter}"

        # Join lines back into a single string
        modified_s = '\n'.join(lines)

        # Write back to file
        with open(r'umat_directory\kMaterialParam.f', 'w') as f:
            f.write(modified_s)

    except FileNotFoundError:
        print("Error: kMaterialParam.f file not found")
        raise
    except Exception as e:
        print(f"Error modifying file: {e}")
        raise