from odbAccess import openOdb
from abaqusConstants import *


def extractData(odbFile):
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