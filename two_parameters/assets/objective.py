from assets.error import *
from assets.model import runModel
from assets.data import extractData

def objectiveFunction(parameter, realData, umatFile, inputFile):
    """
    Objective function for PSO optimization
    Returns residual and errors
    """
    try:
        modify(parameter)
        odbFile = runModel(umatFile, inputFile)
        data = extractData(odbFile)
        residual, rf_error, u3_error = calculateResiduals(data, realData)
        print(f"Parameter: ({parameter[0]:.4f}, {parameter[1]:.4f}), Peak: {data[0]:.4f}, Pile Up: {data[1]:.4f}, Residual: {residual:.4f}")
        print(f"    rf_error: {abs(rf_error):.4e}, u3_error: {abs(u3_error):.4e}")
        return residual, abs(rf_error), abs(u3_error)
    except Exception as e:
        print(f"Error in objective function with parameter {parameter}: {e}")
        return float('inf'), float('inf'), float('inf')
