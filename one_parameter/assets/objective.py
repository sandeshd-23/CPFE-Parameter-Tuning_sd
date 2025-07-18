from assets.error import *
from assets.model import runModel
from assets.data import extractData

def objectiveFunction(parameter, realPeak, umatFile, inputFile):
    """
    Objective function for PSO optimization
    Returns residual between simulated and true data
    """
    try:
        modify(parameter)
        odbFile = runModel(umatFile, inputFile)
        peak = extractData(odbFile)
        residual = calculateResiduals(peak, realPeak)
        print(f"Parameter: {parameter:.4f}, Peak: {peak:.4f}, Residual: {residual:.4f}")
        return residual
    except Exception as e:
        print(f"Error in objective function with parameter {parameter}: {e}")
        return float('inf')  # Return large value for failed evaluations