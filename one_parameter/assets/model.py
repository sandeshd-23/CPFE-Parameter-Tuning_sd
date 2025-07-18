import subprocess
import time
import os
from assets.file_completion import *
from assets.abaqus_completion import *

def runModel(umatFile, inputFile) -> str:
    """Run Abaqus simulation and return ODB file path"""
    odb_file = f"{inputFile}.odb"
    
    # Remove existing ODB file to avoid confusion
    if os.path.exists(odb_file):
        try:
            os.remove(odb_file)
            print(f"Removed existing {odb_file}")
        except:
            print(f"Warning: Could not remove existing {odb_file}")
    
    try:
        print(f"Starting Abaqus simulation: {inputFile}")
        
        # Start Abaqus job
        process = subprocess.Popen(
            f"abaqus job={inputFile} user={umatFile} cpu = 10",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        # Wait for the process to complete
        print("Waiting for Abaqus process to finish...")
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print(f"Abaqus process returned error code {process.returncode}")
            print(f"Error output: {stderr.decode()}")
            raise subprocess.CalledProcessError(process.returncode, f"abaqus job={inputFile} user={umatFile}")
        
        # Additional wait for job completion using status files
        print("Abaqus process finished. Checking job completion...")
        if not wait_for_abaqus_completion(inputFile, timeout=8000):  # 100 minutes timeout
            raise RuntimeError("Abaqus job did not complete successfully")
        
        # Wait for ODB file to be completely written
        print("Waiting for ODB file to be completely written...")
        if not wait_for_file_completion(odb_file, timeout=8000):  # 100 minutes timeout for file completion
            raise RuntimeError(f"ODB file {odb_file} was not completed within timeout")
        
        # Additional delay to ensure file is fully closed after final write
        print("Waiting additional 30 seconds to ensure file is fully closed...")
        time.sleep(30)
        
        print(f"ODB file ready: {odb_file}")
        return odb_file
        
    except subprocess.CalledProcessError as e:
        print(f"Error running Abaqus: {e}")
        raise
    except Exception as e:
        print(f"Unexpected error: {e}")
        raise
    