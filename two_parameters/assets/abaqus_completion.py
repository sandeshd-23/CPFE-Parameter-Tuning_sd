import os
import time

def wait_for_abaqus_completion(job_name, timeout=8000):
    """
    Wait for Abaqus job to complete by monitoring status files
    
    Args:
        job_name: Name of the Abaqus job
        timeout: Maximum time to wait in seconds (default 100 minutes)
    
    Returns:
        True if job completed successfully, False otherwise
    """
    status_file = f"{job_name}.sta"
    msg_file = f"{job_name}.msg"
    
    start_time = time.time()
    print(f"Monitoring Abaqus job completion (timeout: {timeout/60:.0f} minutes)...")
    
    while time.time() - start_time < timeout:
        elapsed_minutes = (time.time() - start_time) / 60
        
        # Check if status file exists and contains completion message
        if os.path.exists(status_file):
            try:
                with open(status_file, 'r') as f:
                    content = f.read()
                    if "COMPLETED" in content.upper():
                        print(f"[{elapsed_minutes:.1f} min] Abaqus job completed successfully")
                        return True
                    elif "ABORTED" in content.upper() or "ERROR" in content.upper():
                        print(f"[{elapsed_minutes:.1f} min] Abaqus job failed - check .msg file for details")
                        return False
            except:
                pass
        
        # Also check message file for completion
        if os.path.exists(msg_file):
            try:
                with open(msg_file, 'r') as f:
                    content = f.read()
                    if "Analysis completed" in content or "Job completed" in content:
                        print(f"[{elapsed_minutes:.1f} min] Abaqus analysis completed")
                        return True
            except:
                pass
        
        # Print status every 5 minutes
        if int(elapsed_minutes) % 5 == 0 and int(elapsed_minutes) > 0 and (time.time() - start_time) % 300 < 30:
            print(f"[{elapsed_minutes:.0f} min] Still running - checking for completion...")
        
        time.sleep(30)  # Check every 30 seconds
    
    print(f"Timeout: Abaqus job did not complete within {timeout/60:.0f} minutes")
    return False