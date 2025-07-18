import os
import time

def wait_for_file_completion(filepath, timeout=8000, check_interval=30):
    """
    Wait for a file to be completely written by monitoring its size stability
    For files that update periodically (every 1.5-2 minutes), we need longer stability checks
    
    Args:
        filepath: Path to the file to monitor
        timeout: Maximum time to wait in seconds (default 100 minutes)
        check_interval: Time between size checks in seconds (30s for periodic updates)
    
    Returns:
        True if file is stable, False if timeout reached
    """
    if not os.path.exists(filepath):
        print(f"Waiting for {filepath} to be created...")
        start_time = time.time()
        while not os.path.exists(filepath):
            if time.time() - start_time > timeout:
                print(f"Timeout: {filepath} was not created within {timeout} seconds")
                return False
            time.sleep(10)  # Check every 10 seconds for file creation
    
    print(f"File {filepath} detected. Monitoring for completion...")
    print("Note: File updates every 1.5-2 minutes during simulation")
    
    # Monitor file size stability - need longer period due to periodic updates
    previous_size = -1
    stable_count = 0
    start_time = time.time()
    last_update_time = time.time()
    
    while time.time() - start_time < timeout:
        try:
            current_size = os.path.getsize(filepath)
            current_time = time.time()
            elapsed_minutes = (current_time - start_time) / 60
            
            if current_size != previous_size:
                # File size changed - simulation still running
                stable_count = 0
                last_update_time = current_time
                print(f"[{elapsed_minutes:.1f} min] File updated - Size: {current_size:,} bytes")
            else:
                # File size unchanged
                stable_count += 1
                time_since_update = (current_time - last_update_time) / 60
                
                # Since file updates every 1.5-2 minutes, wait at least 4 minutes of no changes
                # before considering it complete (2x the update interval for safety)
                if stable_count >= 8 and time_since_update >= 4.0:  # 8 checks * 30s = 4 minutes
                    print(f"[{elapsed_minutes:.1f} min] File appears complete - No updates for {time_since_update:.1f} minutes")
                    print(f"Final size: {current_size:,} bytes")
                    return True
                elif stable_count % 4 == 0:  # Print status every 2 minutes
                    print(f"[{elapsed_minutes:.1f} min] No update for {time_since_update:.1f} minutes - still monitoring...")
            
            previous_size = current_size
            time.sleep(check_interval)
            
        except OSError:
            # File might be locked or still being written
            stable_count = 0
            time.sleep(check_interval)
    
    print(f"Timeout: File {filepath} monitoring exceeded {timeout/60:.1f} minutes")
    return False