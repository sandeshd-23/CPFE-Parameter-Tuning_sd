import random
import sys
import numpy as np
from assets.objective import objectiveFunction

# Redirect all output to a file
class Tee:
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()

# Open log file and redirect stdout/stderr
log_file = open('output.txt', 'w')
sys.stdout = Tee(sys.stdout, log_file)
sys.stderr = Tee(sys.stderr, log_file)


#-------------------------  Modify Input Data Here  -----------------------------------

# Global input file and real experimental data
inputFile = r"unimp_whole_amp4_ID100" # .inp file
realPeak = 5  # TODO: Set this to actual peak from experimental data
realPileUp = 0.25  # TODO: Set this to actual peak from experimental data
realData = realPeak, realPileUp
umatFile = r'umat_directory\umat.for'
# .for file template

#-------------------------  Modify initial guess Here  --------------------------------

# Initial parameter guess (if no guess then comment out line 35) and change line 376
initial_peak_parameter = 360
initial_pileUp_parameter = 0.01
initial_guess = initial_peak_parameter, initial_pileUp_parameter

#-------------------------  Modify PSO parameters Here  -------------------------------

# PSO hyperparameters
n_particles = 5 # int(input("Enter number of particles to test"))
iterations = 4 # int(input("Enter number of iterations"))
# Time complexity = O(n_particles*iterations + n_particles)*O(runModel)
# print("Enter w, c1, c2:")
w = 0.5
c1 = 1.5
c2 = 1.5

#--------------------------------------------------------------------------------------

# --- Particle Swarm Optimization core ---
def apply_bounds_with_reflection(position, bounds):
    """Apply bounds with boundary reflection to prevent particle sticking"""
    PEAK_MIN, PEAK_MAX = bounds[0]
    PILEUP_MIN, PILEUP_MAX = bounds[1]
    
    peak, pileup = position
    
    # Reflect off boundaries instead of simple clamping
    if peak < PEAK_MIN:
        peak = PEAK_MIN + abs(peak - PEAK_MIN)
        if peak > PEAK_MAX:  # Handle case where reflection goes beyond upper bound
            peak = PEAK_MAX
    elif peak > PEAK_MAX:
        peak = PEAK_MAX - abs(peak - PEAK_MAX)
        if peak < PEAK_MIN:  # Handle case where reflection goes beyond lower bound
            peak = PEAK_MIN
            
    if pileup < PILEUP_MIN:
        pileup = PILEUP_MIN + abs(pileup - PILEUP_MIN)
        if pileup > PILEUP_MAX:
            pileup = PILEUP_MAX
    elif pileup > PILEUP_MAX:
        pileup = PILEUP_MAX - abs(pileup - PILEUP_MAX)
        if pileup < PILEUP_MIN:
            pileup = PILEUP_MIN
            
    return (peak, pileup)

def pso(n_particles, iterations, w, c1, c2, initial_guess=None):
    """
    Particle Swarm Optimization algorithm for 2D parameters (tuples).
    Improved with velocity clamping and boundary reflection.
    """
    print("Starting PSO optimization...")

    # Define consistent bounds
    PEAK_MIN, PEAK_MAX = 50, 500
    PILEUP_MIN, PILEUP_MAX = 0.001, 10
    bounds = [(PEAK_MIN, PEAK_MAX), (PILEUP_MIN, PILEUP_MAX)]
    
    # Define maximum velocity (10% of parameter range)
    max_velocity = [(PEAK_MAX - PEAK_MIN) * 0.1, (PILEUP_MAX - PILEUP_MIN) * 0.1]
    
    # Initialize particles and velocities
    particles = []
    velocities = []

    if initial_guess is not None:
        print(f"Using initial guess: {initial_guess}")
        peak_init, pileup_init = initial_guess
        for i in range(n_particles):
            if i == 0:
                # Ensure initial guess is within bounds
                peak_clipped = np.clip(peak_init, PEAK_MIN, PEAK_MAX)
                pileup_clipped = np.clip(pileup_init, PILEUP_MIN, PILEUP_MAX)
                particles.append((peak_clipped, pileup_clipped))
            else:
                # Add random variation around initial guess
                peak_var = random.uniform(-0.5 * peak_init, 1.5 * peak_init)
                pileup_var = random.uniform(-0.5 * pileup_init, 10 * pileup_init)
                peak_new = np.clip(peak_init + peak_var, PEAK_MIN, PEAK_MAX)
                pileup_new = np.clip(pileup_init + pileup_var, PILEUP_MIN, PILEUP_MAX)
                particles.append((peak_new, pileup_new))
    else:
        for _ in range(n_particles):
            particles.append((
                random.uniform(PEAK_MIN, PEAK_MAX),     # Peak parameter range
                random.uniform(PILEUP_MIN, PILEUP_MAX)   # Pile-up parameter range
            ))

    # Initialize velocities with small random values
    velocities = [(random.uniform(-max_velocity[0]*0.1, max_velocity[0]*0.1),
                   random.uniform(-max_velocity[1]*0.1, max_velocity[1]*0.1)) 
                  for _ in range(n_particles)]
    
    pbest = particles[:]
    pbest_scores = []

    print("Evaluating initial particles...")
    for i, p in enumerate(particles):
        try:
            score, rf_err, u3_err = objectiveFunction(p, realData, umatFile, inputFile)
        except Exception as e:
            print(f"Error in objective function for initial particle {i}: {e}")
            score, rf_err, u3_err = float('inf'), float('inf'), float('inf')

        # Early stopping check
        if rf_err <= 1e-2 and u3_err <= 1e-2:
            print("\nEarly stopping: Errors are within tolerance.")
            print(f"\nOptimization completed. Best parameters: {p}, Best score: {score:.6f}")
            return p

        pbest_scores.append(score)
        print(f"Particle {i+1}/{n_particles}: Score = {score:.6f}")

    # Find initial global best
    gbest_index = pbest_scores.index(min(pbest_scores))
    gbest = pbest[gbest_index]
    gbest_score = pbest_scores[gbest_index]

    print(f"Initial global best: {gbest} with score: {gbest_score:.6f}")

    # Main PSO iteration loop
    for iteration in range(iterations):
        print(f"\n--- Iteration {iteration + 1}/{iterations} ---")

        for i in range(n_particles):
            r1, r2 = random.random(), random.random()
            v_old = velocities[i]
            p_curr = particles[i]
            p_best = pbest[i]

            # Update velocity for each component
            new_velocity = tuple(
                w * v_old[j] +
                c1 * r1 * (p_best[j] - p_curr[j]) +
                c2 * r2 * (gbest[j] - p_curr[j])
                for j in range(2)
            )
            
            # Clamp velocity to prevent excessive movement
            clamped_velocity = tuple(
                np.clip(new_velocity[j], -max_velocity[j], max_velocity[j])
                for j in range(2)
            )
            
            velocities[i] = clamped_velocity

            # Update position
            new_position = tuple(
                p_curr[j] + clamped_velocity[j] for j in range(2)
            )
            
            # Apply bounds with reflection to prevent boundary sticking
            bounded_position = apply_bounds_with_reflection(new_position, bounds)
            particles[i] = bounded_position

            # Evaluate new particle
            try:
                score, rf_err, u3_err = objectiveFunction(bounded_position, realData, umatFile, inputFile)
            except Exception as e:
                print(f"Error in objective function for particle {i}: {e}")
                score, rf_err, u3_err = float('inf'), float('inf'), float('inf')

            # Early stopping check
            if rf_err <= 1e-2 and u3_err <= 1e-2:
                print("\nEarly stopping: Errors are within tolerance.")
                print(f"\nOptimization completed. Best parameters: {bounded_position}, Best score: {score:.6f}")
                return bounded_position

            # Update personal best
            if score < pbest_scores[i]:
                pbest[i] = bounded_position
                pbest_scores[i] = score
                print(f"Particle {i+1}: New personal best = {score:.6f}")

        # Update global best
        current_best_index = pbest_scores.index(min(pbest_scores))
        current_best_score = pbest_scores[current_best_index]

        if current_best_score < gbest_score:
            gbest = pbest[current_best_index]
            gbest_score = current_best_score
            print(f"New global best: {gbest} with score: {gbest_score:.6f}")

        print(f"Current best: {gbest} | Score: {gbest_score:.6f}")
        
        # Optional: Add convergence check
        if iteration > 0:
            velocity_magnitudes = [np.sqrt(v[0]**2 + v[1]**2) for v in velocities]
            avg_velocity = np.mean(velocity_magnitudes)
            print(f"Average velocity magnitude: {avg_velocity:.6f}")
            
            # Stop if particles have converged (very low velocities)
            if avg_velocity < 1e-6:
                print("Early stopping: Particles have converged (low velocity)")
                break

    print(f"\nOptimization completed. Best parameters: {gbest}, Best score: {gbest_score:.6f}")
    return gbest


# --- Main execution ---
if __name__ == "__main__":
    # TODO: Load realData
    realData = (realPeak, realPileUp)

    # Initial guess as a tuple
    initial_parameter = (initial_peak_parameter, initial_pileUp_parameter)

    if realPeak == -1:  # Default value, likely needs to be changed
        print("Warning: realPeak is set to default value (-1). Please verify this is correct.")

    if realPileUp == -1:  # Default value, likely needs to be changed
        print("Warning: realPileUp is set to default value (-1). Please verify this is correct.")
    
    # Run the PSO optimization
    try:
        finalParameter = pso(n_particles, iterations, w, c1, c2, initial_parameter)
        print(f"\nOptimization completed!")
        print(f"Final optimized parameter: ({finalParameter[0]:.6f}, {finalParameter[1]:.6f})")

        print(f"Final parameter applied to kMaterialParam.f and kCRSS.f")
        
    except Exception as e:
        print(f"Optimization failed: {e}")
        
    # Close log file
    log_file.close()
        
    # print(extractData("unimp_whole_amp4_ID100.odb"))