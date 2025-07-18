import random
import sys
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
realPeak = 2  # TODO: Set this to actual peak from experimental data
umatFile = r'umat_directory\umat.for' # .for file template

#-------------------------  Modify initial guess Here  --------------------------------

# Set initial parameter
initial_parameter = 360

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
def pso(n_particles, iterations, w, c1, c2, initial_guess=None):
    """
    Particle Swarm Optimization algorithm
    """
    print("Starting PSO optimization...")
    
    # Initialize particles around the initial guess if provided
    if initial_guess is not None:
        print(f"Using initial guess: {initial_guess}")
        # Create particles around the initial guess with some variation
        search_range = initial_guess * 0.5  # 50% variation around initial guess
        particles = []
        for i in range(n_particles):
            if i == 0:
                # First particle is exactly the initial guess
                particles.append(initial_guess)
            else:
                # Other particles are distributed around the initial guess
                variation = random.uniform(-search_range, search_range)
                particles.append(initial_guess + variation)
    else:
        # Fallback to random initialization if no initial guess
        particles = [random.uniform(-10, 10) for _ in range(n_particles)]
    
    velocities = [0.0 for _ in range(n_particles)]
    pbest = particles[:]
    
    # Evaluate initial particles
    print("Evaluating initial particles...")
    pbest_scores = []
    for i, p in enumerate(particles):
        score = objectiveFunction(p, realPeak, umatFile, inputFile)
        pbest_scores.append(score)
        print(f"Particle {i+1}/{n_particles}: {score:.6f}")

    # Initialize global best
    gbest_index = pbest_scores.index(min(pbest_scores))
    gbest = pbest[gbest_index]
    gbest_score = pbest_scores[gbest_index]
    
    print(f"Initial global best: {gbest:.4f} with score: {gbest_score:.6f}")

    # Main PSO loop
    for iteration in range(iterations):
        print(f"\nIteration {iteration + 1}/{iterations}")
        
        for i in range(n_particles):
            r1 = random.random()
            r2 = random.random()

            # Update velocity
            velocities[i] = (w * velocities[i] +
                             c1 * r1 * (pbest[i] - particles[i]) +
                             c2 * r2 * (gbest - particles[i]))

            # Update position
            particles[i] += velocities[i]

            # Evaluate new position
            score = objectiveFunction(particles[i], realPeak, umatFile, inputFile)
            
            # Update personal best
            if score < pbest_scores[i]:
                pbest[i] = particles[i]
                pbest_scores[i] = score

        # Update global best
        current_best_index = pbest_scores.index(min(pbest_scores))
        current_best_score = pbest_scores[current_best_index]
        
        if current_best_score < gbest_score:
            gbest = pbest[current_best_index]
            gbest_score = current_best_score
            print(f"New global best found: {gbest:.4f} with score: {gbest_score:.6f}")
        
        print(f"Current global best: {gbest:.4f} with score: {gbest_score:.6f}")

    return gbest

# --- Main execution ---
if __name__ == "__main__":
    # TODO: Load realPeak from experimental.csv
    
    if not realPeak:
        print("Warning: realPeak is set to 0. Please set the correct experimental peak value.")
    
    # Run the PSO optimization
    try:
        finalParameter = pso(n_particles, iterations, w, c1, c2, initial_parameter)
        print(f"\nOptimization completed!")
        print(f"Final optimized parameter: {finalParameter:.6f}")
        
        print(f"Final parameter applied to kMaterialParam.f")
        
    except Exception as e:
        print(f"Optimization failed: {e}")
        '''
    print(extractData("unimp_whole_amp3_ID500.odb"))
    #'''