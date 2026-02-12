const scriptData = {
    "scripts/cnsr_diploid_conflictic_selection/demo_simulation.py": `# Simulation Demo Script
import numpy as np

def run_simulation(steps=100):
    """Simple demo of a population simulation step."""
    population = np.random.rand(steps)
    print(f"Simulation finished for {steps} steps.")
    return np.mean(population)

result = run_simulation()
print(f"Average population value: {result}")
`,
    "scripts/cnsr_diploid_conflictic_selection/testing_funs.py": `#
#  Collabrative non-self recognition (CNSR) based self-incompatibility (SI) systems in plants
##
import numpy as np
########################################################################################################

#
# finding probabilities of different females; SC (self), SC (out-cross), SI (always out-crossed)
def diploid_female_probability(input_Var):

    alpha, delta, fi, fhsc, ffsc = input_Var
    fi, fhsc, ffsc = fi/(fi+fhsc+ffsc), fhsc/(fi+fhsc+ffsc), ffsc/(fi+fhsc+ffsc)

    alpha_prime = (alpha/2)/(1 - (alpha/2)) # alpha/2
    tot_dying_fraction = ffsc*alpha*delta + fhsc*alpha_prime*delta
    surviving_fraction = 1 - tot_dying_fraction

    # full self-compatible self-fertilized
    ffsc_self = ffsc * alpha * (1 - delta) / surviving_fraction
    # ffsc_self = 1 * alpha * (1 - delta) / surviving_fraction
    # full self-compatible out-crossed
    ffsc_out = ffsc * (1 - alpha) / surviving_fraction
    # ffsc_out = 1 * (1 - alpha) / surviving_fraction

    # half self-compatible self-fertilized
    fhsc_self = fhsc * alpha_prime * (1 - delta) / surviving_fraction
    # fhsc_self = 1 * alpha_prime * (1 - delta) / surviving_fraction
    # half self-compatible out-crossed
    fhsc_out = fhsc * (1 - alpha_prime) / surviving_fraction
    # fhsc_out = 1 * (1 - alpha_prime) / surviving_fraction

    # self-incompatible out-crossed
    fi_out = fi / surviving_fraction
    # fi_out = 1 / surviving_fraction

    # output = [ffsc_self, ffsc_out, fhsc_self, fhsc_out, fi_out]
    output = [ffsc_self + ffsc_out, fhsc_self + fhsc_out, fi_out]
    output = np.array(output)/sum(output)

    return output*0.5#, surviving_fraction
    
input_Var = [0.6, 0.98, 0.48, 0.16, 0.36]
diploid_pop_size = 500
print('FSI, HSC, FSC:', (diploid_pop_size*np.array(input_Var[2:])).astype(int))
f_outcome = diploid_female_probability(input_Var)
print(sum(f_outcome), f_outcome)

########################################################################################################
`,
    "scripts/cnsr_diploid_phase_transitions/demo_simulation.py": `# Simulation Demo Script
import numpy as np

def run_simulation(steps=100):
    """Simple demo of a population simulation step."""
    population = np.random.rand(steps)
    print(f"Simulation finished for {steps} steps.")
    return np.mean(population)

result = run_simulation()
print(f"Average population value: {result}")
`
};