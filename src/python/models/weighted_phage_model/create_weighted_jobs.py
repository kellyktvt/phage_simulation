import random

# Set random seed for reproducibility
random.seed(42)

# Parameter combinations
fops = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
opt_weight = 1.0
nonopt_weight = 0.25

# Generate 10 random seeds between 1 and 1000
seeds = [random.randint(1, 100) for _ in range(10)]

# Generate the commands
commands = []
for fop in fops:
    for seed in seeds:
        cmd = f"python3 src/python/models/weighted_phage_model/revised_weighted_phage_model.py {fop} {opt_weight} {nonopt_weight} {seed}"
        commands.append(cmd)

# Write to jobs file
with open('phage_simulation/src/python/models/weighted_phage_model/jobs_revised_weighted.txt', 'w') as f:
    for cmd in commands:
        f.write(cmd + '\n')
