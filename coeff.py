from itertools import combinations_with_replacement

# Define atom types and their masses (for reference)
atom_types = {
    1: ("HGP2", 1.008),
    2: ("HGA2", 1.008),
    3: ("CG324", 12.01),
    4: ("CG321", 12.01),
    5: ("OG301", 16.0),
    6: ("NG3P3", 14.01),
    7: ("Pb", 207.2),
    8: ("I", 126.9)
}

# Define Buckingham parameters for specific pairs (A, rho, C)
# Format: (type1, type2): (A, rho, C)
buckingham_params = {
    (7, 7): (77395897.29267, 0.1443838, 0.0),
    (7, 8): (113845.746311, 0.3539107, 0.0),

    (8, 8): (25072.672438, 0.5304387, 766.644),
}
# Define Lennard-Jones parameters for other pairs (epsilon, sigma)
# Format: (type1, type2): (epsilon, sigma)
lj_params = {
    (1, 1): (0.0460, 0.2245),
    (2, 2): (0.0350, 1.3400),
    (3, 3): (0.0550, 2.1750),
    (4, 4): (0.0560, 2.0100),
    (5, 5): (0.1000, 1.6500),
    (6, 6): (0.2000, 1.8500)
}

# Generate all unique pairs (i, j) with i <= j
all_pairs = list(combinations_with_replacement(atom_types.keys(), 2))

# Generate pair_coeff commands
pair_coeff_commands = []

for i, j in all_pairs:
    if 7 in (i, j):
        if 1 in (i, j) or 2 in(i, j):
            pair_coeff_commands.append(f"pair_coeff {i} {j} lj/charmm/coul/long 1.4E-002 2.70999")
        else:
            # Use Buckingham if either atom is Pb or I
            key = tuple(sorted((i, j)))
            if key in buckingham_params:
                A, rho, C = buckingham_params[key]
            else:
                # Default Buckingham parameters if not specified
                A, rho, C = (32690390.937995, 0.1509470, 0.0)
            pair_coeff_commands.append(f"pair_coeff {i} {j} buck/coul/long {A} {rho} {C}")
    elif 8 in (i, j):
        if 1 in (i, j) or 2 in(i, j):
            pair_coeff_commands.append(f"pair_coeff {i} {j} lj/charmm/coul/long 5.74E-002 3.1")
        else:
            # Use Buckingham if either atom is Pb or I
            key = tuple(sorted((i, j)))
            if key in buckingham_params:
                A, rho, C = buckingham_params[key]
            else:
                # Default Buckingham parameters if not specified
                A, rho, C = (112936.7142130, 0.3424260, 0.0)
            pair_coeff_commands.append(f"pair_coeff {i} {j} buck/coul/long {A} {rho} {C}")
    else:
        # Use Lennard-Jones
        key = tuple(sorted((i, j)))
        if key in lj_params:
            epsilon, sigma = lj_params[key]
        else:
            # Use geometric mixing rule for LJ if not explicitly defined
            epsilon_i, sigma_i = lj_params.get((i, i), (0.1, 1.0))
            epsilon_j, sigma_j = lj_params.get((j, j), (0.1, 1.0))
            epsilon = (epsilon_i * epsilon_j) ** 0.5
            sigma = (sigma_i + sigma_j) / 2
        pair_coeff_commands.append(f"pair_coeff {i} {j} lj/charmm/coul/long {epsilon:.4f} {sigma:.4f}")

# Output the commands
with open('paircoef.csv', 'w') as f:
    for cmd in pair_coeff_commands:
        f.write(f'{cmd}\n')

