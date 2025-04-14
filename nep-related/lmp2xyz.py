def parse_masses(data_file):
    """
    Parse the Masses section in the LAMMPS data file to map atom types to atomic masses.
    """
    with open(data_file, 'r') as file:
        lines = file.readlines()

    start = False
    skip_next_line = False
    mass_dict = {}
    for line in lines:
        if 'Masses' in line:
            start = True
            skip_next_line = True  # 设置标志以跳过下一行
            continue
        if start:
            if skip_next_line:  # 跳过“Masses”后的第一行（通常是空行）
                skip_next_line = False
                continue
            if line.strip() == '':  # 遇到空行时结束
                break
            parts = line.split()
            mass_dict[int(parts[0])] = float(parts[1])
    return mass_dict
    
def determine_element(mass):
    """
    Determine the element based on the atomic mass.
    This function needs to be refined based on the specific mass-to-element mapping.
    """
    # This is a simplified mapping, needs to be refined
    if abs(mass - 1.00794) < 0.1:
        return 'H'
    elif abs(mass - 12.011) < 0.1:
        return 'C'
    elif abs(mass - 14.01) < 0.1:
        return 'N'
    elif abs(mass - 16.00) < 0.1:
        return 'O'
    # Add more elements as needed
    else:
        return 'X'  # Unknown element

def convert_to_xyz(data_file, xyz_file):
    """
    Convert LAMMPS data file to XYZ file format.
    """
    mass_dict = parse_masses(data_file)

    with open(data_file, 'r') as file:
        lines = file.readlines()

    atoms = []
    for line in lines:
        if line.startswith('Atoms'):  # Start of Atoms section
            break
    for line in lines[lines.index(line) + 2:]:
        if line.strip() == '':
            break
        parts = line.split()
        atom_type = int(parts[2])
        x, y, z = parts[4:7]
        print(mass_dict)
        element = determine_element(mass_dict[atom_type])
        atoms.append(f"{element} {x} {y} {z}")

    with open(xyz_file, 'w') as file:
        file.write(f"{len(atoms)}\n")
        file.write("Converted from LAMMPS data file\n")
        for atom in atoms:
            file.write(atom + '\n')

# Convert em.data to mol.xyz
convert_to_xyz('em.data', 'mol.xyz')
