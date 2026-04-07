#calculate the total charge of the input xyz file
def calculate_total_charge(xyz_file):
    total_charge = 0
    with open(xyz_file, 'r') as file:
        lines = file.readlines()
        # Skip the first two lines (number of atoms and comment)
        for line in lines[2:]:
            parts = line.split()
            element = parts[0]
            # Define a simple charge dictionary for common elements
            charge_dict = {
                'H': 1,
                'He': 0,
                'Li': 1,
                'Be': 2,
                'B': 3,
                'C': 4,
                'N': -3,
                'O': -2,
                'F': -1,
                'Ne': 0,
                'Na': 1,
                'Mg': 2,
                'Al': 3,
                'Si': 4,
                'P': -3,
                'S': -2,
                'Cl': -1,
                'Ar': 0,
                'Er': 3,
                'Y': 3,
                'F': -1
            }
            total_charge += charge_dict.get(element, 0)  # Default to 0 if element not found
    return total_charge

if __name__ == "__main__":
    xyz_file = 'aimp.xyz'
    total_charge = calculate_total_charge(xyz_file)
    print(f'Total Charge: {total_charge}')