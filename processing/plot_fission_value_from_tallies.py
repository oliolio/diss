import re
import numpy as np
import matplotlib.pyplot as plt

class MeshTallyParser:
    def __init__(self, filename):
        self.filename = filename
        self.mesh_data = {}

    def parse_file(self):
        # Regex patterns to capture mesh index and fission rate
        mesh_pattern = re.compile(r'Mesh Index \((\d+), (\d+)\)')
        rate_pattern = re.compile(r'Fission Rate\s+([\d\.e\+\-]+) \+/- ([\d\.e\+\-]+)')

        with open(self.filename, 'r') as file:
            lines = file.readlines()

        current_mesh = None

        for line in lines:
            mesh_match = mesh_pattern.search(line)
            rate_match = rate_pattern.search(line)

            if mesh_match:
                current_mesh = (int(mesh_match.group(1)), int(mesh_match.group(2)))

            if rate_match and current_mesh:
                fission_rate = float(rate_match.group(1))
                self.mesh_data[current_mesh] = fission_rate
                current_mesh = None  # Reset for the next mesh

    def plot_fission_rates(self):
        # Determine the size of the grid based on the mesh indices
        max_x = max(index[0] for index in self.mesh_data.keys())
        max_y = max(index[1] for index in self.mesh_data.keys())

        # Create a grid and populate it with fission rates
        fission_grid = np.zeros((max_x, max_y))

        for (x, y), rate in self.mesh_data.items():
            fission_grid[x-1, y-1] = np.power(rate, 1)  # Subtract 1 to convert to 0-based index
        # fission_grid[fission_grid==0] = np.nan
        # Plot the fission rates as a heatmap
        plt.imshow(fission_grid, cmap='hot', interpolation='nearest')
        plt.colorbar(label='(Fission Rate)')
        plt.title('Fission Rate Heatmap')
        plt.xlabel('Mesh X Index')
        plt.ylabel('Mesh Y Index')
        plt.show()
        plt.savefig('plotey.png') 

# Example usage
if __name__ == "__main__":
    parser = MeshTallyParser('/home/hpcatta/submit1/75_1.35/tallies.out')
    parser.parse_file()
    parser.plot_fission_rates()
