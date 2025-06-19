import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# Define the path to the results directory
results_dir = "results"

# Load solution and error data from the results directory
solution_df = pd.read_csv(os.path.join(results_dir, "solution.csv"))
error_df = pd.read_csv(os.path.join(results_dir, "error.csv"))

# Extract and sort unique x and y values
x_vals = np.sort(solution_df['x'].unique())
y_vals = np.sort(solution_df['y'].unique())
nx = len(x_vals)
ny = len(y_vals)

# Reshape u_numeric and error fields into 2D arrays for plotting
u_numeric = solution_df['u_numeric'].values.reshape((ny, nx))
error = error_df['error'].values.reshape((ny, nx))

# Create a meshgrid for 3D plotting
X, Y = np.meshgrid(x_vals, y_vals)

# Plot the numerical solution and error field
fig = plt.figure(figsize=(12, 5))

# Subplot for numerical solution
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(X, Y, u_numeric, cmap='viridis')
ax1.set_title("Numerical Solution (u)")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("u")

# Subplot for error field
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(X, Y, error, cmap='inferno')
ax2.set_title("Error Field (u - u_exact)")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("Error")

plt.tight_layout()

# Save the combined plot
plot_file = os.path.join(results_dir, "solution_and_error_surface.png")
plt.savefig(plot_file, dpi=300)
print(f"Plot saved to {plot_file}")

plt.show()
