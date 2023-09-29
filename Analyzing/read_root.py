import ROOT
import numpy as np
import matplotlib.pyplot as plt

# Open the ROOT file for reading
root_file = ROOT.TFile("/home/teikiet/Data/kiet.root", "READ")

keys = root_file.GetListOfKeys()

# Loop through the keys and identify TGraph objects
G_name = []
for key in keys:
    obj = key.ReadObj()
    if isinstance(obj, ROOT.TGraph):
        graph_name = key.GetName()
        G_name.append(graph_name)

def read_root_file(graph_name):
    # Get the TGraph by name
    graph = root_file.Get(graph_name)

    # Check if the TGraph exists
    if graph:
        # Access information about the TGraph
        n_points = graph.GetN()

        # Convert LowLevelViews to Python lists
        x_values = np.frombuffer(graph.GetX(), dtype=np.float64)
        y_values = np.frombuffer(graph.GetY(), dtype=np.float64)
        return graph_name, n_points, x_values, y_values
    else:
        print("TGraph", graph_name, "not found in the ROOT file.")
        
graph_name, n_points, x_values, y_values = read_root_file(G_name[2])
# Close the ROOT file when done
root_file.Close()

print(graph_name, n_points)
plt.plot(x_values, y_values, "-", label = f"{graph_name}, {n_points}")
plt.savefig("test.pdf")


