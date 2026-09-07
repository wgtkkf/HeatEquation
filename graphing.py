# =============================================
# This code was generated supported by Gemini =
# =============================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

class Graph:
    def __init__(self, x1, y1):
        self.x1 = x1
        self.y1 = y1

    def display(self):
        csfont = {'fontname':'DejaVu Sans'} 
        plt.figure() 
        plt.grid(True, which="both", ls="-", alpha=0.5) # Enhanced grid for log scale

        # Plot 1: Analytical (Python)
        plt.plot(self.x1, self.y1, 'blue', label="Analytical solution")         

        ## graph information    
        plt.xlabel(r'Distance of y-axis [m]', **csfont)            
        plt.ylabel(r'Temperature [$^\circ\text{C}$]', **csfont)

        # font for legend
        font = font_manager.FontProperties(family='DejaVu Sans',
                                           weight='bold',
                                           style='normal', size=10)
        plt.legend(loc='lower right', prop=font) 

        # Log-Log Style
        #plt.xscale('log')
        #plt.yscale('log')        

        #
        plt.xlim(0, 0.05) # Sets x-axis from 0 to 100
        plt.ylim(0, 100)  # Useful for log scales

        # plot options
        plt.xticks([0, 0.01, 0.02, 0.03, 0.04, 0.05], **csfont)
        plt.yticks([0, 50, 100], **csfont)

        # graph save & display
        plt.savefig("comparison.png", dpi=300) 
        plt.show()                  
        return 0

# ==========================================
# DATA LOADING LOGIC (EXTERNAL)
# ==========================================

# 2. Load Python data (Dynamic row length)
py_data = np.loadtxt("./excel/analytical.txt", skiprows=1)
x_py = py_data[:, 0]
y_py = py_data[:, 1]

# 3. Create the instance with BOTH datasets
my_graph = Graph(x_py, y_py)

# 4. Generate the single combined plot
my_graph.display()