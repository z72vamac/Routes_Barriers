The folder contains 3 files:

- barriers2.csv
- barriers3.csv
- circles.csv

Each barrier file is labeled with a number representing the example used in the manuscript. 

To obtain a feasible problem, it is mandatory to take a circle instance with a barrier instance.

Here, an example Python code to load correctly the instance example of the TSPHN is provided:

# Segments loading
segments = np.genfromtxt('./instances_example/barriers2.csv', delimiter = ',')

barriers = []
for segment in segments:
    barriers.append([[segment[0], segment[1]], [segment[2], segment[3]]])

# Neighbourhoods loading
circles = np.genfromtxt('./instances_example/circles.csv', delimiter = ',')

neighbourhoods = [Circle(center = [centerx, centery], radii = radius) for centerx, centery, radius in circles]