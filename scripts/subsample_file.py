import random
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

# Parameters
chunk_size = 50
total_lines = 1000

# Read the file
with open(infile) as f:
    lines = f.readlines()

# # Shuffle the lines
# random.shuffle(lines)

# Select random chunks of 50 lines at a time
selected_lines = []
for i in range(0, total_lines, chunk_size):
    selected_lines.extend(lines[i:i+chunk_size])

# Save the output
with open(outfile, "w") as f:
    f.writelines(selected_lines)