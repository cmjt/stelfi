import re
import pandas as pd

# Initialize an empty list to store the data
data = []

# Read the file
with open('ratesep_vector_output.txt', 'r') as file:
    lines = file.readlines()

current_time = None
current_loci = None
current_bvn_value = None
current_transformed_bvn_value = None

count = 0 
for i, line in enumerate(lines):
    line = line.strip()
    if line.startswith("TIME:"):
        current_time = float(line.split(":")[1].strip())
        count += 1
    elif line.startswith("loci:"):
        loci_value1 = float(line.split(":")[1].strip())
        loci_value2 = float(lines[i + 1].strip())
        current_loci = (loci_value1, loci_value2)
    elif line.startswith("bvn_value:"):
        current_bvn_value = float(line.split(":")[1].strip())
    elif line.startswith("transformed bvn_value:"):
        current_transformed_bvn_value = float(line.split(":")[1].strip())
    elif line == "====":
        if current_time is not None and current_loci is not None:
            data.append({
                'time': current_time,
                'loci': current_loci,
                'bvn_value': current_bvn_value,
                'transformed_bvn_value': current_transformed_bvn_value
            })
        # Reset current values
        current_loci = None
        current_bvn_value = None
        current_transformed_bvn_value = None
    
    if count == 2:
        break

# Convert data to a Pandas DataFrame
df = pd.DataFrame(data)

# Sort the data by time
df = df.sort_values(by='time')

# Group by loci and check for increases in transformed bvn_value
increases = []

for loci, group in df.groupby('loci'):
    transformed_values = group['transformed_bvn_value'].values
    if any(y > x for x, y in zip(transformed_values, transformed_values[1:])):
        increases.append(loci)

# Output the loci with increases in transformed bvn_value over time
print("Loci with increases in transformed bvn_value over time:")
for loci in increases:
    print(loci)
