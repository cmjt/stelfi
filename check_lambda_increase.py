import re

def parse_simulation_output(filename):
    lambda_values = []
    time_values = []

    with open(filename, 'r') as file:
        for line in file:
            if "lambdaXsasep" in line:
                if "SUM" not in line:

                    values = [float(x) for x in line[14:].split(', ')]
                    lambda_values.append(values)

            if "Confirmed time:" in line:
                time_values.append(float(line[16:]))


    return lambda_values, time_values

def check_for_increasing_indices(lambda_values, time_values):
    increasing_indices = set()
    # increasing_indices = list()
    
    # Compare values at each index across consecutive iterations
    for i in range(len(lambda_values) - 1):
        current_values = lambda_values[i]
        next_values = lambda_values[i + 1]
        
        for idx in range(min(len(current_values), len(next_values))):
            if next_values[idx] > current_values[idx] and current_values[idx] != 0:
                increasing_indices.add(idx)
                # increasing_indices.add((idx, time_values[i]))
                # increasing_indices.append((idx, time_values[i]))

    
    return increasing_indices

def main():
    filename = 'output.txt'
    lambda_values, time_values = parse_simulation_output(filename)
    increasing_indices = check_for_increasing_indices(lambda_values, time_values)
    
    print("Indices where lambdaXsasep elements increase in the next iteration:")
    print(increasing_indices)
    # print(increasing_indices[0:30])

if __name__ == "__main__":
    main()