import argparse
import matplotlib.pyplot as plt
import numpy as np
import gzip

def read_data(file_path):
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as f:
            data = [tuple(map(int, line.split()))[::-1] for line in f]  # Swap columns
    else:
        with open(file_path, 'r') as f:
            data = [tuple(map(int, line.split()))[::-1] for line in f]  # Swap columns
    return data

def plot_insert_sizes(data, output_file, min_x=None, max_x=None, min_y=None, max_y=None):
    # Filter data based on min_x and max_x
    if min_x is not None or max_x is not None:
        data = [(length, freq) for length, freq in data if (min_x is None or length >= min_x) and (max_x is None or length <= max_x)]
    
    if not data:
        print("No data points within the specified x range.")
        return
    
    lengths, frequencies = zip(*data)
    
    min_x_data, max_x_data = min(lengths), max(lengths)
    min_y_data, max_y_data = min(frequencies), max(frequencies)
    
    plt.figure(figsize=(10, 6))
    plt.bar(lengths, frequencies, color='darkblue', edgecolor='darkblue', width=0.8)
    plt.xlabel("Insert Size", fontsize=14, fontweight='bold')
    plt.ylabel("Frequency", fontsize=14, fontweight='bold')
    plt.title("Insert Size Distribution", fontsize=16, fontweight='bold')
    
    plt.xlim(min_x if min_x is not None else min_x_data, max_x if max_x is not None else max_x_data)
    plt.ylim(min_y if min_y is not None else min_y_data, max_y if max_y is not None else max_y_data)
    
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")
    print(f"X-axis range: {min_x_data} - {max_x_data}")
    print(f"Y-axis range: {min_y_data} - {max_y_data}")

def main():
    parser = argparse.ArgumentParser(description="Plot insert size distribution for NGS data.")
    parser.add_argument("-i", "--input", required=True, help="Input file containing insert size data (text or gzipped text)")
    parser.add_argument("-o", "--output", required=True, help="Output file name for the plot (e.g., output.png)")
    parser.add_argument("--min_x", type=int, default=None, help="Minimum value for x-axis")
    parser.add_argument("--max_x", type=int, default=None, help="Maximum value for x-axis")
    parser.add_argument("--min_y", type=int, default=None, help="Minimum value for y-axis")
    parser.add_argument("--max_y", type=int, default=None, help="Maximum value for y-axis")
    
    args = parser.parse_args()
    
    data = read_data(args.input)
    plot_insert_sizes(data, args.output, args.min_x, args.max_x, args.min_y, args.max_y)

if __name__ == "__main__":
    main()
