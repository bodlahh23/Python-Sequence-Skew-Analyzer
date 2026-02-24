import matplotlib.pyplot as plt

# -----------------------
# Example DNA sequence
# -----------------------
dna_sequence = (
    "ATGCGTACGTAGCTAGCTAGCTAGCTAGCGTACGTAGCTAGCTAGCATGCTAGCTAGCTA"
)

# -----------------------
# Function: compute GC content, GC skew, and AT skew in sliding windows
# -----------------------
def compute_skews(sequence, window_size=10):
    gc_content_values = []
    gc_skew_values = []
    at_skew_values = []

    for i in range(0, len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        G = window.count('G')
        C = window.count('C')
        A = window.count('A')
        T = window.count('T')

        # GC content
        gc_content = ((G + C) / window_size) * 100
        gc_content_values.append(gc_content)

        # GC skew = (G - C) / (G + C), handle divide by zero
        if (G + C) != 0:
            gc_skew = (G - C) / (G + C)
        else:
            gc_skew = 0
        gc_skew_values.append(gc_skew)

        # AT skew = (A - T) / (A + T), handle divide by zero
        if (A + T) != 0:
            at_skew = (A - T) / (A + T)
        else:
            at_skew = 0
        at_skew_values.append(at_skew)

    return gc_content_values, gc_skew_values, at_skew_values

# -----------------------
# Run computation
# -----------------------
window_size = 10
gc_content, gc_skew, at_skew = compute_skews(dna_sequence, window_size)

# -----------------------
# Plot results
# -----------------------
positions = range(len(gc_content))

plt.figure(figsize=(12,4))
plt.plot(positions, gc_content, label='GC Content (%)', color='green', marker='o')
plt.plot(positions, gc_skew, label='GC Skew', color='blue', marker='x')
plt.plot(positions, at_skew, label='AT Skew', color='red', marker='s')
plt.title(f'DNA Sequence Analysis (Window Size={window_size})')
plt.xlabel('Window Start Position')
plt.ylabel('Value')
plt.legend()
plt.grid(True)
plt.show()