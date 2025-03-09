import matplotlib.pyplot as plt

# Data values
thresholds = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]

male_spearman = [0.734, 0.723, 0.688, 0.659, 0.628, 0.601, 0.587, 0.583, 0.574]
male_pearson = [0.575, 0.553, 0.522, 0.499, 0.482, 0.472, 0.462, 0.463, 0.460]

female_spearman = [-0.520, -0.499, -0.493, -0.456, -0.420, -0.373, -0.349, -0.288, -0.219]
female_pearson = [-0.535, -0.524, -0.519, -0.496, -0.471, -0.440, -0.422, -0.380, -0.334]

def plot_with_regions_and_labels(ax, x, y, title, ylabel, color, label, is_male):
    # Define correlation regions
    for (low, high, shade_color, text) in [
        (0.0, 0.4, 'lightgray', "Low or no correlation"),
        (0.4, 0.6, 'lightblue', "Moderately strong correlation"),
        (0.6, 0.8, 'lightgreen', "Very strong correlation"),
        (0.8, 1.0, 'lightcoral', "Near-perfect correlation"),
        (-1.0, -0.8, 'lightcoral', "Near-perfect negative correlation"),
        (-0.8, -0.6, 'lightgreen', "Very strong negative correlation"),
        (-0.6, -0.4, 'lightblue', "Moderately strong negative correlation"),
        (-0.4, 0.0, 'lightgray', "Low or no negative correlation")
    ]:
        if is_male:  # Only plot positive ranges for male
            if low >= 0:
                ax.axhspan(low, high, color=shade_color, alpha=0.3)
                ax.text(1.05, (low + high) / 2, text, va='center', fontsize=8, bbox=dict(facecolor='white', alpha=0.7))
        else:  # Only plot negative ranges for female
            if high <= 0:
                ax.axhspan(low, high, color=shade_color, alpha=0.3)
                ax.text(1.05, (low + high) / 2, text, va='center', fontsize=8, bbox=dict(facecolor='white', alpha=0.7))

    ax.plot(x, y, marker='o', linestyle='-', color=color, label=label)
    ax.set_title(title)
    ax.set_xlabel("LFC cut-off")
    ax.set_ylabel(ylabel)
    
    if is_male:
        ax.set_ylim(0, 1)  # Male data range
    else:
        ax.set_ylim(-1, 0)  # Female data range
    
    ax.grid(True)
    ax.legend()

# Create subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Generate graphs
plot_with_regions_and_labels(axes[0, 0], thresholds, male_spearman, "Male - Spearman Correlation", "Spearman r-value", 'b', 'Spearman', is_male=True)
plot_with_regions_and_labels(axes[0, 1], thresholds, male_pearson, "Male - Pearson Correlation", "Pearson r-value", 'r', 'Pearson', is_male=True)
plot_with_regions_and_labels(axes[1, 0], thresholds, female_spearman, "Female - Spearman Correlation", "Spearman r-value", 'g', 'Spearman', is_male=False)
plot_with_regions_and_labels(axes[1, 1], thresholds, female_pearson, "Female - Pearson Correlation", "Pearson r-value", 'm', 'Pearson', is_male=False)

# Adjust layout and show plot
plt.tight_layout()
plt.show()
