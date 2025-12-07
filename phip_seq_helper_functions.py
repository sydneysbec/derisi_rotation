def plot_per_peptide_comparison(df, peptide, diseased_columns, healthy_columns):
    """
    Create per-peptide comparison plot between diseased and healthy sample columns.

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with samples as columns and peptides as rows
    peptide : str
        Peptide index/name in df
    diseased_columns : list of str
        List of column names corresponding to diseased samples
    healthy_columns : list of str
        List of column names corresponding to healthy samples
    """
    import numpy as np
    import matplotlib.pyplot as plt

    try:
        # Checks that all columns exist for peptide
        df.loc[peptide, diseased_columns]
        df.loc[peptide, healthy_columns]
    except KeyError:
        raise ValueError(f"Peptide '{peptide}' or some column names not found in DataFrame.")

    all_groups = [('Diseased', diseased_columns, '#d62728'),
                  ('Healthy', healthy_columns, '#1f77b4')]  # red, blue

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    # Strip plot
    ax1 = axes[0]
    xticks = []
    xticklabels = []

    for i, (label, columns, color) in enumerate(all_groups):
        y_vals = df.loc[peptide, columns]
        # y_vals may be 1d np array or pd.Series, but always treat as array of values
        y_vals = np.array(y_vals).flatten()
        if len(y_vals) > 0:
            x_vals = i + np.random.normal(0, 0.04, size=len(y_vals))
            ax1.scatter(
                x_vals, y_vals, alpha=0.6, s=80, color=color,
                label=f'{label} (n={len(y_vals)})', edgecolors='black', linewidths=0.5
            )
            mean_val = np.mean(y_vals)
            ax1.hlines(mean_val, i - 0.2, i + 0.2, colors='black', linewidths=3)
            ax1.text(i + 0.25, mean_val, f'{mean_val:.2f}', fontsize=9, va='center', fontweight='bold')
            xticks.append(i)
            xticklabels.append(label)

    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xticklabels, rotation=45, ha='right')
    ax1.set_ylabel('Log2 FC (over beads)', fontsize=12, fontweight='bold')
    ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # Boxplot
    ax2 = axes[1]
    box_data = []
    labels = []
    bcolors = []
    for label, columns, color in all_groups:
        y_vals = df.loc[peptide, columns]
        y_vals = np.array(y_vals).flatten()
        if len(y_vals) > 0:
            box_data.append(y_vals)
            labels.append(label)
            bcolors.append(color)
    bp = ax2.boxplot(box_data, labels=labels, patch_artist=True)
    for patch, color in zip(bp['boxes'], bcolors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax2.set_ylabel('Log2 FC (over beads)', fontsize=12, fontweight='bold')
    ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax2.grid(axis='y', alpha=0.3)
    plt.suptitle(f'{str(peptide)[:80]}', fontsize=10, fontweight='bold')
    plt.tight_layout()
    plt.show()
