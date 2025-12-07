def plot_per_peptide_comparison(df, peptide, diseased_columns, healthy_columns, clean_label=None):
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
    clean_label : str, optional
        Pre-computed clean label for the peptide.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

    # Get clean label for title
    if clean_label is None:
        if 'clean_label' in df.columns:
            clean_label = df.loc[peptide, 'clean_label']
        else:
            clean_label = peptide[:60] + '...' if len(peptide) > 60 else peptide

    # Color palette
    all_groups = [('ROHHAD', diseased_columns, '#E74C3C'),
                  ('Healthy', healthy_columns, '#3498DB')]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor('white')
    
    # Left panel: Strip plot with jitter
    ax1 = axes[0]
    group_data = {}

    for i, (label, columns, color) in enumerate(all_groups):
        y_vals = np.array(df.loc[peptide, columns], dtype=float).flatten()
        group_data[label] = y_vals
        
        if len(y_vals) > 0:
            x_vals = i + np.random.normal(0, 0.08, size=len(y_vals))
            ax1.scatter(x_vals, y_vals, alpha=0.7, s=100, color=color,
                       label=f'{label} (n={len(y_vals)})', edgecolors='white', linewidths=1)
            
            mean_val = np.mean(y_vals)
            std_val = np.std(y_vals, ddof=1) if len(y_vals) > 1 else 0
            ax1.hlines(mean_val, i - 0.25, i + 0.25, colors='#2C3E50', linewidths=3)
            if std_val > 0:
                ax1.errorbar(i, mean_val, yerr=std_val, color='#2C3E50', capsize=8, capthick=2, linewidth=2)
            ax1.text(i + 0.35, mean_val, f'{mean_val:.2f}±{std_val:.2f}', 
                    fontsize=10, va='center', fontweight='bold', color='#2C3E50')

    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(['ROHHAD', 'Healthy'], fontsize=12, fontweight='bold')
    ax1.set_ylabel('Log2 Fold Change (vs MockIP)', fontsize=12, fontweight='bold', color='#2C3E50')
    ax1.axhline(0, color='#7F8C8D', linestyle='--', linewidth=1.5, alpha=0.7)
    ax1.legend(loc='upper right', frameon=True, fancybox=True)
    ax1.yaxis.grid(True, linestyle='--', alpha=0.3, color='#BDC3C7')
    ax1.set_axisbelow(True)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_title('Individual Sample Values', fontsize=11, fontweight='bold', color='#2C3E50')

    # Right panel: Boxplot only (robust to edge cases)
    ax2 = axes[1]
    box_data = []
    labels = []
    bcolors = []
    
    for label, columns, color in all_groups:
        y_vals = np.array(df.loc[peptide, columns], dtype=float).flatten()
        if len(y_vals) > 0:
            box_data.append(y_vals)
            labels.append(label)
            bcolors.append(color)
    
    if box_data:
        bp = ax2.boxplot(box_data, positions=range(len(box_data)), widths=0.5, 
                         patch_artist=True, showfliers=True)
        for patch, color in zip(bp['boxes'], bcolors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        for median in bp['medians']:
            median.set_color('#2C3E50')
            median.set_linewidth(2)
    
    ax2.set_xticks(range(len(labels)))
    ax2.set_xticklabels(labels, fontsize=12, fontweight='bold')
    ax2.set_ylabel('Log2 Fold Change (vs MockIP)', fontsize=12, fontweight='bold', color='#2C3E50')
    ax2.axhline(0, color='#7F8C8D', linestyle='--', linewidth=1.5, alpha=0.7)
    ax2.yaxis.grid(True, linestyle='--', alpha=0.3, color='#BDC3C7')
    ax2.set_axisbelow(True)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_title('Distribution Comparison', fontsize=11, fontweight='bold', color='#2C3E50')

    # P-value if both groups have enough data
    if len(group_data.get('ROHHAD', [])) >= 2 and len(group_data.get('Healthy', [])) >= 2:
        _, pval = stats.ttest_ind(group_data['ROHHAD'], group_data['Healthy'], equal_var=False)
        sig_text = f"p = {pval:.2e}" if pval < 0.001 else f"p = {pval:.4f}"
        if pval < 0.001:
            sig_text += " ***"
        elif pval < 0.01:
            sig_text += " **"
        elif pval < 0.05:
            sig_text += " *"
        ax2.text(0.5, 0.95, sig_text, transform=ax2.transAxes, fontsize=11, 
                fontweight='bold', color='#8E44AD', ha='center', va='top')

    plt.suptitle(clean_label, fontsize=13, fontweight='bold', color='#1A252F', y=1.02)
    plt.tight_layout()
    plt.show()
