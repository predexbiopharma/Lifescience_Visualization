import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

ADSL_PATH = 'ADSL.csv'
DOSES     = {1: '100mg', 2: '200mg', 3: '400mg', 4: '800mg'}

cmap = LinearSegmentedColormap.from_list(
    'orr', ['#2166ac','#92c5de','#f7f7f7','#f4a582','#d6604d','#b2182b'], N=256
)


def build_df(adsl):
    cohorts    = sorted(adsl['COHORT'].dropna().unique(), reverse=True)
    tumortypes = sorted(adsl['TUMORTYPE'].dropna().unique())
    records    = []
    for c in cohorts:
        for t in tumortypes:
            sub = adsl[(adsl['COHORT'] == c) & (adsl['TUMORTYPE'] == t)]
            n   = len(sub)
            cr  = (sub['BESTRSPC'] == 'CR').sum() if n > 0 else 0
            pr  = (sub['BESTRSPC'] == 'PR').sum() if n > 0 else 0
            sd  = (sub['BESTRSPC'] == 'SD').sum() if n > 0 else 0
            pd_ = (sub['BESTRSPC'] == 'PD').sum() if n > 0 else 0
            orr = round((cr + pr) / n * 100, 1) if n > 0 else np.nan
            records.append({'cohort': int(c), 'tumor': t, 'n': n, 'orr': orr,
                            'cr': cr, 'pr': pr, 'sd': sd, 'pd': pd_})
    return pd.DataFrame(records), cohorts, tumortypes


def draw_heatmap(hm, cohorts, tumortypes, filename):
    nx = len(tumortypes)
    ny = len(cohorts)

    fig, ax = plt.subplots(figsize=(12, 5.5))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    for xi, t in enumerate(tumortypes):
        for yi, c in enumerate(cohorts):
            row      = hm[(hm['cohort'] == int(c)) & (hm['tumor'] == t)].iloc[0]
            n        = row['n']

            if n == 0:
                ax.add_patch(plt.Rectangle((xi - 0.5, yi - 0.5), 1, 1,
                                           facecolor='#e8e8e8', edgecolor='white',
                                           linewidth=2, zorder=2))
                ax.text(xi, yi, 'n/a', ha='center', va='center',
                        fontsize=9, color='#aaaaaa', style='italic', zorder=3)
                continue

            orr_val  = row['orr']
            norm_val = orr_val / 100.0
            color    = cmap(norm_val)

            ax.add_patch(plt.Rectangle((xi - 0.5, yi - 0.5), 1, 1,
                                       facecolor=color, edgecolor='white',
                                       linewidth=2.5, zorder=2))

            txt_color = 'white' if (norm_val > 0.65 or norm_val < 0.25) else '#1a1a1a'

            ax.text(xi, yi + 0.12, f'{orr_val:.0f}%',
                    ha='center', va='center',
                    fontsize=13, fontweight='bold',
                    color=txt_color, zorder=4)

            breakdown = f"CR:{row['cr']}  PR:{row['pr']}  SD:{row['sd']}  PD:{row['pd']}"
            ax.text(xi, yi - 0.18, breakdown,
                    ha='center', va='center',
                    fontsize=6.5, color=txt_color, alpha=0.85,
                    family='monospace', zorder=4)

            ax.text(xi, yi - 0.36, f'n={int(n)}',
                    ha='center', va='center',
                    fontsize=7.5, color=txt_color, alpha=0.7, zorder=4)

    ax.set_xlim(-0.5, nx - 0.5)
    ax.set_ylim(-0.5, ny - 0.5)
    ax.set_xticks(range(nx))
    ax.set_xticklabels(tumortypes, fontsize=11, fontweight='bold', color='#2c3e50')
    ax.set_yticks(range(ny))
    ax.set_yticklabels([f'Cohort {int(c)}  ({DOSES[int(c)]})' for c in cohorts],
                       fontsize=10, color='#2c3e50')
    ax.tick_params(length=0, pad=10)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xlabel('Tumor Type', fontsize=12, fontweight='bold',
                  color='#2c3e50', labelpad=12)
    ax.set_ylabel('Dose Cohort', fontsize=12, fontweight='bold',
                  color='#2c3e50', labelpad=12)
    ax.set_title('Objective Response Rate (ORR) by Cohort and Tumor Type',
                 fontsize=13, fontweight='bold', color='#1a1a2e', pad=16, loc='left')
    ax.text(0, 1.015,
            'Color intensity reflects ORR  ·  CR/PR/SD/PD counts per cell',
            transform=ax.transAxes, fontsize=8, color='#7f8c8d', style='italic')

    sm   = ScalarMappable(cmap=cmap, norm=Normalize(vmin=0, vmax=100))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.85, pad=0.025, aspect=20)
    cbar.set_label('ORR (%)', fontsize=10, color='#2c3e50', labelpad=8)
    cbar.set_ticks([0, 25, 50, 75, 100])
    cbar.ax.tick_params(labelsize=9, colors='#2c3e50', length=3)
    cbar.outline.set_edgecolor('#cccccc')

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'Saved: {filename}')


def main():
    adsl              = pd.read_csv(ADSL_PATH)
    hm, cohorts, tumortypes = build_df(adsl)
    draw_heatmap(hm, cohorts, tumortypes, 'response_heatmap.png')


if __name__ == '__main__':
    main()
