import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ADSL_PATH   = 'ADSL.csv'
RESP_COLORS = {'CR':'#2ecc71','PR':'#3498db','SD':'#f39c12','PD':'#c41e3a','NE':'#95a5a6'}
RESP_ORDER  = ['CR','PR','SD','PD','NE']
DOSES       = {1:'100mg', 2:'200mg', 3:'400mg', 4:'800mg'}


def calc_pcts(subset):
    n = len(subset)
    if n == 0:
        return {r: 0 for r in RESP_ORDER}, n
    counts = subset['BESTRSPC'].fillna('NE').value_counts()
    return {r: counts.get(r, 0) / n * 100 for r in RESP_ORDER}, n


def draw_bor(adsl, filename):
    cohorts = sorted(adsl['COHORT'].dropna().unique())
    total   = len(adsl)

    fig = plt.figure(figsize=(16, 9), facecolor='white')
    ax1 = fig.add_axes([0.06, 0.22, 0.42, 0.68])
    ax2 = fig.add_axes([0.55, 0.22, 0.42, 0.68])

    for ax in [ax1, ax2]:
        ax.set_facecolor('white')
        ax.tick_params(length=0, colors='#2c3e50')
        ax.grid(axis='y', color='#e8e8e8', linewidth=0.8, zorder=0)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        ax.spines['left'].set_color('#cccccc')
        ax.spines['bottom'].set_color('#cccccc')

    # ── Panel 1: stacked bar by cohort ──
    x       = np.arange(len(cohorts))
    bar_w   = 0.55
    bottoms = np.zeros(len(cohorts))

    for resp in RESP_ORDER:
        vals = [calc_pcts(adsl[adsl['COHORT'] == c])[0][resp] for c in cohorts]
        ax1.bar(x, vals, bar_w, bottom=bottoms,
                color=RESP_COLORS[resp], edgecolor='white',
                linewidth=1.5, zorder=3, label=resp)
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v > 5:
                ax1.text(xi, b + v / 2, f'{v:.0f}%',
                         ha='center', va='center',
                         fontsize=8.5, fontweight='bold', color='white', zorder=4)
        bottoms += np.array(vals)

    for xi, c in enumerate(cohorts):
        _, n = calc_pcts(adsl[adsl['COHORT'] == c])
        ax1.text(xi, 102, f'n={n}', ha='center', va='bottom',
                 fontsize=9, color='#2c3e50', fontweight='bold')

    ax1.set_xticks(x)
    ax1.set_xticklabels([f'Cohort {int(c)}\n({DOSES[int(c)]})' for c in cohorts],
                        fontsize=10, color='#2c3e50')
    ax1.set_ylim(0, 115)
    ax1.set_ylabel('Patients (%)', fontsize=11, fontweight='bold', color='#2c3e50')
    ax1.set_title('Best Overall Response by Cohort', fontsize=12,
                  fontweight='bold', color='#1a1a2e', pad=14, loc='left')

    legend_elements = [mpatches.Patch(facecolor=RESP_COLORS[r], edgecolor='white', label=r)
                       for r in RESP_ORDER]
    ax1.legend(handles=legend_elements, loc='upper left',
               bbox_to_anchor=(1.02, 1.0), bbox_transform=ax1.transAxes,
               fontsize=9, frameon=True, fancybox=False,
               edgecolor='#cccccc', title='Response', title_fontsize=9)

    # ── Panel 2: grouped bar — exclude NE if all zero ──
    active_resp = [r for r in RESP_ORDER
                   if any(calc_pcts(adsl[adsl['COHORT'] == c])[0][r] > 0 for c in cohorts)]

    n_resp     = len(active_resp)
    group_w    = 0.70
    bar_w2     = group_w / len(cohorts)
    coh_colors = ['#e74c3c', '#2ecc71', '#3498db', '#9b59b6']

    for ci, c in enumerate(cohorts):
        pcts, _ = calc_pcts(adsl[adsl['COHORT'] == c])
        vals    = [pcts[r] for r in active_resp]
        xi      = np.arange(n_resp) + ci * bar_w2 - group_w / 2 + bar_w2 / 2
        ax2.bar(xi, vals, bar_w2,
                color=coh_colors[ci], edgecolor='white',
                linewidth=1.2, zorder=3, alpha=0.88,
                label=f'Cohort {int(c)} ({DOSES[int(c)]})')

    ax2.set_xticks(range(n_resp))
    ax2.set_xticklabels(active_resp, fontsize=11, fontweight='bold', color='#2c3e50')
    ax2.set_ylabel('Patients (%)', fontsize=11, fontweight='bold', color='#2c3e50')
    ax2.set_title('Response Distribution Across Cohorts', fontsize=12,
                  fontweight='bold', color='#1a1a2e', pad=14, loc='left')
    ax2.set_ylim(0, 75)
    ax2.legend(fontsize=9, frameon=True, fancybox=False,
               edgecolor='#cccccc', loc='upper right')

    # ── Summary box: bottom center, outside plots ──
    n_orr = adsl['BESTRSPC'].isin(['CR', 'PR']).sum()
    n_cr  = (adsl['BESTRSPC'] == 'CR').sum()
    n_pr  = (adsl['BESTRSPC'] == 'PR').sum()
    n_sd  = (adsl['BESTRSPC'] == 'SD').sum()
    n_pd  = (adsl['BESTRSPC'] == 'PD').sum()

    summary = (f"  Total  {total}    "
               f"ORR  {n_orr} ({n_orr/total*100:.0f}%)    "
               f"CR  {n_cr} ({n_cr/total*100:.0f}%)    "
               f"PR  {n_pr} ({n_pr/total*100:.0f}%)    "
               f"SD  {n_sd} ({n_sd/total*100:.0f}%)    "
               f"PD  {n_pd} ({n_pd/total*100:.0f}%)  ")

    fig.text(0.50, 0.08, summary,
             ha='center', va='center',
             fontsize=10, fontweight='bold',
             color='white', family='monospace',
             bbox=dict(boxstyle='round,pad=0.7',
                       facecolor='#2c3e50',
                       edgecolor='#1a1a2e',
                       linewidth=2))

    fig.suptitle(f'Best Overall Response — Synthetic Oncology Phase I (N={total})',
                 fontsize=13, fontweight='bold', color='#1a1a2e', y=0.97)

    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'Saved: {filename}')


def main():
    adsl = pd.read_csv(ADSL_PATH)
    draw_bor(adsl, 'bor_chart.png')


if __name__ == '__main__':
    main()
