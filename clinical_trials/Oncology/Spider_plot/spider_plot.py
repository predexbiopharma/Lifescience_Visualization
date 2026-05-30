import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ADSL_PATH = 'ADSL.csv'
ADTR_PATH = 'ADTR.csv'


def build_df(adsl, adtr):
    rows = []
    for _, s in adsl.iterrows():
        uid = s['USUBJID']

        bl           = adtr[(adtr['USUBJID'] == uid) & (adtr['AVISIT'] == 'BASELINE')]
        baseline_sum = pd.to_numeric(bl['AVAL'].iloc[0], errors='coerce') if len(bl) > 0 else np.nan
        baseline_dt  = bl['ADT'].iloc[0] if len(bl) > 0 else pd.NaT

        if pd.isna(baseline_sum) or baseline_sum == 0 or pd.isna(baseline_dt):
            continue

        rows.append({'usubjid': uid, 'cohort': s.get('COHORT', ''), 'months': 0.0, 'pct_change': 0.0})

        post = adtr[(adtr['USUBJID'] == uid) & (adtr['AVISIT'] != 'BASELINE')]
        for _, r in post.iterrows():
            current_sum = pd.to_numeric(r['AVAL'], errors='coerce')
            current_dt  = r['ADT']
            if pd.isna(current_sum) or pd.isna(current_dt):
                continue
            mo  = ((current_dt - baseline_dt).days + 1) / 30.4
            pct = ((current_sum - baseline_sum) / baseline_sum) * 100
            rows.append({'usubjid': uid, 'cohort': s.get('COHORT', ''), 'months': mo, 'pct_change': pct})

    return pd.DataFrame(rows)


def draw_spider(data, title, filename):
    pids  = data['usubjid'].unique()
    y_top = max(100, min(data['pct_change'].max(), 200))
    x_max = int(np.ceil(data['months'].max() / 3) * 3)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, x_max)
    ax.set_ylim(-100, y_top)
    ax.set_xlabel('Months from Treatment Start', fontsize=11, fontweight='bold')
    ax.set_ylabel('Change from baseline in target lesion (%)', fontsize=11, fontweight='bold')
    ax.set_title(title, fontsize=12, fontweight='bold', pad=14)

    y_ticks = list(range(-100, int(y_top) + 1, 20))
    if int(y_top) not in y_ticks:
        y_ticks.append(int(y_top))
    ax.set_xticks(range(0, x_max + 1, 3))
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f'{y}%' for y in y_ticks])

    ax.axhline(y=0,   color='black', linewidth=1.5, zorder=3)
    ax.axhline(y=-30, color='black', linewidth=1.5, linestyle='--', zorder=3)
    ax.axhline(y=20,  color='black', linewidth=1.5, linestyle='--', zorder=3)
    ax.text(x_max * 0.98,  22, '+20%: PD threshold', ha='right', va='bottom', fontsize=9)
    ax.text(x_max * 0.98, -32, '-30%: PR threshold', ha='right', va='top',   fontsize=9)
    ax.grid(color='lightgray', linestyle='-', linewidth=0.5, alpha=0.3)

    for pid in pids:
        pt = data[data['usubjid'] == pid].sort_values('months').copy()
        pt['pct_change'] = pt['pct_change'].clip(upper=y_top)
        ax.plot(pt['months'], pt['pct_change'], color='black', linewidth=1.2)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {filename}')


def main():
    adsl = pd.read_csv(ADSL_PATH, parse_dates=['TRTSDT', 'TRTEDT', 'CUTDTC'])
    adtr = pd.read_csv(ADTR_PATH, parse_dates=['ADT'])

    df = build_df(adsl, adtr)
    print(f'Spider patients: {df["usubjid"].nunique()}')

    draw_spider(df,
                f'Spider Plot — All Cohorts (N={df["usubjid"].nunique()})',
                'spider_all.png')

    for cohort in sorted(df['cohort'].dropna().unique()):
        sub = df[df['cohort'] == cohort]
        draw_spider(sub,
                    f'Spider Plot — Cohort {int(cohort)} (N={sub["usubjid"].nunique()})',
                    f'spider_cohort{int(cohort)}.png')


if __name__ == '__main__':
    main()
