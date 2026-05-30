import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Patch
import sys

ADSL_PATH = 'ADSL.csv'
ADRS_PATH = 'ADRS.csv'

COLOR_TRT     = '#1f4788'
COLOR_OFF     = '#808080'
RESP_COLORS   = {'CR': '#00CC00', 'PR': '#5a9bd4', 'SD': '#f0c929', 'PD': '#c41e3a'}
RESP_MARKERS  = {'CR': 'o', 'PR': 's', 'SD': '^', 'PD': 'D'}
RESP_SIZES    = {'CR': 9, 'PR': 8, 'SD': 9, 'PD': 8}
MIN_BAR       = 0.3
SYMBOL_OFFSET = 1.2


def months(start, end):
    return ((end - start).days + 1) / 30.4


def build_df(adsl, adrs, cutoff):
    rows = []
    for _, s in adsl.iterrows():
        uid      = s['USUBJID']
        trtsdt   = s['TRTSDT']
        trtedt   = s['TRTEDT']
        still_on = s['EOSSTT'] == 'ONGOING'
        trt_mo   = months(trtsdt, cutoff if still_on else trtedt)

        off_end    = None
        eosdtc_val = str(s.get('EOSDTC', '')).strip()
        ltfu_val   = str(s.get('LTFU', 'N')).strip().upper()
        if not still_on:
            if eosdtc_val not in ['', 'nan', 'NaT']:
                eos_mo = months(trtsdt, pd.to_datetime(eosdtc_val))
                if eos_mo > trt_mo:
                    off_end = eos_mo
            elif ltfu_val == 'Y':
                cutoff_mo = months(trtsdt, cutoff)
                if cutoff_mo > trt_mo:
                    off_end = cutoff_mo

        pt_rs = adrs[(adrs['USUBJID'] == uid) & (adrs['PARAMCD'] == 'OVRLRESP')]
        resps = [{'months': months(trtsdt, r['ADT']), 'response': r['AVALC']}
                 for _, r in pt_rs.iterrows()]

        rows.append({
            'usubjid':    uid,
            'cohort':     s.get('COHORT', ''),
            'trt_mo':     trt_mo,
            'still_on':   still_on,
            'dcsreas':    str(s.get('DCSREAS', '')).strip(),
            'liver_mets': s.get('LIVERMETS', 'N') == 'Y',
            'best_resp':  s.get('BESTRSPC', 'NE'),
            'responses':  resps,
            'off_end':    off_end,
        })

    df = pd.DataFrame(rows).sort_values('trt_mo', ascending=True).reset_index(drop=True)
    df['y'] = range(len(df) - 1, -1, -1)
    return df


def draw_swimmer(df, title, filename):
    n   = len(df)
    fig, ax = plt.subplots(figsize=(16, max(8, n * 0.42)))
    fig.subplots_adjust(right=0.72)

    max_mo = df['trt_mo'].max()
    if df['off_end'].notna().any():
        max_mo = max(max_mo, df['off_end'].dropna().max())
    max_mo += 3

    ax.set_xlim(-5, max_mo)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_xlabel('Months from Treatment Start', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=13, fontweight='bold', pad=14)
    ax.set_xticks(range(0, int(np.ceil(max_mo)) + 1, 3))
    ax.set_yticks(df['y'])
    ax.set_yticklabels(df['usubjid'], fontsize=7, fontweight='bold')
    ax.grid(axis='x', color='gray', linestyle='--', linewidth=0.5, alpha=0.3)

    for _, r in df.iterrows():
        ax.add_patch(Rectangle((0, r['y'] - 0.3), max(r['trt_mo'], MIN_BAR), 0.6,
                                facecolor=COLOR_TRT, edgecolor='none'))
        if pd.notna(r['off_end']) and not r['still_on']:
            ax.add_patch(Rectangle((r['trt_mo'], r['y'] - 0.3),
                                    r['off_end'] - r['trt_mo'], 0.6,
                                    facecolor=COLOR_OFF, edgecolor='none'))

    for _, r in df.iterrows():
        for resp in r['responses']:
            rc = resp['response']
            if rc in RESP_COLORS:
                ax.plot(resp['months'], r['y'],
                        marker=RESP_MARKERS[rc], markersize=RESP_SIZES[rc],
                        markerfacecolor=RESP_COLORS[rc], markeredgecolor=RESP_COLORS[rc],
                        markeredgewidth=1.2, zorder=4)

    for _, r in df.iterrows():
        candidates = [r['trt_mo']]
        if pd.notna(r['off_end']): candidates.append(r['off_end'])
        mo_list = [rr['months'] for rr in r['responses']]
        if mo_list: candidates.append(max(mo_list))
        sx  = max(candidates) + SYMBOL_OFFSET
        dcs = r['dcsreas']

        if r['still_on']:
            ax.plot(sx, r['y'], marker='>', markersize=14, color=COLOR_TRT,
                    markerfacecolor=COLOR_TRT, markeredgewidth=0, zorder=5)
        elif dcs == 'Death':
            ax.plot(sx, r['y'], marker='x', markersize=10, markeredgewidth=2.5,
                    color='black', linestyle='None', zorder=6)
        elif dcs in ['Withdrawal of Consent', 'Lost to Follow-up']:
            ax.plot(sx, r['y'], marker='x', markersize=10, markeredgewidth=2.5,
                    color='#c41e3a', linestyle='None', zorder=6)
        elif dcs != '':
            ax.plot(sx, r['y'], marker='x', markersize=10, markeredgewidth=2.5,
                    color='black', linestyle='None', zorder=6)

    ax.text(-3.5, n + 0.3, 'Best\nResp',  ha='center', va='bottom', fontsize=8, fontweight='bold')
    ax.text(-1.5, n + 0.3, 'Liver\nMets', ha='center', va='bottom', fontsize=8, fontweight='bold')
    for _, r in df.iterrows():
        ax.text(-3.5, r['y'], r['best_resp'], ha='center', va='center',
                fontsize=7, fontweight='bold', color=RESP_COLORS.get(r['best_resp'], 'gray'))
        ax.text(-1.5, r['y'], 'Y' if r['liver_mets'] else 'N',
                ha='center', va='center', fontsize=7)

    legend_elements = [
        Patch(color=COLOR_TRT, label='On Treatment'),
        Patch(color=COLOR_OFF, label='Off Treatment / On Study Follow-up'),
        plt.Line2D([0],[0], marker='>', color='w', markerfacecolor=COLOR_TRT,
                   markersize=12, label='Still on Treatment', linestyle='None'),
        plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=RESP_COLORS['CR'],
                   markersize=9, label='CR', linestyle='None'),
        plt.Line2D([0],[0], marker='s', color='w', markerfacecolor=RESP_COLORS['PR'],
                   markersize=8, label='PR', linestyle='None'),
        plt.Line2D([0],[0], marker='^', color='w', markerfacecolor=RESP_COLORS['SD'],
                   markersize=9, label='SD', linestyle='None'),
        plt.Line2D([0],[0], marker='D', color='w', markerfacecolor=RESP_COLORS['PD'],
                   markersize=8, label='PD', linestyle='None'),
        plt.Line2D([0],[0], marker='x', color='black', markersize=10,
                   markeredgewidth=2.5, label='Death / Early Termination', linestyle='None'),
        plt.Line2D([0],[0], marker='x', color='#c41e3a', markersize=10,
                   markeredgewidth=2.5, label='Withdrawal / Lost to Follow-up', linestyle='None'),
    ]
    legend = ax.legend(handles=legend_elements, loc='upper left',
                       bbox_to_anchor=(1.02, 1.0), bbox_transform=ax.transAxes,
                       fontsize=9, title='Legend', title_fontsize=10,
                       frameon=True, fancybox=False, borderpad=1, labelspacing=1.0)
    legend.get_frame().set_linewidth(1.5)
    legend.get_frame().set_edgecolor('black')

    n_ongoing = int(df['still_on'].sum())
    summary   = f"Summary\nOngoing   {n_ongoing:2d}\nOff Trt   {n - n_ongoing:2d}\nTotal     {n:2d}"
    ax.annotate(summary, xy=(1.02, 0.0), xycoords='axes fraction',
                fontsize=10, fontweight='bold', color='white', family='monospace',
                va='bottom', ha='left',
                bbox=dict(boxstyle='round,pad=0.6', facecolor='#2b5f8a',
                          edgecolor='black', linewidth=2))

    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {filename}')


def main():
    adsl = pd.read_csv(ADSL_PATH, parse_dates=['TRTSDT', 'TRTEDT', 'CUTDTC'])
    adrs = pd.read_csv(ADRS_PATH, parse_dates=['ADT'])
    cutoff = adsl['CUTDTC'].max()

    df_all = build_df(adsl, adrs, cutoff)
    draw_swimmer(df_all,
                 f'Swimmer Plot — Synthetic Oncology Phase I (All Cohorts, N={len(df_all)})',
                 'swimmer_plot_all.png')

    for cohort in sorted(adsl['COHORT'].dropna().unique()):
        subset = adsl[adsl['COHORT'] == cohort].reset_index(drop=True)
        df_c   = build_df(subset, adrs, cutoff)
        draw_swimmer(df_c,
                     f'Swimmer Plot — Synthetic Oncology Phase I (Cohort {int(cohort)}, N={len(df_c)})',
                     f'swimmer_plot_cohort{int(cohort)}.png')


if __name__ == '__main__':
    main()
