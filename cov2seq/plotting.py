import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatch

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

font = {'size'   : 12}
matplotlib.rc('font', **font)

def avgcov(cov, window):
    avgs =  np.nanmean(np.pad(cov.astype(float), 
                              (0, window - cov.size%window), 
                              mode='constant', constant_values=np.NaN
                             ).reshape(-1, window), axis=1)
    centers = np.array([i*window + 0.5*window for i in range(len(avgs))])
    return centers, avgs

def make_regions_visible(regions, visibility_threshold=20):
    '''converts regions to non-overlapping, 0-based regions with a minimal width to ensure their visibility when plotted'''
    regions['start'] -= 1 # make 0-based
    regions.loc[regions.width < visibility_threshold, 'start'] += (regions.loc[regions.width < visibility_threshold, 'width'] - visibility_threshold).abs() // 2
    regions.loc[regions.width < visibility_threshold, 'end'] += (regions.loc[regions.width < visibility_threshold, 'width'] - visibility_threshold).abs() // 2
    no_overlaps = []
    for i,(start,end,width) in regions.iterrows():
        if no_overlaps:
            prev_start,prev_end,prev_width = no_overlaps[-1]
            if (start - prev_end) < visibility_threshold:
                no_overlaps[-1] = (prev_start, end, prev_width + width)
                continue
        no_overlaps.append((start,end,width))
    return pd.DataFrame(no_overlaps, columns=['start', 'end', 'width'])

def create_summary_plot(sample_schemes, cov_primertrimmed, cov_illumina, cov_sanger, cov_pools, cov_limit_regions, cov_low_regions,
                        snv_info, reference, reference_genes, amplicons, final, window=20, snv_annotation_layers=None, 
                        snv_annotation_offset=0, savepath=None):
    if snv_annotation_layers is None:
        snv_count = len(snv_info.loc[pd.notnull(snv_info[('longshot', 'cov')]) | (snv_info[('final', 'decision')] == 'confirmed')].index.drop_duplicates())
        snv_annotation_layers = int(np.round(snv_count/3 + 0.5))
    fig = plt.figure(figsize=(15, 8+len(sample_schemes)*0.2))
    gs = gridspec.GridSpec(4,1,hspace=0.0,height_ratios=[1,0.1*len(sample_schemes),0.05,0.04*snv_annotation_layers])
    ax1 = plt.Subplot(fig, gs[0])
    ax2 = plt.Subplot(fig, gs[1])
    ax3 = plt.Subplot(fig, gs[2])
    ax4 = plt.Subplot(fig, gs[3])
    ax1.set_xlim(0, len(reference.seq))
    ax1_ylim = (0,max(cov_primertrimmed)*1.1)
    ax1.set_ylim(ax1_ylim)
    ax1.yaxis.grid(alpha=0.2)
    ax1.xaxis.tick_top()
    ax1.set_xticks(np.arange(0,len(reference.seq),5000))
    ax1.set_xticks(np.arange(0,len(reference.seq),1000), minor=True)
    pool_colors = {}
    last_color = 1
    for scheme in cov_pools:
        pool_colors[scheme] = {}
        for pool in cov_pools[scheme]:
            pool_colors[scheme][pool] = "C{}".format(last_color)
            last_color += 1
    # coverage plots
    if np.any(cov_primertrimmed):
        X, cov_primertrimmed_vals = avgcov(cov_primertrimmed, window)
        ax1.plot(X, cov_primertrimmed_vals, linewidth=0.5, color="black", label="total Nanopore coverage")
    if np.any(cov_sanger):
        X, cov_sanger_vals = avgcov(cov_sanger, window)
        ax1.plot(X, cov_sanger_vals*100, linewidth=0.5, color="green", linestyle='dotted', label="approx. Sanger coverage (factor 100)")
    if np.any(cov_illumina):
        X, cov_illumina_vals = avgcov(cov_illumina, window)
        ax1.plot(X, cov_illumina_vals, linewidth=0.5, color="blue", linestyle='dashed', label="Illumina coverage")
    for scheme in cov_pools:
        for pool in cov_pools[scheme]:
            X, cov_pool_vals = avgcov(cov_pools[scheme][pool], window)
            ax1.plot(X, cov_pool_vals, linewidth=0.3, color=pool_colors[scheme][pool], label=pool)
    ax1.legend(loc="upper right", bbox_to_anchor=(1.0,0.95))
    ax1.set_ylabel('coverage ({}nt window)'.format(window))
    # highlight low coverage regions
    for i,(start,end,width) in make_regions_visible(cov_low_regions).iterrows():
        ax1.axvspan(start, end, facecolor='orange', alpha=0.2, linewidth=None)
    for i,(start,end,width) in make_regions_visible(cov_limit_regions).iterrows():
        if start == 0:
            halign = 'left'
        elif end == len(reference.seq):
            halign = 'right'
        else:
            halign = 'center'
        ax1.axvspan(start, end, facecolor='red', alpha=0.2, linewidth=None)
        ax1.text(np.mean([start, end]), ax1_ylim[1]*0.99, "{}".format(width), 
                 horizontalalignment=halign, verticalalignment="top")
    # plot amplicons
    # 'name', pool', pstart', 'start', 'end', 'pend', 'mingc'
    ax2.set_yticks([])
    pool_loc = {}
    for scheme in sample_schemes:
        for i,amplicon in amplicons[scheme].iterrows():
            pool = amplicon['pool']
            if pool not in pool_loc:
                pool_loc[pool] = len(pool_loc)
            rl = mpatch.Rectangle((amplicon['pstart'],pool_loc[pool]), 
                                  amplicon['start']-amplicon['pstart'], 1, 
                                  facecolor="grey",
                                  edgecolor=None)
            r  = mpatch.Rectangle((amplicon['start'],pool_loc[pool]), 
                                  amplicon['end']-amplicon['start'], 1, 
                                  #facecolor="C{}".format(pool_loc[pool]+1),
                                  facecolor=pool_colors[scheme][pool],
                                  edgecolor=None)
            rr = mpatch.Rectangle((amplicon['end'],pool_loc[pool]), 
                                  amplicon['pend']-amplicon['end'], 1, 
                                  facecolor="grey",
                                  edgecolor=None)
            ax2.annotate(amplicon['name'].split("_")[-1], 
                         (amplicon['start']+(amplicon['end']-amplicon['start'])/2, pool_loc[pool]+0.5), 
                         color='w', weight='bold', 
                         fontsize=10, ha='center', va='center')
            #rr = mpatch.Rectangle((2,2), 8, 2)
            ax2.add_patch(r)
            ax2.add_patch(rl)
            ax2.add_patch(rr)
    ax2.set_ylim((0,len(pool_loc)))
    
    # plot genes
    ax3.set_yticks([])
    for i,(gene, strand, start, end, protein_id, product, note) in reference_genes.iterrows():
        r  = mpatch.Rectangle((start,0), 
                               end-start, 1, 
                               facecolor="C{}".format(3 + i%2),
                               edgecolor=None)
        ax3.annotate(gene, (start+(end-start)/2, 0.5), color='w', weight='bold', 
                     fontsize=10, ha='center', va='center')
        ax3.add_patch(r)
    ax3.set_ylim((0,1))
    # show SNVs
    i = snv_annotation_offset
    prev_index = None
    for index, v in snv_info.iterrows():
        if index == prev_index:
            continue
        prev_index = index
        site = v[('medaka variant', 'site')]
        alt = v[('medaka variant', 'alt')]
        ref = v[('medaka variant', 'ref')]
        scov = v[('longshot', 'cov')]
        passed = v[('ARTIC', 'snv_filter')]
        decision = v[('final', 'decision')]
        clades_ = v[('nextstrain', 'clades')]

        #clades_ = v[('nextstrain', 'clades')]
        if np.isnan(scov) and decision != 'confirmed':
            continue
        height = snv_annotation_layers - 1 - i % snv_annotation_layers
        color = "green" if passed == True else "red"
        ax4.plot([site, site], [snv_annotation_layers+1,height+1], color=color, linewidth=1.)
        txtcolor = "black"
        if decision == True:
            txtcolor = 'green'
        elif clades_:
            #txtcolor = clades_df.loc[(clades_df['site'] == site) & (clades_df['alt'] == alt),
            #                         'color'].to_list()[0]
            txtcolor = "blue"
        ax4.text(site, height, index, 
                 horizontalalignment='left', verticalalignment="bottom", color=txtcolor)
        i += 1

    ax4.set_ylim((0,snv_annotation_layers+1))
    ax4.set_yticks([])
    ax4.spines['bottom'].set_visible(False)
    ax4.spines['left'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlim(0, len(reference.seq))
        if ax != ax1:
            ax.set_xticks([])
        fig.add_subplot(ax)
    
    if savepath:
        plt.savefig(savepath)
    else:
        plt.show()
    return