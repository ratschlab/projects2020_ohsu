import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from config import *
from matplotlib.ticker import ScalarFormatter
import numpy as np

def plot_venn(count1, count2, count_int, isDownload, title = 'Comparison of unique kmers counts between whole data sets',nameDownload='picture', isPNG=True, isPDF = True):
    plt.figure(figsize=(50,50))
    total = count1 + count2 - count_int
    out = venn2(
        subsets = (count1 - count_int, count2 - count_int, count_int),
        set_labels = ('BRCA(ETH)', 'BRCA(OHSU)'),
        set_colors=[ETH_COLOR,OHSU_COLOR]
    )
    for text in out.set_labels:
        text.set_fontsize(TEXT_SIZE)
    for text in out.subset_labels:
        text.set_fontsize(TEXT_SIZE)
    plt.title(title, fontsize=TEXT_SIZE)
    #save_plot(title)
    if isDownload:
        if isPNG:
            plt.savefig(nameDownload+'.png', bbox_inches='tight',dpi=300)
        if isPDF:
            plt.savefig(nameDownload+'.pdf', bbox_inches='tight',dpi=300)
    plt.show()
    return plt


def plot_venn_perc(count1, count2, count_int, isDownload, title = 'Comparison of unique kmers counts between whole data sets',nameDownload='picture', isPNG=True, isPDF = True):
    plt.figure(figsize=(50,50))
    total = count1 + count2 - count_int
    out = venn2(
        subsets = (count1 - count_int, count2 - count_int, count_int),
        set_labels = ('BRCA(ETH)', 'BRCA(OHSU)'),
        set_colors=[ETH_COLOR,OHSU_COLOR],
        subset_label_formatter=lambda x: f"{(x/total):0.01%}"
    )
    for text in out.set_labels:
        text.set_fontsize(TEXT_SIZE)
    for text in out.subset_labels:
        text.set_fontsize(TEXT_SIZE)
    plt.title(title, fontsize=TEXT_SIZE)
    if isDownload:
        if isPNG:
            plt.savefig(nameDownload+'.png', bbox_inches='tight',dpi=300)
        if isPDF:
            plt.savefig(nameDownload+'.pdf', bbox_inches='tight',dpi=300)
    plt.show()
    return plt

def barplot_column_common_unique(count1, count2, count_inter,  isDownload, ylabel = 'Unique kmers counts per all genes',
                                title = 'Comparison of unique kmers counts between BRCA(ETH) and BRCA(OHSU) whole data sets',nameDownload='picture', isPNG=True, isPDF = True):
    data = [count1 - count_inter, count2 - count_inter, count_inter]
    labels = ['ETH\OHSU', 'OHSU\ETH',
              'OHSU & ETH']
    plt.bar(labels, data, color=[ETH_COLOR, OHSU_COLOR, OHSU_ETH_COLOR], alpha=0.5)
    # Add labels to the bars
    for i, v in enumerate(data):
        plt.text(i, v, str(v), ha='center', va='bottom')
    plt.ylabel(ylabel)
    plt.title(title)
    if isDownload:
        if isPNG:
            plt.savefig(nameDownload+'.png', bbox_inches='tight',dpi=300)
        if isPDF:
            plt.savefig(nameDownload+'.pdf', bbox_inches='tight',dpi=300)
    plt.show()


def plotting_filtering_barplot(fffb,v,name,bottom,bar_position,sample,filfor,x_secondary,salt,path,path_pdf,coor=True,ySalt=''):
    fig, ax = plt.subplots(nrows=1,ncols=1)
    bw=0.25
    p1 = ax.bar(np.arange(len(fffb))-bw,v[0],bw,label=name[0],bottom=bottom, edgecolor='white',color=COLORS_COOR['size_ohsu\eth_coor'] if coor else COLORS_COOR['size_ohsu\eth'],zorder=2)
    p1 = ax.bar(np.arange(len(fffb)),v[1],bw,label=name[1],bottom=bottom, edgecolor='white',color=COLORS_COOR['size_intersection_coor'] if coor else COLORS_COOR['size_intersection'] ,zorder=3)
    p1 = ax.bar(np.arange(len(fffb))+bw,v[2],bw,label=name[2],bottom=bottom, edgecolor='white',color=COLORS_COOR['size_eth\ohsu_coor'] if coor else COLORS_COOR['size_eth\ohsu']  ,zorder=4)
        
    SALT = f"{PLOT_SORT_BY}"
    plt.suptitle('-'.join([sample[0:4], sample[4:6], sample[6:10], sample[10:13], sample[13:16], sample[16:20], sample[20:22]])+'\n'+salt,size=15,x=0.35,y=1.05)
    ax.set_xticks([pos+bw/2 for pos in bar_position])
    ax.set_xticklabels(filfor,rotation=90,ha='center', fontsize=8)
    ax.tick_params(labelsize=8)
    ax.set_xlabel('Filter foreground',size=10)
    type=' junctions (coordinates)' if coor else ' kmers'
    ax.set_ylabel(f'Number of{type}{ySalt}',size=10)
    ax.legend(fontsize=8,loc='center left',bbox_to_anchor=(1,0.5))
    ax.grid(axis='y',zorder=0)
    # Add second axis
    ax_sec = ax.secondary_xaxis('top')
    ax_sec.set_xticks([pos+bw/2 for pos in bar_position])
    ax_sec.set_xticklabels(x_secondary,rotation=90,ha='center',fontsize=8)
    ax_sec.set_xlabel('Filter background',size=10)
    # Save plot
    
    plt.tight_layout()
    
    plt.savefig(path, dpi=300,bbox_inches='tight')
    plt.savefig(path_pdf, dpi=300,bbox_inches='tight')
    plt.show()


def plotting_nf_barplot(axis,fffb2,bottoml,sample,path_sample_nf,path_sample_nf_pdf,salt='',coor=True):
    fig, ax = plt.subplots(nrows=1,ncols=1)
    for data, data_count in axis.items():
            p = ax.bar(fffb2,data_count,width=0.1,label=data,bottom=bottoml, edgecolor='white',color=COLORS[data],zorder=3)
            # Use to add number of k-mer to center bars (or only intersection)
            ax.bar_label(p,label_type='center', size =10,zorder=6)
            bottoml+=data_count
        
    y_formatter = ScalarFormatter()
    y_formatter.set_scientific(False)  # Отключение научной нотации
    plt.gca().yaxis.set_major_formatter(y_formatter)
    
    plt.suptitle('-'.join([sample[0:4], sample[4:6], sample[6:10], sample[10:13], sample[13:16], sample[16:20], sample[20:22]]),size=15,x=0.4,y=1.05)
    type=' junctions (coordinates)' if coor else ' kmers'
  
    ax.set_title('Non-filtered',size=10)
    ax.set_ylabel(f'Number of{type}{salt}',size=10)
    ax.legend(fontsize=10,loc='center left',bbox_to_anchor=(1,0.5))
    ax.set_xlim(-0.08,0.08)
    ax.grid(axis='y',zorder=0)
    ax.tick_params(labelsize=10)
    
    plt.tight_layout()
    plt.savefig(path_sample_nf, dpi=300,bbox_inches='tight')
    plt.savefig(path_sample_nf_pdf, dpi=300,bbox_inches='tight')
    plt.show()