import pandas as pd
import glob
import os 
import numpy as np 
from collections import defaultdict




def reader_assign_conf_pep(path, FDR_threshold, col_seq, col_qval):
    print(f'Reading {path}')
    if os.path.isfile(path):
        df = pd.read_csv(path, sep = '\t')
        tot_peptides = len(df[col_seq].unique())
        print(f'With Shape: {df.shape[0]}')
        print(f'With unique peptides: {tot_peptides}')
        assert('sequence' in df.columns)
        df_filtered = df.loc[df[col_qval] < FDR_threshold]
        print(f'Number of validated psm: {df_filtered.shape}')
        peptides = set(df_filtered[col_seq])
        val = len(peptides)
        if tot_peptides:
            val_rate = np.round(val / tot_peptides * 100 , 2)
        else:
            val_rate = 0.0
        print(f'Number of validated unique peptides: {val}')
        print(f'Validation Rate: {val_rate } percent')
        return val, val_rate, peptides, df_filtered
    else:
        return 0, 0.0, set(), None
    
    
def plot_text(Y, T, position='top', color='black', font=None):
    if max(Y) > 0:
        font['color'] = color
        Y = np.array(Y)
        T = np.array(T)
        Y[np.where(Y > 0 )[0]]
        T[np.where(Y > 0 )[0]]
        change_val = [i for i in np.arange(len(Y) - 1) if Y[i] != Y[i - 1]]    
        weighted = [change_val[i] + (change_val[i+1] - change_val[i]) / 2 for i, x in enumerate(change_val[:-1])]
        X = [np.floor(change_val[i] + (change_val[i+1] - change_val[i]) / 2) for i, x in enumerate(change_val[:-1])]
        Y = Y[np.array(change_val[:-1])]
        T = T[np.array(change_val[:-1])]
        p_prev = 0 
        percent_diff = 20
        min_x = min(X)

        previous_plot = 0 
        for x, y, p in zip(X, Y, T):
            if position == 'bottom':
                delta = - (y/5.5)
            elif position == 'top': 
                delta = (y/12)
            if (p > p_prev + (p_prev/percent_diff)) or (p < p_prev - (p_prev/percent_diff)):
                previous_plot += 1
            if previous_plot == 2 or x == min_x: # delay the plotting
                if y != 0 : # Because log scale
                    plt.text(x - 0.5 , y + delta , p, ha='left', **font)
                    previous_plot = 0

            p_prev = p 

        
def plot_text_all(X, Y, T):
    for x, y, p in zip(X, Y, T):
        plt.text(x, y, p)
        #plt.text(x - 0.5 , y + (y/10), p)
        
        
def plot_intersection_bars(back_ticks, front_ticks, ticks_fontsize, axislabels_fontsize, 
                          legend_fontsize, axes_fontsize, run_type, 
                           serie_index, serie_intersection, serie_eth, serie_ohsu,
                           y_label, save, plot_dir, base_plot, name_plot):

    text_font = {'size':'12', 'weight':'normal'}

    xupper_axis_label = 'GTEX (reads, samples)'
    xlower_axis_label = f'{run_type.upper()} (min, reads, samples)'


    colorgrid = 'grey'
    alpha_grid = 0.3
    marker_size = 10
    mew = 4
    color1 = 'gold'
    color2 = 'crimson'
    color3 = 'mediumblue'
    color4 = 'tomato'
    colorgrid = 'grey'
    width = 0.4

    fig, ax1 = plt.subplots(figsize=(15, 6))
    ax2 =  ax1.secondary_xaxis('top')   
    plt.grid(b=True, axis = 'both', which='major', color=colorgrid, linestyle='-', alpha=alpha_grid)
    plt.grid(b=False, axis = 'both', which='minor', color=colorgrid, linestyle='--', alpha=alpha_grid)


    plt.bar(serie_index, serie_intersection, width=width, 
            color=color1, label='Intersection size')
    plt.plot(serie_index, serie_eth, alpha=1, color=color3,
             linestyle = 'None', markerfacecolor='None', marker="_", markersize=marker_size, markeredgewidth=mew,
             label = 'Total set size Graph Pipeline')
    plt.plot(serie_index, serie_ohsu, alpha=1, color=color2,
             linestyle = 'None', markerfacecolor='None', marker="_", markersize=marker_size, markeredgewidth=mew,
             label = 'Total set size Junction Pipeline')

    plot_text(serie_ohsu, serie_ohsu, 'top', color=color2, font=text_font)
    plot_text(serie_eth, serie_eth, 'top', color=color3, font=text_font)
    plot_text(serie_intersection, serie_intersection, color=color4, font=text_font)

    #plt.yscale('log')
    max_scale = np.max([serie_ohsu.values, serie_eth.values])

    ax1.set_xticks(serie_index, 
               labels = front_ticks,
               rotation = 90, 
               ha = 'center', 
               fontsize=ticks_fontsize)

    ax2.set_xticks(serie_index, 
               labels = back_ticks,
               rotation = 90, 
               ha = 'center', 
               fontsize=ticks_fontsize)


    plt.legend(fontsize=legend_fontsize)
    plt.ylabel(y_label, fontsize=axes_fontsize)
    ax2.set_xlabel(xupper_axis_label, fontsize=axes_fontsize)
    ax1.set_xlabel(xlower_axis_label, fontsize=axes_fontsize)

    save_path = os.path.join(plot_dir, f'{base_plot}_{name_plot}.pdf')
    print("saving path is to {}".format(save_path))
    if save:
        print("Saving!")
        plt.savefig(save_path, bbox_inches='tight')

    plt.show()
