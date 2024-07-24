import pandas as pd
import glob
import os 
import numpy as np 
from collections import defaultdict
import matplotlib.pyplot as plt
import math

def reader_assign_conf_pep(path, FDR_threshold, col_seq, col_qval, input_trypPep):
    print(f'Reading {path}')
    if os.path.isfile(path) and os.path.isfile(input_trypPep):
        # Read
        df = pd.read_csv(path, sep = '\t')
        df_all_pep = pd.read_csv(input_trypPep, sep = '\t')
        tryptic_peptides = set(df_all_pep['sequence'].unique())
        
        # Total
        tot_peptides = len(df_all_pep['sequence'].unique())
        
        print(f'Shape of FDR file: {df.shape[0]}')
        print(f'Total input tryptic junction peptides: {tot_peptides}')
        assert('sequence' in df.columns)
        
        # Validated peptides
        df_filtered = df.loc[df[col_qval] < FDR_threshold]
        peptides = set(df_filtered[col_seq])
        val = len(peptides)

        # Validation rate
        if tot_peptides:
            val_rate = np.round(val / tot_peptides * 100 , 2)
        else:
            val_rate = 0.0
        
        print(f'Number of validated unique peptides: {val}')
        print(f'Validation Rate: {val_rate } percent')
        return val, val_rate, peptides, df_filtered, tryptic_peptides
    else:
        return 0, 0.0, set(), None, set()
    
    
def sort_filters(df, order_background, order_foreground):
    def prepare_backticks(a, b):
        return  [f'({i}, {j})' for i, j in zip(a.values, b.values)]
    def prepare_frontticks(a, b, c):
        return  [f'({i}, {j}, {k})' for i, j, k in zip(a.values, b.values, c.values)]

    ## Separate filter values
    df['foreground_pattern']  = [i[0:3] for i in df['filter_'].values]
    df['background_pattern']  = [i[3:5] for i in df['filter_'].values]


    for i, level in enumerate(order_background):
        df.loc[df['background_pattern'] == level, 'leniency_background'] = i
    for i, level in enumerate(order_foreground):
        df.loc[df['foreground_pattern'] == level, 'leniency_foreground'] = i
    
    df = df.sort_values(['leniency_foreground', 'leniency_background'])
    df['index'] = np.arange(len(df))
    
    df['filter_foreground_target']  = [i[0] for i in df['filter_'].values]
    df['filter_foreground_reads'] =  [i[1] for i in df['filter_'].values]
    df['filter_foreground_samples'] = [i[2] for i in df['filter_'].values]
    df['filter_background_reads'] = [i[3] for i in df['filter_'].values]
    df['filter_background_samples'] = [i[4] for i in df['filter_'].values]
    
    df['filter_background'] = prepare_backticks(df['filter_background_reads'], df['filter_background_samples'])
    df['filter_foreground'] = prepare_frontticks(df['filter_foreground_target'], df['filter_foreground_reads'], df['filter_foreground_samples'])

    return df


def calculate_mean_std(df, group_cols, target_cols, run_type_plot_dir, sample, decimals=1):
    # Restrict table to cohort
    all_samples = [s for s in run_type_plot_dir 
                   if run_type_plot_dir[sample] == run_type_plot_dir[s]]
    df = df.set_index('sample').loc[all_samples].reset_index()
    
    
    # Compute means and std
    df_means = df.groupby(group_cols)[target_cols].mean().reset_index().rename({col : 'mean_' + col 
                                                                     for col in target_cols}, axis = 1)
    df_std = df.groupby(group_cols)[target_cols].std().reset_index().rename({col : 'std_' + col 
                                                                     for col in target_cols}, axis = 1)

    df = df.merge(df_means, on = group_cols).merge(df_std, on = group_cols)
    
    for col in df.columns:
        if ('mean_' in col):
            df[col] = [np.round(i, decimals) for i in df[col]]
            if decimals == 0:
                df[col] = [int(i) for i in df[col]]
    return df
    
    
class plotting_parameters():

    def __init__(self, ticks_fontsize, axislabels_fontsize, legend_fontsize, axes_fontsize):
        self.ticks_fontsize = ticks_fontsize
        self.axislabels_fontsize = axislabels_fontsize
        self.legend_fontsize = legend_fontsize
        self.axes_fontsize = axes_fontsize
        self.log_scale = False
    
    def edit_scale(self, is_log_scale):
        self.log_scale = is_log_scale
        
    def add_saving_instructions(self, save, run_type_plot_dir, sample_plot_dir):
        self.save = save
        self.run_type_plot_dir = run_type_plot_dir #TODO simplify
        self.sample_plot_dir = sample_plot_dir  #TODO simplify
        
    def add_ticks(self, back_ticks, front_ticks):
        self.back_ticks = back_ticks
        self.front_ticks = front_ticks
        
    def add_y_label(self, y_label):
        self.y_label = y_label
        
    def add_x_label(self, xupper_axis_label, xlower_axis_label):
        self.xupper_axis_label = xupper_axis_label
        self.xlower_axis_label = xlower_axis_label
        
    def add_saving_path(self, plot_dir, base_plot, name_plot):
        self.plot_dir = plot_dir
        self.base_plot = base_plot
        self.name_plot = name_plot
    
    def add_plotting_data(self, data_both, data_eth, data_ohsu, 
                          serie_index, serie_intersection, serie_eth, serie_ohsu):
        self.data_both = data_both
        self.data_eth = data_eth
        self.data_ohsu = data_ohsu
        self.serie_index = serie_index
        self.serie_intersection = serie_intersection
        self.serie_eth = serie_eth
        self.serie_ohsu = serie_ohsu

    def add_labels(self, intersection_label, eth_label, ohsu_label):
        self.intersection_label = intersection_label
        self.eth_label = eth_label
        self.ohsu_label = ohsu_label
    
    def add_color_options(self, color1, color2, color3, color4, colorgrid):
        self.color1 = color1
        self.color2 = color2
        self.color3 = color3
        self.color4 = color4
        self.colorgrid = colorgrid
    
    def edit_marker(self, marker_type, marker_size, markeredgewidth):
        self.marker_type = marker_type
        self.marker_size = marker_size
        self.markeredgewidth =  markeredgewidth
        
    def add_title(self, title):
        self.title = title
        

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
        if len(change_val) > 1:
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
        
        
def plot_text_dev(Y, T, position='top', color='black', font=None):
    def x_centering(y):
        if isinstance(y, float) or isinstance(y, np.float64):
            delta_x = - 0.5
        elif isinstance(y, int) or isinstance(y, np.int64):
            digits = int(math.log10(y)) + 1 if y else 1
            delta_x = - (0.25 + (digits / 8)) if digits > 1 else - 0.25 
        else:
            assert(False)
        return delta_x
    
    def y_shift(Y, position):
        if position == 'bottom':
            delta_y = - (Y[x] /5.5)
        elif position == 'top': 
            delta_y = Y[x] / 12
        return delta_y
    
    start = 0 
    step = 3
    #if max(Y) > 0:
    font['color'] = color
    Y = np.array(Y)
    T = np.array(T)
    for x in np.arange(start, len(Y), step):
        delta_x = x_centering(Y[x])
        # TESTING print(x, Y[x], delta_x)
        delta_y = y_shift(Y, position)
        plt.text(x + delta_x , Y[x] + delta_y , T[x], ha='left', **font)
        #print('next') 
        
def print_statistics(serie, label):
    stat_text(serie, label)
    
    
def print_ratio(serie_JP, serie_GP, label='ratio JP/GP'):
    serie = np.round(np.divide(list(serie_JP), list(serie_GP)), 2)
    stat_text(serie, label)
    
def stat_text(serie, label):
      print(f'Stats {label}', 
          f'/ min: {np.min(serie)}', 
          f'/ max: {np.round(np.max(serie), 2)}', 
          f'/ mean: {np.round(np.mean(serie), 2)}', 
          f'/ median: {np.median(serie)}', 
          f'/ non_zero: {sum(serie > 0 )}/{len(serie)}')
        
def plot_intersection_bars(param):
    # Get series 
    if param.serie_intersection is not None:
         intersection = param.data_both[param.serie_intersection]

    index = np.arange(len(param.back_ticks))

    eth = param.data_eth[param.serie_eth]
    ohsu = param.data_ohsu[param.serie_ohsu]

    text_font = {'size':'12', 'weight':'bold'}
    alpha_grid = 0.3
    width = 0.4

    fig, ax1 = plt.subplots(figsize=(15, 6))
    ax2 =  ax1.secondary_xaxis('top')   
    plt.grid(b=True, axis = 'both', which='major', color=param.colorgrid, linestyle='-', alpha=alpha_grid)
    plt.grid(b=False, axis = 'both', which='minor', color=param.colorgrid, linestyle='--', alpha=alpha_grid)
    
    print(param.y_label)
    if param.serie_intersection is not None:
        plt.bar(index, intersection, width=width, 
                color=param.color1, label=param.intersection_label)
        print_statistics(intersection, param.intersection_label)
    plt.plot(index, eth, alpha=1, color=param.color3,
             linestyle = 'None', markerfacecolor='None', marker=param.marker_type,
             markersize=param.marker_size, markeredgewidth=param.markeredgewidth,
             label = param.eth_label) 
    print_statistics(eth, param.eth_label)
    plt.plot(index, ohsu, alpha=1, color=param.color2,
             linestyle = 'None', markerfacecolor='None', marker=param.marker_type, 
             markersize=param.marker_size, markeredgewidth=param.markeredgewidth,
             label = param.ohsu_label)
    print_statistics(ohsu, param.ohsu_label)
    
    print_ratio(ohsu, eth)
    print_ratio(intersection, ohsu, 'stats inter/OHSU')
    print_ratio(intersection, eth, 'stats inter/ETH')

    

    lower_bound = param.data_eth[param.serie_eth] - param.data_eth[param.serie_eth.replace('mean', 'std')]
    upper_bound = param.data_eth[param.serie_eth] + param.data_eth[param.serie_eth.replace('mean', 'std')]
    plt.fill_between(index, 
                     lower_bound.apply(lambda x: max(x, 0)), upper_bound, 
                     alpha=0.35, color=param.color3) 
    lower_bound = param.data_ohsu[param.serie_ohsu] - param.data_ohsu[param.serie_ohsu.replace('mean', 'std')]
    upper_bound = param.data_ohsu[param.serie_ohsu] + param.data_ohsu[param.serie_ohsu.replace('mean', 'std')]
    plt.fill_between(index, 
                     lower_bound.apply(lambda x: max(x, 0)), upper_bound,
                     alpha=0.35, color=param.color2)
    
    plot_text_dev(ohsu, ohsu, 'top', color=param.color2, font=text_font)
    plot_text_dev(eth, eth, 'top', color=param.color3, font=text_font)
    if (param.serie_intersection) is not None and (param.log_scale): # Skip the intersection size if not log scale
        plot_text_dev(intersection, intersection, color=param.color4, font=text_font)


    if param.log_scale:
        plt.yscale('symlog')
        
    max_scale = np.max([ohsu.values, eth.values])
    min_scale = np.min([ohsu.values, eth.values])
    
    # TESTING
#     plt.plot(np.arange(0,len(index), 3), [np.mean([max_scale, min_scale]) ] * len(np.arange(0,len(index), 3)), linestyle='None', marker='|', markersize=15, color='maroon', markeredgewidth=3)
    

        
    ax1.set_xticks(index, 
               labels = param.front_ticks,
               rotation = 90, 
               ha = 'center', 
               fontsize=param.ticks_fontsize)

    ax2.set_xticks(index, 
               labels = param.back_ticks,
               rotation = 90, 
               ha = 'center', 
               fontsize=param.ticks_fontsize)


    plt.legend(fontsize=param.legend_fontsize)
    plt.ylabel(param.y_label, fontsize=param.axes_fontsize)
    ax2.set_xlabel(param.xupper_axis_label, fontsize=param.axes_fontsize)
    ax1.set_xlabel(param.xlower_axis_label, fontsize=param.axes_fontsize)
    plt.title(param.title)

    save_path = os.path.join(param.plot_dir, f'{param.base_plot}_{param.name_plot}.pdf')
    print("\n saving path is to {}".format(save_path))
    if param.save:
        print("Saving!")
        plt.savefig(save_path, bbox_inches='tight')

    plt.show()
