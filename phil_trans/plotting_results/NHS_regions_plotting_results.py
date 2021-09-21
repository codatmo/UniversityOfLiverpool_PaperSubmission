import os
import pandas as pd
import numpy as np
import argparse
import re
import glob
import json
from datetime import datetime, timedelta, date
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.dates as mdates


def plot_daily_deaths(deaths_np, data, sdate, edate, actual_deaths, place_name, plotting_name, plotting_source, beta, NEES, lb_std, ub_std):
    
    
    sdate_deaths = date(2020,3,23)
    edate_deaths = date(2021,1,6)
    date_list = pd.date_range(sdate,edate-timedelta(days=1),freq='d')[-200:]
    mean_deaths = np.mean(deaths_np, axis=0)[-200:]

    y_max = np.max(np.concatenate([ub_std[-200:],data['deaths'][-193:]])) + 20
    y_min = np.min(np.concatenate([lb_std[-200:],data['deaths'][-193:]])) - 10.0
    deaths_summary_fig = plt.figure('deaths_summary', figsize=(15.0, 10.0))
    fig_title = r'$\sigma_\beta$ = {}'.format(beta) + ', NEES = {}'.format(round(NEES, 3))
    plt.title(fig_title, fontsize=50)
    plt.xlabel('Date', fontsize=35)
    plt.ylabel('Daily Deaths', fontsize=35)
    #plt.ylim((y_min, y_max))
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d-%b-%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=30))
    plt.plot(date_list[-200:-7], data['deaths'][-193:], 'k.', markersize=19)

    plt.plot(date_list[-7:], actual_deaths, 'k.', color='g', markersize=19)
    plt.plot(date_list, mean_deaths, color="#D55E00")
    plt.vlines(date_list[-7], y_min, y_max, color='r', linestyles='dashed')
    plt.fill_between(date_list, lb_std[-200:], ub_std[-200:], alpha=0.25, color="#D55E00")
    plt.gcf().autofmt_xdate()
    plt.grid()
    plt.savefig(place_name + '.png', dpi=300)
    plt.clf()
    
    

NHS_names = ['East_of_England', 'London', 'Midlands', 'North_East_and_Yorkshire', 'North_West', 'South_East', 'South_West']
plotting_NHS_names = ['East of England', 'London', 'Midlands', 'North East and Yorkshire', 'North West', 'South East', 'South West']

actual_deaths_ = './nhs_regions_daily_deaths_.csv'

actual_deaths_df = pd.read_csv(actual_deaths_)
actual_deaths_df = actual_deaths_df[actual_deaths_df.columns[-27:]]

sources_list = ['deaths']
plotting_sources_list = ['deaths']

NEES_list = []
RMSE_list = []
MAE_list = []
        
data_folders = './example/'
       
for idx, (name, plotting_name) in enumerate(zip(NHS_names, plotting_NHS_names)):
          #  print(name)
            
    path_ = data_folders + name
            
    json_path = path_ + '/' + name +'.json'
        
    with open(json_path) as json_file:
       data = json.load(json_file)
            
            
    samples_df_list = []
    sample_csv_files = glob.glob(path_+'/output/*.csv')
    for filename in sample_csv_files:
        samples_df_list.append(pd.read_csv(filename, comment='#'))
    samples_df = pd.concat(samples_df_list)
    
    column_names_deaths = ['pred_deaths.' + str(t) for t in range(1, int(data['deaths_length']+8))]
    deaths_df = samples_df.loc[:, column_names_deaths]
    deaths_np = deaths_df.to_numpy()
          
    mean = np.mean(deaths_np, axis=0)
    mean_pred = mean[-7:]
        
    # 95% and 5% confidence intervals
    #median = np.median(deaths_np, axis=0)
    #ub_quan = np.quantile(deaths_np, 0.95, axis=0)
    #lb_quan = np.quantile(deaths_np, 0.05, axis=0)
    
    num_of_samples = len(deaths_np)
    num_of_deaths = len(mean_pred)
   
    sdate_deaths = date(2020,3,23)
    edate_deaths = date(2021,1,6)
  
        
    actual_deaths = [int(x) for x in actual_deaths_df.iloc[idx].tolist()[-7:]]
            
    first_ = deaths_np.T - (mean.reshape((len(mean), 1)) * np.ones(num_of_samples))
    second_ = (deaths_np.T - (mean.reshape((len(mean), 1)) * np.ones(num_of_samples))).T
    third_ = np.dot(first_, second_)
    Q_ = (1/(num_of_samples-1))*third_

    Q_diag_ = Q_.diagonal()
          
    ub_std = mean + np.sqrt(Q_diag_)
    lb_std = mean - np.sqrt(Q_diag_)
    
    summed_arr = 0
            
    actual_deaths = [int(x) for x in actual_deaths_df.iloc[idx].tolist()[-7:]]
        
    Q_diag_pred = Q_diag_[-7:]
            
    for idx, (death_, mean_sample_) in enumerate(zip(actual_deaths, mean_pred)):   
        summed_arr += (mean_sample_ - death_) * 1/Q_diag_pred[idx] * (mean_sample_ - death_)
              
    NEES = 1/len(mean_pred) * summed_arr
            
    rmse = np.linalg.norm(mean_pred - actual_deaths) / np.sqrt(len(mean_pred))
    RMSE_list.append(rmse)
    NEES_list.append(NEES)
    plot_daily_deaths(deaths_np, data, sdate_deaths, edate_deaths, 
                          actual_deaths, name, plotting_name, '111 Calls and Deaths', '0.01', NEES, lb_std, ub_std)
            
print(NEES_list)
print(RMSE_list)
        
        

