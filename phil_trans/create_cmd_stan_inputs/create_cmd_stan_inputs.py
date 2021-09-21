import json
import io
import requests
import pandas as pd
import numpy as np
import datetime
from datetime import date, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import random
from scipy.special import expit, logit

init_list = []
NHS_names = ['East_of_England', 'London', 'Midlands', 
               'North_East_and_Yorkshire', 'North_West', 
               'South_East', 'South_West']

# Choose how many warm-up and samples you want to run with

warm_up_names = ['warm_up_256']
warm_up_num = ['256']

samples_names = ['samples_256']
samples_num = ['256']

# Change directory to location of Stan model
        
for wu_dir, wu_num in zip(warm_up_names, warm_up_num):
    for samples_dir, samp_num in zip(samples_names, samples_num):
        for place in NHS_names:
            for x in range(1,7):
                iu = '/mnt/data2/users/bighyp/us_states/src/models/predictions/one_feed_new/test_model/model_25 sample num_samples={} num_warmup={} adapt delta=0.9 algorithm=hmc engine=nuts max_depth=12 random seed=42 id={} data file=/mnt/data2/users/bighyp/us_states/src/daily_deaths/multiple_uk/phil_trans/_0.25/jsons/{}/{}.json init=/mnt/data2/users/bighyp/us_states/src/daily_deaths/multiple_uk/phil_trans/_0.25/jsons/{}/inits/inits{}.json output file=/mnt/data2/users/bighyp/us_states/src/daily_deaths/multiple_uk/phil_trans/_0.25/jsons/{}/output/samples{}.csv'.format(samp_num, wu_num, x, place, place, place, x, place, x)
                init_list.append(iu)


textfile = open('{}.txt'.format(date.today()), 'w') 
data = '{}'.format('\n'.join(init_list))
textfile.write(data)

#samp_num, wu_num, x, wu_dir, samples_dir, place, place, 