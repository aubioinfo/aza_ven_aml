###############################################################################################
# Aim: Score vs. VEN/AZA response probability
# Description: Relationship between RF8 score vs. VEN/AZA response probability
###############################################################################################

import sys
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns

import sklearn.neighbors._base
sys.modules['sklearn.neighbors.base'] = sklearn.neighbors._base

plt.rcParams.update({'font.size': 9})
plt.rcParams["font.family"] = "Arial"
palette = sns.color_palette("deep")


if __name__ == "__main__":

    bs_number = 1000  # bootstrapping
    random.seed(1)

    bin_size = 0.1

    start_time = time.time()
    print('Raw data read in ...')

    # Read RF8 prediction data
    fnIn = '01.RF8_SurvivalAnalysis.xlsx' 
    y_pred_RF8 = []
    y_true = []
    output_curve_fn = 'RF8_vs_ORR' + '.png'
    
    # Reading only one sheet, since you have only one sheet in the file
    data = pd.read_excel(fnIn, sheet_name=0, header=0, index_col=0)
    
    # Extract relevant columns for prediction and true response
    y_pred_RF8 = data['RF8.prob.CR'].tolist()  # RF8 predicted response probability
    y_true = data['Response'].tolist()         # True response (0 or 1)

    y_true = np.array(y_true)
    y_pred_RF8 = np.array(y_pred_RF8)
    score_list_RF8 = np.arange(0.0, 1.01, 0.01)
    RF8_num = len(score_list_RF8)

    RF8high_ORR_list = [[] for _ in range(RF8_num)]
    RF8low_ORR_list = [[] for _ in range(RF8_num)]
    RF8_ORR_list = [[] for _ in range(RF8_num)]
    RF8_patientNUM_list = [[] for _ in range(RF8_num)]
    sampleNUM = len(y_true)
    idx_list = range(sampleNUM)
    print('Sample num:', sampleNUM)

    # Bootstrap resampling and ORR calculation
    for bs in range(bs_number):
        idx_resampled = random.choices(idx_list, k=sampleNUM)
        y_true_resampled = y_true[idx_resampled]
        y_pred_RF8_resampled = y_pred_RF8[idx_resampled]

        for score_i in range(len(score_list_RF8)):
            score = score_list_RF8[score_i]
            idx_high_interval = y_pred_RF8_resampled >= score
            y_true_high = y_true_resampled[idx_high_interval]
            Rhigh_num = sum(y_true_high)
            tot_high_num = len(y_true_high)
            patientRatio_temp = sum(y_pred_RF8_resampled < score) / sampleNUM
            RF8_patientNUM_list[score_i].append(patientRatio_temp)
            
            if not tot_high_num:
                RF8high_ORR_list[score_i].append(RF8high_ORR_list[score_i-1][-1])
            else:
                ORRhigh_temp = Rhigh_num / tot_high_num
                RF8high_ORR_list[score_i].append(ORRhigh_temp)

            idx_low_interval = y_pred_RF8_resampled < score
            y_true_low = y_true_resampled[idx_low_interval]
            Rlow_num = sum(y_true_low)
            tot_low_num = len(y_true_low)
            
            if not tot_low_num:
                RF8low_ORR_list[score_i].append(0)
            else:
                ORRlow_temp = Rlow_num / tot_low_num
                RF8low_ORR_list[score_i].append(ORRlow_temp)

            if sum(y_pred_RF8_resampled <= score + bin_size / 2) < 0.01 * len(y_pred_RF8_resampled):
                idx_interval = []
            elif sum(y_pred_RF8_resampled > score - bin_size / 2) < 0.01 * len(y_pred_RF8_resampled):
                idx_interval = (y_pred_RF8_resampled > score - bin_size / 2)
            else:
                idx_interval = (y_pred_RF8_resampled <= score + bin_size / 2) & (y_pred_RF8_resampled > score - bin_size / 2)
                
            y_true_temp = y_true_resampled[idx_interval]
            R_num = sum(y_true_temp)
            tot_num = len(y_true_temp)
            
            if not tot_num:
                RF8_ORR_list[score_i].append(0)
            else:
                ORR_temp = R_num / tot_num
                RF8_ORR_list[score_i].append(ORR_temp)
            
            if sum(y_pred_RF8_resampled > score - bin_size / 2) < 0.01 * len(y_pred_RF8_resampled):
                break

    # Remove empty elements for high scores
    for i in range(len(RF8high_ORR_list)):
        if len(RF8high_ORR_list[i]) == 0:
            break
    RF8high_ORR_list = RF8high_ORR_list[0:i]
    RF8low_ORR_list = RF8low_ORR_list[0:i]
    RF8_ORR_list = RF8_ORR_list[0:i]
    RF8_patientNUM_list = RF8_patientNUM_list[0:i]
    score_list_RF8 = score_list_RF8[0:i]

    # Compute mean and confidence intervals
    RF8high_ORR_mean = [np.mean(c) for c in RF8high_ORR_list]
    RF8high_ORR_05 = [np.quantile(c, 0.05) for c in RF8high_ORR_list]
    RF8high_ORR_95 = [np.quantile(c, 0.95) for c in RF8high_ORR_list]
    RF8low_ORR_mean = [np.mean(c) for c in RF8low_ORR_list]
    RF8low_ORR_05 = [np.quantile(c, 0.05) for c in RF8low_ORR_list]
    RF8low_ORR_95 = [np.quantile(c, 0.95) for c in RF8low_ORR_list]
    RF8low_patientRatio_mean = [np.mean(c) for c in RF8_patientNUM_list]
    RF8_ORR_mean = [np.mean(c) for c in RF8_ORR_list]
    RF8_ORR_05 = [np.quantile(c, 0.05) for c in RF8_ORR_list]
    RF8_ORR_95 = [np.quantile(c, 0.95) for c in RF8_ORR_list]
    RF8_patientRatio_mean = [np.mean(c) for c in RF8_patientNUM_list]

    print('RF8 response odds:')
    for i in range(len(RF8high_ORR_95)):
        print(score_list_RF8[i], RF8high_ORR_mean[i], RF8low_ORR_mean[i], RF8_ORR_mean[i], RF8low_patientRatio_mean[i])

    # Save results to CSV files
    df = pd.DataFrame({'RF8_score': score_list_RF8, 'Prob_mean': RF8_ORR_mean, 'Prob_lower': RF8_ORR_05, 'Prob_upper': RF8_ORR_95})
    df.to_csv('source_data_forplot.csv', index=False)

    # Plotting Score-Prob curve
    fig1, axes = plt.subplots(1, 1, figsize=(6.5, 2.8))
    fig1.subplots_adjust(left=0.1, bottom=0.15, right=0.98, top=0.96)
    axes.plot(score_list_RF8, RF8_ORR_mean, '-', color='r')
    axes.fill_between(score_list_RF8, RF8_ORR_05, RF8_ORR_95, facecolor='r', alpha=0.25)
    axes.set_ylabel("Response probability (%)", color="k")
    axes.set_xlabel('RF8')  # RF8 score

    axes.set_ylim([-0.02, 1.02])
    axes.set_yticks([0, 0.25, 0.5, 0.75, 1])
    axes.set_yticklabels([0, 25, 50, 75, 100])
    axes.spines['right'].set_visible(False)
