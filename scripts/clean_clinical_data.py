#!/usr/bin/env python
# coding: utf-8
 
import pandas as pd
from datetime import datetime
import argparse
 
def main(infile):
 
    df = pd.read_csv(infile, delimiter=',', encoding='latin-1')
    print(df.loc)
 
    noted_ids = df["1.2|GLASS-NL ID"]
    df = df.loc[noted_ids.isnull() == False]
    noted_ids = df["1.2|GLASS-NL ID"]
    cent_vumc = df["CENT#VUmc"]
    df['Abbrev'] = ''
    df.loc[df["CENT#VUmc"] == 1, 'Abbrev'] = 'AUMC'
    df.loc[df["CENT#ErasmusMC"] == 1, 'Abbrev'] = 'EMCR'
    df.loc[df["CENT#UMCU"] == 1, 'Abbrev'] = 'UMCU'
 
    ids = noted_ids.apply(round)
    ids = ids.apply("{:03d}".format)
 
 
    df['GLASS_ID'] = ids
    df['INSTITUTE'] = df['Abbrev']
    df['FULL_ID'] = 'GLNL_' + df['INSTITUTE'] + '_' + df['GLASS_ID']
    new_df = df[['GLASS_ID', 'INSTITUTE', 'FULL_ID', 'Record Id']]
 
    print(f"Number of unique patients: {len(new_df['GLASS_ID'].unique())}")
    print(f"Number of records: {len(new_df['GLASS_ID'])}")
 
    operation_date_names = ["DATE_1stSURG",
                                                "DATE_2ndSURG",
                               "DATE_3ndSURG_1",
                               "DATE_SURG4_",
                               "DATE_SURG5",
                               "DATE_SURG5_1",
                                                "DATE_SURG5_1_1"]
 
    treatment_checks = ["treatment_diagnosis",
                    "treatment_1stprog",
                    "TREAT_PROG3",
                    "TREAT_PROG4",
                    "TREAT_PROG5",
                    "TREAT_PROG5_1",
                    "TREAT_PROG5_1_1"]
 
    radioth_names = ["RT_START",
                "3.56|Start radiation therapy",       
                "4.39|Start Radiotherapy",
        "5.39|Start Radiotherapy",
        "6.36|Start Radiotherapy",
        "7.36|Start Radiotherapy",
        "8.36|Start Radiotherapy"]
 
    chemo_names = ["2.55|Start Chemotherapy",
                "3.60|Start Chemotherapy",
                "4.42|Start Chemotherapy",
                "5.43|Start Chemotherapy",
                "6.40|Start Chemotherapy",
                "7.40|Start Chemotherapy",
                "8.40|Start Chemotherapy"]
 
    prog_names = ["DATE_DIAGN",
                "DATE_PROGR",
                "PROG_2ND_1",
                "PROG_4",
                "PROG_5",
                "PROG_5_1",
                "PROG_5_1_1"]
 
    kps_score_prog_names = ["KPS_DIAGN",
                            "KPS_PROG2",
                            "KPS_PROG3",
                            "KPS_PROG4",
                            "KPS_PROG5",
                            "KPS_PROG5_2",
                            "KPS_PROG5_2_1"]
 
    who_score_prog_names = ["WHODIAGN",
                            "WHO_PROG2",
                            "WHO_PROG3",
                            "WHO_PROG4",
                            "WHO_PROG5",
                            "WHO_PROG5_1",
                            "WHO_PROG5_1_1"]
 
    kps_score_prog_names_after = ["KPS_SURG1",
                            "KPS_SURG2",
                            "KPS_PROG3_1",
                            "KPS_PROG4_1",
                            "KPS_PROG5_1",
                            "KPS_PROG5_1_1",
                            "KPS_PROG5_1_1_1"]
 
    who_score_prog_names_after = ["WHO_SURG1",
                            "WHO_SURG2",
                            "WHO_PROG3_1",
                            "WHO_PROG4_1",
                            "WHO_PROG5_2",
                            "WHO_PROG5_2_1",
                            "WHO_PROG5_2_1_1"]
 
    chemo_type_names = ["CHEMO_TYPE",
                    "CHEMO_TYPE_2nd_1",
                    "CHEMO_TYPE_2nd_2",
                    "CHEMO_TYPE_4",
                    "CHEMO_TYPE_5",
                    "CHEMO_TYPE_5_1",
                    "CHEMO_TYPE_5_1_1"]
 
    chemo_stop_names = ["PRESTOP_CHEMO1",
                    "PRESTOP_CHEMO_SURG2",
                    "PRESTOP_CHEMO3",
                    "PRESTOP_CHEMO4",
                    "PRESTOP_CHEMO5",
                    "PRESTOP_CHEMO5_1",
                    "PRESTOP_CHEMO5_1_1"]
 
    extent_resection_names = ["BR_SURG_RES_TP",
                            "BR_SURG_RES_TPx_2nd_1",
                            "BR_SURG_RES_TP_1",
                            "BR_SURG_RES_TP_4",
                            "BR_SURG_RES_TP_5",
                            "BR_SURG_RES_TP_5_1",
                            "BR_SURG_RES_TP_5_1_1"]
 
 
    def get_date_info(value, alternative):
        try:
            trdate = datetime.strptime(value, '%d-%m-%Y')
            date_found = True
        except TypeError:
            trdate = alternative
            date_found = False
        return trdate, date_found
 
    def get_chemo_data(row, i):
        chemo_type = row[chemo_type_names[i]]
        type_str = ""
        if chemo_type == 0:
            type_str = "Temozolomide"
        elif chemo_type == 1:
            type_str = "Lomustine"
        elif chemo_type == 2:
            type_str = "PCV"
        elif chemo_type == 3:
            type_str = "Experimental"
           
        chemo_stop = row[chemo_stop_names[i]]
        stop_str = ""
        if chemo_stop == 1:
            stop_str = "Yes"
        elif chemo_stop == 2:
            stop_str = "No"
        elif chemo_stop == 3:
            stop_str = "Unknown"
        return {
            'Chemotherapy type': type_str,
            'Chemotherapy stopped prematurely': stop_str
        }
 
    def get_surgery_data(row, i):
        sur_type = row[extent_resection_names[i]]
        type_str = ""
        if sur_type == 0:
            type_str = "Biopsy"
        elif sur_type == 1:
            type_str = "Partial resection"
        elif sur_type == 2:
            type_str = "Complete resection"
        elif sur_type == 3:
            type_str = "Resection"
 
        return {
            'Extent of resection': type_str
        }
 
    def get_kps_score(value):
        if value == 110:
            return 'Unknown'
        else:
            return value
       
    def get_who_score(value):
        if value == 6:
            return 'Unknown'
        else:
            return value
       
    def get_cause_of_death(value):
        if value == 0:
            return 'Tumor'
        if value == 1:
            return 'Treatment related'
        if value == 2:
            return 'Other'
       
 
 
    df_total = pd.DataFrame()
    for index, row in df.iterrows():
   
        for i, dname in enumerate(prog_names):
            value = row[dname]
            if not pd.isnull(value):
                dttime = datetime.strptime(row[dname], '%d-%m-%Y')
                birthdate = datetime.strptime(row['BIRTH'], '%d-%m-%Y')
               
                if i == 0:
                    nm = 'Diagnosis'
                else:
                    nm = 'Progression'
                df_total = df_total.append({
                    'FULL_ID': row['FULL_ID'],
                    'RecordID': row['Record Id'],
                    'date': dttime,
                    'prog_count': i,
                    'event': nm,
                    'date found': True,
                    'date_of_birth': birthdate,
                    'age': dttime.year - birthdate.year - ((dttime.month, birthdate.day) < (dttime.month, birthdate.day)),
                    'kps score before treatment': get_kps_score(row[kps_score_prog_names[i]]),
                    'who score before treatment': get_who_score(row[who_score_prog_names[i]]),
                    'kps score after treatment': get_kps_score(row[kps_score_prog_names_after[i]]),
                    'who score after treatment': get_who_score(row[who_score_prog_names_after[i]])
                }, ignore_index=True)
               
            treatment_sur = treatment_checks[i] + "#Surgery"
            date = row[operation_date_names[i]]
            trdate, date_found = get_date_info(date, dttime)
            if row[treatment_sur] == 1 or date_found:
                date = row[operation_date_names[i]]
                trdate, date_found = get_date_info(date, dttime)
                data = {
                    'FULL_ID': row['FULL_ID'],
                    'RecordID': row['Record Id'],
                    'date': trdate,
                    'prog_count': i,
                    'event': 'Surgery',
                    'date found': date_found,
                    'kps score before treatment': get_kps_score(row[kps_score_prog_names[i]]),
                    'who score before treatment': get_who_score(row[who_score_prog_names[i]]),
                    'kps score after treatment': get_kps_score(row[kps_score_prog_names_after[i]]),
                    'who score after treatment': get_who_score(row[who_score_prog_names_after[i]])
                }
               
                data.update(get_surgery_data(row, i))
                df_total = df_total.append(data, ignore_index=True)
               
            treatment_rt = treatment_checks[i] + "#Radiation"
            date = row[radioth_names[i]]
            trdate, date_found = get_date_info(date, dttime)
            if row[treatment_rt] == 1 or date_found:
                date = row[radioth_names[i]]
                trdate, date_found = get_date_info(date, dttime)
                df_total = df_total.append({
                    'FULL_ID': row['FULL_ID'],
                    'RecordID': row['Record Id'],
                    'date': trdate,
                    'prog_count': i,
                    'event': 'Radiotherapy',
                    'date found': date_found,
                    'kps score before treatment': get_kps_score(row[kps_score_prog_names[i]]),
                    'who score before treatment': get_who_score(row[who_score_prog_names[i]]),
                    'kps score after treatment': get_kps_score(row[kps_score_prog_names_after[i]]),
                    'who score after treatment': get_who_score(row[who_score_prog_names_after[i]])
                }, ignore_index=True)
               
            treatment_chemo = treatment_checks[i] + "#Chemotherapy"
            date = row[chemo_names[i]]
            trdate, date_found = get_date_info(date, dttime)
            if row[treatment_chemo] == 1 or date_found:
                date = row[chemo_names[i]]
                trdate, date_found = get_date_info(date, dttime)
                data = {
                    'FULL_ID': row['FULL_ID'],
                    'RecordID': row['Record Id'],
                    'date': trdate,
                    'prog_count': i,
                    'event': 'Chemotherapy',
                    'date found': date_found,
                   'kps score before treatment': get_kps_score(row[kps_score_prog_names[i]]),
                    'who score before treatment': get_who_score(row[who_score_prog_names[i]]),
                    'kps score after treatment': get_kps_score(row[kps_score_prog_names_after[i]]),
                    'who score after treatment': get_who_score(row[who_score_prog_names_after[i]]),
                   
                }
                data.update(get_chemo_data(row, i))
 
                df_total = df_total.append(data, ignore_index=True)
               
            
                
        death_date, death_found = get_date_info(row['DATE_DEATH'], row['DATE_DEATH'])
        last_followup, followup_found = get_date_info(row['LAST_SEEN'], row['LAST_SEEN'])
        if death_found:
            data = {
                    'FULL_ID': row['FULL_ID'],
                    'RecordID': row['Record Id'],
                    'date': death_date,
                    'event': 'Death',
                    'cause_of_death': get_cause_of_death(row['COD_TP'])
            }
            df_total = df_total.append(data, ignore_index=True)
        if followup_found:
            data = {
                    'FULL_ID': row['FULL_ID'],
                    'RecordID': row['Record Id'],
                    'date': last_followup,
                    'event': 'Last followup',
            }
            df_total = df_total.append(data, ignore_index=True)
       
    ### Get surgery counts
    ### NB: assumes the data is already correctly sorted
 
    df_total['is_surgery'] = df_total['event'] == 'Surgery'
    df_total = df_total.sort_values(['RecordID', 'date'], ascending = [True, True])
 
    df_total['surgery_count'] = df_total.groupby('RecordID')['is_surgery'].cumsum(axis = 0)
    return df_total
 
if __name__ == '__main__':
 
    parser = argparse.ArgumentParser(description='Coregister and resample for segmentation.')
    parser.add_argument('-i', help='Path to Castor output file', required=True)
    parser.add_argument('-o', help='Path to write output .csv file', required=True)
    args = parser.parse_args()
 
    result_dataframe = main(args.i)
 
    result_dataframe.to_csv(args.o)
