#!/usr/bin/env python


header = None

def clean_line(params):
    for i in range(len(params)):
        params[i] = params[i].strip("\ufeff")
        params[i] = params[i].strip('"')
        params[i] = params[i].strip()
    
    return params


def find_key(params, key):
    for i in range(len(params)):
        if params[i] == key:
            return i
    
    print(params)
    
    print("Error, could not find key: " + key)
    import sys
    sys.exit(1)


with open("/home/youri/Downloads/Glioma_Longitudinal_AnalySiS_export_20230508.csv", "r") as fh:
    for line in fh:
        if not header:
            header = clean_line(line.split(";"))
            
            key_pid = find_key(header, "Participant Id")
            key_gender = find_key(header, "Gender")
            
            #print("castor_id\tgender")
        else:
            params = clean_line(line.split(";"))
            
            pid = params[key_pid]
            gender = {'0':'male','1':'female'}[params[key_gender]]
            
            print("UPDATE patients SET")
            print("    sex = '" + gender + "',")
            print("    notes = COALESCE(notes || ' ' , '') || '[2023-05-08 gender extracted from castor: ' || COALESCE(sex , 'NA') || ' -> "+gender+"]'")
            print("WHERE")
            print("    castor_id = '"+pid+"' LIMIT 1;\n\n")



