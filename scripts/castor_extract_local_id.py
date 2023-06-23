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
            key_local_id = find_key(header, "LOC_ID")
            
            print("castor_id\tlocal_id")
        else:
            params = clean_line(line.split(";"))
            
            castor_id = params[key_pid]
            local_id = params[key_local_id].lstrip('0')
            
            if len(local_id.strip()) >= 6:
              #print(pid + "\t" + local_id)
              print("UPDATE patients SET")
              print("    castor_id = '" + castor_id + "',")
              print("    notes = COALESCE(notes || ' ' , '') || '[2023-05-08 castor id added -> "+castor_id+"]'")
              print("WHERE")
              print("    castor_id IS NULL AND")
              print("    local_id = '"+local_id+"';")
              print("\n\n")

#UPDATE patients SET castor_id = 1337
#WHERE castor_id IS NULL AND local_id = 8264462;


