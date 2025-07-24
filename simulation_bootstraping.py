import numpy as np
import csv
import os 
import subprocess
import datetime

def linspace_circ(clump_count, num, h_width=5):
    return np.array([np.linspace(h_width/(2*clump_count), h_width - (h_width/(2*clump_count)), clump_count)[i%clump_count] for i in range(num)])

def create_temp_csv(fish_params, fish_counts, f_ic, s_ic, temp_path = "temp_ics.csv"):
    ic_csv_arr = [f_ic]
    params=[]
    for param_index, fish_count in enumerate(fish_counts):
        for _ in range(fish_count):
            params.append(fish_params[param_index])
    params = np.transpose(np.array(params))
    for p in params:
        ic_csv_arr.append(p)
    ic_csv_arr = np.transpose(ic_csv_arr)
    np.random.shuffle(ic_csv_arr)
    ic_csv_arr = np.transpose(ic_csv_arr)


    with open(temp_path, mode="w") as f:
        wtr = csv.writer(f)
        for row in ic_csv_arr:
            wtr.writerow(row)
        wtr.writerow(s_ic)

    current_path = os.getcwd()
    return current_path + f"/{temp_path}"

def mp_run(fish_params, fish_counts, f_ic, s_ic, out_path, only_final=True):
    temp_csv_path = str(datetime.datetime.now()).replace(' ','_').replace(":",".") +".csv"
    csv_path = create_temp_csv(fish_params, fish_counts, f_ic, s_ic,temp_path = temp_csv_path)
    if only_final:
        pros = subprocess.Popen(f"cd sim/target/release && ./seeking_service_main {out_path} {csv_path} true", shell=True, stdout=subprocess.DEVNULL)
        # pros = subprocess.Popen(f"cd sp/target/release && ./fp {out_path} {csv_path} true")
    else:
        pros = subprocess.Popen(f"cd sim/target/release && ./seeking_service_main {out_path} {csv_path}", shell=True, stdout=subprocess.DEVNULL)
    pros.wait()
    os.system(f"rm {csv_path}")

if __name__ == "__main__":
    f_ic = linspace_circ(1000,1000,5)
    s_ic = linspace_circ(25,25,5)
    fp = [[0.2,0,0],[0,0,0.2]]
    fcs = [500,500]
    os.getcwd()
    save_dir = f"{os.getcwd()}/trial_runs/multi_wide_bootstrap_test/"
    mp_run(fp,fcs,f_ic,s_ic, save_dir, only_final=True)
