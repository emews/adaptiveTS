import csv


def format(val):
    try:
        return int(val)
    except:
        return val


def csv2cmd_lines(fname):
    rows = []
    with open(fname) as fin:
        reader = csv.reader(fin)
        header = next(reader)
        for items in reader:
            row = {header[i]: format(item) for i, item in enumerate(items)}
            # print(row, flush=True)
            cmd_line = f"--exp_seed={row['exp_seed']} --init_npar={row['init_npar']} --nrep={row['nrep']} --sim_budget={row['sim_budget']} --grid_npar={row['grid_npar']} --nTS_samp={row['nTS_samp']} --exp_id={row['exp_id']}"
            rows.append(f"{row['exp_id']}|{cmd_line}")
    
    return rows
