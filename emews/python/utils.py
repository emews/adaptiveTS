import csv
import json
from datetime import datetime


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


def format_parameters(param_str: str, tick0_date: str = '2020-03-02', school_closure_date: str = '2020-03-17',
                      stay_at_home_date: str = '2020-03-21'):
    """Formats a single set of parameters from the specified JSON string."""

    # print("Parameter Parse Start {}".format(datetime.now().strftime("%%FT%%T%%z")))
    p_map = json.loads(param_str)

    tick0_date = datetime.strptime(tick0_date, '%Y-%m-%d')
    school_closure_date = datetime.strptime(school_closure_date, '%Y-%m-%d')
    stay_at_home_date = datetime.strptime(stay_at_home_date, '%Y-%m-%d')
    feb1_date = datetime.strptime('2020-02-01', '%Y-%m-%d')
    stop_date = datetime.strptime('2020-06-04', '%Y-%m-%d')

    school_closure_tick = (school_closure_date - tick0_date).days * 24
    stay_at_home_tick = (stay_at_home_date - tick0_date).days * 24
    # how many previous to starting date
    feb1_tick = (tick0_date - feb1_date).days * 24
    stop_at_tick = (stop_date - tick0_date).days * 24

    try:
        # scheduled initial infections
        infected_count = int(round(p_map['infected.count']))
        init_exposure_tick = p_map['initial.exposure.tick']
        init_exposure_tick = int(round(init_exposure_tick / 24.0)) * 24
        p_map['wild.infection.count.1'] = '{}|{}'.format(init_exposure_tick, infected_count)

        if 'start.home.isolation.delay.max' in p_map:
            p_map['start.home.isolation.delay.max'] = int(round(p_map['start.home.isolation.delay.max']))

        # Assumes V3 sythn pop
        p_map['intervention.ptype.all_but_school'] = 'workplace,other household,restaurant,supermarket,store,convenience store,place of worship,recreation,nursing home,library,hospital,museum'

        p_map['intervention.close_school'] = "{},0,school".format(school_closure_tick)
        sah_p = float(p_map['stay.at.home.probability'])
        p_map['intervention.stay_at_home_1'] = "{},{},all_but_school".format(stay_at_home_tick, 1 - sah_p)
        sba_p = p_map['stoe.behavioral.adjustment.probability']
        p_map['stoe.behavioral.adjustment.lt18.1'] = "{},{}".format(stay_at_home_tick, sba_p)
        p_map['stoe.behavioral.adjustment.18_40.1'] = "{},{}".format(stay_at_home_tick, sba_p)
        p_map['stoe.behavioral.adjustment.40_60.1'] = "{},{}".format(stay_at_home_tick, sba_p)
        p_map['stoe.behavioral.adjustment.gte60.1'] = "{},{}".format(stay_at_home_tick, sba_p)
        p_map['stop.at'] = stop_at_tick
        # puts peak on feb 1
        p_map['seasonality.peak'] = feb1_tick
    except KeyError:
        print(p_map, flush=True)
        raise

    return p_map

    # print("Parameter Parse End {}".format(datetime.now().strftime("%%FT%%T%%z")))

