import files;
import string;
import sys;
import io;
import python;
import location;
import unix;
import emews;

import EQSQL;
import covid_model;

// deletes the specified directory
app (void o) rm_dir(string dirname) {
  "rm" "-rf" dirname;
}

// deletes the specified directories
app (void o) rm_dirs(file dirnames[]) {
  "rm" "-rf" dirnames;
}

app (void o) rm(file fname) {
  "rm" fname;
}

app (void o) rm(string fname) {
  "rm" fname;
}

string emews_root = getenv("EMEWS_PROJECT_ROOT");
string turbine_output = getenv("TURBINE_OUTPUT");
int resident_work_rank = string2int(getenv("RESIDENT_WORK_RANK"));
int procs_per_run = toint(getenv("PROCS_PER_RUN"));

int TASK_TYPE = string2int(argv("task_type", "0"));
int BATCH_SIZE = string2int(argv("batch_size"));
int BATCH_THRESHOLD = string2int(argv("batch_threshold", "1"));
string WORKER_POOL_ID = argv("worker_pool_id", "default");

string model_props = argv("model_props");

string tick0_date = argv("tick0_date", "2020-03-02");
string school_closure_date = argv("school_closure_date", "2020-03-17");
string stay_at_home_date = argv("stay_at_home_date", "2020-03-21");

string param_code_t = """
from utils import format_parameters
import json

param_str = '%s'
tick0_date = '%s'
school_closure_date = '%s'
stay_at_home_date = '%s'
params = format_parameters(param_str, tick0_date, school_closure_date, stay_at_home_date)

params['output.directory'] = '%s'
params['global.random.seed'] = params["seed"]

if 'instance' in params:
    del params['instance']
if 'replicates' in params:
    del params['replicates']
if "seed" in params:
    del params["seed"]

input_str  = '{}'.format(json.dumps(params))
""";


// // app function used to run the task
// app (file out, file err) run_task_app(file shfile, string task_payload, string output_file) {
//     "bash" shfile task_payload output_file emews_root @stdout=out @stderr=err;
// }

// (string result) run_obj(string task_payload, string instance_id) {
//     string tmp_dir = turbine_output + "/tmp";
//     file out <tmp_dir + "/" + instance_id+"_out.txt">;
//     file err <tmp_dir + "/" + instance_id+"_err.txt">;
//     string output_file = "%s/output_%s.csv" % (tmp_dir, instance_id);
//     (out, err) = run_task_app(model_sh, task_payload, output_file) =>
//     result = trim(read(input(output_file))) =>
//     rm(out) =>
//     rm(err) =>
//     rm(output_file);
// }

(string result) run_obj(string task_payload, string instance_id) {
  string instance = "%s/instances/instance_%s" % (turbine_output, instance_id);
  mkdir(instance) => {
      string p_code = param_code_t % (task_payload, tick0_date, school_closure_date,
                                      stay_at_home_date, instance);
      string json_str = python_persist(p_code, "input_str");
      // printf(json_str);
      @par=procs_per_run covid_model_run(model_props, json_str) =>
      // TODO: instance + counts file, if that's the appropriate one
      result = instance + "/output/counts_r1.csv";
  }
}

(string obj_result) run_task(int task_id, string task_payload) {
    string instance_id = "%i" % (task_id);
    obj_result = run_obj(task_payload, instance_id);
}


run(message msgs[]) {
  // printf("MSGS SIZE: %d", size(msgs));
  foreach msg, i in msgs {
    // printf("Running: %s", msg.payload);
    string result_payload = run_task(msg.eq_task_id, msg.payload);
    eq_task_report(msg.eq_task_id, TASK_TYPE, result_payload);
  }
}


(void v) loop(location querier_loc) {
  for (boolean b = true;
       b;
       b=c)
  {
    message msgs[] = eq_batch_task_query(querier_loc);
    boolean c;
    if (msgs[0].msg_type == "status") {
      if (msgs[0].payload == "EQ_STOP") {
        printf("loop.swift: STOP") =>
          v = propagate() =>
          c = false;
      } else {
        // sleep to give time for Python etc.
        // to flush messages
        sleep(5);
        printf("loop.swift: got %s: exiting!", msgs[0].payload) =>
        v = propagate() =>
        c = false;
      }
    } else {
      run(msgs);
      c = true;
    }
  }
}

(void o) start() {
  location querier_loc = locationFromRank(resident_work_rank);
  eq_init_batch_querier(querier_loc, WORKER_POOL_ID, BATCH_SIZE, BATCH_THRESHOLD, TASK_TYPE) =>
  loop(querier_loc) => {
    eq_stop_batch_querier(querier_loc);
    o = propagate();
  }
}

start() => printf("worker pool: normal exit.");