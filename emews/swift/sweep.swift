import io;
import sys;
import files;
import string;
import python;
import stats;
import R;

string emews_root = getenv("EMEWS_PROJECT_ROOT");
string turbine_output = getenv("TURBINE_OUTPUT");
string r_path = "%s/../script/R" % emews_root;
string r_arg =  "--r_path=%s" % r_path;
string exp_path_arg = "--exp_path=%s/results" % turbine_output;
string tmp_path = "%s/tmp" % turbine_output;
file run_sh = input(("%s/scripts/run.sh" % emews_root));

string upf = argv("f");

string parse_csv = """
import utils

csv_file = '%s'
lines = utils.csv2cmd_lines(csv_file)
joined_lines = '!'.join(lines)
""";

app (file out, file err) run(string input_line)
{
    "bash" run_sh r_path input_line @stdout=out @stderr=err;
}

// call this to create any required directories
app (void o) mkdir(string dirname) {
    "mkdir" "-p" dirname;
}

(string result) stage(string line) {
    string parts[] = split(line, "|");
    string exp_id = parts[0];
    string cmd_line = parts[1];
    full_line = "%s %s %s" % (cmd_line, exp_path_arg, r_arg);
    file out <tmp_path + "/" + exp_id + "_out.txt">;
    file err <tmp_path + "/" + exp_id + "_err.txt">;
    (out, err) = run(full_line) => 
    result = "done";
    
}


// Iterate over each line in the upf file, passing each line 
// to the model script to run
main() {
    // printf("workers: %d", turbine_workers());
    string code = parse_csv % upf;
    string joined_lines = python_persist(code, "joined_lines");    
    string upf_lines[] = split(joined_lines, "!");

    string results[];
    foreach line, i in upf_lines {
        results[i] = stage(line);
    }
}
    