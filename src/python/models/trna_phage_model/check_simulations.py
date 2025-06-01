import pandas as pd
import os
import glob
from collections import defaultdict

# Configuration for both models
CONFIGS = {
    'weighted': {
        'opt': 1.0,
        'nonopt': 0.25
    },
    'dynamic': {
        'pref': 0.7,
        'charge': 0.5,
        'ribospeed': 1.5,
        'trna': 50000
    }
}

def main():
    rerun_commands = defaultdict(list)
    base_dir = os.path.join('data', 'simulation', 'phage')
    
    # Check weighted model directory
    weighted_dir = f"revised_weighted_opt{CONFIGS['weighted']['opt']:.1f}_nonopt{CONFIGS['weighted']['nonopt']}"
    weighted_path = os.path.join(base_dir, weighted_dir)
    
    # Check dynamic model directory
    dynamic_dir = (f"revised_dynamic_pref{CONFIGS['dynamic']['pref']:.1f}_"
                  f"charge{CONFIGS['dynamic']['charge']}_"
                  f"ribospeed{CONFIGS['dynamic']['ribospeed']}")
    dynamic_path = os.path.join(base_dir, dynamic_dir)
    
    for model_type, dir_path in [('weighted', weighted_path), ('dynamic', dynamic_path)]:
        if not os.path.exists(dir_path):
            print(f"Warning: Directory {dir_path} not found")
            continue
            
        for filepath in glob.glob(os.path.join(dir_path, "*.tsv")):
            try:
                filename = os.path.basename(filepath)
                fop_str = filename.split('_fop')[1].split('_')[0]
                
                # Read the file and check max time
                df = pd.read_csv(filepath, sep='\t')
                max_time = round(df['time'].max())
                
                if max_time < 1200:
                    fop = float(fop_str)
                    
                    if model_type == 'weighted':
                        # Extract seed for weighted model
                        seed = int(filename.split('nonopt')[1].split('_')[1])
                        command = (f"python3 src/python/models/weighted_model/weighted_model.py "
                                 f"{fop} {CONFIGS['weighted']['opt']} {CONFIGS['weighted']['nonopt']} {seed}")
                        screen_name = f"sim_weighted_opt{int(CONFIGS['weighted']['opt']*10)}_fop{fop:.1f}_seed{seed}"
                    else:
                        # Extract seed for dynamic model
                        seed = int(filename.split('pref')[1].split('_')[1])
                        command = (f"python3 src/python/models/trna_phage_model/trna_phage_model.py "
                                 f"{fop} {CONFIGS['dynamic']['charge']} {CONFIGS['dynamic']['pref']} "
                                 f"{seed} {CONFIGS['dynamic']['ribospeed']} {CONFIGS['dynamic']['trna']}")
                        screen_name = f"sim_dynamic_pref{int(CONFIGS['dynamic']['pref']*10)}_fop{fop:.1f}_seed{seed}"
                    
                    screen_name = screen_name.replace('.', '_')
                    rerun_commands[os.path.basename(dir_path)].append((screen_name, command, filename, max_time))
            
            except Exception as e:
                print(f"Error processing {filepath}: {str(e)}")
    
    # Write rerun script
    script_path = os.path.join(os.path.dirname(__file__), 'rerun_short_sims.sh')
    with open(script_path, 'w') as f:
        f.write("#!/bin/bash\n\n")
        
        f.write("# To kill all simulation screens:\n")
        f.write("#   screen -ls | grep sim_ | cut -d. -f1 | awk '{print $1}' | xargs -I {} screen -X -S {} quit\n\n")
        
        for directory in sorted(rerun_commands.keys()):
            f.write(f"# Directory: {directory}\n")
            for screen_name, cmd, filename, time in sorted(rerun_commands[directory]):
                f.write(f"# File: {filename} (reached {time}s)\n")
                f.write(f"screen -dmS {screen_name} bash -c \"nice -n 30 {cmd} && exit\"\n")
                f.write(f"echo \"Started {screen_name}\"\n\n")
            f.write("\n")
    
    os.chmod(script_path, 0o755)
    print(f"Created rerun_short_sims.sh with {sum(len(cmds) for cmds in rerun_commands.values())} commands")

if __name__ == "__main__":
    main()
