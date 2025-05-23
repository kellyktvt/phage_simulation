import pandas as pd
import os
import glob
from collections import defaultdict

# Configuration
CONFIGS = [
    {'pref': 0.5, 'charge': 25.0, 'ribospeed': 75.0, 'trna': 1000},
    {'pref': 0.5, 'charge': 2.5, 'ribospeed': 7.5, 'trna': 10000},
    {'pref': 0.5, 'charge': 0.25, 'ribospeed': 0.75, 'trna': 100000},
    {'pref': 0.7, 'charge': 25.0, 'ribospeed': 75.0, 'trna': 1000},
    {'pref': 0.7, 'charge': 2.5, 'ribospeed': 7.5, 'trna': 10000},
    {'pref': 0.7, 'charge': 0.25, 'ribospeed': 0.75, 'trna': 100000},
]

def format_number(value):
    """Format numbers with appropriate decimal places"""
    if value < 1:
        return f"{value:.2f}"
    return f"{value:.1f}"

def is_one_decimal_fop(fop_str):
    """Check if FOP is in the set: 0.0, 0.1, 0.2, ..., 1.0"""
    try:
        fop = float(fop_str)
        valid_fops = [i/10 for i in range(11)]  # [0.0, 0.1, 0.2, ..., 1.0]
        return any(abs(fop - valid) < 0.001 for valid in valid_fops)
    except:
        return False

def main():
    rerun_commands = defaultdict(list)
    base_dir = os.path.join('data', 'simulation', 'phage')
    
    for config in CONFIGS:
        pref = f"{config['pref']:.1f}"
        charge = format_number(config['charge'])
        ribospeed = format_number(config['ribospeed'])
        
        directory = f"revised_dynamic_pref{pref}_charge{charge}_ribospeed{ribospeed}"
        dir_path = os.path.join(base_dir, directory)
        
        if not os.path.exists(dir_path):
            print(f"Warning: Directory {dir_path} not found")
            continue
            
        for filepath in glob.glob(os.path.join(dir_path, "*.tsv")):
            try:
                filename = os.path.basename(filepath)
                fop_str = filename.split('_fop')[1].split('_')[0]
                
                if not is_one_decimal_fop(fop_str):
                    continue
                
                max_time = round(pd.read_csv(filepath, sep='\t')['time'].max())
                
                if max_time < 1000:
                    fop = float(fop_str)
                    seed = int(filename.split('pref')[1].split('_')[1])
                    
                    command = (f"python3 src/python/models/trna_phage_model/trna_phage_model.py "
                             f"{fop} {config['charge']} {config['pref']} "
                             f"{seed} {config['ribospeed']} {config['trna']}")
                    
                    screen_name = f"sim_pref{int(float(pref)*10)}_fop{fop:.1f}_seed{seed}_trna{config['trna']}"
                    screen_name = screen_name.replace('.', '_')
                    
                    rerun_commands[directory].append((screen_name, command, filename, max_time))
            except Exception as e:
                print(f"Error processing {filepath}: {str(e)}")
    
    script_path = os.path.join(os.path.dirname(__file__), 'rerun_short_sims.sh')
    with open(script_path, 'w') as f:
        f.write("#!/bin/bash\n\n")
        
        # Add usage instructions as comments
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
