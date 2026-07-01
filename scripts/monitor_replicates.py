#!/usr/bin/env python3
import os
import sys
import time
import re
import datetime
import subprocess

# Styling
RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
CYAN = "\033[36m"
BLUE = "\033[34m"
WHITE = "\033[37m"
BG_BLUE = "\033[44m"

def get_gpu_usage():
    try:
        cmd = ["nvidia-smi", "--query-gpu=utilization.gpu,temperature.gpu,memory.used,memory.total", "--format=csv,noheader,nounits"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            parts = [x.strip() for x in res.stdout.strip().split(',')]
            return {
                "util": parts[0],
                "temp": parts[1],
                "used": parts[2],
                "total": parts[3]
            }
    except:
        pass
    return None

def parse_replicate_progress(progress_file):
    if not os.path.exists(progress_file):
        return {"status": "Pending", "progress": 0.0, "current_ns": 0.0, "total_ns": 0.0, "speed": 0.0, "eta": "N/A"}
        
    try:
        with open(progress_file, "r") as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
            
        if not lines:
            return {"status": "Starting", "progress": 0.0, "current_ns": 0.0, "total_ns": 0.0, "speed": 0.0, "eta": "N/A"}
            
        is_completed = any("completed successfully" in l for l in lines)
        
        # Look for progress lines
        prog_match = None
        for line in reversed(lines):
            prog_match = re.search(r"Progress:\s+([\d\.]+)/([\d\.]+)\s+ns.*Speed:\s+([\d\.]+)\s+ns/day.*ETA:\s+([^\s\|]+)", line)
            if prog_match:
                break
                
        if prog_match:
            curr_ns = float(prog_match.group(1))
            total_ns = float(prog_match.group(2))
            speed = float(prog_match.group(3))
            eta = prog_match.group(4)
            pct = (curr_ns / total_ns) * 100.0 if total_ns > 0 else 0.0
            
            return {
                "status": "Completed" if is_completed else "Running",
                "progress": pct,
                "current_ns": curr_ns,
                "total_ns": total_ns,
                "speed": speed,
                "eta": "Done" if is_completed else eta
            }
            
        # Search for initialization
        for line in reversed(lines):
            if "NPT Equilibration starting" in line:
                return {"status": "NPT Equilibration", "progress": 0.0, "current_ns": 0.0, "total_ns": 0.0, "speed": 0.0, "eta": "Equilibrating..."}
            if "NVT Equilibration starting" in line:
                return {"status": "NVT Equilibration", "progress": 0.0, "current_ns": 0.0, "total_ns": 0.0, "speed": 0.0, "eta": "Equilibrating..."}
            if "Energy Minimization starting" in line:
                return {"status": "Energy Minimization", "progress": 0.0, "current_ns": 0.0, "total_ns": 0.0, "speed": 0.0, "eta": "Minimizing..."}
                
        return {"status": "Starting", "progress": 0.0, "current_ns": 0.0, "total_ns": 0.0, "speed": 0.0, "eta": "N/A"}
    except Exception as e:
        return {"status": f"Error: {e}", "progress": 0.0, "current_ns": 0.0, "total_ns": 0.0, "speed": 0.0, "eta": "N/A"}

def draw_bar(pct, width=20):
    filled = int((pct / 100.0) * width)
    filled = min(max(filled, 0), width)
    return "█" * filled + "░" * (width - filled)

def monitor_dir(base_dir):
    # Find all progress files matching md_replicate_seed_*_progress.log
    # or progress log files in base_dir
    files = []
    if os.path.exists(base_dir):
        for f in os.listdir(base_dir):
            if f.endswith("_progress.log") or (f.startswith("md_progress_seed_") and f.endswith(".log")):
                files.append(os.path.join(base_dir, f))
    return sorted(files)

def main():
    parser = argparse.ArgumentParser(description="Replicate MD Simulations Progress Monitor")
    parser.add_argument("--dir", default="results/step9_md", help="Directory where replicate runs are stored")
    args = parser.parse_args()
    
    try:
        while True:
            # Clear screen
            print("\033[2J\033[H", end="")
            
            print(f"{BOLD}{BG_BLUE}{WHITE}  REPLICATE MD SIMULATIONS REAL-TIME MONITOR  {RESET}")
            print(f"{DIM}Refreshes every 5 seconds. Press Ctrl+C to exit.{RESET}\n")
            
            # GPU status
            gpu = get_gpu_usage()
            if gpu:
                print(f"  {BOLD}GPU Load:{RESET} {gpu['util']}%  |  {BOLD}Temp:{RESET} {gpu['temp']}°C  |  {BOLD}VRAM:{RESET} {float(gpu['used'])/1024:.2f} / {float(gpu['total'])/1024:.2f} GB\n")
            
            progress_files = monitor_dir(args.dir)
            
            if not progress_files:
                # Fallback check for any log files in subdirectories
                if os.path.exists(args.dir):
                    for root, dirs, files in os.walk(args.dir):
                        for f in files:
                            if f.endswith("_progress.log") or "progress" in f.lower() and f.endswith(".log"):
                                progress_files.append(os.path.join(root, f))
                                
            progress_files = sorted(list(set(progress_files)))
            
            if not progress_files:
                print(f"  {DIM}No active replicate progress files found in {args.dir}.{RESET}")
                print(f"  {DIM}Ensure scripts are running and outputting logs to results/step9_md/.{RESET}")
            else:
                print(f"  {BOLD}{'Replicate Log File':<35}{'Status':<18}{'Progress Bar':<22}{'Time (ns)':<15}{'Speed':<12}{'ETA':<12}{RESET}")
                print(f"  " + "-"*105)
                for pf in progress_files:
                    name = os.path.basename(pf)
                    stats = parse_replicate_progress(pf)
                    
                    status_color = GREEN if stats["status"] == "Completed" else (YELLOW if stats["status"] == "Running" else CYAN)
                    status_str = f"{status_color}{stats['status']}{RESET}"
                    
                    bar_str = f"[{draw_bar(stats['progress'])}]"
                    time_str = f"{stats['current_ns']:.1f} / {stats['total_ns']:.1f} ns"
                    speed_str = f"{stats['speed']:.1f} ns/d" if stats['speed'] > 0 else "N/A"
                    
                    print(f"  {name:<35}{status_str:<18}{bar_str:<22}{time_str:<15}{speed_str:<12}{stats['eta']:<12}")
                    
            time.sleep(5.0)
    except KeyboardInterrupt:
        print(f"\n{BOLD}Exiting monitor. Replicates continue running in the background.{RESET}\n")

if __name__ == "__main__":
    main()
