#!/usr/bin/env python3
import os
import sys
import time
import subprocess
import json

# Terminal Colors & Styling
CLEAR_SCREEN = "\033[2J\033[H"
RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"
UNDERLINE = "\033[4m"

# Foregrounds
BLACK = "\033[30m"
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
MAGENTA = "\033[35m"
CYAN = "\033[36m"
WHITE = "\033[37m"

# Backgrounds
BG_CYAN = "\033[46m"
BG_BLUE = "\033[44m"

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MD_LOG = os.path.join(BASE_DIR, "md_apo_log.csv")
PANGENOME_DIR = os.path.join(BASE_DIR, "pangenome_expansion_highcpu")
GENOMES_DIR = os.path.join(PANGENOME_DIR, "genomes")
PROKKA_DIR = os.path.join(PANGENOME_DIR, "prokka_out")
CHECKPOINT_FILE = os.path.join(PANGENOME_DIR, "pangenome_checkpoint.json")
FASTTREE_JSON = os.path.join(PANGENOME_DIR, "fasttree_progress.json")

TOTAL_MD_STEPS = 50000000
TOTAL_GENOMES = 2000

def get_cpu_usage(interval=0.2):
    try:
        with open('/proc/stat', 'r') as f:
            line1 = f.readline()
        time.sleep(interval)
        with open('/proc/stat', 'r') as f:
            line2 = f.readline()
        
        def parse_stat(line):
            fields = [float(x) for x in line.strip().split()[1:]]
            idle = fields[3] + fields[4] # idle + iowait
            total = sum(fields)
            return idle, total
        
        idle1, total1 = parse_stat(line1)
        idle2, total2 = parse_stat(line2)
        
        diff_total = total2 - total1
        diff_idle = idle2 - idle1
        if diff_total > 0:
            return 100.0 * (1.0 - diff_idle / diff_total)
    except:
        pass
    return 0.0

def get_ram_usage():
    try:
        with open('/proc/meminfo', 'r') as f:
            lines = f.readlines()
        mem_total = 0
        mem_avail = 0
        for line in lines:
            if line.startswith('MemTotal:'):
                mem_total = int(line.split()[1]) # in kB
            elif line.startswith('MemAvailable:'):
                mem_avail = int(line.split()[1]) # in kB
        if mem_total > 0:
            mem_used = mem_total - mem_avail
            return mem_used / 1024 / 1024, mem_total / 1024 / 1024 # in GB
    except:
        pass
    return 0.0, 0.0

def get_gpu_usage():
    try:
        cmd = ["nvidia-smi", "--query-gpu=utilization.gpu,temperature.gpu,memory.used,memory.total", "--format=csv,noheader,nounits"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            parts = [x.strip() for x in res.stdout.strip().split(',')]
            if len(parts) >= 4:
                return float(parts[0]), float(parts[1]), float(parts[2])/1024.0, float(parts[3])/1024.0
    except:
        pass
    return None

def is_analysis_active(out_dir):
    try:
        cmd = ["ps", "-eo", "args"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            for line in res.stdout.splitlines():
                if "analyze_md_trajectory.py" in line and os.path.basename(out_dir) in line:
                    return True
    except:
        pass
    return False

def get_md_target_status(name, out_dir, progress_file):
    log_file = os.path.join(out_dir, "md_apo_log.csv")
    TOTAL_STEPS = 50000000
    
    if not os.path.exists(out_dir):
        return {"status": "Pending", "progress": 0.0, "current_ns": 0.0, "total_ns": 100.0, "speed": 0.0, "temp": 0.0, "eta": "Waiting..."}
        
    is_completed = False
    if os.path.exists(progress_file):
        try:
            with open(progress_file, 'r') as f:
                content = f.read()
                if "MD Simulation completed successfully." in content:
                    parts = content.split("MD Simulation completed successfully.")
                    last_part = parts[-1]
                    if "Resuming simulation" not in last_part and "Progress:" not in last_part:
                        is_completed = True
        except:
            pass

    temp = 0.0
    if os.path.exists(log_file):
        try:
            with open(log_file, 'r') as f:
                lines = [line.strip() for line in f.readlines() if line.strip() and not line.startswith('#')]
            if lines:
                last_line = lines[-1]
                parts = last_line.split(',')
                if len(parts) >= 7:
                    step = int(parts[0].replace('"', ''))
                    temp = float(parts[4])
                    if step >= TOTAL_STEPS:
                        is_completed = True
        except:
            pass

    if is_completed:
        timeseries_csv = os.path.join(out_dir, "md_analysis_timeseries.csv")
        if os.path.exists(timeseries_csv) and os.path.getsize(timeseries_csv) > 0:
            return {"status": "Completed", "progress": 100.0, "current_ns": 100.0, "total_ns": 100.0, "speed": 0.0, "temp": temp, "eta": "Done"}
        elif is_analysis_active(out_dir):
            return {"status": "Analyzing Trajectory", "progress": 100.0, "current_ns": 100.0, "total_ns": 100.0, "speed": 0.0, "temp": temp, "eta": "Analyzing..."}
        else:
            return {"status": "Awaiting Analysis", "progress": 100.0, "current_ns": 100.0, "total_ns": 100.0, "speed": 0.0, "temp": temp, "eta": "Queued..."}

    if os.path.exists(log_file):
        try:
            with open(log_file, 'r') as f:
                lines = [line.strip() for line in f.readlines() if line.strip() and not line.startswith('#')]
            if lines:
                last_line = lines[-1]
                parts = last_line.split(',')
                if len(parts) >= 7:
                    step = int(parts[0].replace('"', ''))
                    temp = float(parts[4])
                    speed = float(parts[6])
                    
                    pct = (step / TOTAL_STEPS) * 100.0
                    curr_ns = step * 2.0 * 1e-6
                    
                    eta_str = "Calculating..."
                    if speed > 0:
                        remaining_ns = 100.0 - curr_ns
                        eta_sec = (remaining_ns / speed) * 86400.0
                        eta_str = format_eta(eta_sec)
                        
                    return {
                        "status": "Running",
                        "progress": pct,
                        "current_ns": curr_ns,
                        "total_ns": 100.0,
                        "speed": speed,
                        "temp": temp,
                        "eta": eta_str
                    }
        except Exception as e:
            pass

    if os.path.exists(progress_file):
        try:
            with open(progress_file, 'r') as f:
                lines = [line.strip() for line in f.readlines() if line.strip()]
            if lines:
                last_line = lines[-1]
                status = "Preparing"
                if "NPT Equilibration starting" in last_line:
                    status = "NPT Equilibration"
                elif "NVT Equilibration starting" in last_line:
                    status = "NVT Equilibration"
                elif "Energy Minimization starting" in last_line:
                    status = "Energy Minimization"
                elif "Solvation completed" in last_line or "Solvating" in last_line:
                    status = "System Solvation"
                elif "System Preparation starting" in last_line:
                    status = "System Preparation"
                elif "Resuming simulation" in last_line:
                    status = "Starting Production"
                elif "Equilibration completed" in last_line:
                    status = "Equilibration Done"
                    
                return {"status": status, "progress": 0.0, "current_ns": 0.0, "total_ns": 100.0, "speed": 0.0, "temp": 0.0, "eta": "Equilibrating..."}
        except:
            pass
            
    return {"status": "Pending", "progress": 0.0, "current_ns": 0.0, "total_ns": 100.0, "speed": 0.0, "temp": 0.0, "eta": "Waiting..."}

def get_trajectory_stats(out_dir):
    csv_file = os.path.join(out_dir, "md_analysis_timeseries.csv")
    if os.path.exists(csv_file):
        try:
            with open(csv_file, 'r') as f:
                lines = f.readlines()
            if len(lines) > 1:
                header = lines[0].strip().split(',')
                rmsd_idx = -1
                time_idx = -1
                for idx, col in enumerate(header):
                    if 'rmsd' in col.lower():
                        rmsd_idx = idx
                    elif 'time' in col.lower():
                        time_idx = idx
                if rmsd_idx != -1:
                    rmsds = []
                    times = []
                    for line in lines[1:]:
                        parts = line.strip().split(',')
                        if len(parts) > rmsd_idx:
                            try:
                                rmsds.append(float(parts[rmsd_idx]))
                                if time_idx != -1:
                                    times.append(float(parts[time_idx]))
                            except:
                                pass
                    if rmsds:
                        import numpy as np
                        mean_rmsd = np.mean(rmsds)
                        std_rmsd = np.std(rmsds)
                        total_time = times[-1] if times else len(rmsds) * 0.01
                        return f"{total_time:.2f} ns finished; backbone RMSD = {mean_rmsd:.2f} ± {std_rmsd:.2f} Å"
        except:
            pass
    return None

def count_pangenome_progress():
    # Count downloaded genomes
    downloaded = 0
    if os.path.exists(GENOMES_DIR):
        downloaded = len([f for f in os.listdir(GENOMES_DIR) if f.endswith(".fna")])
        
    # Count annotated genomes
    annotated = 0
    if os.path.exists(PROKKA_DIR):
        for root, dirs, files in os.walk(PROKKA_DIR):
            for f in files:
                if f.endswith(".gff") and os.path.getsize(os.path.join(root, f)) > 0:
                    annotated += 1
                    
    # Read checkpoint for state status
    state_desc = "Initializing/Downloading"
    if os.path.exists(CHECKPOINT_FILE):
        try:
            with open(CHECKPOINT_FILE, 'r') as f:
                state = json.load(f)
                if state.get("packaged_completed"):
                    state_desc = "Packaging Completed!"
                elif state.get("fasttree_completed"):
                    state_desc = "Reconstructing Core Phylogeny (FastTree)"
                elif state.get("panaroo_completed"):
                    state_desc = "Running Core Phylogeny (FastTree)"
                elif len(state.get("completed_annotations", [])) >= TOTAL_GENOMES:
                    state_desc = "Running Pangenome Analysis (Panaroo)"
                elif len(state.get("completed_downloads", [])) >= TOTAL_GENOMES:
                    state_desc = "Annotating Genomes (Prokka)"
        except:
            pass
            
    return downloaded, annotated, state_desc

def draw_progress_bar(val, total, width=30, color=CYAN):
    ratio = val / total if total > 0 else 0.0
    ratio = min(max(ratio, 0.0), 1.0)
    filled_width = int(ratio * width)
    bar = color + "█" * filled_width + RESET + DIM + "░" * (width - filled_width) + RESET
    pct = f"{ratio * 100.0:.1f}%"
    return f"[{bar}] {BOLD}{pct}{RESET}"

def format_eta(seconds):
    if seconds is None or seconds < 0:
        return "Calculating..."
    if seconds > 86400:
        return f"{seconds / 86400:.1f} days"
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = int(seconds % 60)
    return f"{h:02d}:{m:02d}:{s:02d}"

STATE_TRACKER = {
    'start_time': None,
    'initial_downloaded': None,
    'initial_annotated': None
}

def display_dashboard():
    global STATE_TRACKER
    # Gather System Info
    cpu_usage = get_cpu_usage(0.1)
    ram_used, ram_total = get_ram_usage()
    gpu_info = get_gpu_usage()
    
    # Gather MD Info
    # parsed dynamically inside Panel 3
    
    # Gather Pangenome Info
    pg_downloaded, pg_annotated, pg_state = count_pangenome_progress()
    
    # Initialize trackers if None
    if STATE_TRACKER['start_time'] is None:
        STATE_TRACKER['start_time'] = time.time()
    if STATE_TRACKER['initial_downloaded'] is None:
        STATE_TRACKER['initial_downloaded'] = pg_downloaded
    if STATE_TRACKER['initial_annotated'] is None:
        STATE_TRACKER['initial_annotated'] = pg_annotated
        
    elapsed = time.time() - STATE_TRACKER['start_time']
    
    # Calculate Download ETA & Speed
    download_eta = "Calculating..."
    download_speed_str = ""
    downloaded_diff = pg_downloaded - STATE_TRACKER['initial_downloaded']
    if downloaded_diff > 0 and elapsed > 5:
        speed = downloaded_diff / elapsed  # files/sec
        download_speed_str = f" ({speed:.2f} genomes/s)"
        remaining = TOTAL_GENOMES - pg_downloaded
        if speed > 0:
            download_eta = format_eta(remaining / speed)
            
    # Calculate Annotation ETA & Speed
    annotation_eta = "Calculating..."
    annotation_speed_str = ""
    annotated_diff = pg_annotated - STATE_TRACKER['initial_annotated']
    if annotated_diff > 0 and elapsed > 5:
        speed = annotated_diff / elapsed  # files/sec
        annotation_speed_str = f" ({speed*60:.1f} genomes/min)"
        remaining = TOTAL_GENOMES - pg_annotated
        if speed > 0:
            annotation_eta = format_eta(remaining / speed)
            
    # Render screen
    os.system("clear")
    print(f"{BOLD}{BG_BLUE}{WHITE}  WORKSTATION PIPELINE DASHBOARD (Real-Time)  {RESET}")
    print(f"{DIM}Press Ctrl+C to exit dashboard.{RESET}")
    print()
    
    # Panel 1: System Metrics
    print(f"{BOLD}{BLUE}┌── System Resources ────────────────────────────────────────────────────────┐{RESET}")
    print(f"│  {BOLD}CPU Usage:{RESET} {cpu_usage:5.1f}% ({os.cpu_count()} threads total)")
    print(f"│  {BOLD}RAM Usage:{RESET} {ram_used:5.2f} GB / {ram_total:5.2f} GB")
    if gpu_info:
        gpu_util, gpu_temp, gpu_vram_used, gpu_vram_total = gpu_info
        print(f"│  {BOLD}GPU Util :{RESET} {gpu_util:5.1f}% ({CYAN}RTX 5070 Ti{RESET})")
        print(f"│  {BOLD}GPU Temp :{RESET} {gpu_temp:5.1f}°C  |  {BOLD}VRAM Usage:{RESET} {gpu_vram_used:.2f} GB / {gpu_vram_total:.2f} GB")
    else:
        print(f"│  {BOLD}GPU Util :{RESET} Not detected / Offline")
    print(f"{BLUE}└────────────────────────────────────────────────────────────────────────────┘{RESET}")
    print()
    
    # Panel 2: Pangenome Expansion (CPU-bound)
    print(f"{BOLD}{MAGENTA}┌── Process: Pangenome Expansion ────────────────────────────────────────────┐{RESET}")
    print(f"│  {BOLD}Pipeline Status :{RESET} {BOLD}{CYAN}{pg_state}{RESET}")
    print(f"│  {BOLD}Downloads       :{RESET} {pg_downloaded:,} / {TOTAL_GENOMES:,} genomes{download_speed_str}")
    print(f"│  {BOLD}Download Prog   :{RESET} {draw_progress_bar(pg_downloaded, TOTAL_GENOMES, color=MAGENTA)}  |  {BOLD}ETA:{RESET} {YELLOW}{download_eta}{RESET}")
    print(f"│  {BOLD}Annotations     :{RESET} {pg_annotated:,} / {TOTAL_GENOMES:,} annotated (Prokka){annotation_speed_str}")
    print(f"│  {BOLD}Annotation Prog :{RESET} {draw_progress_bar(pg_annotated, TOTAL_GENOMES, color=BLUE)}  |  {BOLD}ETA:{RESET} {YELLOW}{annotation_eta}{RESET}")
    
    # Try reading FastTree progress
    fasttree_status = None
    if os.path.exists(FASTTREE_JSON):
        try:
            with open(FASTTREE_JSON, 'r') as f:
                fasttree_status = json.load(f)
        except:
            pass
            
    if fasttree_status and fasttree_status.get("stage") != "completed":
        stage = fasttree_status.get("stage", "Unknown")
        pct = fasttree_status.get("progress_percent", 0.0)
        est_mem = fasttree_status.get("current_est_memory_gb", 0.0)
        eta_sec = fasttree_status.get("eta_seconds", -1.0)
        msg = fasttree_status.get("message", "")
        
        eta_str = format_eta(eta_sec) if eta_sec > 0 else "Calculating..."
        
        print(f"├────────────────────────────────────────────────────────────────────────────┤")
        print(f"│  {BOLD}{UNDERLINE}FastTree Reconstruct Progress{RESET}")
        print(f"│  {BOLD}Active Stage    :{RESET} {YELLOW}{stage}{RESET}")
        print(f"│  {BOLD}Stage Progress  :{RESET} {draw_progress_bar(pct, 100.0, color=CYAN)}")
        print(f"│  {BOLD}Est. Memory     :{RESET} {CYAN}{est_mem:.2f} GB{RESET} / stage peak")
        print(f"│  {BOLD}Stage ETA       :{RESET} {BOLD}{YELLOW}{eta_str}{RESET}")
        print(f"│  {BOLD}Current Status  :{RESET} {DIM}{msg[:60]}{RESET}")
        
    print(f"{MAGENTA}└────────────────────────────────────────────────────────────────────────────┘{RESET}")
    print()
    
    # Panel 3: Molecular Dynamics (MD) Simulations (GPU-bound)
    print(f"{BOLD}{CYAN}┌── Process: Molecular Dynamics Simulations (GPU-bound) ─────────────────────┐{RESET}")
    
    # Live info for GNAT_02769
    gnat_stats = get_trajectory_stats(BASE_DIR)
    gnat_info = f" ({gnat_stats})" if gnat_stats else " (127.78 ns finished; backbone RMSD = 2.73 ± 0.20 Å)"
    print(f"│  {BOLD}1. GNAT_02769{RESET}          : {GREEN}Completed{gnat_info}{RESET}")
    
    # Live info for Enterohemolysin_00061
    status_00061 = get_md_target_status("Enterohemolysin_00061", os.path.join(BASE_DIR, "md_00061"), os.path.join(BASE_DIR, "00061_md_progress.txt"))
    status_str_00061 = f"{BOLD}{YELLOW}{status_00061['status']}{RESET}"
    if status_00061['status'] == "Running":
        status_str_00061 = f"{BOLD}{GREEN}Running{RESET} ({status_00061['current_ns']:.2f} / {status_00061['total_ns']:.2f} ns)"
    elif status_00061['status'] == "Completed":
        stats_00061 = get_trajectory_stats(os.path.join(BASE_DIR, "md_00061"))
        info_00061 = f" ({stats_00061})" if stats_00061 else ""
        status_str_00061 = f"{BOLD}{GREEN}Completed{info_00061}{RESET}"
        
    print(f"│  {BOLD}2. Enterohemolysin_00061{RESET}: {status_str_00061}")
    if status_00061['status'] in ["Running"]:
        print(f"│     Prog: {draw_progress_bar(status_00061['progress'], 100.0, color=GREEN)}  |  Speed: {YELLOW}{status_00061['speed']:.2f} ns/day{RESET}  |  ETA: {YELLOW}{status_00061['eta']}{RESET}")
    elif status_00061['status'] not in ["Pending", "Completed"]:
        print(f"│     State: {YELLOW}{status_00061['status']}{RESET}  |  {status_00061['eta']}")
        
    # Live info for O-antigen Polymerase_03161
    status_03161 = get_md_target_status("O-antigen Polymerase_03161", os.path.join(BASE_DIR, "md_03161"), os.path.join(BASE_DIR, "03161_md_progress.txt"))
    status_str_03161 = f"{BOLD}{YELLOW}{status_03161['status']}{RESET}"
    if status_03161['status'] == "Running":
        status_str_03161 = f"{BOLD}{GREEN}Running{RESET} ({status_03161['current_ns']:.2f} / {status_03161['total_ns']:.2f} ns)"
    elif status_03161['status'] == "Completed":
        stats_03161 = get_trajectory_stats(os.path.join(BASE_DIR, "md_03161"))
        info_03161 = f" ({stats_03161})" if stats_03161 else ""
        status_str_03161 = f"{BOLD}{GREEN}Completed{info_03161}{RESET}"
        
    print(f"│  {BOLD}3. O-antigen Poly_03161{RESET} : {status_str_03161}")
    if status_03161['status'] in ["Running"]:
        print(f"│     Prog: {draw_progress_bar(status_03161['progress'], 100.0, color=GREEN)}  |  Speed: {YELLOW}{status_03161['speed']:.2f} ns/day{RESET}  |  ETA: {YELLOW}{status_03161['eta']}{RESET}")
    elif status_03161['status'] not in ["Pending", "Completed"]:
        print(f"│     State: {YELLOW}{status_03161['status']}{RESET}  |  {status_03161['eta']}")
        
    print(f"{CYAN}└────────────────────────────────────────────────────────────────────────────┘{RESET}")
    print()

def main():
    try:
        # Loop every 2 seconds to refresh progress
        while True:
            display_dashboard()
            time.sleep(2.0)
    except KeyboardInterrupt:
        print(f"\n{BOLD}Exiting dashboard. Process continues in background.{RESET}\n")

if __name__ == "__main__":
    main()
