#!/usr/bin/env python3
import os
import sys
import time
import subprocess
import json
import re
import argparse
import shutil

# Terminal Colors & Styling
RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"
UNDERLINE = "\033[4m"

# Foregrounds
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
MAGENTA = "\033[35m"
CYAN = "\033[36m"
WHITE = "\033[37m"

# Backgrounds
BG_MAGENTA = "\033[45m"
BG_BLUE = "\033[44m"

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PANGENOME_DIR = os.path.join(BASE_DIR, "pangenome_expansion_highcpu")
PANAROO_OUT = os.path.join(PANGENOME_DIR, "panaroo_out")
CHECKPOINT_FILE = os.path.join(PANGENOME_DIR, "pangenome_checkpoint.json")
LOG_FILE = os.path.join(PANGENOME_DIR, "pangenome_run.log")

def get_log_file_path():
    default_log = os.path.join(PANGENOME_DIR, "pangenome_run.log")
    if os.path.exists(default_log) and os.path.getsize(default_log) > 100:
        return default_log
    try:
        brain_dir = "/home/hp/.gemini/antigravity-ide/brain"
        if os.path.exists(brain_dir):
            subdirs = [os.path.join(brain_dir, d) for d in os.listdir(brain_dir)]
            subdirs = [d for d in subdirs if os.path.isdir(d)]
            if subdirs:
                latest_conv = max(subdirs, key=os.path.getmtime)
                tasks_dir = os.path.join(latest_conv, ".system_generated", "tasks")
                if os.path.exists(tasks_dir):
                    log_files = [os.path.join(tasks_dir, f) for f in os.listdir(tasks_dir) if f.startswith("task-") and f.endswith(".log")]
                    if log_files:
                        return max(log_files, key=os.path.getmtime)
    except:
        pass
    return default_log

TOTAL_GENOMES = 2000

def get_running_total_genomes():
    try:
        cmd = ["ps", "-eo", "args"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            for line in res.stdout.splitlines():
                if "run_pangenome_expansion_highcpu.py" in line:
                    match = re.search(r"--total-genomes\s+(\d+)", line)
                    if match:
                        return int(match.group(1))
    except:
        pass
    return None

def format_size(bytes_val):
    if bytes_val < 1024:
        return f"{bytes_val:.0f} B"
    elif bytes_val < 1024 * 1024:
        return f"{bytes_val/1024:.2f} KB"
    elif bytes_val < 1024 * 1024 * 1024:
        return f"{bytes_val/1024/1024:.2f} MB"
    else:
        return f"{bytes_val/1024/1024/1024:.2f} GB"

def get_cpu_usage(interval=0.1):
    try:
        with open('/proc/stat', 'r') as f:
            line1 = f.readline()
        time.sleep(interval)
        with open('/proc/stat', 'r') as f:
            line2 = f.readline()
        
        def parse_stat(line):
            fields = [float(x) for x in line.strip().split()[1:]]
            idle = fields[3] + fields[4]
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
                mem_total = int(line.split()[1])
            elif line.startswith('MemAvailable:'):
                mem_avail = int(line.split()[1])
        if mem_total > 0:
            mem_used = mem_total - mem_avail
            return mem_used / 1024 / 1024, mem_total / 1024 / 1024
    except:
        pass
    return 0.0, 0.0

def get_disk_usage(path=BASE_DIR):
    try:
        total, used, free = shutil.disk_usage(path)
        return total / (1024**3), used / (1024**3), free / (1024**3)
    except:
        return 0.0, 0.0, 0.0

def get_gpu_usage():
    try:
        cmd = ["nvidia-smi", "--query-gpu=utilization.gpu,temperature.gpu,memory.used,memory.total,name", "--format=csv,noheader,nounits"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            parts = [x.strip() for x in res.stdout.strip().split(',')]
            if len(parts) >= 5:
                return float(parts[0]), float(parts[1]), float(parts[2])/1024.0, float(parts[3])/1024.0, parts[4]
    except:
        pass
    return None

def get_pipeline_processes():
    active = []
    try:
        cmd = ["ps", "-eo", "pid,ppid,pcpu,pmem,comm,args"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            lines = res.stdout.strip().split('\n')
            for line in lines[1:]:
                parts = line.strip().split(None, 5)
                if len(parts) >= 6:
                    pid, ppid, pcpu, pmem, comm, args = parts
                    args_lower = args.lower()
                    comm_lower = comm.lower()
                    
                    keywords = ["panaroo", "fasttree", "cd-hit", "mcl", "mafft", "blast", "prank", "clustal", "python"]
                    matched = None
                    for kw in keywords:
                        if kw in comm_lower or kw in args_lower:
                            if kw == "python" and "run_pangenome_expansion" not in args_lower:
                                continue
                            matched = kw
                            break
                    if matched:
                        active.append({
                            'pid': pid,
                            'name': comm,
                            'cpu': float(pcpu),
                            'mem': float(pmem),
                            'args': args,
                            'type': matched
                        })
    except:
        pass
    return active

class FileTracker:
    def __init__(self, filename):
        self.filename = filename
        self.path = os.path.join(PANAROO_OUT, filename)
        self.last_size = os.path.getsize(self.path) if os.path.exists(self.path) else 0
        self.last_time = time.time()
        self.speed = 0.0
        
    def update(self):
        now = time.time()
        if os.path.exists(self.path):
            size = os.path.getsize(self.path)
        else:
            size = 0
            
        dt = now - self.last_time
        if dt > 0:
            self.speed = (size - self.last_size) / dt
            
        self.last_size = size
        self.last_time = now
        
    def get_status_str(self):
        if not os.path.exists(self.path):
            return f"{DIM}Not created yet{RESET}"
        
        size_str = format_size(self.last_size)
        if self.speed > 1024:
            speed_str = f" {GREEN}(+{format_size(self.speed)}/s){RESET}"
        elif self.speed > 0:
            speed_str = f" {GREEN}(+{self.speed:.0f} B/s){RESET}"
        else:
            speed_str = f" {DIM}(static){RESET}"
            
        return f"{BOLD}{size_str}{RESET}{speed_str}"

def read_last_log_lines(n=8):
    log_path = get_log_file_path()
    if not os.path.exists(log_path):
        return ["Log file not created yet."]
    try:
        with open(log_path, 'r') as f:
            lines = f.readlines()
            return [line.strip() for line in lines[-n:]]
    except Exception as e:
        return [f"Error reading log: {e}"]

def get_panaroo_tqdm_progress():
    log_path = get_log_file_path()
    if not os.path.exists(log_path):
        return None
    try:
        with open(log_path, 'r') as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            seek_pos = max(0, size - 4096)
            f.seek(seek_pos)
            chunk = f.read()
            
        # Parse typical tqdm format: " 40%|████      | 100/250 [12:34<18:51,  4.53s/it]"
        pattern = r"([0-9]+)%\s*\|.*\|\s*([0-9]+)/([0-9]+)\s*\[([0-9:]+)<([0-9:]+)"
        matches = list(re.finditer(pattern, chunk))
        if matches:
            last_match = matches[-1]
            pct = int(last_match.group(1))
            curr = int(last_match.group(2))
            total = int(last_match.group(3))
            elapsed = last_match.group(4)
            remaining = last_match.group(5)
            return pct, curr, total, elapsed, remaining
    except:
        pass
    return None

def parse_tqdm_time(time_str):
    try:
        parts = time_str.split(':')
        if len(parts) == 2:
            return int(parts[0]) * 60 + int(parts[1])
        elif len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + int(parts[2])
    except:
        pass
    return None

def format_eta(seconds):
    if seconds is None or seconds < 0:
        return "Calculating..."
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = int(seconds % 60)
    if h > 0:
        return f"{h}h {m}m"
    elif m > 0:
        return f"{m}m {s}s"
    else:
        return f"{s}s"

def get_overall_progress_stats():
    # Count completed genomes and annotations
    genomes_dir = os.path.join(PANGENOME_DIR, "genomes")
    prokka_dir = os.path.join(PANGENOME_DIR, "prokka_out")
    
    downloaded = 0
    if os.path.exists(genomes_dir):
        downloaded = len([f for f in os.listdir(genomes_dir) if f.endswith(".fna")])
        
    annotated = 0
    if os.path.exists(prokka_dir):
        for root, dirs, files in os.walk(prokka_dir):
            for f in files:
                if f.endswith(".gff") and os.path.getsize(os.path.join(root, f)) > 0:
                    annotated += 1
                    
    # Read checkpoint
    ncbi_done = False
    panaroo_done = False
    fasttree_done = False
    packaged_done = False
    
    if os.path.exists(CHECKPOINT_FILE):
        try:
            with open(CHECKPOINT_FILE, 'r') as f:
                state = json.load(f)
                ncbi_done = state.get("ncbi_query_completed", False)
                panaroo_done = state.get("panaroo_completed", False)
                fasttree_done = state.get("fasttree_completed", False)
                packaged_done = state.get("packaged_completed", False)
        except:
            pass
            
    processes = get_pipeline_processes()
    is_running = len(processes) > 0
    
    # Identify status string
    if packaged_done:
        pipeline_status = "Pipeline Completed & Zipped!"
    elif fasttree_done:
        pipeline_status = "Packaging Results..."
    elif panaroo_done:
        pipeline_status = "Running Core Phylogeny (FastTree)" if is_running else "FastTree Paused"
    elif annotated >= TOTAL_GENOMES + 1:
        pipeline_status = "Running Pangenome Analysis (Panaroo)" if is_running else "Panaroo Paused"
    elif downloaded >= TOTAL_GENOMES:
        pipeline_status = "Annotating Genomes (Prokka)" if is_running else "Prokka Paused"
    elif ncbi_done:
        pipeline_status = "Downloading Genomes" if is_running else "Downloads Paused"
    else:
        pipeline_status = "Querying NCBI" if is_running else "NCBI Query Paused"
        
    # Calculate weighted overall progress percentage
    pct = 0.0
    if ncbi_done:
        pct += 2.0
    pct += (min(downloaded, TOTAL_GENOMES) / TOTAL_GENOMES) * 8.0
    pct += (min(annotated, TOTAL_GENOMES + 1) / (TOTAL_GENOMES + 1)) * 40.0
    
    if panaroo_done:
        pct += 40.0
    else:
        prog = get_panaroo_tqdm_progress()
        if prog:
            pct += (prog[0] / 100.0) * 40.0
            
    if fasttree_done:
        pct += 8.0
    else:
        if any(p['type'] == 'fasttree' for p in processes):
            pct += 4.0  # partial FastTree progress
            
    if packaged_done:
        pct += 2.0
        
    pct = min(pct, 100.0)
    
    # Calculate ETAs
    stage_eta = "N/A"
    overall_eta = "N/A"
    
    if is_running:
        if not ncbi_done:
            stage_eta = "~5 mins"
            overall_eta = "~6.5 hours"
        elif downloaded < TOTAL_GENOMES:
            stage_eta = "Calculating..."
            overall_eta = "Calculating..."
        elif annotated < TOTAL_GENOMES + 1:
            stage_eta = "Calculating..."
            overall_eta = "Calculating..."
        elif not panaroo_done:
            prog = get_panaroo_tqdm_progress()
            if prog:
                _, _, _, _, remaining_str = prog
                stage_eta = remaining_str
                rem_sec = parse_tqdm_time(remaining_str)
                if rem_sec:
                    # remaining Panaroo + 2.0h (FastTree) + 5m (Packaging)
                    overall_eta = format_eta(rem_sec + 2 * 3600 + 300)
                else:
                    overall_eta = f"Panaroo ({remaining_str}) + ~2h"
            else:
                last_lines = read_last_log_lines(15)
                if any("Sanitizing for Panaroo" in line for line in last_lines) and not any("Executing Panaroo" in line for line in last_lines):
                    stage_eta = "~3 mins (Sanitizing)"
                    overall_eta = "~6.0 hours"
                else:
                    stage_eta = "CD-HIT (~2-3 hours)"
                    overall_eta = "~6.0 hours"
        elif not fasttree_done:
            stage_eta = "~2.0 hours"
            overall_eta = "~2.0 hours"
        elif not packaged_done:
            stage_eta = "~5 mins"
            overall_eta = "~5 mins"
        else:
            stage_eta = "0s"
            overall_eta = "0s"
    else:
        if packaged_done:
            stage_eta = "0s (Complete)"
            overall_eta = "0s (Complete)"
        else:
            stage_eta = "Paused"
            overall_eta = "Paused"
            
    return downloaded, annotated, pipeline_status, pct, stage_eta, overall_eta, is_running

def render_dashboard(trackers):
    # Update trackers
    for t in trackers.values():
        t.update()
        
    cpu_usage = get_cpu_usage(0.1)
    ram_used, ram_total = get_ram_usage()
    disk_total, disk_used, disk_free = get_disk_usage()
    gpu_info = get_gpu_usage()
    processes = get_pipeline_processes()
    
    downloaded, annotated, pipeline_status, overall_pct, stage_eta, overall_eta, is_running = get_overall_progress_stats()
    
    os.system("clear")
    print(f"{BOLD}{BG_MAGENTA}{WHITE}  PANAROO & FASTTREE REAL-TIME MONITOR  {RESET}")
    print(f"{DIM}Refreshes every 2s. Press Ctrl+C to exit.{RESET}\n")
    
    # Panel 1: System Resources (including Disk and GPU)
    print(f"{BOLD}{CYAN}┌── System Resources ────────────────────────────────────────────────────────┐{RESET}")
    print(f"│  {BOLD}CPU Usage:{RESET} {cpu_usage:5.1f}% ({os.cpu_count()} threads total)")
    print(f"│  {BOLD}RAM Usage:{RESET} {ram_used:5.2f} GB / {ram_total:5.2f} GB")
    
    # Disk Usage
    disk_free_pct = (disk_free / disk_total) * 100.0 if disk_total > 0 else 0
    disk_alert = RED if disk_free < 10.0 else (YELLOW if disk_free < 25.0 else RESET)
    print(f"│  {BOLD}Disk space:{RESET} {disk_used:.1f} GB used / {disk_alert}{disk_free:.1f} GB free{RESET} ({disk_free_pct:.1f}% free on /media/hp/Data)")
    
    # GPU Usage
    if gpu_info:
        gpu_util, gpu_temp, gpu_vram_used, gpu_vram_total, gpu_name = gpu_info
        print(f"│  {BOLD}GPU Util :{RESET} {gpu_util:5.1f}% ({CYAN}{gpu_name}{RESET})")
        print(f"│  {BOLD}GPU Temp :{RESET} {gpu_temp:5.1f}°C  |  {BOLD}VRAM Usage:{RESET} {gpu_vram_used:.2f} GB / {gpu_vram_total:.2f} GB")
    else:
        print(f"│  {BOLD}GPU Util :{RESET} Not detected / Offline")
    print(f"{CYAN}└────────────────────────────────────────────────────────────────────────────┘{RESET}\n")
    
    # Panel 2: Pangenome Progress & ETAs
    print(f"{BOLD}{MAGENTA}┌── Pipeline Progress: {pipeline_status} ──────────────────────────────────────┐{RESET}")
    
    # Draw overall weighted progress bar
    filled_width = int((overall_pct / 100.0) * 30)
    bar = MAGENTA + "█" * filled_width + RESET + DIM + "░" * (30 - filled_width) + RESET
    print(f"│  {BOLD}Overall Progress:{RESET} [{bar}] {BOLD}{overall_pct:.1f}%{RESET}")
    
    # Pipeline stats and ETAs
    print(f"│  {BOLD}Stage Completed :{RESET} NCBI Query (100%), Genomes Download ({downloaded}/{TOTAL_GENOMES}), Annotations ({annotated}/{TOTAL_GENOMES + 1})")
    
    # Dynamic display during Panaroo
    panaroo_prog = ""
    prog = get_panaroo_tqdm_progress()
    if prog:
        pct, curr, total_jobs, elapsed, remaining = prog
        sub_filled = int((pct / 100.0) * 20)
        sub_bar = "█" * sub_filled + "░" * (20 - sub_filled)
        panaroo_prog = f"[{sub_bar}] {pct}% ({curr}/{total_jobs} jobs) | Elapsed: {elapsed}"
    else:
        last_lines = read_last_log_lines(15)
        if any("Sanitizing for Panaroo" in line for line in last_lines) and not any("Executing Panaroo" in line for line in last_lines):
            panaroo_prog = f"Sanitizing {TOTAL_GENOMES} GFF files... (approx. 1-2 mins)"
        elif is_running and "Panaroo" in pipeline_status:
            panaroo_prog = "Initializing CD-HIT / Graph building..."
        else:
            panaroo_prog = "Waiting to resume / paused..."
            
    print(f"│  {BOLD}Panaroo Progress:{RESET} {CYAN}{panaroo_prog}{RESET}")
    print(f"│")
    print(f"│  {BOLD}Stage ETA       :{RESET} {BOLD}{YELLOW}{stage_eta}{RESET}")
    print(f"│  {BOLD}Pipeline ETA    :{RESET} {BOLD}{GREEN}{overall_eta}{RESET}")
    print(f"│")
    
    # Tracked Files
    print(f"│  {BOLD}Pangenome Output Files Status:{RESET}")
    print(f"│    - DNA Sequences (combined_DNA_CDS.fasta)     : {trackers['dna'].get_status_str()}")
    print(f"│    - Protein Sequences (combined_protein.fasta) : {trackers['protein'].get_status_str()}")
    print(f"│    - Gene Data Table (gene_data.csv)            : {trackers['gene_csv'].get_status_str()}")
    print(f"│    - Core Gene Alignment (core_gene_alignment.aln): {trackers['alignment'].get_status_str()}")
    print(f"│    - Reconstructed Tree (core_phylogeny.tree)     : {trackers['tree'].get_status_str()}")
    print(f"{MAGENTA}└────────────────────────────────────────────────────────────────────────────┘{RESET}\n")
    
    # Panel 3: Active Processes
    print(f"{BOLD}{YELLOW}┌── Active Pipeline Processes ────────────────────────────────────────────────┐{RESET}")
    if not processes:
        print(f"│  {DIM}No active Panaroo, FastTree, CD-HIT, MCL, or MAFFT processes found.{RESET}")
    else:
        print(f"│  {BOLD}{'PID':<8}{'Command':<18}{'CPU %':<8}{'Mem %':<8}{'Purpose/Args':<32}{RESET}")
        for proc in sorted(processes, key=lambda x: x['cpu'], reverse=True)[:8]:
            # Label tools nicely
            purpose = proc['name']
            if proc['type'] == 'python' and 'run_pangenome_expansion' in proc['args']:
                purpose = "Pipeline Coordinator"
            elif proc['type'] == 'cd-hit':
                purpose = "CD-HIT Clustering"
            elif proc['type'] == 'mcl':
                purpose = "MCL Graph Partitioning"
            elif proc['type'] == 'mafft':
                purpose = "MAFFT Alignment"
            elif proc['type'] == 'fasttree':
                purpose = "FastTree Tree Construction"
            elif proc['type'] == 'panaroo':
                purpose = "Panaroo Graph Builder"
            print(f"│  {proc['pid']:<8}{proc['name']:<18}{proc['cpu']:<8.1f}{proc['mem']:<8.1f}{purpose:<32}")
    print(f"{YELLOW}└────────────────────────────────────────────────────────────────────────────┘{RESET}\n")
    
    # Panel 4: Recent Log Lines
    log_path = get_log_file_path()
    log_name = os.path.basename(log_path)
    print(f"{BOLD}{WHITE}┌── Recent Pipeline Log Output ({log_name}) ──────────────────────────┐{RESET}")
    log_lines = read_last_log_lines(8)
    for line in log_lines:
        disp_line = line if len(line) < 74 else line[:71] + "..."
        print(f"│  {disp_line}")
    print(f"{WHITE}└────────────────────────────────────────────────────────────────────────────┘{RESET}")

def print_summary_snapshot():
    cpu_usage = get_cpu_usage(0.1)
    ram_used, ram_total = get_ram_usage()
    disk_total, disk_used, disk_free = get_disk_usage()
    gpu_info = get_gpu_usage()
    
    downloaded, annotated, pipeline_status, overall_pct, stage_eta, overall_eta, is_running = get_overall_progress_stats()
    
    # File sizes
    def get_file_size_str(filename):
        path = os.path.join(PANAROO_OUT, filename)
        if os.path.exists(path):
            return format_size(os.path.getsize(path))
        return "Not created yet"
        
    print("======================================================================")
    print("           PANGENOME PIPELINE MONITOR - STATUS SUMMARY REPORT         ")
    print("======================================================================")
    print(f"Timestamp        : {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    print(f"Pipeline Status  : {pipeline_status} ({'Running' if is_running else 'Paused/Stopped'})")
    print(f"Overall Progress : {overall_pct:.1f}%")
    print(f"Current Stage ETA: {stage_eta}")
    print(f"Pipeline ETA     : {overall_eta}")
    print("----------------------------------------------------------------------")
    print("System Resources:")
    print(f"  CPU Usage: {cpu_usage:.1f}% ({os.cpu_count()} threads total)")
    print(f"  RAM Usage: {ram_used:.2f} GB / {ram_total:.2f} GB")
    disk_free_pct = (disk_free / disk_total) * 100.0 if disk_total > 0 else 0
    print(f"  Disk Usage: {disk_used:.1f} GB used / {disk_free:.1f} GB free ({disk_free_pct:.1f}% free on /media/hp/Data)")
    if gpu_info:
        gpu_util, gpu_temp, gpu_vram_used, gpu_vram_total, gpu_name = gpu_info
        print(f"  GPU Usage: {gpu_util:.1f}% utilization ({gpu_name}) | Temp: {gpu_temp:.0f}°C | VRAM: {gpu_vram_used:.2f} GB / {gpu_vram_total:.2f} GB")
    else:
        print("  GPU Usage: Not detected or Offline")
    print("----------------------------------------------------------------------")
    print("Tracked Output Files:")
    print(f"  - DNA Sequences (combined_DNA_CDS.fasta)     : {get_file_size_str('combined_DNA_CDS.fasta')}")
    print(f"  - Protein Sequences (combined_protein.fasta) : {get_file_size_str('combined_protein_CDS.fasta')}")
    print(f"  - Gene Data Table (gene_data.csv)            : {get_file_size_str('gene_data.csv')}")
    print(f"  - Core Gene Alignment (core_gene_alignment.aln): {get_file_size_str('core_gene_alignment.aln')}")
    print(f"  - Reconstructed Tree (core_phylogeny.tree)     : {get_file_size_str('core_phylogeny.tree')}")
    print("======================================================================")

def main():
    global TOTAL_GENOMES
    parser = argparse.ArgumentParser(description="High-Performance Pangenome Pipeline Monitor")
    parser.add_argument("--summary", action="store_true", help="Print a single summary snapshot of progress and resources, then exit")
    parser.add_argument("--total-genomes", type=int, default=None, help="The target number of genomes (overrides auto-detection)")
    args = parser.parse_args()
    
    if args.total_genomes:
        TOTAL_GENOMES = args.total_genomes
    else:
        detected = get_running_total_genomes()
        if detected:
            TOTAL_GENOMES = detected
            
    if args.summary:
        print_summary_snapshot()
        sys.exit(0)
        
    trackers = {
        'dna': FileTracker('combined_DNA_CDS.fasta'),
        'protein': FileTracker('combined_protein_CDS.fasta'),
        'gene_csv': FileTracker('gene_data.csv'),
        'alignment': FileTracker('core_gene_alignment.aln'),
        'tree': FileTracker('core_phylogeny.tree')
    }
    
    try:
        while True:
            render_dashboard(trackers)
            time.sleep(2.0)
    except KeyboardInterrupt:
        print(f"\n{BOLD}Exiting monitor. Pipeline continues running in the background.{RESET}\n")

if __name__ == "__main__":
    main()
