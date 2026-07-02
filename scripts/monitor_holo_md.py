#!/usr/bin/env python3
import os
import sys
import time
import subprocess
import re
import datetime

# Terminal Colors
RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
MAGENTA = "\033[35m"
CYAN = "\033[36m"
WHITE = "\033[37m"
BG_BLUE = "\033[44m"

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
HOLO_DIR = os.path.join(BASE_DIR, "AMR_Work", "results", "md_simulation", "holo_gnat")
PROGRESS_FILE = os.path.join(HOLO_DIR, "md_holo_progress.log")

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
        mem_total, mem_avail = 0, 0
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

def get_gpu_usage():
    try:
        cmd = ["nvidia-smi", "--query-gpu=utilization.gpu,temperature.gpu,memory.used,memory.total", "--format=csv,noheader,nounits"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            parts = [x.strip() for x in res.stdout.strip().split(',')]
            if len(parts) >= 4:
                return {
                    "util": float(parts[0]),
                    "temp": float(parts[1]),
                    "vram_used": float(parts[2]),
                    "vram_total": float(parts[3])
                }
    except:
        pass
    return None

def parse_progress_file():
    if not os.path.exists(PROGRESS_FILE):
        return {"status": "Pending", "progress": 0.0, "current_ns": 0.0, "total_ns": 20.0, "speed": 0.0, "eta": "Waiting..."}
        
    try:
        with open(PROGRESS_FILE, "r") as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
            
        if not lines:
            return {"status": "Starting", "progress": 0.0, "current_ns": 0.0, "total_ns": 20.0, "speed": 0.0, "eta": "Initializing..."}
            
        # Check if completed
        is_completed = any("Holo-State MD Simulation completed successfully." in l for l in lines)
        
        # Look for a progress line by scanning backwards
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
            pct = (curr_ns / total_ns) * 100.0
            
            return {
                "status": "Completed" if is_completed else "Running",
                "progress": pct,
                "current_ns": curr_ns,
                "total_ns": total_ns,
                "speed": speed,
                "eta": "Done" if is_completed else eta
            }
            
        # If no progress line, search for initialization states
        status = "Preparing"
        for line in reversed(lines):
            if "NPT Equilibration starting" in line:
                status = "NPT Equilibration"
                break
            elif "NVT Equilibration starting" in line:
                status = "NVT Equilibration"
                break
            elif "Energy Minimization starting" in line:
                status = "Energy Minimization"
                break
            elif "Solvation completed" in line or "Solvating" in line:
                status = "System Solvation"
                break
            elif "complex solvation starting" in line:
                status = "Complex Prep"
                break
                
        return {"status": status, "progress": 0.0, "current_ns": 0.0, "total_ns": 20.0, "speed": 0.0, "eta": "Equilibrating..."}
        
    except Exception as e:
        return {"status": f"Error: {e}", "progress": 0.0, "current_ns": 0.0, "total_ns": 20.0, "speed": 0.0, "eta": "Unknown"}

def draw_progress_bar(val, total, width=40, color=GREEN):
    ratio = val / total if total > 0 else 0.0
    ratio = min(max(ratio, 0.0), 1.0)
    filled_width = int(ratio * width)
    bar = color + "█" * filled_width + RESET + DIM + "░" * (width - filled_width) + RESET
    pct = f"{ratio * 100.0:.2f}%"
    return f"[{bar}] {BOLD}{pct}{RESET}"

def display_dashboard():
    cpu_usage = get_cpu_usage(0.1)
    ram_used, ram_total = get_ram_usage()
    gpu = get_gpu_usage()
    md = parse_progress_file()
    
    # Clear screen
    print("\033[2J\033[H", end="")
    
    # Title
    print(f"{BOLD}{BG_BLUE}{WHITE}  EXPERIMENT 3: GNAT_KA27 HOLO-STATE MD MONITOR  {RESET}")
    print(f"{DIM}Press Ctrl+C to exit monitor. Process will continue in the background.{RESET}\n")
    
    # 1. Hardware Status Panel
    print(f"{BOLD}{BLUE}┌── Hardware & Resource Usage ───────────────────────────────────────────────┐{RESET}")
    print(f"│  {BOLD}CPU Usage:{RESET} {cpu_usage:5.1f}% ({os.cpu_count()} threads total)")
    print(f"│  {BOLD}RAM Usage:{RESET} {ram_used:5.2f} GB / {ram_total:5.2f} GB")
    if gpu:
        print(f"│  {BOLD}GPU Util :{RESET} {gpu['util']:5.1f}%  |  {BOLD}GPU Temp:{RESET} {gpu['temp']:4.1f}°C")
        print(f"│  {BOLD}VRAM     :{RESET} {gpu['vram_used']/1024.0:5.2f} GB / {gpu['vram_total']/1024.0:5.2f} GB")
    else:
        print(f"│  {BOLD}GPU Util :{RESET} Not detected / Offline")
    print(f"{BLUE}└────────────────────────────────────────────────────────────────────────────┘{RESET}")
    print()
    
    # 2. Simulation Status Panel
    print(f"{BOLD}{CYAN}┌── Simulation Progress (20.0 ns Holo-State MD) ─────────────────────────────┐{RESET}")
    print(f"│  {BOLD}Status    :{RESET} {BOLD}{GREEN if md['status'] in ['Running','Completed'] else YELLOW}{md['status']}{RESET}")
    print(f"│  {BOLD}Time      :{RESET} {md['current_ns']:.2f} ns / {md['total_ns']:.2f} ns")
    print(f"│  {BOLD}Progress  :{RESET} {draw_progress_bar(md['progress'], 100.0, color=CYAN)}")
    print(f"│  {BOLD}Speed     :{RESET} {YELLOW}{md['speed']:.2f} ns/day{RESET}")
    print(f"│  {BOLD}ETA       :{RESET} {BOLD}{YELLOW}{md['eta']}{RESET}")
    print(f"{CYAN}└────────────────────────────────────────────────────────────────────────────┘{RESET}")
    
    # Update a summary file for tracking
    summary_file = os.path.join(HOLO_DIR, "holo_md_summary.json")
    try:
        import json
        with open(summary_file, "w") as sf:
            json.dump({
                "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "status": md["status"],
                "progress_percent": md["progress"],
                "current_ns": md["current_ns"],
                "total_ns": md["total_ns"],
                "speed_ns_day": md["speed"],
                "eta": md["eta"],
                "gpu_util": gpu["util"] if gpu else 0.0,
                "gpu_temp": gpu["temp"] if gpu else 0.0,
                "vram_gb": (gpu["vram_used"]/1024.0) if gpu else 0.0,
                "cpu_usage": cpu_usage,
                "ram_gb": ram_used
            }, sf, indent=4)
    except:
        pass

def main():
    try:
        while True:
            display_dashboard()
            time.sleep(2.0)
    except KeyboardInterrupt:
        print(f"\n{BOLD}Exiting monitor. Simulation continues in the background.{RESET}\n")

if __name__ == "__main__":
    main()
