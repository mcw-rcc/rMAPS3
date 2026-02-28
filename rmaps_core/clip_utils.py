"""
Shared utilities for CLIP-seq analysis and RNA map generation.

This module provides common functions used across all CLIP-seq event types (A3SS, A5SS, ES/SE, MXE, RI).
"""

from __future__ import annotations

import os
import shutil
import subprocess
import logging
from pathlib import Path
from typing import Dict, Tuple, Optional


def setup_logging(output_path: Path, version: str = "3.0.0") -> None:
    """
    Configure logging for CLIP-seq analysis.
    
    Args:
        output_path: Output directory for log files
        version: Version string for the log filename
    """
    log_file = output_path / f"log.CLIPSeq{version}.txt"
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(message)s',
        filename=str(log_file),
        filemode='w'
    )
    logging.debug('Version: %s', version)


def run_perl_command(cmd: str) -> Tuple[int, str]:
    """
    Run a Perl command and return status and output.
    
    Python 3 replacement for commands.getstatusoutput().
    
    Args:
        cmd: Command string to execute
        
    Returns:
        Tuple of (return_code, output_string)
    """
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            check=False
        )
        return result.returncode, result.stdout + result.stderr
    except Exception as e:
        return -1, str(e)


def copy_file(src: str, dest: str) -> Tuple[int, str]:
    """
    Copy a file from source to destination.
    
    Python 3 replacement for 'cp' command.
    
    Args:
        src: Source file path
        dest: Destination file path
        
    Returns:
        Tuple of (return_code, output_message)
    """
    try:
        shutil.copy2(src, dest)
        return 0, f"Successfully copied {src} to {dest}"
    except Exception as e:
        return 1, f"Error copying file: {str(e)}"


def parse_rMATS_line(
    line: str,
    sig_fdr: float,
    sig_delta_psi: float,
    bg_fdr: float,
    min_psi: float,
    max_psi: float,
    coord_indices: Tuple[int, ...]
) -> Optional[Tuple[str, str, str]]:
    """
    Parse a single rMATS line and categorize it as up/down/background.
    
    Args:
        line: Tab-delimited rMATS line
        sig_fdr: FDR cutoff for significant events
        sig_delta_psi: Delta PSI cutoff for significant events
        bg_fdr: FDR cutoff for background events
        min_psi: Minimum PSI threshold for background
        max_psi: Maximum PSI threshold for background
        coord_indices: Tuple of column indices for coordinates (event-specific)
        
    Returns:
        Tuple of (category, key, coord_string) where category is 'up', 'down', or 'bg'
        Returns None if event doesn't meet any criteria
    """
    ele = line.strip().split('\t')
    if len(ele) < 11:
        return None
        
    # Extract coordinates based on event-specific indices
    try:
        coords = [ele[i] for i in coord_indices]
        key = ':'.join(coords[:4])  # chr:strand:start:end typically
        coord_string = '\t'.join(coords)
    except IndexError:
        return None
    
    try:
        fdr = float(ele[-4])
    except (ValueError, IndexError):
        return None
        
    # Check for significant events
    if fdr < sig_fdr:
        try:
            delta_psi = float(ele[-1])
            if delta_psi >= sig_delta_psi:
                return ('up', key, coord_string)
            elif delta_psi <= -sig_delta_psi:
                return ('down', key, coord_string)
        except (ValueError, IndexError):
            pass
    
    # Check for background events
    elif fdr > bg_fdr:
        # Calculate average PSI for both samples
        try:
            psi_1_str = ele[-3].replace('"', '').split(',')
            psi_1_values = [float(p) for p in psi_1_str if p != "NA"]
            psi_1 = sum(psi_1_values) / len(psi_1_values) if psi_1_values else 0.0
            
            psi_2_str = ele[-2].replace('"', '').split(',')
            psi_2_values = [float(p) for p in psi_2_str if p != "NA"]
            psi_2 = sum(psi_2_values) / len(psi_2_values) if psi_2_values else 0.0
            
            if min(psi_1, psi_2) < min_psi and max(psi_1, psi_2) > max_psi:
                return ('bg', key, coord_string)
        except (ValueError, IndexError):
            pass
    
    return None


def deduplicate_events(
    up: Dict[str, list],
    down: Dict[str, list],
    bg: Dict[str, list]
) -> Tuple[Dict[str, list], Dict[str, list], Dict[str, list]]:
    """
    Remove events that appear in multiple categories.
    
    Marks events that appear in more than one dictionary so they can be excluded.
    
    Args:
        up: Dictionary of upregulated events {key: [count, coord_string]}
        down: Dictionary of downregulated events
        bg: Dictionary of background events
        
    Returns:
        Tuple of (up, down, bg) dictionaries with counts updated
    """
    # Mark duplicates across categories
    for key in up:
        if key in down:
            down[key][0] += 1
            up[key][0] += 1
        if key in bg:
            bg[key][0] += 1
            up[key][0] += 1
    
    for key in down:
        if key in bg:
            bg[key][0] += 1
            down[key][0] += 1
    
    return up, down, bg


def write_coordinate_files(
    up: Dict[str, list],
    down: Dict[str, list],
    bg: Dict[str, list],
    exon_path: Path,
    header: str
) -> Tuple[int, int, int]:
    """
    Write coordinate files for up/down/background events.
    
    Only writes events with count == 1 (i.e., unique to one category).
    
    Args:
        up: Dictionary of upregulated events
        down: Dictionary of downregulated events
        bg: Dictionary of background events
        exon_path: Directory to write coordinate files
        header: Header line for coordinate files
        
    Returns:
        Tuple of (num_up, num_down, num_bg) - counts of written events
    """
    nu = nd = nb = 0
    
    with open(exon_path / 'up.coord.txt', 'w') as u_file:
        u_file.write(header + '\n')
        for key in up:
            if up[key][0] == 1:
                nu += 1
                u_file.write(up[key][1] + '\n')
    
    with open(exon_path / 'dn.coord.txt', 'w') as d_file:
        d_file.write(header + '\n')
        for key in down:
            if down[key][0] == 1:
                nd += 1
                d_file.write(down[key][1] + '\n')
    
    with open(exon_path / 'bg.coord.txt', 'w') as b_file:
        b_file.write(header + '\n')
        for key in bg:
            if bg[key][0] == 1:
                nb += 1
                b_file.write(bg[key][1] + '\n')
    
    return nu, nd, nb


def process_rmats_file(
    rmats_path: str,
    exon_path: Path,
    header: str,
    coord_indices: Tuple[int, ...],
    sig_fdr: float,
    sig_delta_psi: float,
    bg_fdr: float = 0.3,
    min_psi: float = 0.95,
    max_psi: float = 0.05
) -> Tuple[int, int, int]:
    """
    Process rMATS file and generate coordinate files.
    
    Args:
        rmats_path: Path to rMATS input file
        exon_path: Output directory for coordinate files
        header: Header line for output files
        coord_indices: Event-specific coordinate column indices
        sig_fdr: FDR cutoff for significant events
        sig_delta_psi: Delta PSI cutoff for significant events
        bg_fdr: FDR cutoff for background events
        min_psi: Minimum PSI for background
        max_psi: Maximum PSI for background
        
    Returns:
        Tuple of (num_up, num_down, num_bg)
    """
    up: Dict[str, list] = {}
    down: Dict[str, list] = {}
    bg: Dict[str, list] = {}
    
    logging.debug("Making input files from rMATS output")
    
    with open(rmats_path, 'r') as r_file:
        # Skip header
        next(r_file)
        
        for line in r_file:
            result = parse_rMATS_line(
                line, sig_fdr, sig_delta_psi, bg_fdr, min_psi, max_psi, coord_indices
            )
            if result:
                category, key, coord_string = result
                if category == 'up':
                    up[key] = [1, coord_string]
                elif category == 'down':
                    down[key] = [1, coord_string]
                elif category == 'bg':
                    bg[key] = [1, coord_string]
    
    logging.debug(
        "Done populating initial dictionaries with possible duplicates"
    )
    logging.debug(
        "Number of up, down, and background exons are: %d, %d, %d",
        len(up), len(down), len(bg)
    )
    
    logging.debug("Removing exons included in more than one dictionaries..")
    up, down, bg = deduplicate_events(up, down, bg)
    
    nu, nd, nb = write_coordinate_files(up, down, bg, exon_path, header)
    
    logging.debug("Number of upregulated (high in sample_1) exons: %d", nu)
    logging.debug("Number of downregulated (high in sample_2) exons: %d", nd)
    logging.debug("Number of background exons: %d", nb)
    
    return nu, nd, nb


def convert_miso_to_rmats(
    miso_path: str,
    output_path: Path,
    bin_path: Path,
    event_type: str
) -> str:
    """
    Convert MISO format to rMATS format using Perl converter.
    
    Args:
        miso_path: Path to MISO input file
        output_path: Directory for converted output
        bin_path: Directory containing conversion scripts
        event_type: Event type (se, a3ss, a5ss, ri, mxe)
        
    Returns:
        Path to converted rMATS file
        
    Raises:
        RuntimeError: If conversion fails
    """
    rmats_output = output_path / 'converted.rMATS.txt'
    
    event_key = event_type.lower()
    converters = {
        "se": "miso2rMATS.SE.pl",
        "a3ss": "miso2rMATS.A3SS.pl",
        "a5ss": "miso2rMATS.A5SS.pl",
        "ri": "miso2rMATS.RI.pl",
        "mxe": "miso2rMATS.MXE.pl",
    }
    if event_key not in converters:
        raise RuntimeError(f"Unsupported MISO conversion event type: {event_type}")

    converter = bin_path / converters[event_key]
    
    conv_cmd = f'perl "{converter}" 1 10 "{miso_path}" "{rmats_output}"'
    status, output = run_perl_command(conv_cmd)
    
    logging.debug("Conversion from miso to rMATS is done. Status: %s", status)
    
    if status != 0:
        logging.debug("error in converting miso file")
        logging.debug("error detail: %s", output)
        raise RuntimeError(f"MISO conversion failed with status {status}")
    
    logging.debug(output)
    return str(rmats_output)
