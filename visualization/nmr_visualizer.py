# visualization/nmr_visualizer.py

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from collections import defaultdict


class NMRVisualizer:
    
    def __init__(self):
        # Cores para cada sistema de spin (máximo 20)
        self.system_colors = [
            '#FF0000', '#00FF00', '#0000FF', '#FFA500', '#800080',  # red, green, blue, orange, purple
            '#FFC0CB', '#008080', '#FF00FF', '#4B0082', '#00FFFF',  # pink, teal, magenta, indigo, cyan
            '#FFD700', '#A0522D', '#7FFF00', '#DC143C', '#00FA9A',  # gold, brown, chartreuse, crimson, medium spring green
            '#8B008B', '#FF8C00', '#9932CC', '#FF1493', '#1E90FF'   # dark magenta, dark orange, dark orchid, deep pink, dodger blue
        ]

    def plot_spin_systems_on_tocsy(self, spin_systems, tocsy_peaks, title="TOCSY Spin Systems"):
        """Plot TOCSY spectrum with each spin system in a different color"""
        
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Build peak index for ALL systems
        system_peaks = defaultdict(list)
        assigned_peak_ids = set()
        
        for idx, ss in enumerate(spin_systems):
            proton_shifts = set()
            for node in ss.graph.nodes:
                if node.startswith('HN:'):
                    proton_shifts.add(round(float(node.split(':')[1]), 3))
                elif ':spin:' in node:
                    proton_shifts.add(round(float(node.split(':')[-1]), 3))
            
            for peak in tocsy_peaks:
                w1 = round(float(peak.w1), 3)
                w2 = round(float(peak.w2), 3)
                if w1 in proton_shifts and w2 in proton_shifts:
                    system_peaks[idx].append(peak)
                    assigned_peak_ids.add(id(peak))
        
        # Background unassigned peaks
        for peak in tocsy_peaks:
            if id(peak) not in assigned_peak_ids:
                w1, w2 = float(peak.w1), float(peak.w2)
                ax.plot(w1, w2, 'o', color='lightgray', markersize=3, alpha=0.3)
                ax.plot(w2, w1, 'o', color='lightgray', markersize=3, alpha=0.3)
        
        # Plot ALL systems
        for idx, ss in enumerate(spin_systems):
            color = self.system_colors[idx % len(self.system_colors)]
            
            for peak in system_peaks[idx]:
                w1, w2 = float(peak.w1), float(peak.w2)
                ax.plot(w1, w2, 'o', color=color, markersize=5, alpha=0.7)
                ax.plot(w2, w1, 'o', color=color, markersize=5, alpha=0.7)
            
            # Label
            hn_shift = None
            ha_shift = None
            for node in ss.graph.nodes:
                if node.startswith('HN:'):
                    hn_shift = float(node.split(':')[1])
                elif ':spin:' in node:
                    shift = float(node.split(':')[-1])
                    if 3.5 <= shift <= 5.5:
                        ha_shift = shift
            
            if hn_shift and ha_shift:
                residue = getattr(ss, 'assigned_residue', '?')
                if residue and residue != "UNKNOWN":
                    label = f"S{idx+1}:{residue}"
                else:
                    label = f"S{idx+1}"
                
                ax.annotate(label, xy=(hn_shift, ha_shift), xytext=(5, 5),
                        textcoords='offset points', fontsize=8,
                        color=color, alpha=0.8, fontweight='bold')
        
        ax.set_xlim(10, 0)
        ax.set_ylim(10, 0)
        ax.set_xlabel('1H (ppm)')
        ax.set_ylabel('1H (ppm)')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
    
    def plot_spin_systems_on_hsqc(self, spin_systems, hsqc_peaks, title="Assigned HSQC Spin Systems"):
        """Plot HSQC spectrum with each spin system in a different color"""
        
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Plot all HSQC peaks as gray background
        for peak in hsqc_peaks:
            h = float(peak.w1)
            c = float(peak.w2)
            ax.plot(h, c, 'o', color='lightgray', markersize=4, alpha=0.5)
        
        # Plot each spin system in its color
        for idx, ss in enumerate(spin_systems):
            if not hasattr(ss, 'assigned_residue'):
                continue
            
            color = self.system_colors[idx % len(self.system_colors)]
            
            # Find spins with carbon annotation
            for node in ss.graph.nodes:
                if ':spin:' in node and 'carbon_shift' in ss.graph.nodes[node]:
                    h_shift = float(node.split(':')[-1])
                    c_shift = ss.graph.nodes[node]['carbon_shift']
                    
                    # Check if it's an aliphatic spin (not HN)
                    if h_shift < 6.0:  # Aliphatic region
                        ax.plot(h_shift, c_shift, 'o', color=color, markersize=8, alpha=0.8)
                        
                        # Add label for the atom type if available
                        if 'carbon_type' in ss.graph.nodes[node]:
                            c_type = ss.graph.nodes[node]['carbon_type']
                            label = f"{ss.assigned_residue}-{c_type}"
                        else:
                            label = ss.assigned_residue
                        
                        ax.annotate(label, xy=(h_shift, c_shift), xytext=(3, 3),
                                   textcoords='offset points', fontsize=7,
                                   color=color, alpha=0.7)
        
        ax.set_xlim(10, 0)  # Reverse axis for NMR
        ax.set_ylim(180, 10)  # Carbon range
        ax.set_xlabel('1H (ppm)')
        ax.set_ylabel('13C (ppm)')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()

    def plot_spin_systems_summary(self, spin_systems, tocsy_peaks, hsqc_peaks, save_path=None):
        """Create a summary figure with TOCSY and HSQC side by side"""
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # === Build TOCSY peak index for ALL spin systems ===
        system_peaks = defaultdict(list)
        assigned_peak_ids = set()
        
        for idx, ss in enumerate(spin_systems):
            # Coleta todos os shifts de próton deste sistema (mesmo não assinalado)
            proton_shifts = set()
            for node in ss.graph.nodes:
                if node.startswith('HN:'):
                    proton_shifts.add(round(float(node.split(':')[1]), 3))
                elif ':spin:' in node:
                    proton_shifts.add(round(float(node.split(':')[-1]), 3))
            
            # Encontra picos TOCSY
            for peak in tocsy_peaks:
                w1 = round(float(peak.w1), 3)
                w2 = round(float(peak.w2), 3)
                
                if w1 in proton_shifts and w2 in proton_shifts:
                    system_peaks[idx].append(peak)
                    assigned_peak_ids.add(id(peak))
            
            # Se não encontrou picos baseado em shifts, tenta usar ss.peaks (fallback)
            if not system_peaks[idx] and hasattr(ss, 'peaks'):
                for peak in ss.peaks:
                    system_peaks[idx].append(peak)
                    assigned_peak_ids.add(id(peak))
        
        # === TOCSY panel ===
        # Background unassigned peaks (picos que não pertencem a nenhum sistema)
        for peak in tocsy_peaks:
            if id(peak) not in assigned_peak_ids:
                w1, w2 = float(peak.w1), float(peak.w2)
                ax1.plot(w1, w2, 'o', color='lightgray', markersize=3, alpha=0.3)
                ax1.plot(w2, w1, 'o', color='lightgray', markersize=3, alpha=0.3)
        
        # Plot ALL spin systems in their colors (even unassigned)
        for idx, ss in enumerate(spin_systems):
            if idx >= len(self.system_colors):
                color = self.system_colors[idx % len(self.system_colors)]
            else:
                color = self.system_colors[idx]
            
            # Plot peaks
            for peak in system_peaks[idx]:
                w1, w2 = float(peak.w1), float(peak.w2)
                ax1.plot(w1, w2, 'o', color=color, markersize=5, alpha=0.7)
                ax1.plot(w2, w1, 'o', color=color, markersize=5, alpha=0.7)
            
            # Label HN-HA with system number and residue (if assigned)
            hn_shift = None
            ha_shift = None
            for node in ss.graph.nodes:
                if node.startswith('HN:'):
                    hn_shift = float(node.split(':')[1])
                elif ':spin:' in node:
                    shift = float(node.split(':')[-1])
                    if 3.5 <= shift <= 5.5:  # HA region
                        ha_shift = shift
            
            if hn_shift and ha_shift:
                residue = getattr(ss, 'assigned_residue', '?')
                if residue and residue != "UNKNOWN":
                    label = f"S{idx+1}:{residue}"
                else:
                    label = f"S{idx+1}"
                
                ax1.annotate(label, xy=(hn_shift, ha_shift), xytext=(3, 3),
                            textcoords='offset points', fontsize=7,
                            color=color, alpha=0.8, fontweight='bold')
        
        ax1.set_xlim(10, 0)
        ax1.set_ylim(10, 0)
        ax1.set_xlabel('1H (ppm)')
        ax1.set_ylabel('1H (ppm)')
        ax1.set_title('TOCSY - All Spin Systems (Colored by System)')
        ax1.grid(True, alpha=0.3)
        
        # === HSQC panel ===
        # Background all peaks
        for peak in hsqc_peaks:
            h = float(peak.w1)
            c = float(peak.w2)
            ax2.plot(h, c, 'o', color='lightgray', markersize=3, alpha=0.3)
        
        # Plot HSQC annotations for ALL spin systems
        for idx, ss in enumerate(spin_systems):
            if idx >= len(self.system_colors):
                color = self.system_colors[idx % len(self.system_colors)]
            else:
                color = self.system_colors[idx]
            
            for node in ss.graph.nodes:
                if ':spin:' in node and 'carbon_shift' in ss.graph.nodes[node]:
                    h_shift = float(node.split(':')[-1])
                    c_shift = ss.graph.nodes[node]['carbon_shift']
                    
                    if h_shift < 6.0:  # Aliphatic only
                        ax2.plot(h_shift, c_shift, 'o', color=color, markersize=6, alpha=0.8)
                        
                        residue = getattr(ss, 'assigned_residue', '?')
                        if residue and residue != "UNKNOWN":
                            label = f"S{idx+1}"
                        else:
                            label = f"S{idx+1}"
                        
                        ax2.annotate(label, xy=(h_shift, c_shift), xytext=(3, 3),
                                    textcoords='offset points', fontsize=6,
                                    color=color, alpha=0.7)
        
        ax2.set_xlim(10, 0)
        ax2.set_ylim(180, 10)
        ax2.set_xlabel('1H (ppm)')
        ax2.set_ylabel('13C (ppm)')
        ax2.set_title('HSQC - Carbon Assignments (Colored by System)')
        ax2.grid(True, alpha=0.3)
        
        plt.suptitle('NMR Assignment Summary', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"✅ Figure saved to: {save_path}")
        
        plt.show()
    
    def plot_single_spin_system(self, spin_system, tocsy_peaks, hsqc_peaks, system_id):
        """Plot a single spin system with its assignments"""
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        color = self.system_colors[system_id % len(self.system_colors)]
        
        # === TOCSY ===
        for peak in tocsy_peaks:
            w1, w2 = float(peak.w1), float(peak.w2)
            ax1.plot(w1, w2, 'o', color='lightgray', markersize=3, alpha=0.3)
            ax1.plot(w2, w1, 'o', color='lightgray', markersize=3, alpha=0.3)
        
        for peak in spin_system.peaks:
            w1, w2 = float(peak.w1), float(peak.w2)
            ax1.plot(w1, w2, 'o', color=color, markersize=6, alpha=0.8)
            ax1.plot(w2, w1, 'o', color=color, markersize=6, alpha=0.8)
        
        # Label HN-HA
        hn_shift = None
        ha_shift = None
        for node in spin_system.graph.nodes:
            if node.startswith('HN:'):
                hn_shift = float(node.split(':')[1])
            elif ':spin:' in node:
                shift = float(node.split(':')[-1])
                if 3.5 <= shift <= 5.5:
                    ha_shift = shift
        
        if hn_shift and ha_shift:
            residue = getattr(spin_system, 'assigned_residue', 'Unknown')
            ax1.annotate(f"{residue}", xy=(hn_shift, ha_shift), xytext=(5, 5),
                        textcoords='offset points', fontsize=10,
                        color=color, fontweight='bold')
        
        ax1.set_xlim(10, 0)
        ax1.set_ylim(10, 0)
        ax1.set_xlabel('1H (ppm)')
        ax1.set_ylabel('1H (ppm)')
        ax1.set_title(f'Spin System {system_id+1} - TOCSY')
        ax1.grid(True, alpha=0.3)
        
        # === HSQC ===
        for peak in hsqc_peaks:
            h = float(peak.w1)
            c = float(peak.w2)
            ax2.plot(h, c, 'o', color='lightgray', markersize=3, alpha=0.3)
        
        for node in spin_system.graph.nodes:
            if ':spin:' in node and 'carbon_shift' in spin_system.graph.nodes[node]:
                h_shift = float(node.split(':')[-1])
                c_shift = spin_system.graph.nodes[node]['carbon_shift']
                
                if h_shift < 6.0:
                    ax2.plot(h_shift, c_shift, 'o', color=color, markersize=8, alpha=0.8)
                    
                    if 'carbon_type' in spin_system.graph.nodes[node]:
                        c_type = spin_system.graph.nodes[node]['carbon_type']
                        label = c_type[:4]
                        ax2.annotate(label, xy=(h_shift, c_shift), xytext=(3, 3),
                                    textcoords='offset points', fontsize=7,
                                    color=color, alpha=0.7)
        
        ax2.set_xlim(10, 0)
        ax2.set_ylim(180, 10)
        ax2.set_xlabel('1H (ppm)')
        ax2.set_ylabel('13C (ppm)')
        ax2.set_title(f'Spin System {system_id+1} - HSQC')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
