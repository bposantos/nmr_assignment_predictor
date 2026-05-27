# graph/graph_builder.py

import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Set, Optional, Tuple

from core.peak import Peak
from core.spin import Spin, SpinSystemStructure
from core.spin_system import SpinSystem
from core.node_helpers import is_hn_node, is_spin_node, extract_hn_shift, extract_spin_shift, extract_hn_from_spin_node


@dataclass
class SpinCluster:
    """Representa um cluster LOCAL de picos para um mesmo HN"""
    shift: float  # shift médio
    hn_shift: float  # HN associado
    members: List[float] = field(default_factory=list)
    peaks: List[Peak] = field(default_factory=list)
    
    @property
    def n(self) -> int:
        return len(self.members)
    
    @property
    def variance(self) -> float:
        if self.n < 2:
            return 0.0
        mean = self.shift
        return sum((m - mean) ** 2 for m in self.members) / self.n


class GraphBuilder:

    def __init__(self, proton_tolerance=0.03, intensity_threshold=500):
        self.proton_tolerance = proton_tolerance
        self.intensity_threshold = intensity_threshold
        
        # Regiões espectrais
        self.HN_REGION = (7.0, 9.5)
        self.ALIPHATIC_REGION = (0.5, 5.5)
        
        # HSQC index for fast lookup
        self.hsqc_index = []  # Lista de dicts com H, C, peak

    def find_existing_cluster(self, shift, clusters, tolerance=0.04):
        """Encontra um cluster existente dentro da tolerância"""
        for c in clusters:
            if abs(c - shift) < tolerance:
                return c
        return None
    
    def build_noesy_graph(self, noesy_peaks):
        G = nx.Graph()
        for peak in noesy_peaks:
            p1 = round(float(peak.w1), 3)
            p2 = round(float(peak.w2), 3)
            G.add_node(p1)
            G.add_node(p2)
            G.add_edge(p1, p2, peak=peak, intensity=peak.intensity, spectrum="NOESY")
        return G

    def build_hsqc_graph(self, hsqc_peaks):
        """Build HSQC graph and populate index for fast lookup"""
        G = nx.Graph()
        self.hsqc_index = []  # Reset index
        
        for peak in hsqc_peaks:
            proton = round(float(peak.w1), 3)
            carbon = round(float(peak.w2), 2)
            
            proton_float = float(proton)
            if self.HN_REGION[0] <= proton_float <= self.HN_REGION[1]:
                proton_node = f"HN:{proton}"
            else:
                proton_node = f"H:{proton}"
            
            carbon_node = f"C:{carbon}"
            G.add_node(proton_node)
            G.add_node(carbon_node)
            G.add_edge(proton_node, carbon_node, peak=peak, intensity=peak.intensity, spectrum="HSQC")
            
            # Add to index for fast lookup
            self.hsqc_index.append({
                "H": proton_float,
                "C": carbon,
                "peak": peak,
                "intensity": peak.intensity
            })
        
        # Sort by intensity for better match selection
        self.hsqc_index.sort(key=lambda x: x["intensity"], reverse=True)
        
        return G

    def build_tocsy_graph(self, tocsy_peaks):
        """
        Constrói grafo TOCSY com clustering LOCAL (apenas HNs).
        Spins alifáticos são criados como nós locais associados a cada HN.
        """
        G = nx.Graph()
        
        # Filtra ruído por intensidade
        filtered_peaks = [p for p in tocsy_peaks if p.intensity > self.intensity_threshold]
        
        # Clusteriza APENAS HNs (global)
        hn_clusters = []
        hn_shift_map = {}
        
        for peak in filtered_peaks:
            w1 = float(peak.w1)
            w2 = float(peak.w2)
            
            if self.HN_REGION[0] <= w1 <= self.HN_REGION[1]:
                existing = self.find_existing_cluster(w1, hn_clusters, self.proton_tolerance)
                if existing is None:
                    hn_clusters.append(w1)
                    hn_shift_map[w1] = w1
                else:
                    hn_shift_map[w1] = existing
            
            if self.HN_REGION[0] <= w2 <= self.HN_REGION[1]:
                existing = self.find_existing_cluster(w2, hn_clusters, self.proton_tolerance)
                if existing is None:
                    hn_clusters.append(w2)
                    hn_shift_map[w2] = w2
                else:
                    hn_shift_map[w2] = existing
        
        print(f"  HN clusters: {len(hn_clusters)}")
        
        # Agrupa picos por HN clusterizado
        hn_groups = defaultdict(list)
        
        for peak in filtered_peaks:
            w1 = float(peak.w1)
            w2 = float(peak.w2)
            
            if self.HN_REGION[0] <= w1 <= self.HN_REGION[1]:
                hn_shift = hn_shift_map.get(w1, w1)
                other_shift = w2
                hn_groups[hn_shift].append((other_shift, peak))
            
            elif self.HN_REGION[0] <= w2 <= self.HN_REGION[1]:
                hn_shift = hn_shift_map.get(w2, w2)
                other_shift = w1
                hn_groups[hn_shift].append((other_shift, peak))
        
        # Constrói o grafo com nós LOCAIS para spins alifáticos
        for hn_shift, correlations in hn_groups.items():
            hn_float = float(hn_shift)
            if not (self.HN_REGION[0] <= hn_float <= self.HN_REGION[1]):
                continue
            
            hn_node = f"HN:{hn_shift:.3f}"
            G.add_node(hn_node, atom_type='HN', shift=hn_float)
            
            # Clusteriza spins alifáticos LOCALMENTE (apenas para este HN)
            local_spin_clusters = []
            local_spin_map = {}
            
            # Primeiro, clusteriza os shifts deste HN específico
            for other_shift, peak in correlations:
                existing = self.find_existing_cluster(
                    other_shift, local_spin_clusters, self.proton_tolerance
                )
                if existing is None:
                    local_spin_clusters.append(other_shift)
                    local_spin_map[other_shift] = other_shift
                else:
                    local_spin_map[other_shift] = existing
            
            # Cria nós para os clusters locais
            for cluster_shift in local_spin_clusters:
                # Nó LOCAL: contém o HN shift para evitar colisões
                spin_node = f"{hn_shift:.3f}:spin:{cluster_shift:.3f}"
                G.add_node(
                    spin_node, 
                    atom_type='aliphatic', 
                    shift=cluster_shift,
                    hn_anchor=hn_float
                )
                
                # Adiciona edge HN-spin
                G.add_edge(hn_node, spin_node)
            
            # Agora adiciona edges entre spins alifáticos (HA-HB, etc.)
            # Isso é feito baseado em picos que conectam dois spins alifáticos
            # Para cada par de spins no mesmo HN que têm correlação direta
            spin_nodes = [f"{hn_shift:.3f}:spin:{cs:.3f}" for cs in local_spin_clusters]
            
            # Conecta spins que aparecem juntos em picos
            for i, (other_shift1, peak1) in enumerate(correlations):
                for other_shift2, peak2 in correlations[i+1:]:
                    # Verifica se há correlação direta entre estes dois spins
                    # (simplificado: se aparecem no mesmo pico ou em picos relacionados)
                    # Para TOCSY, spins alifáticos correlacionados aparecem em picos cruzados
                    pass  # Implementação mais simples: deixar apenas HN-spin por enquanto
        
        # Remove nós isolados
        isolated = [node for node in G.nodes if G.degree[node] == 0]
        G.remove_nodes_from(isolated)
        
        hn_nodes = [n for n in G.nodes if n.startswith('HN:')]
        spin_nodes = [n for n in G.nodes if ':spin:' in n]
        
        print(f"TOCSY graph: {len(G.nodes)} nodes, {len(G.edges)} edges")
        print(f"  HN nodes: {len(hn_nodes)}")
        print(f"  Spin nodes: {len(spin_nodes)}")

        #print(f"\n  DEBUG: Todos os HN nodes no grafo:")
        hn_nodes = [n for n in G.nodes if n.startswith('HN:')]
        for node in sorted(hn_nodes, key=lambda x: float(x.split(':')[1])):
            hn_shift = float(node.split(':')[1])
            neighbors = list(G.neighbors(node))
            n_spins = len([n for n in neighbors if ':spin:' in n])
            #print(f"    HN={hn_shift:.3f}: {n_spins} spins alifáticos, {len(neighbors)} total vizinhos")
        
        return G  # Retorna None para cluster_objects (não usado mais)

    def build_spin_system_structure(self, spin_system):
        """
        Build a complete spin system structure from the spin system.
        """
        from core.spin import SpinSystemStructure, Spin
        
        # Cria nova estrutura
        structure = SpinSystemStructure(hn_shift=spin_system.hn_shift)

        #print(f"\n  DEBUG: Building structure for HN={spin_system.hn_shift:.3f}")
        print(f"    Spin nodes: {len(spin_system.spin_nodes)}")        
        
        # Converte spin_nodes para objetos Spin
        for spin_node in spin_system.spin_nodes:
            proton_shift = getattr(spin_node, 'proton_shift', None)
            carbon_shift = getattr(spin_node, 'carbon_shift', None)

            # DEBUG: Mostrar cada par HC
            #if carbon_shift:
                #print(f"      HC pair: H={proton_shift:.3f}, C={carbon_shift:.1f}")

            carbon_type = getattr(spin_node, 'carbon_type', None)
            
            # Se não tem carbon_type mas tem carbon_shift, infere
            # CORREÇÃO: Usar ranges mais precisos
            if carbon_shift and not carbon_type:
                if carbon_shift < 22:  # Methyl reduzido
                    carbon_type = 'methyl'
                elif 22 <= carbon_shift < 50:  # Methylene
                    carbon_type = 'methylene'
                elif 50 <= carbon_shift < 100:  # Methine
                    carbon_type = 'methine'
                elif carbon_shift > 100:  # Aromatic
                    carbon_type = 'aromatic'
            
            # Cria objeto Spin
            spin = Spin(
                proton_shift=proton_shift,
                carbon_shift=carbon_shift,
                carbon_type=carbon_type,
                multiplicity=2 if carbon_type == 'methylene' else 1,
                source_peaks=[]
            )
            
            structure.spins.append(spin)
        
        # O __post_init__ já chama group_by_carbon, então não precisamos chamar novamente
        # Mas forçamos a chamada para garantir que carbon_groups está populado
        if len(structure.spins) > 0:
            structure.group_by_carbon()
        
        print(f"    Structure: {structure.n_spins} spins, {structure.n_carbon_annotated} with carbon annotation, "
            f"{len(structure.carbon_groups)} carbon groups, fingerprint={structure.fingerprint}")
        
        # ============================================================
        # CORREÇÃO CRÍTICA: Propagar fingerprint e carbon_groups de volta para o spin_system
        # ============================================================
        spin_system.structural_fingerprint = structure.fingerprint.copy()
        spin_system.carbon_groups = structure.carbon_groups
        spin_system.fingerprint = structure.fingerprint.copy()  # Para compatibilidade
        
        # Debug para verificar propagação
        #print(f"    → Propagated fingerprint to spin_system: {spin_system.fingerprint}")
        
        return structure

    def _find_ca_group(self, structure):
        """Encontra o grupo que provavelmente é o carbono alfa"""
        best_ca = None
        best_score = 0.0
        
        for group in structure.carbon_groups:
            score = 0.0
            
            # CA típico: 50-65 ppm
            if 50 <= group.carbon_shift <= 65:
                score += 0.4
            
            # Próton HA típico: 3.5-5.0
            if 3.5 <= group.average_proton <= 5.0:
                score += 0.4
            
            # Deve ser CH (metionina) para CA
            if group.carbon_type == 'methine':
                score += 0.2
            
            if score > best_score:
                best_score = score
                best_ca = group
        
        return best_ca

    def extract_spin_systems(self, tocsy_graph, hsqc_graph=None, min_spins=2):
        """
        Extrai spin systems de forma SIMPLES: cada HN com seus vizinhos diretos.
        """
        from core.spin_system import SpinSystem as GraphSpinSystem
        
        spin_systems = []
        
        # Encontra todos os nós HN
        hn_nodes = [n for n in tocsy_graph.nodes if n.startswith('HN:')]
        
        print(f"  Building spin systems from {len(hn_nodes)} HN nodes...")
        
        all_hn_shifts = sorted([float(n.split(':')[1]) for n in hn_nodes])
        print(f"  Todos os HNs detectados: {[f'{s:.3f}' for s in all_hn_shifts]}")
        
        system_counter = 0
        
        for hn_node in hn_nodes:
            hn_shift = float(hn_node.split(':')[1])
            
            # Filtra HNs aromáticos (shift < 7.5)
            if hn_shift < 7.5:
                print(f"  Ignorando HN={hn_shift:.3f} (possível sinal aromático)")
                continue
            
            if not (7.0 <= hn_shift <= 9.5):
                continue
            
            system_counter += 1
            
            # Cria subgrafo para este componente
            component = set([hn_node])
            
            # Adiciona vizinhos diretos do HN
            # CORREÇÃO: Verifica se o vizinho contém ':spin:' (não começa com)
            for neighbor in tocsy_graph.neighbors(hn_node):
                if ':spin:' in neighbor:  # <--- CORREÇÃO AQUI
                    component.add(neighbor)
            
            # Cria subgrafo
            subgraph = tocsy_graph.subgraph(component).copy()
            
            # Cria sistema de spin
            ss = GraphSpinSystem(
                hn_shift=hn_shift,
                system_id=system_counter,
                graph=subgraph,
                peaks=[]
            )
            ss.hn_node = hn_node
            ss.spin_nodes = []
            
            # Adiciona spins ao sistema
            for node in subgraph.nodes:
                if node.startswith('HN:'):
                    continue
                
                # Extrai shift do spin (formato: "7.769:spin:4.229")
                try:
                    parts = node.split(':')
                    if len(parts) >= 3 and parts[1] == 'spin':
                        spin_shift = float(parts[2])
                    else:
                        continue
                except:
                    continue
                
                # Cria objeto spin simples
                class SpinNode:
                    pass
                
                spin_node = SpinNode()
                spin_node.node_id = node
                spin_node.proton_shift = spin_shift
                spin_node.carbon_shift = None
                spin_node.carbon_type = None
                
                ss.spin_nodes.append(spin_node)
            
            print(f"  HN={hn_shift:.3f}: {len(ss.spin_nodes)} spins")
            
            if len(ss.spin_nodes) >= min_spins - 1:
                spin_systems.append(ss)
                print(f"  ✓ Sistema {system_counter}: HN={hn_shift:.3f}, {len(ss.spin_nodes)} spins alifáticos")
            else:
                print(f"  Rejeitando HN={hn_shift:.3f}: apenas {len(ss.spin_nodes)} spins (min={min_spins-1})")
        
        print(f"\nSpin systems: {len(spin_systems)} sistemas alifáticos válidos")
        return spin_systems

    def annotate_spin_system_with_hsqc(self, spin_system, tolerance=0.04):
        """
        Annotate each spin in the system with its carbon shift from HSQC.
        """
        if not hasattr(self, 'hsqc_peaks') or not self.hsqc_peaks:
            print("  Warning: No HSQC peaks available")
            return spin_system
        
        annotated_count = 0
        
        for spin_node in spin_system.spin_nodes:
            proton_shift = getattr(spin_node, 'proton_shift', None)
            if not proton_shift:
                continue
            
            best_match = None
            best_diff = float('inf')
            
            for hsqc_peak in self.hsqc_peaks:
                h_shift = float(hsqc_peak.w1) if 'w1' in hsqc_peak.__dict__ else float(hsqc_peak.w2)
                c_shift = float(hsqc_peak.w2) if 'w1' in hsqc_peak.__dict__ else float(hsqc_peak.w1)
                
                diff = abs(proton_shift - h_shift)
                if diff < tolerance and diff < best_diff:
                    best_diff = diff
                    best_match = (c_shift, hsqc_peak)
            
            if best_match:
                c_shift, _ = best_match
                spin_node.carbon_shift = c_shift
                
                # Inferir carbon_type do shift (mais confiável que proton count)
                if c_shift < 25:
                    spin_node.carbon_type = 'methyl'
                elif 25 <= c_shift < 55:
                    spin_node.carbon_type = 'methylene'
                elif 55 <= c_shift < 100:
                    spin_node.carbon_type = 'methine'
                elif c_shift > 100:
                    spin_node.carbon_type = 'aromatic'
                else:
                    spin_node.carbon_type = 'unknown'
                    
                annotated_count += 1
        
        print(f"  Annotated {annotated_count} spins with HSQC data")
        return spin_system

    def visualize_spin_systems(self, spin_systems):
        if not spin_systems:
            print("No spin systems to visualize")
            return
            
        n_systems = len(spin_systems)
        cols = min(4, n_systems)
        rows = (n_systems + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
        if n_systems == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
        
        for i, system in enumerate(spin_systems):
            ax = axes[i]
            G = system.graph
            pos = nx.spring_layout(G, seed=42)
            
            node_colors = []
            node_labels = {}
            
            for node in G.nodes:
                if node.startswith('HN:'):
                    node_colors.append('red')
                    shift = node.split(':')[1]
                    node_labels[node] = f"HN:{shift}"
                elif ':spin:' in node:
                    node_colors.append('lightblue')
                    # Extrai shift do spin
                    parts = node.split(':')
                    if len(parts) >= 3:
                        shift = parts[2]
                        node_labels[node] = f"spin:{shift}"
                    else:
                        node_labels[node] = node
                else:
                    node_colors.append('gray')
                    node_labels[node] = node
            
            nx.draw(G, pos, ax=ax, with_labels=True, labels=node_labels,
                   node_color=node_colors, node_size=500, font_size=7)
            
            hn_nodes = [n for n in G.nodes if n.startswith('HN:')]
            hn_info = f" HN={hn_nodes[0].split(':')[1]}ppm" if hn_nodes else ""
            ax.set_title(f"Spin System {system.system_id}{hn_info}")
        
        for j in range(i + 1, len(axes)):
            axes[j].axis("off")
        
        plt.tight_layout()
        plt.show()

    def _calculate_contextual_fingerprint(self, structure, stage='pre_noesy'):
        """
        Calcula fingerprint considerando apenas átomos esperados no estágio
        """
        from database.statistical_ranges import get_expected_atoms, ATOM_TO_TYPE
        
        if not structure.assignment:
            return structure.fingerprint.copy()
        
        expected_atoms = get_expected_atoms(structure.assignment, stage)
        
        # Mapeia átomos esperados para tipos
        expected_types = {}
        for atom in expected_atoms:
            atom_type = ATOM_TO_TYPE.get(atom, 'unknown')
            expected_types[atom_type] = expected_types.get(atom_type, 0) + 1
        
        # Calcula fingerprint baseado apenas no que era esperado
        contextual_fp = {'methyl': 0, 'methylene': 0, 'methine': 0, 'aromatic': 0}
        
        for spin in structure.spins:
            # Mapeia spin para tipo (simplificado)
            spin_type = self._infer_spin_type(spin)
            if spin_type in expected_types:
                contextual_fp[spin_type] = contextual_fp.get(spin_type, 0) + 1
        
        return contextual_fp

    def _infer_spin_type(self, spin):
        """Infere tipo estrutural do spin baseado no shift de carbono"""
        if spin.c_shift:
            if spin.c_shift < 25:
                return 'methyl'
            elif 25 <= spin.c_shift < 45:
                return 'methylene'
            elif 45 <= spin.c_shift < 65:
                return 'methine'
            elif spin.c_shift > 100:
                return 'aromatic'
        
        # Fallback baseado no nome
        if 'methyl' in spin.name.lower() or 'ch3' in spin.name.lower():
            return 'methyl'
        elif 'ch2' in spin.name.lower():
            return 'methylene'
        elif 'ch' in spin.name.lower():
            return 'methine'
        
        return 'unknown'

    def visualize_graph(self, graph, title="Graph", show_cluster_info=False):
        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(graph, seed=42)
        
        node_colors = []
        node_labels = {}
        
        for node in graph.nodes:
            if node.startswith('HN:'):
                node_colors.append('red')
                node_labels[node] = node
            elif ':spin:' in node:
                node_colors.append('lightblue')
                parts = node.split(':')
                if len(parts) >= 3:
                    shift = parts[2]
                    node_labels[node] = f"spin:{shift}"
                else:
                    node_labels[node] = node
            else:
                node_colors.append('gray')
                node_labels[node] = node
        
        nx.draw(graph, pos, with_labels=True, labels=node_labels,
               node_color=node_colors, node_size=700, edge_color='gray', font_size=8)
        plt.title(title)
        plt.tight_layout()
        plt.show()

    def validate_spin_system(self, structure):
        """
        Valida quimicamente um spin system antes de aceitá-lo.
        Retorna False se houver inconsistências graves.
        """
        if not structure.carbon_groups:
            return True
        
        # THR: se tem methyl (~20) e CA (~58), precisa ter CB oxigenado (60-75)
        if (structure.fingerprint.get('methyl', 0) >= 1 and 
            structure.fingerprint.get('methine', 0) >= 1):
            
            # Procura por carbonos oxigenados
            has_oxygenated = False
            for group in structure.carbon_groups:
                if 60 <= group.carbon_shift <= 75:
                    has_oxygenated = True
                    break
            
            # Se tem methyl + methine mas não tem oxigenado, pode ser ALA ou outro
            # Apenas avisa, não rejeita
            if not has_oxygenated:
                print(f"    Warning: THR-like pattern but no oxygenated carbon found")
        
        # GLY: não deve ter carbonos > 55 ppm
        if structure.fingerprint.get('methylene', 0) == 1:
            for group in structure.carbon_groups:
                if group.carbon_shift > 55:
                    print(f"    Warning: GLY-like but has carbon at {group.carbon_shift:.1f} ppm")
        
        return True

    def validate_connectivity(self, structure):
        """
        Verifica se os carbon groups fazem sentido quimicamente.
        Ex: THR não deve ter methylene alifático no lugar do CB.
        """
        if not structure.carbon_groups:
            return True
        
        # Ordena carbon groups por shift
        sorted_groups = sorted(structure.carbon_groups, key=lambda g: g.carbon_shift)
        
        # Verifica gaps grandes (possível fragmentação)
        for i in range(len(sorted_groups) - 1):
            diff = sorted_groups[i+1].carbon_shift - sorted_groups[i].carbon_shift
            if diff > 30:  # Gap grande
                print(f"    Warning: Large gap ({diff:.1f} ppm) between carbons")
                print(f"      {sorted_groups[i].carbon_shift:.1f} -> {sorted_groups[i+1].carbon_shift:.1f}")
        
        return True
