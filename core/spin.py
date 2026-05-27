# core/spin.py
from dataclasses import dataclass, field
from typing import List, Optional, Dict
from .carbon_environment import CarbonEnvironment, CarbonClustering
import numpy as np

@dataclass
class CarbonGroup:
    def __init__(self, carbon_shift: float, proton_shifts: List[float], peaks: List[Dict] = None):
        self.carbon_shift = carbon_shift
        self.proton_shifts = proton_shifts
        self.peaks = peaks or []
        
    @property
    def average_proton_shift(self) -> Optional[float]:
        """Average proton shift for this carbon group."""
        if not self.proton_shifts:
            return None
        return np.mean(self.proton_shifts)

    @property
    def carbon_type(self) -> str:
        """
        Classifica carbono baseado em multiplicidade observada, NÃO em shift.
        Shift é usado apenas como fallback.
        """
        from .carbon_environment import CarbonEnvironment
        
        n_protons = len(self.proton_shifts)
        c_shift = self.carbon_shift
        h_avg = self.average_proton_shift
        
        # REGRA ESPECIAL: 2 prótons no mesmo carbono = methylene (a menos que seja oxigenado)
        if n_protons == 2:
            # Exceção: se for oxigenado (SER/THR), mantém como methine
            if 55 <= c_shift <= 75 and h_avg and h_avg > 3.5:
                return 'methine'
            return 'methylene'
        
        # Para outros casos, usa a classificação normal
        return CarbonEnvironment.classify_carbon(c_shift, h_avg)
    
    @property
    def proton_count(self) -> int:
        """Keep for compatibility, but carbon_type is now primary."""
        return len(self.proton_shifts)
   

@dataclass
class Spin:
    """Representa um spin com seus shifts de próton e carbono"""
    proton_shift: float
    carbon_shift: Optional[float] = None
    carbon_type: Optional[str] = None  # methyl, methylene, methine, aromatic
    multiplicity: int = 1
    source_peaks: List = field(default_factory=list)
    
    def __post_init__(self):
        if self.source_peaks is None:
            self.source_peaks = []
    
    @property
    def is_aromatic(self) -> bool:
        """Verifica se é aromático baseado no carbono"""
        return self.carbon_shift is not None and self.carbon_shift > 100
    
    @property
    def fingerprint(self) -> Dict[str, int]:
        """Retorna fingerprint estrutural baseado no carbono"""
        return {
            'methyl': 1 if self.carbon_type == 'methyl' else 0,
            'methylene': 1 if self.carbon_type == 'methylene' else 0,
            'methine': 1 if self.carbon_type == 'methine' else 0,
            'aromatic': 1 if self.is_aromatic else 0
        }

class SpinSystem:
    """Sistema de spin com agrupamento por carbono"""
    
    def __init__(self, hn_shift=None):
        self.hn_shift = hn_shift
        self.spins = []  # Lista de Spin objects (mantido para compatibilidade)
        self.carbon_groups = []  # Lista de CarbonGroup objects
        self.fingerprint = {'methyl': 0, 'methylene': 0, 'methine': 0, 'aromatic': 0}

    def group_by_carbon(self, tolerance=0.2):
        """
        Agrupa spins por carbono usando tolerância, não rounding.
        tolerance: ppm de tolerância para considerar mesmo carbono
        """
        if not hasattr(self, 'spins') or not self.spins:
            return
        
        # Coleta todos os spins com carbono
        carbon_spins = [(s.proton_shift, s.carbon_shift, s) 
                        for s in self.spins if s.carbon_shift is not None]
        
        if not carbon_spins:
            return
        
        # Ordena por carbon shift
        carbon_spins.sort(key=lambda x: x[1])
        
        # Agrupa com tolerância
        groups = []
        current_group = [carbon_spins[0]]
        
        for ps, cs, spin in carbon_spins[1:]:
            if abs(cs - current_group[-1][1]) <= tolerance:
                current_group.append((ps, cs, spin))
            else:
                groups.append(current_group)
                current_group = [(ps, cs, spin)]
        
        if current_group:
            groups.append(current_group)
        
        # Cria CarbonGroups
        self.carbon_groups = []
        for group in groups:
            if group:
                carbon_shift = np.mean([g[1] for g in group])
                proton_shifts = [g[0] for g in group]
                spins = [g[2] for g in group]
                
                carbon_group = CarbonGroup(
                    carbon_shift=carbon_shift,
                    proton_shifts=proton_shifts,
                    peaks=spins
                )
                self.carbon_groups.append(carbon_group)
        
        # Atualiza fingerprint
        self._update_fingerprint_from_groups()

    @property
    def n_unique_carbons(self) -> int:
        """Número de carbonos únicos no sistema"""
        return len(self.carbon_groups)

    def get_carbon_group_by_shift(self, shift, tolerance=0.5):
        """Encontra grupo de carbono por shift"""
        for group in self.carbon_groups:
            if abs(group.carbon_shift - shift) <= tolerance:
                return group
        return None

    def get_structural_fingerprint(self) -> Dict[str, int]:
        """Retorna fingerprint estrutural baseado em carbon groups"""
        return self.fingerprint.copy()
    
    def _update_fingerprint_from_groups(self):
        """Atualiza fingerprint baseado em carbon groups"""
        self.fingerprint = {'methyl': 0, 'methylene': 0, 'methine': 0, 'aromatic': 0}
        
        for group in self.carbon_groups:
            carbon_type = group.carbon_type
            if carbon_type in self.fingerprint:
                self.fingerprint[carbon_type] += 1

    def validate_fingerprint(self) -> bool:
        """
        Valida se o fingerprint é quimicamente plausível.
        Retorna False se detectar inconsistências.
        """
        total_carbons = sum(self.fingerprint.values())
        
        # Caso extremo: muitos methyls
        if self.fingerprint.get('methyl', 0) > 3:
            print(f"  Warning: {self.fingerprint.get('methyl', 0)} methyls é improvável")
            return False
        
        # Caso extremo: muitos methines sem methylene
        if (self.fingerprint.get('methine', 0) > 3 and 
            self.fingerprint.get('methylene', 0) == 0):
            print(f"  Warning: {self.fingerprint.get('methine', 0)} methines sem methylene")
            return False
        
        # GLY não pode ter methyl
        if (self.fingerprint.get('methyl', 0) > 0 and 
            self.fingerprint.get('methylene', 0) == 1 and
            self.fingerprint.get('methine', 0) == 0):
            print(f"  Warning: GLY-like fingerprint com methyl")
            return False
        
        return True

    def debug_carbon_groups(self):
        """Debug: mostra detalhes dos carbon groups"""
        #print(f"\n  Carbon Groups for HN={self.hn_shift}:")
        #for i, group in enumerate(self.carbon_groups):
        #    print(f"    Group {i}: C={group.carbon_shift:.2f}, "
        #        f"H={[f'{h:.3f}' for h in group.proton_shifts]}, "
        #        f"type={group.carbon_type}, "
        #        f"n_protons={group.proton_count}")
        #print(f"    Fingerprint: {self.fingerprint}")

@dataclass
class SpinSystemStructure:
    """Estrutura completa de um sistema de spin com anotação HSQC"""
    hn_shift: Optional[float] = None
    spins: List[Spin] = field(default_factory=list)
    carbon_groups: List['CarbonGroup'] = field(default_factory=list)
    ca_group: Optional['CarbonGroup'] = None

    
    def __post_init__(self):
        """Agrupa por carbono após inicialização"""
        if self.spins:
            self.group_by_carbon()

    @property
    def has_oxygenated_carbon(self) -> bool:
        """Verifica se há carbonos oxigenados no sistema"""
        for group in self.carbon_groups:
            if group.carbon_type == 'methine' and 60 <= group.carbon_shift <= 75:
                return True
            if group.carbon_type == 'methylene' and group.average_proton_shift and group.average_proton_shift > 3.5:
                return True
        return False
    
    @property
    def fingerprint(self) -> Dict[str, int]:
        """Fingerprint baseado em carbon groups"""
        result = {'methyl': 0, 'methylene': 0, 'methine': 0, 'aromatic': 0}
        for group in self.carbon_groups:
            carbon_type = group.carbon_type
            if carbon_type in result:
                result[carbon_type] += 1
        return result    
    
    def group_by_carbon(self, tolerance=0.3):
        """Agrupa spins por carbono usando tolerância"""
        if not self.spins:
            return
        
        from .carbon_environment import CarbonClustering
        
        # Coleta todos os pares H-C
        hc_pairs = []
        for spin in self.spins:
            if spin.carbon_shift is not None:
                hc_pairs.append({
                    'proton_shift': spin.proton_shift,
                    'carbon_shift': spin.carbon_shift,
                    'spin': spin
                })
        
        if not hc_pairs:
            return
        
        # DEBUG: Mostrar antes do clustering
        #print(f"    HC pairs before clustering:")
        #for p in hc_pairs:
        #    print(f"      H={p['proton_shift']:.3f}, C={p['carbon_shift']:.1f}")
        
        # CLUSTERING
        clusters = CarbonClustering.cluster_carbons(hc_pairs)
        
        # DEBUG: Mostrar clusters
        #print(f"    Clusters formed:")
        #for i, cluster in enumerate(clusters):
        #    shifts = [p['carbon_shift'] for p in cluster]
        #    print(f"      Cluster {i}: C={shifts}")
        
        self.carbon_groups = []
        for cluster in clusters:
            if cluster:
                carbon_shift = np.mean([p['carbon_shift'] for p in cluster])
                proton_shifts = [p['proton_shift'] for p in cluster]
                spins = [p['spin'] for p in cluster]
                
                print(f"      Creating CarbonGroup: C={carbon_shift:.1f}, H={proton_shifts}")
                
                group = CarbonGroup(
                    carbon_shift=carbon_shift,
                    proton_shifts=proton_shifts,
                    peaks=spins
                )
                self.carbon_groups.append(group)
    
    def _update_fingerprint(self):
        """Atualiza o fingerprint baseado nos carbon_groups"""
        # Este método pode ser chamado diretamente ou via property
        pass
    
    
    def _find_ca_group(self):
        """Encontra o grupo que provavelmente é o carbono alfa"""
        best_ca = None
        best_score = 0.0
        
        for group in self.carbon_groups:
            score = 0.0
            
            # CA típico: 50-65 ppm
            if 50 <= group.carbon_shift <= 65:
                score += 0.4
            
            # Próton HA típico: 3.5-5.0
            if 3.5 <= group.average_proton_shift <= 5.0:
                score += 0.4
            
            # Deve ser CH (methine) para CA
            if group.carbon_type == 'methine':
                score += 0.2
            
            if score > best_score:
                best_score = score
                best_ca = group
        
        self.ca_group = best_ca
    
    @property
    def fingerprint(self) -> Dict[str, int]:
        """Fingerprint baseado em carbon groups (não em spins individuais)"""
        result = {'methyl': 0, 'methylene': 0, 'methine': 0, 'aromatic': 0}
        for group in self.carbon_groups:
            carbon_type = group.carbon_type
            if carbon_type in result:
                result[carbon_type] += 1
        return result
    
    @property
    def n_spins(self) -> int:
        return len(self.spins)
    
    @property
    def n_carbon_annotated(self) -> int:
        """Número de spins com anotação de carbono"""
        return sum(1 for s in self.spins if s.carbon_shift is not None)
    
    @property
    def n_unique_carbons(self) -> int:
        """Número de carbonos únicos no sistema"""
        return len(self.carbon_groups)
    
    @property
    def carbon_annotation_ratio(self) -> float:
        """Proporção de spins com anotação de carbono"""
        if self.n_spins == 0:
            return 0.0
        return self.n_carbon_annotated / self.n_spins
    
    def get_hc_pairs(self) -> List[tuple]:
        """Retorna pares (H, C) para todos os spins com carbono"""
        return [(s.proton_shift, s.carbon_shift) for s in self.spins if s.carbon_shift is not None]
    
    def get_spins_by_type(self, carbon_type: str) -> List[Spin]:
        """Retorna spins de um tipo específico"""
        return [s for s in self.spins if s.carbon_type == carbon_type]
    
    def get_spins_with_carbon(self) -> List[Spin]:
        """Retorna apenas spins que têm anotação de carbono"""
        return [s for s in self.spins if s.carbon_shift is not None]
    
    def get_carbon_group_by_shift(self, shift, tolerance=0.5):
        """Encontra grupo de carbono por shift"""
        for group in self.carbon_groups:
            if abs(group.carbon_shift - shift) <= tolerance:
                return group
        return None
    
