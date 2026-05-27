# core/structural_profile.py

from typing import List, Optional, Dict
from dataclasses import dataclass, field
from enum import Enum


class MethineType(Enum):
    """Tipos de methine para melhor discriminação"""
    ALPHA = "alpha"           # CA: 50-65 ppm, H 3.5-5.0
    OXYGENATED = "oxygenated" # SER/THR: 60-75 ppm, H 3.5-4.5
    ALIPHATIC = "aliphatic"   # cadeia lateral: 30-50 ppm, H 1.5-2.5
    BETA = "beta"             # beta CH: 35-50 ppm


@dataclass
class StructuralProfile:
    """Perfil estrutural detalhado de um spin system"""
    
    # Contagens básicas
    n_carbons: int = 0
    n_methyl: int = 0
    n_methylene: int = 0
    n_methine: int = 0
    n_aromatic: int = 0
    
    # Tipos de methylene (NOVO!)
    n_alpha_methylene: int = 0      # GLY: CA como CH2 (40-50 ppm)
    n_oxygenated_methylene: int = 0 # SER: CB como CH2OH (60-65 ppm)
    n_aliphatic_methylene: int = 0  # LYS/LEU/ILE: cadeia lateral
    
    # Tipos de methine
    n_alpha_methine: int = 0        # CA: 50-65 ppm
    n_beta_methine: int = 0         # beta CH: 35-50 ppm
    n_sidechain_methine: int = 0    # outros: 25-35 ppm
    n_oxygenated_methine: int = 0   # THR: CB como CHOH (60-75 ppm)
    
    # Características especiais
    has_alpha_carbon: bool = False
    has_oxygenated_carbon: bool = False
    has_branching: bool = False
    
    # Shifts específicos
    alpha_shift: Optional[float] = None
    beta_shift: Optional[float] = None
    methyl_shifts: List[float] = field(default_factory=list)
    _carbon_shifts: List[float] = field(default_factory=list, repr=False)
    _methylene_shifts: List[float] = field(default_factory=list, repr=False)
    _methine_shifts: List[float] = field(default_factory=list, repr=False)
    
    @classmethod
    def from_carbon_groups(cls, carbon_groups: List['CarbonGroup']) -> 'StructuralProfile':
        profile = cls()
        
        for group in carbon_groups:
            profile._carbon_shifts.append(group.carbon_shift)
            profile.n_carbons += 1
            c_type = group.carbon_type
            if c_type == 'methylene':
                profile._methylene_shifts.append(group.carbon_shift)
            elif c_type == 'methine':
                profile._methine_shifts.append(group.carbon_shift)            
            c_shift = group.carbon_shift
            h_shift = group.average_proton_shift
            n_protons = group.proton_count
            
            if c_type == 'methyl':
                profile.n_methyl += 1
                profile.methyl_shifts.append(c_shift)
                
            elif c_type == 'methylene':
                profile.n_methylene += 1
                
                # GLY: alpha methylene (40-50 ppm) - NÃO depende de n_protons!
                if 40 <= c_shift <= 50:
                    profile.n_alpha_methylene += 1
                    profile.has_alpha_carbon = True
                    if profile.alpha_shift is None:
                        profile.alpha_shift = c_shift
                
                # SER: oxygenated methylene (55-68 ppm)
                elif 55 <= c_shift <= 70 and h_shift and h_shift > 3.5:
                    profile.n_oxygenated_methylene += 1
                    profile.has_oxygenated_carbon = True
                
                # Aliphatic methylene (20-55 ppm)
                else:
                    profile.n_aliphatic_methylene += 1

            elif c_type == 'methine':
                profile.n_methine += 1
                
                # CORREÇÃO: Usar fronteira mais segura (65 ppm)
                if 65 <= c_shift <= 80 and h_shift and 3.5 <= h_shift <= 4.8:
                    profile.n_oxygenated_methine += 1
                    profile.has_oxygenated_carbon = True
                
                # Alpha methine (55-65 ppm)
                elif 55 <= c_shift <= 65:
                    profile.n_alpha_methine += 1
                    profile.has_alpha_carbon = True
                    if profile.alpha_shift is None:
                        profile.alpha_shift = c_shift
                
                # Beta methine (35-55 ppm)
                elif 35 <= c_shift < 55:
                    profile.n_beta_methine += 1
                    if profile.beta_shift is None:
                        profile.beta_shift = c_shift
                
                # Sidechain methine (20-35 ppm)
                elif 20 <= c_shift < 35:
                    profile.n_sidechain_methine += 1
                    
            elif c_type == 'aromatic':
                profile.n_aromatic += 1
        
        # CORREÇÃO 3: Branching baseado apenas em methyls
        if profile.n_methyl >= 2:
            profile.has_branching = True
        
        return profile
    
    def to_fingerprint(self) -> Dict[str, int]:
        """Converte para fingerprint simples (compatibilidade)"""
        return {
            'methyl': self.n_methyl,
            'methylene': self.n_methylene,
            'methine': self.n_methine,
            'aromatic': self.n_aromatic
        }

    def score_similarity(self, expected_profile: 'StructuralProfile', 
                        residue_name: str = None) -> float:
        """
        Score com penalidades escaladas por completude.
        """
        score = 0.0
        
        # Calcula completude do sistema observado
        total_expected = (expected_profile.n_methine + 
                        expected_profile.n_methylene + 
                        expected_profile.n_methyl)
        completeness = min(1.0, self.n_carbons / max(total_expected, 1))
        
        # ============================================================
        # HARD INCOMPATIBILITIES (mata o candidato)
        # ============================================================
        
        # CORREÇÃO: GLY - regras ABSOLUTAS e claras
        if residue_name == 'GLY':
            # GLY NÃO PODE ter methyl
            if self.n_methyl > 0:
                return 0.0
            
            # GLY NÃO PODE ter methine (qualquer methine)
            if self.n_methine > 0:
                return 0.0
            
            # GLY NÃO PODE ter methylene não-alfa (CH2 em posições erradas)
            # GLY só pode ter alpha methylene (CA em 40-50 ppm)
            if self.n_methylene > 0 and self.n_alpha_methylene == 0:
                return 0.0
            
            # GLY NÃO PODE ter carbonos oxigenados ( > 55 ppm exceto CA)
            for c_shift in self.methyl_shifts:
                if c_shift > 55:
                    return 0.0
        
        # SER não pode ter alpha methylene (é GLY)
        if residue_name == 'SER' and self.n_alpha_methylene > 0:
            return 0.0

        # SER não pode ter methyl
        if residue_name == 'SER' and self.n_methyl > 0:
            return 0.0
        
        # THR não pode ter methylene
        if residue_name == 'THR' and self.n_methylene > 0:
            # Se o methylene for CB de THR (45-55 ppm)
            is_thr_cb = (self.n_methylene == 1 and 
                        self.beta_shift is not None and 
                        45 <= self.beta_shift <= 55)
            if not is_thr_cb:
                return 0.0
        
        # ALA não pode ter beta methine
        if residue_name == 'ALA' and self.n_beta_methine > 0:
            return 0.0

        # ALA não pode ter methylene
        if residue_name == 'ALA' and self.n_methylene > 0:
            return 0.0
        
        # LYS não pode ter methyl
        if residue_name == 'LYS' and self.n_methyl > 0:
            return 0.0
        
        # ============================================================
        # PENALIDADES FORTES (escaladas por completude)
        # ============================================================
        
        # Ausência de methyl quando esperado (escalado)
        if expected_profile.n_methyl > 0 and self.n_methyl == 0:
            score -= 0.30 * completeness
        
        # Ausência de alpha carbon
        if expected_profile.has_alpha_carbon and not self.has_alpha_carbon:
            score -= 0.25 * completeness
        
        # ============================================================
        # BÔNUS FORTES (features raras/discriminativas)
        # ============================================================
        
        # Aromaticidade (MUITO forte)
        if expected_profile.n_aromatic > 0:
            if self.n_aromatic > 0:
                score += 0.40
            else:
                # Penalidade forte para PHE/TYR sem aromatico
                if residue_name in ['PHE', 'TYR', 'TRP']:
                    score -= 0.30
        
        # SER: oxygenated methylene
        if expected_profile.n_oxygenated_methylene > 0:
            if self.n_oxygenated_methylene > 0:
                score += 0.35
            else:
                score -= 0.20
        
        # THR: oxygenated methine
        if expected_profile.n_oxygenated_methine > 0:
            if self.n_oxygenated_methine > 0:
                score += 0.35
            else:
                score -= 0.20
        
        # GLY: alpha methylene
        if expected_profile.n_alpha_methylene > 0:
            if self.n_alpha_methylene > 0:
                score += 0.35
            else:
                score -= 0.20
        
        # LEU: branching + sidechain methine
        if residue_name == 'LEU':
            if self.n_methyl >= 2 and self.n_sidechain_methine >= 1:
                score += 0.25
        
        # ILE: branching + beta methine
        elif residue_name == 'ILE':
            if self.n_methyl >= 2 and self.n_beta_methine >= 1:
                score += 0.25
        
        # ============================================================
        # BÔNUS FRACOS (features comuns, escalados)
        # ============================================================
        
        # Alpha methine (comum, pouco discriminativo)
        if expected_profile.n_alpha_methine > 0 and self.n_alpha_methine > 0:
            score += 0.05
        
        # Methyl presence
        if expected_profile.n_methyl > 0 and self.n_methyl > 0:
            score += 0.05 * completeness
        
        # Methylene count match
        if expected_profile.n_methylene > 0:
            match = min(self.n_methylene, expected_profile.n_methylene)
            score += (match / expected_profile.n_methylene) * 0.05 * completeness
        
        # Methine count match
        if expected_profile.n_methine > 0:
            match = min(self.n_methine, expected_profile.n_methine)
            score += (match / expected_profile.n_methine) * 0.05 * completeness
        
        # ============================================================
        # BÔNUS DE COBERTURA (proporcional)
        # ============================================================
        
        # Bônus baseado em quantos carbonos observamos
        coverage_bonus = completeness * 0.10
        score += coverage_bonus

        # ============================================================
        # BÔNUS ESPECÍFICOS
        # ============================================================
        
        # PHE: bônus alto para sistemas parciais
        if residue_name == 'PHE':
            # NÃO zere o score! Use o score já acumulado como base
            
            # PHE NÃO pode ter mais que 1 methylene
            if self.n_methylene > 1:
                return 0.0
            
            # PHE NÃO pode ter methyl
            if self.n_methyl > 0:
                return 0.0
            
            # CA é obrigatório
            if self.n_alpha_methine >= 1:
                score += 0.40  # ADICIONA ao score existente
            else:
                return 0.0
            
            # Verificar CB
            has_cb = False
            if self.beta_shift is not None and 35 <= self.beta_shift <= 45:
                has_cb = True
            elif self.n_methylene == 1:
                if hasattr(self, '_methylene_shifts'):
                    for c_shift in self._methylene_shifts:
                        if 35 <= c_shift <= 45:
                            has_cb = True
                            break
            
            if has_cb:
                score += 0.30
            else:
                score -= 0.20
            
            # Bônus para padrão clássico
            if self.n_methylene == 1 and self.n_alpha_methine == 1 and has_cb:
                score += 0.20
        
        # SER: precisa de carbonos oxigenados
        elif residue_name == 'SER':
            # SER deve ter pelo menos um carbono oxigenado (CA ou CB)
            if self.n_alpha_methine >= 1:
                score += 0.30
            if self.n_oxygenated_methylene >= 1 or self.n_oxygenated_methine >= 1:
                score += 0.40
            # Penaliza se não tiver nenhum oxigenado
            if self.n_methyl == 0:
                score += 0.20
        
        # GLY: detecção por dois prótons no mesmo carbono
        elif residue_name == 'GLY':
            # GLY tem alpha methylene (CA como CH2)
            if self.n_alpha_methylene >= 1:
                score += 0.60
            # Dois prótons no mesmo carbono é diagnóstico
            if self.n_methylene >= 1 and self.n_methine == 0:
                score += 0.40
            # Caso específico: S10 tem 2 H no mesmo C (~46.6)
            if self.n_methine == 1 and self.n_methylene == 0:
                # Verifica se é um carbono com 2 prótons (GLY)
                # Isso precisa ser detectado upstream
                score += 0.30
        
        # LYS: múltiplos methylenes
        elif residue_name == 'LYS':
            if self.n_methylene >= 3:
                score += 0.50
            elif self.n_methylene >= 2:
                score += 0.30
            else:
                score -= 0.30  # Penaliza se poucos methylenes
        
        # ILE: bônus para o padrão 2 methyls + methine (CORREÇÃO para S11)
        elif residue_name == 'ILE':
            if self.n_methyl >= 2:
                score += 0.40
            if self.n_beta_methine >= 1:
                score += 0.30
            # Padrão clássico: 2 methyls + 1 methine
            if self.n_methyl >= 2 and self.n_methine >= 1:
                score += 0.20
        
        # LEU: sistemas parciais (CORREÇÃO para S5)
        elif residue_name == 'LEU':
            if self.n_methyl >= 1:
                score += 0.30
            if self.n_methylene >= 2:
                score += 0.30
            # LEU parcial pode ter só CA
            if self.n_alpha_methine >= 1:
                score += 0.20
        
        # THR: methyl + CA + CB oxigenado
        elif residue_name == 'THR':
            if self.n_methyl >= 1:
                score += 0.30
            if self.n_alpha_methine >= 1:
                score += 0.30
            if self.n_oxygenated_methine >= 1:
                score += 0.20
        
        # ALA: simples
        elif residue_name == 'ALA':
            if self.n_methyl == 1 and self.n_methine == 1:
                score += 0.60
        
        return max(0.0, min(1.0, score))

# Perfis esperados para cada resíduo - SEM n_methyl_bearing_methine
# core/structural_profile.py

EXPECTED_PROFILES = {
    'ALA': StructuralProfile(
        n_methyl=1, n_methylene=0, n_methine=1,
        n_alpha_methine=1,
        has_alpha_carbon=True
    ),
    'GLY': StructuralProfile(
        n_methyl=0, n_methylene=1, n_methine=0,
        n_alpha_methylene=1,  # CORRIGIDO: GLY tem alpha methylene
        has_alpha_carbon=True
    ),
    'ILE': StructuralProfile(
        n_methyl=2, n_methylene=1, n_methine=2,
        n_alpha_methine=1,
        n_beta_methine=1,
        has_alpha_carbon=True,
        has_branching=True
    ),
    'LEU': StructuralProfile(
        n_methyl=2, n_methylene=1, n_methine=2,
        n_alpha_methine=1,
        n_sidechain_methine=1,  # CG é methine alifático
        has_alpha_carbon=True,
        has_branching=True
    ),
    'LYS': StructuralProfile(
        n_methyl=0, n_methylene=4, n_methine=1,
        n_alpha_methine=1,
        n_aliphatic_methylene=3,  # CB, CG, CD são alifáticos
        has_alpha_carbon=True
    ),
    'PHE': StructuralProfile(
        n_methyl=0, n_methylene=1, n_methine=1, n_aromatic=5,
        n_alpha_methine=1,
        has_alpha_carbon=True
    ),
    'SER': StructuralProfile(
        n_methyl=0, n_methylene=1, n_methine=1,
        n_alpha_methine=1,
        n_oxygenated_methylene=1,  # CORRIGIDO: SER tem oxygenated methylene, não methine
        has_alpha_carbon=True,
        has_oxygenated_carbon=True
    ),
    'THR': StructuralProfile(
        n_methyl=1, n_methylene=0, n_methine=2,
        n_alpha_methine=1,
        n_oxygenated_methine=1,
        has_alpha_carbon=True,
        has_oxygenated_carbon=True
    ),
}

# Ranges para compatibilidade
@dataclass
class ExpectedRange:
    """Define ranges esperados para sistemas parciais"""
    min_methyl: int = 0
    max_methyl: int = 0
    min_methylene: int = 0
    max_methylene: int = 0
    min_methine: int = 0
    max_methine: int = 0
    min_aromatic: int = 0
    max_aromatic: int = 0
    
    def is_compatible(self, observed: StructuralProfile) -> bool:
        """Verifica compatibilidade com ranges"""
        
        if observed.n_methyl < self.min_methyl:
            return False
        if observed.n_methyl > self.max_methyl:
            return False
        
        if observed.n_methylene < self.min_methylene:
            return False
        if observed.n_methylene > self.max_methylene:
            return False
        
        if observed.n_methine < self.min_methine:
            return False
        if observed.n_methine > self.max_methine:
            return False
        
        if observed.n_aromatic < self.min_aromatic:
            return False
        if observed.n_aromatic > self.max_aromatic:
            return False
        
        return True


EXPECTED_RANGES = {
    'ALA': ExpectedRange(min_methyl=1, max_methyl=1, min_methine=1, max_methine=1),
    'GLY': ExpectedRange(min_methylene=1, max_methylene=2),
    'ILE': ExpectedRange(min_methyl=1, max_methyl=2, min_methine=1, max_methine=2),
    'LEU': ExpectedRange(min_methyl=1, max_methyl=2, min_methylene=1, max_methylene=2, min_methine=1, max_methine=2),
    'LYS': ExpectedRange(min_methylene=2, max_methylene=4, min_methine=1, max_methine=1),
    'PHE': ExpectedRange(min_methine=1, max_methine=2, min_aromatic=0, max_aromatic=5),
    'SER': ExpectedRange(min_methylene=1, max_methylene=1, min_methine=1, max_methine=1),
    'THR': ExpectedRange(min_methyl=1, max_methyl=1, min_methine=2, max_methine=2),
}