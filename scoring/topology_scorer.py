# scoring/topology_scorer.py
"""Scoring baseado em ambientes de carbono e topologia"""

from typing import List
from core.spin import CarbonGroup
from core.spin import CarbonEnvironment

class TopologyScorer:
    """Scorer que usa ambientes de carbono para reconhecimento de padrões"""
    
    def __init__(self, database):
        self.database = database
    
    def score_residue(self, residue_type: str, carbon_groups: List[CarbonGroup]) -> float:
        """Calculate topology-based score with penalties."""
        # Your existing base scoring
        base_score = self._calculate_base_score(residue_type, carbon_groups)
        
        # Add environment-specific rewards
        env_score = self._score_environments(carbon_groups, residue_type)
        base_score += env_score
        
        # Apply penalties for missing expected features
        penalized_score = self._apply_penalties(residue_type, carbon_groups, base_score)
        
        return penalized_score

    def _apply_penalties(self, residue_type: str, carbon_groups: List[CarbonGroup], score: float) -> float:
        """Apply strong penalties for missing expected groups."""
        
        methyl_count = sum(1 for g in carbon_groups if g.carbon_type == 'methyl')
        
        # LEU: needs 2 methyls
        if residue_type == 'LEU':
            if methyl_count < 2:
                score *= 0.3  # Strong penalty
                
        # ILE: needs 2 methyls with specific shifts
        elif residue_type == 'ILE':
            if methyl_count < 2:
                score *= 0.4
            # Check specific methyl shifts
            methyl_shifts = [g.carbon_shift for g in carbon_groups if g.carbon_type == 'methyl']
            if not any(11 <= s <= 18 for s in methyl_shifts):
                score *= 0.5
                
        # LYS: should have NO methyls, but many methylenes
        elif residue_type == 'LYS':
            if methyl_count > 0:
                score *= 0.2  # Very strong penalty
            # Reward multiple methylenes
            methylene_count = sum(1 for g in carbon_groups if g.carbon_type == 'methylene')
            if methylene_count >= 3:
                score *= 1.5
                
        # THR: specific beta carbon pattern
        elif residue_type == 'THR':
            # Should have exactly one methyl and one methine (beta carbon)
            if methyl_count != 1:
                score *= 0.5
                
        return score
    
    def _score_environments(self, carbon_groups: List[CarbonGroup], residue_type: str) -> float:
        """Score based on specific chemical environments."""
        total = 0.0
        
        for group in carbon_groups:
            # EPSILON carbon in LYS (CE)
            if 38 <= group.carbon_shift <= 42:
                if residue_type == 'LYS':
                    total += 2.0
                    
            # Beta carbon in THR (CB)  
            elif 45 <= group.carbon_shift <= 55:
                if residue_type == 'THR' and group.average_proton_shift and group.average_proton_shift > 3.5:
                    total += 2.0
                    
            # Gamma methyl in ILE (CG2)
            elif 11 <= group.carbon_shift <= 16:
                if residue_type == 'ILE':
                    total += 2.0
                    
            # Alpha carbon in GLY (CA)
            elif 40 <= group.carbon_shift <= 45:
                if residue_type == 'GLY' and group.proton_count == 2:
                    total += 2.0
        
        return total

    def _score_by_rules(self, residue, environments, carbon_types, 
                        methyl_count, methylene_count, methine_count, structure):
        """Aplica regras específicas para cada resíduo"""
        
        # PHE e TYR
        if residue in ['PHE', 'TYR']:
            score = 0.0
            # Precisa de alpha e aromatic_beta
            if CarbonEnvironment.ALPHA in environments:
                score += 0.3
            if CarbonEnvironment.AROMATIC_BETA in environments:
                score += 0.4
            # Pode ter mais um CH
            if methine_count >= 2:
                score += 0.1
            return score + 0.2  # Bônus base
        
        # SER
        if residue == 'SER':
            score = 0.0
            if CarbonEnvironment.ALPHA in environments:
                score += 0.3
            if CarbonEnvironment.OXYGENATED in environments:
                score += 0.4
            if methine_count >= 2:
                score += 0.1
            return score + 0.2
        
        # THR
        if residue == 'THR':
            score = 0.0
            if CarbonEnvironment.ALPHA in environments:
                score += 0.3
            if CarbonEnvironment.OXYGENATED in environments:
                score += 0.3
            if methyl_count >= 1:
                score += 0.2
            if methine_count >= 2:
                score += 0.1
            return score + 0.1
        
        # ALA
        if residue == 'ALA':
            score = 0.0
            if CarbonEnvironment.ALPHA in environments:
                score += 0.4
            if methyl_count >= 1:
                score += 0.4
            # ALA NÃO pode ter mais de 2 grupos
            if len(structure.carbon_groups) > 2:
                score -= 0.3
            return score + 0.2
        
        # GLY
        if residue == 'GLY':
            score = 0.0
            # GLY tem CA em 43-46
            for group in structure.carbon_groups:
                if group.environment == CarbonEnvironment.ALPHA:
                    if 43 <= group.carbon_shift <= 47:
                        score += 0.6
                    else:
                        score += 0.2
            # GLY não pode ter methyl ou methylene além do alpha
            if methyl_count > 0 or methylene_count > 1:
                score -= 0.3
            return score + 0.2
        
        # LEU
        if residue == 'LEU':
            score = 0.0
            if CarbonEnvironment.ALPHA in environments:
                score += 0.2
            # LEU tem CG característico ~41
            for group in structure.carbon_groups:
                if group.environment == CarbonEnvironment.GAMMA:
                    if 40 <= group.carbon_shift <= 44:
                        score += 0.4
                    else:
                        score += 0.2
            if methyl_count >= 1:
                score += 0.2
            return score + 0.2
        
        # ILE
        if residue == 'ILE':
            score = 0.0
            if CarbonEnvironment.ALPHA in environments:
                score += 0.2
            # ILE tem metilas assimétricas
            methyl_shifts = [g.carbon_shift for g in structure.carbon_groups 
                           if g.carbon_type == 'methyl']
            if len(methyl_shifts) >= 2:
                if max(methyl_shifts) - min(methyl_shifts) > 3:
                    score += 0.4
                else:
                    score += 0.2
            elif methyl_count >= 1:
                score += 0.2
            if methine_count >= 2:
                score += 0.1
            return score + 0.2
        
        # VAL
        if residue == 'VAL':
            score = 0.0
            if CarbonEnvironment.ALPHA in environments:
                score += 0.2
            # VAL tem duas metilas equivalentes
            if methyl_count >= 2:
                score += 0.4
            elif methyl_count >= 1:
                score += 0.2
            if methine_count >= 2:
                score += 0.2
            return score + 0.2
        
        # LYS
        if residue == 'LYS':
            score = 0.0
            if CarbonEnvironment.ALPHA in environments:
                score += 0.2
            # LYS tem epsilon característico
            for group in structure.carbon_groups:
                if group.environment == CarbonEnvironment.EPSILON:
                    if 38 <= group.carbon_shift <= 42:
                        score += 0.4
                    else:
                        score += 0.2
            # Cadeia longa (múltiplos CH2)
            if methylene_count >= 3:
                score += 0.2
            elif methylene_count >= 2:
                score += 0.1
            return score + 0.2
        
        return 0.0
