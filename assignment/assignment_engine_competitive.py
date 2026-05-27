# assignment/assignment_engine_competitive.py
from assignment.confidence_estimator import ConfidenceEstimator
from assignment.ambiguity_resolver import AmbiguityResolver
from scoring.contextual_scorer import ContextualScorer
from scoring.aromatic_handler import AromaticHandler
from core.assignment_result import AssignmentResult
from collections import defaultdict
import os
import math

class AssignmentEngineCompetitive:
    
    def __init__(self, chemical_database, sequence_file=None):
        self.database = chemical_database
        self.confidence_estimator = ConfidenceEstimator()
        self.ambiguity_resolver = AmbiguityResolver()
        
        # CORREÇÃO: Carrega composição da sequência e guarda como atributo
        self.expected_composition = self._load_sequence_composition(sequence_file)
        self.total_residues = sum(self.expected_composition.values())
        self.contextual_scorer = ContextualScorer(chemical_database, stage='pre_noesy')
        
        print(f"Loaded sequence composition: {dict(self.expected_composition)}")
        print(f"Total residues: {self.total_residues}")

    def _load_sequence_composition(self, sequence_file):
        """Lê sequence.txt e retorna dicionário {residuo: count}"""
        composition = {}
        
        if not sequence_file or not os.path.exists(sequence_file):
            print(f"Warning: Sequence file {sequence_file} not found. Using default limits.")
            return self._get_default_composition()
        
        valid_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                         'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                         'THR', 'TRP', 'TYR', 'VAL'}
        
        with open(sequence_file, 'r') as f:
            for line in f:
                residue = line.strip().upper()
                if residue and residue in valid_residues:
                    composition[residue] = composition.get(residue, 0) + 1
        
        return composition
    
    def _get_default_composition(self):
        """Composição padrão caso não encontre o arquivo"""
        return {
            'PHE': 3, 'ILE': 1, 'GLY': 2, 'LEU': 3, 
            'LYS': 1, 'THR': 1, 'ALA': 1, 'SER': 1
        }
    
    # ============================================================
    # VERSÃO ÚNICA E CORRIGIDA DE _original_scoring
    # ============================================================
    
    def _original_scoring(self, structure, candidate, usage_counts):
        """
        Score combinado usando o perfil estrutural diretamente.
        """
        from core.structural_profile import StructuralProfile, EXPECTED_PROFILES
        
        # Obtém o perfil observado
        if hasattr(structure, 'carbon_groups'):
            observed_profile = StructuralProfile.from_carbon_groups(structure.carbon_groups)
        else:
            observed_profile = StructuralProfile()
        
        # Obtém o perfil esperado
        expected_profile = EXPECTED_PROFILES.get(candidate)
        if not expected_profile:
            return 0.0
        
        # ============================================================
        # USA O SCORE_SIMILARITY DIRETAMENTE
        # ============================================================
        profile_score = observed_profile.score_similarity(expected_profile, candidate)
        
        # ============================================================
        # APENAS COMPOSITION PENALTY É ADICIONADO
        # ============================================================
        comp_penalty = self._calculate_composition_penalty(candidate, usage_counts)
        
        # Score final: profile_score reduzido pela composition penalty
        final_score = profile_score - comp_penalty * 0.3
        
        # Bônus para sistemas parciais
        if candidate == 'PHE' and observed_profile.n_methine >= 1:
            final_score = min(1.0, final_score + 0.1)
        
        return max(0.0, min(1.0, final_score))

    def _get_constraint_dict(self, template):
        """Converte template para dicionário {atom_name: constraint}"""
        if isinstance(template, dict):
            return template
        
        constraint_dict = {}
        
        if hasattr(template, 'atoms'):
            atoms = template.atoms
        elif isinstance(template, list):
            atoms = template
        else:
            return {}
        
        for constraint in atoms:
            name = getattr(constraint, 'name', None)
            if name:
                constraint_dict[name] = constraint
        
        return constraint_dict

    def _score_candidate(self, structure, candidate, usage_counts):
        """
        Calcula score para um candidato usando scoring contextual
        """
        original_score = self._original_scoring(structure, candidate, usage_counts)
        
        # DEBUG
        if candidate == 'PHE':
            print(f"      PHE original_score: {original_score:.3f}")
        
        fp_score = self.contextual_scorer.calculate_fingerprint_score(
            structure.spins, candidate
        )
        
        if candidate == 'PHE':
            print(f"      PHE fp_score: {fp_score:.3f}")
        
        coverage_penalty = self.contextual_scorer.calculate_coverage_penalty(
            structure.spins, candidate
        )
        
        if candidate == 'PHE':
            print(f"      PHE coverage_penalty: {coverage_penalty:.3f}")
        
        combined_score = (original_score * 0.5 + 
                        fp_score * 0.3 + 
                        (1.0 - coverage_penalty) * 0.2)
        
        if candidate == 'PHE':
            print(f"      PHE combined_score before aromatic: {combined_score:.3f}")
        
        combined_score = AromaticHandler.adjust_aromatic_score(
            combined_score, structure, candidate
        )
        
        if candidate == 'PHE':
            print(f"      PHE final_score: {combined_score:.3f}")
        
        return min(1.0, max(0.0, combined_score))
    
    # ============================================================
    # VERSÃO ÚNICA DE _calculate_composition_penalty
    # ============================================================
    
    def _calculate_composition_penalty(self, candidate, usage_counts):
        """Penalidade baseada na composição da sequência - VERSÃO ÚNICA"""
        max_allowed = self.expected_composition.get(candidate, 0)
        
        if max_allowed == 0:
            return 0.8  # Penalidade alta para resíduos não esperados
        
        already_used = usage_counts.get(candidate, 0)
        
        if already_used >= max_allowed:
            return 0.5  # Penalidade por exceder limite
        
        # Penalidade suave baseada na disponibilidade
        remaining = max_allowed - already_used
        availability = remaining / max_allowed
        
        return 0.2 * (1.0 - availability)
    
    # ============================================================
    # CORREÇÃO: _calculate_chemical_score com ranges flexíveis
    # ============================================================
    
    def _calculate_chemical_score(self, structure, candidate):
        """Score químico com ranges flexíveis e score gaussiano"""
        if not structure.carbon_groups:
            return 0.0
        
        # Templates corrigidos - ADICIONAR CG para LEU
        expected_shifts = {
            'ALA': {'CA': (50, 55), 'CB': (15, 22)},
            'GLY': {'CA': (40, 45)},
            'ILE': {'CA': (55, 65), 'CB': (35, 45), 'CG2': (11, 16), 'CG1': (25, 35)},
            'LEU': {'CA': (50, 60), 'CB': (38, 45), 'CG': (24, 30), 'CD1': (20, 25), 'CD2': (20, 25)},
            'LYS': {'CA': (52, 58), 'CB': (28, 35), 'CG': (22, 28), 'CD': (25, 32), 'CE': (38, 42)},
            'PHE': {'CA': (55, 62), 'CB': (35, 42)},
            'SER': {'CA': (55, 62), 'CB': (60, 65)},
            'THR': {'CA': (58, 65), 'CB': (45, 55), 'CG2': (15, 25)},
        }
        
        template = expected_shifts.get(candidate, {})
        if not template:
            return 0.0
        
        total_score = 0.0
        n_matches = 0
        
        for group in structure.carbon_groups:
            c = group.carbon_shift
            best_match_score = 0.0
            
            for atom_name, (c_min, c_max) in template.items():
                # Score gaussiano mesmo fora do range
                center = (c_min + c_max) / 2
                sigma = (c_max - c_min) / 2
                
                if sigma > 0:
                    distance = abs(c - center) / sigma
                    # Gaussiana: exp(-distance^2) - cai suavemente
                    import math
                    match_score = math.exp(-distance ** 2)
                else:
                    match_score = 1.0 if c_min <= c <= c_max else 0.0
                
                best_match_score = max(best_match_score, match_score)
            
            if best_match_score > 0:
                total_score += best_match_score
                n_matches += 1
        
        if n_matches == 0:
            return 0.0
        
        base_score = total_score / max(len(structure.carbon_groups), 1)
        
        # Coverage bonus corrigido (só carbonos esperados)
        n_carbons_expected = sum(1 for k in template if k.startswith('C'))
        if n_carbons_expected > 0:
            coverage_bonus = min(0.3, len(structure.carbon_groups) / n_carbons_expected * 0.15)
        else:
            coverage_bonus = 0.0
        
        return min(1.0, base_score + coverage_bonus)
    
    # ============================================================
    # REMOVER ou SIMPLIFICAR connectivity_score
    # ============================================================
    
    def _calculate_connectivity_score(self, structure, candidate):
        """
        Score de conectividade - TORNADO OPCIONAL e conservador
        Só aplica se tiver conectividade experimental real
        """
        # Por enquanto, retorna neutro para não enviesar
        # Só ativar quando tiver conectividade real de NOESY ou TOCSY
        return 0.5  # Neutro
    
    # ============================================================
    # CORREÇÃO: environment_score com ranges ajustados
    # ============================================================
    
    def _calculate_environment_score(self, structure, candidate):
        """Score baseado em ambientes químicos específicos"""
        if not structure.carbon_groups:
            return 0.0
        
        score = 0.0
        
        for group in structure.carbon_groups:
            c = group.carbon_shift
            h = group.average_proton_shift
            
            # THR: beta carbon
            if candidate == 'THR' and 45 <= c <= 55 and h and h > 3.5:
                score += 0.4
            # ALA: methyl
            elif candidate == 'ALA' and 15 <= c <= 22:
                score += 0.4
            # ILE: gamma methyl
            elif candidate == 'ILE' and 11 <= c <= 16:
                score += 0.3
            # LEU: methyl (característico)
            elif candidate == 'LEU' and 20 <= c <= 25:
                score += 0.2
            # GLY: alpha carbon
            elif candidate == 'GLY' and 40 <= c <= 45 and group.proton_count == 2:
                score += 0.5
            # LYS: epsilon carbon
            elif candidate == 'LYS' and 38 <= c <= 42:
                score += 0.3
        
        return min(1.0, score)
    
    # ============================================================
    # CORREÇÃO: fingerprint_score com regras conservadoras
    # ============================================================
    
    def _calculate_fingerprint_score(self, structure, candidate):
        """Fingerprint score com penalidade para aromaticos ausentes"""
        expected_fingerprints = {
            'ALA': {'methyl': 1, 'methylene': 0, 'methine': 1, 'aromatic': 0},
            'GLY': {'methyl': 0, 'methylene': 1, 'methine': 0, 'aromatic': 0},
            'ILE': {'methyl': 2, 'methylene': 0, 'methine': 2, 'aromatic': 0},
            'LEU': {'methyl': 2, 'methylene': 3, 'methine': 1, 'aromatic': 0},
            'LYS': {'methyl': 0, 'methylene': 4, 'methine': 1, 'aromatic': 0},
            'PHE': {'methyl': 0, 'methylene': 0, 'methine': 2, 'aromatic': 5},
            'SER': {'methyl': 0, 'methylene': 1, 'methine': 1, 'aromatic': 0},
            'THR': {'methyl': 1, 'methylene': 0, 'methine': 2, 'aromatic': 0},
        }
        
        expected = expected_fingerprints.get(candidate, {})
        observed = structure.fingerprint
        
        if not expected:
            return 0.0
        
        score = 0.0
        n_groups = 0
        
        for group_type in ['methyl', 'methylene', 'methine']:
            exp = expected.get(group_type, 0)
            obs = observed.get(group_type, 0)
            
            if exp > 0:
                # Penaliza se tiver muito menos que o esperado
                if obs < exp:
                    score += obs / exp
                else:
                    score += 1.0
                n_groups += 1
        
        # Penalidade ESPECIAL para aromaticos ausentes em PHE/TYR/TRP
        if candidate in ['PHE', 'TYR', 'TRP', 'HIS']:
            if observed.get('aromatic', 0) == 0:
                score *= 0.3  # Penalidade forte
        
        return score / max(n_groups, 1)

    def assign_all_systems(self, spin_structures, candidate_generator, max_iterations=3):
        """
        Atribui todos os spin systems com reatribuição competitiva.
        """
        n_systems = len(spin_structures)
        
        assignments = [None] * n_systems
        usage_counts = defaultdict(int)
        
        impossible_systems = set()
        
        for iteration in range(max_iterations):
            print(f"\n{'='*60}")
            print(f"ITERATION {iteration + 1}")
            print(f"{'='*60}")
            
            unassigned_indices = [i for i, a in enumerate(assignments) if a is None 
                                and i not in impossible_systems]
            
            if not unassigned_indices:
                print("All systems assigned!")
                break
            
            print(f"Unassigned systems: {len(unassigned_indices)}")
            
            all_candidates = []
            
            for idx in unassigned_indices:
                structure = spin_structures[idx]
                print(f"\n  System {idx+1} (HN={structure.hn_shift:.3f}):")
                
                candidates = candidate_generator.generate_candidates(structure)
                
                if not candidates:
                    print(f"    No candidates generated - marking as impossible")
                    impossible_systems.add(idx)
                    continue
                
                # Calcula scores para cada candidato
                scored_candidates = []
                for cand in candidates:
                    score = self._score_candidate(structure, cand.residue_name, usage_counts)
                    cand.score = score
                    scored_candidates.append((score, cand))
                
                scored_candidates.sort(reverse=True, key=lambda x: x[0])
                
                print(f"    Top candidates:")
                for i, (score, cand) in enumerate(scored_candidates[:3]):
                    print(f"      {i+1}. {cand.residue_name} (score={score:.3f})")
                
                if scored_candidates:
                    best_score, best_candidate = scored_candidates[0]
                    all_candidates.append((idx, best_candidate, best_score, best_candidate.residue_name))
            
            if not all_candidates:
                print("No candidates generated for any system!")
                break
            
            all_candidates.sort(key=lambda x: x[2], reverse=True)
            
            newly_assigned = 0
            used_in_iteration = defaultdict(int)
            
            for spin_idx, candidate, score, residue in all_candidates:
                max_allowed = self.expected_composition.get(residue, 0)
                current_used = usage_counts.get(residue, 0) + used_in_iteration.get(residue, 0)
                
                if current_used < max_allowed:
                    assignments[spin_idx] = candidate
                    used_in_iteration[residue] += 1
                    newly_assigned += 1
                    print(f"  ✓ System {spin_idx+1} -> {residue} (score={score:.3f})")
                else:
                    print(f"  ✗ System {spin_idx+1} -> {residue} REJECTED (limit {max_allowed} reached)")
            
            for residue, count in used_in_iteration.items():
                usage_counts[residue] += count
            
            print(f"\n  Assigned {newly_assigned} systems in iteration {iteration+1}")
            
            if newly_assigned == 0 and impossible_systems:
                print(f"  No assignments possible - {len(impossible_systems)} systems impossible to assign")
                break
        
        results = []
        for i, candidate in enumerate(assignments):
            if candidate:
                results.append(candidate)
            else:
                from core.assignment_result import AssignmentResult
                results.append(AssignmentResult(
                    residue_name="UNKNOWN",
                    score=0.0,
                    confidence=0.0,
                    matches=[]
                ))
        
        return results

    def print_system_details(self, spin_structure, result, system_id):
        """Imprime detalhes de um sistema de spin com seus sinais atribuídos"""
        
        print(f"\n{'='*60}")
        print(f"SYSTEM {system_id}: HN={spin_structure.hn_shift:.3f} -> {result.residue_name} (score={result.score:.3f}, conf={result.confidence:.3f})")
        print(f"{'='*60}")
        
        spins_by_type = {
            'methyl': [],
            'methylene': [],
            'methine': [],
            'aromatic': [],
            'unknown': []
        }
        
        for spin in spin_structure.spins:
            if spin.carbon_type:
                spins_by_type[spin.carbon_type].append(spin)
            else:
                spins_by_type['unknown'].append(spin)
        
        print(f"\n📊 STRUCTURAL FINGERPRINT: {spin_structure.fingerprint}")
        
        print(f"\n🔬 ASSIGNED SIGNALS:")
        
        if spins_by_type['methyl']:
            print(f"\n  Methyl groups (CH3):")
            for spin in spins_by_type['methyl']:
                c_info = f"C={spin.carbon_shift:.1f}" if spin.carbon_shift else "C=?"
                print(f"    • H={spin.proton_shift:.3f} ppm, {c_info}")
        
        if spins_by_type['methylene']:
            print(f"\n  Methylene groups (CH2):")
            for spin in spins_by_type['methylene']:
                c_info = f"C={spin.carbon_shift:.1f}" if spin.carbon_shift else "C=?"
                print(f"    • H={spin.proton_shift:.3f} ppm, {c_info}")
        
        if spins_by_type['methine']:
            print(f"\n  Methine groups (CH):")
            for spin in spins_by_type['methine']:
                c_info = f"C={spin.carbon_shift:.1f}" if spin.carbon_shift else "C=?"
                print(f"    • H={spin.proton_shift:.3f} ppm, {c_info}")
        
        if spins_by_type['aromatic']:
            print(f"\n  Aromatic protons:")
            for spin in spins_by_type['aromatic']:
                c_info = f"C={spin.carbon_shift:.1f}" if spin.carbon_shift else "C=?"
                print(f"    • H={spin.proton_shift:.3f} ppm, {c_info}")
        
        if spins_by_type['unknown']:
            print(f"\n  Unclassified spins:")
            for spin in spins_by_type['unknown']:
                c_info = f"C={spin.carbon_shift:.1f}" if spin.carbon_shift else "C=?"
                print(f"    • H={spin.proton_shift:.3f} ppm, {c_info}")
        
        if result.residue_name == "UNKNOWN":
            print(f"\n⚠️  NOT ASSIGNED - Possible reasons:")
            print(f"    • No candidate reached min_score_threshold")
            print(f"    • All candidates exceeded sequence limits")
            print(f"    • Low chemical shift matching")
