#graph/candidate_generator_corrected.py
"""Candidate generator corrigido usando os atributos corretos"""

from core.assignment_result import AssignmentResult
from core.residue_candidate import ResidueCandidate
from core.structural_profile import StructuralProfile, EXPECTED_PROFILES


class CandidateGeneratorCorrected:
    def __init__(self, database, sequence_composition=None, max_candidates=5, min_score_threshold=0.05):
        self.database = database
        self.max_candidates = max_candidates
        self.min_score_threshold = min_score_threshold
        self.sequence_composition = sequence_composition or {}

        # ============================================================
        # ADICIONE AQUI OS RANGES ESPECÍFICOS PARA LYS
        # ============================================================
        self.LYS_RANGES = {
            'min_methylene': 2,  # Pode ter apenas CB+CG em sistemas parciais
            'max_methylene': 4,
            'min_methine': 1,
            'max_methine': 1,
        }

        from scoring.topology_scorer import TopologyScorer
        self.topology_scorer = TopologyScorer(database)

    def is_lys_compatible(self, observed: StructuralProfile) -> float:
        """Score específico para LYS com ranges flexíveis"""
        
        if observed.n_methyl > 0:
            return 0.0
        
        if observed.n_methine < 1:
            return 0.0
        
        # Score baseado em methylene
        if observed.n_methylene < 2:
            return 0.3
        elif observed.n_methylene >= 3:
            return 1.0
        else:
            return 0.7

    def generate_candidates(self, spin_system, sequence_composition=None):
        """Gera candidatos com score contínuo de compatibilidade"""
        
        if sequence_composition is None:
            sequence_composition = self.sequence_composition
        
        # Obtém fingerprint e profile
        if hasattr(spin_system, 'fingerprint'):
            fp = spin_system.fingerprint
        elif hasattr(spin_system, 'structural_fingerprint'):
            fp = spin_system.structural_fingerprint
        else:
            fp = {'methyl': 0, 'methylene': 0, 'methine': 0, 'aromatic': 0}
        
        # Obtém profile estrutural
        if hasattr(spin_system, 'carbon_groups'):
            profile = StructuralProfile.from_carbon_groups(spin_system.carbon_groups)
        else:
            profile = StructuralProfile()
        
        #print(f"  Generating candidates for HN={getattr(spin_system, 'hn_shift', '?')}")
        #print(f"    Fingerprint: {fp}")
        #print(f"    Profile: methyl={profile.n_methyl}, methylene={profile.n_methylene}, methine={profile.n_methine}")
        
        candidates = []
        
        for residue_type, max_count in sequence_composition.items():
            if max_count == 0:
                continue
            
            expected_fp = self.get_expected_fingerprint(residue_type)
            expected_profile = EXPECTED_PROFILES.get(residue_type)
            
            if not expected_profile:
                continue
            
            # Score combinado
            score = self.combined_compatibility_score(
                fp, profile, expected_fp, expected_profile, residue_type
            )
            
            if score >= self.min_score_threshold:
                candidate = ResidueCandidate(
                    residue_name=residue_type,
                    score=score,
                    confidence=score * 0.8,
                    matches={}
                )
                candidates.append(candidate)
                print(f"      {residue_type}: score={score:.3f}")
            else:
                print(f"      {residue_type}: score={score:.3f} (below threshold)")
        
        candidates.sort(key=lambda x: x.score, reverse=True)
        print(f"    Generated {len(candidates)} candidates")
        
        return candidates

    def generate_candidates_with_profile(self, spin_system, sequence_composition=None):
        """Gera candidatos usando perfil estrutural detalhado"""
        
        if sequence_composition is None:
            sequence_composition = self.sequence_composition
        
        # Cria perfil estrutural
        if hasattr(spin_system, 'carbon_groups'):
            profile = StructuralProfile.from_carbon_groups(spin_system.carbon_groups)
        else:
            profile = StructuralProfile()
        
        print(f"  Structural Profile for HN={getattr(spin_system, 'hn_shift', '?')}")
        print(f"    n_carbons={profile.n_carbons}, methyl={profile.n_methyl}, "
            f"methylene={profile.n_methylene}, methine={profile.n_methine}")
        print(f"    has_alpha={profile.has_alpha_carbon}, has_oxygenated={profile.has_oxygenated_carbon}")
        
        candidates = []
        
        for residue_type, max_count in sequence_composition.items():
            if max_count == 0:
                continue
            
            expected = EXPECTED_PROFILES.get(residue_type)
            if not expected:
                continue
            
            # Calcula similaridade
            score = profile.score_similarity(expected)
            
            # Compatibilidade básica (methyls não podem exceder muito)
            if profile.n_methyl > expected.n_methyl + 1:
                continue
            
            candidate = ResidueCandidate(
                residue_name=residue_type,
                score=score,
                confidence=score * 0.8,
                matches={}
            )
            candidates.append(candidate)
            print(f"      {residue_type}: similarity={score:.3f}")
        
        candidates.sort(key=lambda x: x.score, reverse=True)
        print(f"    Generated {len(candidates)} candidates")
        
        return candidates

    def _build_fingerprint_from_groups(self, spin_system):
        """Build fingerprint from carbon_groups if available"""
        fp = {'methyl': 0, 'methylene': 0, 'methine': 0, 'aromatic': 0}
        
        if hasattr(spin_system, 'carbon_groups'):
            for group in spin_system.carbon_groups:
                if hasattr(group, 'carbon_type'):
                    carbon_type = group.carbon_type
                    if carbon_type in fp:
                        fp[carbon_type] += 1
        
        return fp

    def debug_candidate_generation(self, spin_system):
        """Debug: mostra por que nenhum candidato é gerado"""
        # Get fingerprint
        if hasattr(spin_system, 'fingerprint'):
            fp = spin_system.fingerprint
        else:
            fp = self._build_fingerprint_from_groups(spin_system)
        
        print(f"\n=== DEBUG: System HN={spin_system.hn_shift}")
        print(f"Fingerprint: {fp}")
        
        # Lista todos os resíduos possíveis
        residues = ['ALA', 'GLY', 'ILE', 'LEU', 'LYS', 'PHE', 'SER', 'THR']
        
        for res in residues:
            expected = self.get_expected_fingerprint(res)
            print(f"  {res}: expected={expected}")
            
            # Verifica compatibilidade
            compatible = self.is_compatible_fingerprint(fp, expected, strict=False)
            print(f"    Compatible: {compatible}")

    def _score_permissive(self, spin_structure, residue_name):
        """Scoring permissivo para sistemas parciais"""
        if hasattr(spin_structure, 'fingerprint'):
            fp = spin_structure.fingerprint
        else:
            fp = self._build_fingerprint_from_groups(spin_structure)
        
        # Scores simplificados baseados apenas em fingerprint
        permissive_scores = {
            'PHE': 0.3 if fp.get('methine', 0) >= 1 and fp.get('methylene', 0) >= 1 else 0.0,
            'TYR': 0.3 if fp.get('methine', 0) >= 1 and fp.get('methylene', 0) >= 1 else 0.0,
            'SER': 0.3 if fp.get('methine', 0) >= 1 and fp.get('methylene', 0) >= 1 else 0.0,
            'THR': 0.3 if fp.get('methine', 0) >= 2 else 0.0,
            'ALA': 0.3 if fp.get('methyl', 0) >= 1 and fp.get('methine', 0) >= 1 else 0.0,
            'LEU': 0.3 if fp.get('methyl', 0) >= 1 else 0.0,
            'ILE': 0.3 if fp.get('methyl', 0) >= 1 else 0.0,
            'VAL': 0.3 if fp.get('methyl', 0) >= 1 else 0.0,
            'LYS': 0.3 if fp.get('methylene', 0) >= 2 else 0.0,
            'GLY': 0.3 if fp.get('methylene', 0) >= 1 else 0.0,
        }
        
        return permissive_scores.get(residue_name, 0.0)
    
    def _find_matches(self, spin_structure, residue_name):
        """Encontra matches específicos entre grupos de carbono e átomos do resíduo"""
        matches = []
        
        # Mapeamento de ambiente para átomo
        env_to_atom = {
            'alpha': 'HA',
            'aromatic_beta': 'HB',
            'oxygenated': 'HB',
            'methyl': 'HG',
            'gamma': 'HG',
            'delta': 'HD',
            'epsilon': 'HE',
            'beta': 'HB'
        }
        
        if hasattr(spin_structure, 'carbon_groups'):
            for group in spin_structure.carbon_groups:
                # Classify environment
                env = CarbonEnvironment.classify(group.carbon_shift, group.average_proton_shift)
                if env and env in env_to_atom:
                    atom_name = env_to_atom[env]
                    matches.append({
                        'atom': atom_name,
                        'H_obs': group.average_proton_shift,
                        'C_obs': group.carbon_shift,
                        'H_score': 0.8,
                        'C_score': 0.8
                    })
        
        return matches
    
    def _calculate_fp_scores(self, fp):
        """Calcula scores baseados no fingerprint estrutural"""
        scores = {
            'GLY': 0.6 if fp.get('methylene', 0) >= 2 else 0.0,
            'SER': 0.5 if fp.get('methine', 0) >= 1 and fp.get('methylene', 0) >= 1 else 0.0,
            'ALA': 0.5 if fp.get('methine', 0) >= 1 and fp.get('methyl', 0) >= 1 else 0.0,
            'PHE': 0.45 if fp.get('methine', 0) >= 1 and fp.get('methylene', 0) >= 1 else 0.0,
            'TYR': 0.45 if fp.get('methine', 0) >= 1 and fp.get('methylene', 0) >= 1 else 0.0,
            'TRP': 0.45 if fp.get('methine', 0) >= 1 and fp.get('methylene', 0) >= 1 else 0.0,
            'HIS': 0.45 if fp.get('methine', 0) >= 1 and fp.get('methylene', 0) >= 1 else 0.0,
            'LEU': 0.5 if fp.get('methyl', 0) >= 2 else 0.3 if fp.get('methyl', 0) >= 1 else 0.0,
            'ILE': 0.5 if fp.get('methyl', 0) >= 2 else 0.3 if fp.get('methyl', 0) >= 1 else 0.0,
            'VAL': 0.5 if fp.get('methyl', 0) >= 2 else 0.3 if fp.get('methyl', 0) >= 1 else 0.0,
            'LYS': 0.4 if fp.get('methylene', 0) >= 3 else 0.2 if fp.get('methylene', 0) >= 2 else 0.0,
            'ARG': 0.4 if fp.get('methylene', 0) >= 3 else 0.2 if fp.get('methylene', 0) >= 2 else 0.0,
            'GLU': 0.3 if fp.get('methylene', 0) >= 2 else 0.0,
            'GLN': 0.3 if fp.get('methylene', 0) >= 2 else 0.0,
            'ASP': 0.3 if fp.get('methylene', 0) >= 1 and fp.get('methine', 0) >= 1 else 0.0,
            'ASN': 0.3 if fp.get('methylene', 0) >= 1 and fp.get('methine', 0) >= 1 else 0.0,
            'THR': 0.35 if fp.get('methine', 0) >= 2 else 0.0,
            'MET': 0.35 if fp.get('methyl', 0) >= 1 else 0.0,
        }
        return scores
    
    def _get_constraint_dict(self, template):
        """Converte template (lista de AtomConstraint) para dicionário"""
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
    
    def _calculate_chemical_scores(self, spin_structure):
        """Calcula scores baseados em matching de shifts químicos"""
        scores = {}
        
        # Para cada spin, tenta encontrar matches com resíduos
        if hasattr(spin_structure, 'spins'):
            for spin in spin_structure.spins:
                h = spin.proton_shift if hasattr(spin, 'proton_shift') else None
                c = spin.carbon_shift if hasattr(spin, 'carbon_shift') else None
                
                if not h and not c:
                    continue
                
                # Testa cada resíduo
                for residue in ['ALA', 'SER', 'GLY', 'PHE', 'TYR', 'LEU', 'ILE', 'VAL', 'LYS', 'THR']:
                    template = self.database.get_residue_template(residue)
                    if not template:
                        continue
                    
                    constraint_dict = self._get_constraint_dict(template)
                    
                    score = 0.0
                    matched = False
                    
                    for atom_name, constraint in constraint_dict.items():
                        # Verifica HA
                        if atom_name in ['HA', 'HA2', 'HA3'] and h:
                            h_range = getattr(constraint, 'h_range', None)
                            if h_range and h_range[0] <= h <= h_range[1]:
                                score += 0.3
                                matched = True
                        
                        # Verifica HB/CB
                        if atom_name in ['HB', 'HB2', 'HB3', 'CB'] and c:
                            c_range = getattr(constraint, 'c_range', None)
                            if c_range and c_range[0] <= c <= c_range[1]:
                                score += 0.3
                                matched = True
                        
                        # Verifica metilas (para LEU/ILE/VAL/ALA)
                        if atom_name in ['HG1', 'HG2', 'HG12', 'HG13', 'HD1', 'HD2'] and c:
                            c_range = getattr(constraint, 'c_range', None)
                            if c_range and c_range[0] <= c <= c_range[1]:
                                score += 0.2
                                matched = True
                    
                    if matched:
                        scores[residue] = scores.get(residue, 0.0) + score
        
        # Normaliza scores (máximo 1.0)
        for residue in scores:
            scores[residue] = min(1.0, scores[residue])
        
        return scores

    def is_compatible_fingerprint(self, observed, expected, residue_name, strict=False):
        """
        Compatibilidade com regras assimétricas e impossibilidades absolutas.
        """
        if strict:
            return observed == expected
        
        # ============================================================
        # REGRAS ESPECÍFICAS POR RESÍDUO (não globais!)
        # ============================================================
        
        # GLY: NÃO pode ter methyl, NÃO pode ter methine
        if expected.get('methylene', 0) == 1 and expected.get('methine', 0) == 0:
            if observed.get('methyl', 0) > 0:
                return False
            if observed.get('methine', 0) > 0:
                return False
        
        # SER: NÃO pode ter methyl
        elif expected.get('methyl', 0) == 0 and expected.get('methylene', 0) == 1 and expected.get('methine', 0) == 1:
            if observed.get('methyl', 0) > 0:
                return False
        
        # THR: NÃO pode ter methylene
        elif expected.get('methyl', 0) == 1 and expected.get('methylene', 0) == 0 and expected.get('methine', 0) == 2:
            if observed.get('methylene', 0) > 0:
                return False
        
        # ALA: NÃO pode ter methylene
        elif expected.get('methyl', 0) == 1 and expected.get('methylene', 0) == 0 and expected.get('methine', 0) == 1:
            if observed.get('methylene', 0) > 0:
                return False
        
        # LYS: NÃO pode ter methyl
        elif expected.get('methyl', 0) == 0 and expected.get('methylene', 0) >= 3:
            if observed.get('methyl', 0) > 0:
                return False
        
        # PHE/TYR/TRP: precisam de pelo menos um sinal aromático ou methine
        elif expected.get('aromatic', 0) > 0:
            if observed.get('aromatic', 0) == 0 and observed.get('methine', 0) < 1:
                return False
        
        # ============================================================
        # TOLERÂNCIA ASSIMÉTRICA (valores padrão)
        # ============================================================
        
        for group_type in ['methyl', 'methylene', 'methine']:
            obs = observed.get(group_type, 0)
            exp = expected.get(group_type, 0)
            
            # Pode faltar (observado < esperado) - OK, sistema parcial
            if obs <= exp:
                continue
            
            # Não pode ter excesso significativo
            if obs > exp + 1:
                return False
        
        return True

    def get_expected_fingerprint(self, residue_type):
        """Retorna fingerprint estrutural esperado para um resíduo"""
        fingerprints = {
            'ALA': {'methyl': 1, 'methylene': 0, 'methine': 1, 'aromatic': 0},
            'GLY': {'methyl': 0, 'methylene': 1, 'methine': 0, 'aromatic': 0},
            'ILE': {'methyl': 2, 'methylene': 0, 'methine': 2, 'aromatic': 0},
            'LEU': {'methyl': 2, 'methylene': 3, 'methine': 1, 'aromatic': 0},
            'LYS': {'methyl': 0, 'methylene': 4, 'methine': 1, 'aromatic': 0},
            'PHE': {'methyl': 0, 'methylene': 0, 'methine': 2, 'aromatic': 5},
            'SER': {'methyl': 0, 'methylene': 1, 'methine': 1, 'aromatic': 0},
            'THR': {'methyl': 1, 'methylene': 0, 'methine': 2, 'aromatic': 0},
        }
        return fingerprints.get(residue_type, {'methyl': 0, 'methylene': 0, 'methine': 0, 'aromatic': 0})

    def structural_penalty_score(self, observed_profile, residue_name):
        """
        Score de penalidades estruturais.
        Retorna 1.0 se compatível, valor menor se houver violações.
        """
        score = 1.0
        
        # GLY: não pode ter methine
        if residue_name == 'GLY' and observed_profile.n_methine > 0:
            return 0.0
        
        # GLY: não pode ter methyl
        if residue_name == 'GLY' and observed_profile.n_methyl > 0:
            return 0.0
        
        # ALA: não pode ter 2 methines
        if residue_name == 'ALA' and observed_profile.n_methine > 1:
            return 0.0
        
        # ALA: não pode ter methylene
        if residue_name == 'ALA' and observed_profile.n_methylene > 0:
            return 0.0
        
        # SER: não pode ter methyl
        if residue_name == 'SER' and observed_profile.n_methyl > 0:
            return 0.0
        
        # LYS: não pode ter methyl
        if residue_name == 'LYS' and observed_profile.n_methyl > 0:
            return 0.0
        
        # THR: precisa ter methyl (penalidade, não fatal)
        if residue_name == 'THR' and observed_profile.n_methyl == 0:
            score -= 0.5
        
        # THR: precisa ter alpha methine
        if residue_name == 'THR' and observed_profile.n_alpha_methine == 0:
            score -= 0.3
        
        # LEU/ILE: precisam ter pelo menos 1 methyl
        if residue_name in ['LEU', 'ILE'] and observed_profile.n_methyl == 0:
            score -= 0.6
        
        return max(0.0, score)

    def chemical_shift_score(self, observed_profile: StructuralProfile, 
                            expected_profile: StructuralProfile) -> float:
        """
        Score baseado em matching de shifts específicos.
        """
        score = 0.0
        n_checks = 0
        
        # Alpha carbon matching
        if expected_profile.alpha_shift and observed_profile.alpha_shift:
            diff = abs(observed_profile.alpha_shift - expected_profile.alpha_shift)
            if diff < 5:
                score += 1.0 - (diff / 10)
            n_checks += 1
        elif expected_profile.alpha_shift:
            n_checks += 1  # Missing alpha carbon penaliza
        
        # Methyl shift matching (para LEU/ILE/VAL/ALA)
        if expected_profile.n_methyl >= 1 and observed_profile.methyl_shifts:
            avg_methyl = sum(observed_profile.methyl_shifts) / len(observed_profile.methyl_shifts)
            if 15 <= avg_methyl <= 22:  # Typical methyl range
                score += 0.5
            n_checks += 1
        
        # Oxygenated carbon (SER/THR)
        if expected_profile.has_oxygenated_carbon:
            if observed_profile.has_oxygenated_carbon:
                score += 0.5
            n_checks += 1
        
        # Branching (LEU/ILE/VAL)
        if expected_profile.has_branching:
            if observed_profile.has_branching:
                score += 0.3
            n_checks += 1
        
        return score / max(n_checks, 1)

    def combined_compatibility_score(self, observed, observed_profile, 
                                    expected, expected_profile,
                                    residue_name):
        """
        Score combinado para compatibilidade.
        """
        
        # 1. Score de perfil estrutural (50%)
        profile_score = observed_profile.score_similarity(expected_profile, residue_name)
        if profile_score <= 0:
            return 0.0
        
        # 2. Score de penalidades estruturais (25%)
        # CORREÇÃO: Usar o nome correto do método
        penalty_score = self.structural_penalty_score(observed_profile, residue_name)
        
        # 3. Score químico (25%)
        chem_score = self.chemical_shift_score(observed_profile, expected_profile)
        
        # Score final
        final_score = (profile_score * 0.50 + 
                    penalty_score * 0.25 +
                    chem_score * 0.25)
        
        return final_score

    def topology_score(self, observed_profile: StructuralProfile,
                        expected_profile: StructuralProfile,
                        residue_name: str) -> float:
        """
        Score baseado em topologia e conectividade.
        Peso maior que chemical shift individual.
        """
        score = 0.0
        n_checks = 0
        
        # ============================================================
        # TOPOLOGIA > CHEMICAL SHIFT
        # ============================================================
        
        # THR: padrão CH3-CA-CB (methyl ligado a beta)
        if residue_name == 'THR':
            # Tem methyl, tem alpha, tem beta
            if observed_profile.n_methyl >= 1 and observed_profile.n_alpha_methine >= 1:
                score += 0.6
            n_checks += 1
            
            # Padrão de shifts: CH3 ~20, CA ~58-65, CB ~45-55
            if observed_profile.methyl_shifts:
                avg_methyl = sum(observed_profile.methyl_shifts) / len(observed_profile.methyl_shifts)
                if 15 <= avg_methyl <= 25:
                    score += 0.2
            if observed_profile.alpha_shift and 55 <= observed_profile.alpha_shift <= 68:
                score += 0.2
            n_checks += 2
        
        # ALA: padrão CH3-CA (methyl ligado diretamente)
        elif residue_name == 'ALA':
            if observed_profile.n_methyl >= 1 and observed_profile.n_alpha_methine >= 1:
                score += 0.8
            n_checks += 1
            
            # Methyl shift típico de ALA: 15-22
            if observed_profile.methyl_shifts:
                avg_methyl = sum(observed_profile.methyl_shifts) / len(observed_profile.methyl_shifts)
                if 15 <= avg_methyl <= 22:
                    score += 0.2
            n_checks += 1
        
        # LEU: padrão ramificado (múltiplos methyls)
        elif residue_name == 'LEU':
            if observed_profile.n_methyl >= 2:
                score += 0.5
            elif observed_profile.n_methyl >= 1:
                score += 0.3
            n_checks += 1
            
            # Múltiplos methylenes indicam LEU
            if observed_profile.n_methylene >= 2:
                score += 0.3
            n_checks += 1
        
        # ILE: padrão ramificado com beta methine
        elif residue_name == 'ILE':
            if observed_profile.n_methyl >= 2:
                score += 0.4
            elif observed_profile.n_methyl >= 1:
                score += 0.2
            n_checks += 1
            
            # ILE tem beta methine característico
            if observed_profile.n_beta_methine >= 1:
                score += 0.4
            n_checks += 1
        
        # LYS: cadeia longa de methylenes
        elif residue_name == 'LYS':
            if observed_profile.n_methylene >= 3:
                score += 0.6
            elif observed_profile.n_methylene >= 2:
                score += 0.3
            n_checks += 1
            
            if observed_profile.n_alpha_methine >= 1:
                score += 0.2
            n_checks += 1
        
        return score / max(n_checks, 1)