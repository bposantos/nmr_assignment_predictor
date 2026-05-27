# scoring/noe_sequential_scorer.py

from collections import defaultdict
from core.noe_hypothesis import NoeHypothesis, SequentialEdge


class NoeSequentialScorer:
    """Scorer para relações sequenciais baseadas em NOESY"""
    
    # Tolerâncias por tipo de próton
    TOLERANCES = {
        'HN': 0.03,
        'HA': 0.04,
        'HB': 0.06,
        'HG': 0.08,
        'HD': 0.08,
    }
    
    def __init__(self, noesy_peaks):
        self.noesy_peaks = noesy_peaks
    
    def generate_hypotheses(self, spin_systems, assignments):
        """Gera hipóteses de NOE entre sistemas de spin"""
        hypotheses = []
        shift_to_system = self._build_shift_map(spin_systems, assignments)
        
        print(f"\n  === Processing NOESY peaks ===")
        print(f"  Total NOESY peaks: {len(self.noesy_peaks)}")
        
        ha_hn_count = 0
        
        for peak in self.noesy_peaks:
            shift1 = float(peak.w1)
            shift2 = float(peak.w2)
            intensity = float(getattr(peak, 'intensity', 1.0))
            
            matches1 = self._find_system_for_shift(shift1, shift_to_system)
            matches2 = self._find_system_for_shift(shift2, shift_to_system)
            
            for sys1, atom1, _ in matches1:
                for sys2, atom2, _ in matches2:
                    if sys1 == sys2:
                        continue
                    
                    # Conta HA-HN
                    if (atom1 == 'HA' and atom2 == 'HN') or (atom1 == 'HN' and atom2 == 'HA'):
                        ha_hn_count += 1
                    
                    hypothesis = NoeHypothesis(
                        donor_system_id=sys1,
                        donor_atom=atom1,
                        acceptor_system_id=sys2,
                        acceptor_atom=atom2,
                        intensity=intensity,
                        peak_shift1=shift1,
                        peak_shift2=shift2
                    )
                    hypotheses.append(hypothesis)
        
        print(f"  HA-HN hypotheses found: {ha_hn_count}")
        return hypotheses
    
    def _build_shift_map(self, spin_systems, assignments):
        """Constrói mapeamento de shifts para sistemas e átomos"""
        shift_map = defaultdict(list)
        
        print(f"\n  === Mapped HA shifts for each system ===")
        
        for i, (ss, ass) in enumerate(zip(spin_systems, assignments)):
            if ass and ass.residue_name != 'UNKNOWN':
                print(f"    S{i+1} ({ass.residue_name}):")
                
                # HN shift
                if ss.hn_shift:
                    rounded = round(ss.hn_shift, 3)
                    shift_map[rounded].append((i, 'HN', ss.hn_shift))
                    print(f"      HN: {ss.hn_shift:.3f}")
                
                # Carbon groups
                if hasattr(ss, 'carbon_groups'):
                    for group in ss.carbon_groups:
                        atom_type = self._determine_atom_type(group)
                        h_shift = group.average_proton_shift
                        if atom_type and h_shift:
                            rounded = round(h_shift, 3)
                            shift_map[rounded].append((i, atom_type, h_shift))
                            if atom_type == 'HA':
                                print(f"      HA: H={h_shift:.3f}, C={group.carbon_shift:.1f}")
                            elif atom_type == 'HB':
                                print(f"      HB: H={h_shift:.3f}, C={group.carbon_shift:.1f}")
        
        return shift_map
    
    def _determine_atom_type(self, group):
        """Determina o tipo de átomo baseado no carbon group"""
        c_type = group.carbon_type
        c_shift = group.carbon_shift
        h_shift = group.average_proton_shift

        # HA: próton na faixa 3.5-5.0 com carbono > 50 ou GLY
        if h_shift and 3.5 <= h_shift <= 5.0:
            if c_shift > 50:
                return 'HA'
            if c_type == 'methylene' and 40 <= c_shift <= 50:
                return 'HA'
        
        # HB: prótons alifáticos típicos
        if h_shift and 1.0 <= h_shift <= 3.5:
            return 'HB'
        
        # HN: prótons amida
        if h_shift and 7.0 <= h_shift <= 9.0:
            return 'HN'
        
        if c_type == 'methyl':
            return 'HB' if c_shift < 25 else 'HG'
        elif c_type == 'methylene':
            if 40 <= c_shift <= 50:
                return 'HA'
            else:
                return 'HB'
        elif c_type == 'methine':
            if 50 <= c_shift <= 65:
                return 'HA'
            else:
                return 'HB'
        return None
    
    def _find_system_for_shift(self, shift, shift_map):
        """Encontra sistema para um shift com tolerância"""
        matches = []
        # Usa tolerância para HN (0.03) como padrão
        tolerance = self.TOLERANCES.get('HN', 0.03)
        
        for tol_shift, systems in shift_map.items():
            if abs(shift - tol_shift) <= tolerance:
                matches.extend(systems)
        return matches
    
    def build_sequential_graph(self, spin_systems, assignments):
        """Constrói grafo sequencial baseado em NOEs"""
        
        print("\n=== DEBUG: Construindo grafo sequencial ===")
        
        # Mostra sistemas com assignments
        print("\n  Systems with assignments:")
        for i, (ss, ass) in enumerate(zip(spin_systems, assignments)):
            if ass and ass.residue_name != 'UNKNOWN':
                print(f"    S{i+1}: {ass.residue_name} (HN={ss.hn_shift:.3f})")
            else:
                print(f"    S{i+1}: UNKNOWN (HN={ss.hn_shift:.3f})")
        
        # Gera hipóteses
        hypotheses = self.generate_hypotheses(spin_systems, assignments)
        
        print(f"\n  Total hypotheses generated: {len(hypotheses)}")
        
        # Conta tipos
        from collections import defaultdict
        by_type = defaultdict(int)
        for h in hypotheses:
            by_type[(h.donor_atom, h.acceptor_atom)] += 1
        
        print(f"  Hypothesis types:")
        for (d, a), count in sorted(by_type.items(), key=lambda x: -x[1])[:10]:
            print(f"    {d} -> {a}: {count}")
        
        # Constrói edges
        edges = {}
        for h in hypotheses:
            key = (h.donor_system_id, h.acceptor_system_id)
            
            if key not in edges:
                from core.noe_hypothesis import SequentialEdge
                edges[key] = SequentialEdge(
                    source_system_id=h.donor_system_id,
                    target_system_id=h.acceptor_system_id
                )
            
            edge = edges[key]
            edge.supporting_peaks.append(h)
            
            if h.donor_atom == 'HN' and h.acceptor_atom == 'HN':
                edge.hn_hn_score = max(edge.hn_hn_score, h.confidence)
            elif h.donor_atom == 'HA' and h.acceptor_atom == 'HN':
                edge.ha_hn_score = max(edge.ha_hn_score, h.confidence)
            elif h.donor_atom == 'HN' and h.acceptor_atom == 'HA':
                edge.ha_hn_score = max(edge.ha_hn_score, h.confidence)  # Bidirecional
            elif h.donor_atom == 'HB' and h.acceptor_atom == 'HN':
                edge.hb_hn_score = max(edge.hb_hn_score, h.confidence)
            elif h.donor_atom == 'HN' and h.acceptor_atom == 'HB':
                edge.hb_hn_score = max(edge.hb_hn_score, h.confidence)  # Bidirecional
        
        print(f"\n  Unique edges: {len(edges)}")
        
        # Mostra edges com HA-HN
        ha_hn_edges = [e for e in edges.values() if e.ha_hn_score > 0]
        print(f"  Edges with HA-HN: {len(ha_hn_edges)}")
        
        return list(edges.values())