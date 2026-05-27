# scoring/aromatic_anchor.py

from collections import defaultdict
from typing import List, Dict, Tuple, Optional


class AromaticAnchor:
    """Ancora sinais aromáticos (HD/HE/HZ) aos sistemas PHE usando NOEs"""
    
    # Faixas típicas para prótons aromáticos de PHE
    AROMATIC_PROTON_RANGES = {
        'HD': (6.8, 7.6),   # delta (orto) - geralmente mais desblindado
        'HE': (7.0, 7.5),   # epsilon (meta)
        'HZ': (7.1, 7.4),   # zeta (para)
    }
    
    # Faixas para carbonos aromáticos
    AROMATIC_CARBON_RANGES = {
        'HD': (131, 132),
        'HE': (130.0, 130.9),
        'HZ': (129.0, 129.9),
    }
    
    def __init__(self, hsqc_peaks, noesy_peaks, spin_systems, assignments):
        self.hsqc_peaks = hsqc_peaks
        self.noesy_peaks = noesy_peaks
        self.spin_systems = spin_systems
        self.assignments = assignments
    
    def extract_aromatic_peaks(self) -> List[Dict]:
        """Extrai todos os picos aromáticos do HSQC (C > 100 ppm)"""
        
        aromatic_peaks = []
        
        for peak in self.hsqc_peaks:
            h_shift = float(peak.w1)
            c_shift = float(peak.w2)
            
            # Verifica se é aromático
            if c_shift > 100:
                aromatic_peaks.append({
                    'proton': h_shift,
                    'carbon': c_shift,
                    'intensity': getattr(peak, 'intensity', 1.0),
                    'used': False
                })
        
        # Ordena por intensidade (mais intensos primeiro)
        aromatic_peaks.sort(key=lambda x: x['intensity'], reverse=True)
        
        print(f"\n  Found {len(aromatic_peaks)} aromatic peaks in HSQC")
        for p in aromatic_peaks[:10]:
            print(f"    H={p['proton']:.3f}, C={p['carbon']:.1f}")
        
        return aromatic_peaks
    
    def find_noes_to_phe(self, phe_system_id: int, phe_ha: float, phe_hb: float) -> List[Dict]:
        """Encontra NOEs entre um sistema PHE e picos aromáticos"""
        
        connections = []
        
        for peak in self.noesy_peaks:
            shift1 = float(peak.w1)
            shift2 = float(peak.w2)
            intensity = float(getattr(peak, 'intensity', 1.0))
            
            # Verifica se o pico envolve HA ou HB do PHE
            is_phe_ha = phe_ha and (abs(shift1 - phe_ha) < 0.03 or abs(shift2 - phe_ha) < 0.03)
            is_phe_hb = phe_hb and (abs(shift1 - phe_hb) < 0.03 or abs(shift2 - phe_hb) < 0.03)
            
            if not (is_phe_ha or is_phe_hb):
                continue
            
            # Encontra o shift aromático (o outro)
            aromatic_shift = None
            donor_atom = None
            
            if is_phe_ha:
                aromatic_shift = shift2 if abs(shift1 - phe_ha) < 0.03 else shift1
                donor_atom = 'HA'
            elif is_phe_hb:
                aromatic_shift = shift2 if abs(shift1 - phe_hb) < 0.03 else shift1
                donor_atom = 'HB'
            
            connections.append({
                'phe_system': phe_system_id,
                'phe_atom': donor_atom,
                'aromatic_proton': aromatic_shift,
                'intensity': intensity,
                'peak': (shift1, shift2)
            })
        
        return connections

    def assign_aromatic_to_phe(self, phe_system_id: int, phe_ha: float, phe_hb: float,
                                    aromatic_peaks: List[Dict]) -> List[Dict]:
        """Atribui picos aromáticos a um sistema PHE baseado em NOEs"""
        
        TOLERANCE = 0.05
        
        noes = self.find_noes_to_phe(phe_system_id, phe_ha, phe_hb)
        
        if not noes:
            return []
        
        # Extrai shifts aromáticos dos NOEs (na faixa 6.5-8.0)
        aromatic_noe_shifts = []
        for noe in noes:
            shift = noe['aromatic_proton']
            if 6.5 <= shift <= 8.0:
                aromatic_noe_shifts.append(shift)
        
        print(f"      Aromatic NOE shifts: {[f'{s:.3f}' for s in aromatic_noe_shifts]}")
        
        assigned = []
        for peak in aromatic_peaks:
            if peak['used']:
                continue
            
            h_shift = peak['proton']
            
            # Procura o NOE mais próximo
            best_match = None
            best_diff = TOLERANCE
            
            for noe_shift in aromatic_noe_shifts:
                diff = abs(h_shift - noe_shift)
                if diff < best_diff:
                    best_diff = diff
                    best_match = noe_shift
            
            if best_match is not None:
                c_shift = peak['carbon']
                proton_name, carbon_name = self._classify_aromatic_atom(c_shift, h_shift)
                
                assigned.append({
                    'proton': h_shift,
                    'proton_name': proton_name,      # ← USAR proton_name
                    'carbon': c_shift,
                    'carbon_name': carbon_name,
                    'atom_type': proton_name,        # ← ADICIONAR para compatibilidade
                    'phe_atom': 'HA/HB',
                    'intensity': peak['intensity'],
                    'noe_intensity': 1.0
                })
                peak['used'] = True
                print(f"      Matched: NOE {best_match:.3f} → HSQC H={h_shift:.3f}, C={c_shift:.1f} (diff={best_diff:.3f})")
        
        return assigned

    def _classify_aromatic_atom(self, carbon: float, proton: float) -> tuple:
        """Classifica o átomo aromático, retorna (proton_name, carbon_name)"""
        
        # Classifica carbono (CD, CE, CZ)
        carbon_name = 'C?'  # fallback
        for name, (c_min, c_max) in self.AROMATIC_CARBON_RANGES.items():
            if c_min <= carbon <= c_max:
                carbon_name = name
                break
        
        # Classifica próton (HD, HE, HZ)
        proton_name = 'H?'  # fallback
        for name, (h_min, h_max) in self.AROMATIC_PROTON_RANGES.items():
            if h_min <= proton <= h_max:
                proton_name = name
                break
        
        return proton_name, carbon_name
    
    def anchor_all_phe_systems(self) -> Dict:
        """Ancora sinais aromáticos a todos os sistemas PHE"""
        
        print("\n" + "="*80)
        print("AROMATIC ANCHORING TO PHE SYSTEMS")
        print("="*80)
        
        # Extrai picos aromáticos do HSQC
        aromatic_peaks = self.extract_aromatic_peaks()
        
        if not aromatic_peaks:
            print("\n  No aromatic peaks found in HSQC!")
            return {}
        
        # Identifica sistemas PHE
        phe_systems = []
        for i, (ss, ass) in enumerate(zip(self.spin_systems, self.assignments)):
            if ass and ass.residue_name == 'PHE':
                ha_shift = None
                hb_shift = None
                
                for group in ss.carbon_groups:
                    if group.carbon_type == 'methine' and 50 <= group.carbon_shift <= 65:
                        ha_shift = group.average_proton_shift
                    elif group.carbon_type == 'methylene' and 35 <= group.carbon_shift <= 45:
                        hb_shift = group.average_proton_shift
                
                phe_systems.append({
                    'id': i,
                    'ha': ha_shift,
                    'hb': hb_shift,
                    'spin_system': ss
                })
        
        print(f"\n  Found {len(phe_systems)} PHE systems")
        
        # Ancora cada PHE
        results = {}
        for phe in phe_systems:
            print(f"\n  Anchoring PHE S{phe['id']+1}: HA={phe['ha']:.3f}, HB={phe['hb']:.3f}")
            
            assigned = self.assign_aromatic_to_phe(
                phe['id'], phe['ha'], phe['hb'], aromatic_peaks
            )
            
            if assigned:
                results[phe['id']] = assigned
                print(f"    Assigned {len(assigned)} aromatic spins:")
                for a in assigned:
                    # Usa proton_name em vez de atom_type
                    atom_name = a.get('proton_name', a.get('atom_type', '?'))
                    print(f"      {atom_name}: H={a['proton']:.3f}, C={a['carbon']:.1f} "
                        f"(connected via {a.get('phe_atom', '?')})")
        
        # Mostra picos não utilizados
        unused = [p for p in aromatic_peaks if not p['used']]
        if unused:
            print(f"\n  Unused aromatic peaks: {len(unused)}")
            for p in unused[:5]:
                print(f"    H={p['proton']:.3f}, C={p['carbon']:.1f}")
        
        return results

    def add_aromatic_to_spin_system(self, phe_system_id: int, aromatic_spins: List[Dict]):
        """Adiciona spins aromáticos ao sistema PHE"""
        
        ss = self.spin_systems[phe_system_id]
        
        from core.spin import Spin
        
        # Verifica se o sistema tem o atributo correto
        if hasattr(ss, 'spin_nodes'):
            spin_container = ss.spin_nodes
        elif hasattr(ss, 'spins'):
            spin_container = ss.spins
        else:
            ss.spin_nodes = []
            spin_container = ss.spin_nodes
        
        for aro in aromatic_spins:
            # Cria um novo spin para o átomo aromático
            spin = Spin(
                proton_shift=aro['proton'],
                carbon_shift=aro['carbon'],
                carbon_type='aromatic',
                multiplicity=1,
                source_peaks=[]
            )
            spin_container.append(spin)
            
            # Usa as chaves corretas
            proton_name = aro.get('proton_name', aro.get('atom_type', 'H?'))
            carbon_name = aro.get('carbon_name', 'C?')
            
            print(f"    Added {proton_name} ({carbon_name}) to PHE S{phe_system_id+1}: "
                f"H={aro['proton']:.3f}, C={aro['carbon']:.1f}")
        
        # Reagrupa os carbonos do sistema
        if hasattr(ss, 'group_by_carbon'):
            ss.group_by_carbon()
        elif hasattr(ss, 'carbon_groups'):
            self._rebuild_carbon_groups(ss)
        
        # Atualiza fingerprint
        if hasattr(ss, '_update_fingerprint_from_groups'):
            ss._update_fingerprint_from_groups()

    def _rebuild_carbon_groups(self, ss):
        """Reconstrói carbon_groups a partir dos spins"""
        from collections import defaultdict
        import numpy as np
        
        carbon_map = defaultdict(list)
        
        # Coleta todos os spins com carbono
        spin_nodes = getattr(ss, 'spin_nodes', getattr(ss, 'spins', []))
        
        for spin in spin_nodes:
            if hasattr(spin, 'carbon_shift') and spin.carbon_shift:
                carbon_map[spin.carbon_shift].append(spin)
        
        # Cria CarbonGroups
        from core.spin import CarbonGroup
        ss.carbon_groups = []
        for c_shift, spins in carbon_map.items():
            proton_shifts = [s.proton_shift for s in spins if hasattr(s, 'proton_shift')]
            group = CarbonGroup(
                carbon_shift=c_shift,
                proton_shifts=proton_shifts,
                peaks=spins
            )
            ss.carbon_groups.append(group)