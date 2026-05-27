# assignment/sequential_validator.py

from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple


@dataclass
class SequentialEvidence:
    """Evidência sequencial entre dois resíduos"""
    i: int
    j: int
    res_i: str
    res_j: str
    hn_hn_score: float = 0.0
    ha_hn_score: float = 0.0
    hb_hn_score: float = 0.0
    
    @property
    def total_score(self) -> float:
        return (self.hn_hn_score * 0.4 + 
                self.ha_hn_score * 0.5 + 
                self.hb_hn_score * 0.1)
    
    @property
    def is_sequential(self) -> bool:
        return abs(self.i - self.j) == 1


class SequentialValidator:
    """Valida assignments sequenciais contra a sequência real"""

    MAX_LONG_RANGE_DISTANCE = 4  # Só mostra conexões até i,i+4
    
    def __init__(self, sequence: List[str], spin_systems, assignments, edges):
        self.sequence = sequence
        self.spin_systems = spin_systems
        self.assignments = assignments
        self.edges = edges
        
        # Mapeia sistema para posição na sequência
        self.system_to_position = self._assign_positions_to_systems()
        
        # Constrói matriz de conexões
        self.connection_matrix = self._build_connection_matrix()
       
    def validate_sequential_assignments(self) -> Dict:
        """Valida assignments sequenciais usando posições mapeadas"""
        
        results = {
            'correct_sequential': [],
            'missing_sequential': [],
            'long_range': [],
            'statistics': {}
        }
        
        # Mapeia edges por par de sistemas
        edge_map = {}
        for edge in self.edges:
            key = (edge.source_system_id, edge.target_system_id)
            edge_map[key] = edge
        
        # Constrói matriz de conexões entre posições
        pos_connections = defaultdict(lambda: defaultdict(float))
        
        for (sys_i, sys_j), edge in edge_map.items():
            if sys_i in self.system_to_position and sys_j in self.system_to_position:
                pos_i = self.system_to_position[sys_i]
                pos_j = self.system_to_position[sys_j]
                
                # Acumula o melhor score para cada par de posições
                current = pos_connections[pos_i][pos_j]
                pos_connections[pos_i][pos_j] = max(current, edge.total_score)
        
        # Identifica conexões sequenciais (i, i+1)
        found_sequential = set()
        for i in range(len(self.sequence) - 1):
            if pos_connections[i][i+1] > 0:
                evidence = SequentialEvidence(
                    i=i, j=i+1,
                    res_i=self.sequence[i], res_j=self.sequence[i+1],
                    hn_hn_confidence=0.0,  # Pode detalhar mais
                    ha_hn_confidence=0.0,
                    hb_hn_confidence=0.0
                )
                # Aqui você pode adicionar os scores específicos
                evidence.hn_hn_confidence = pos_connections[i][i+1]
                evidence.ha_hn_confidence = pos_connections[i][i+1]
                
                results['correct_sequential'].append(evidence)
                found_sequential.add((i, i+1))
            else:
                results['missing_sequential'].append((i, i+1, self.sequence[i], self.sequence[i+1]))
        
        # Identifica conexões de longo alcance
        for i in range(len(self.sequence)):
            for j in range(i+2, len(self.sequence)):
                if pos_connections[i][j] > 0 or pos_connections[j][i] > 0:
                    score = max(pos_connections[i][j], pos_connections[j][i])
                    results['long_range'].append((i, j, self.sequence[i], self.sequence[j], score))
        
        # Estatísticas
        total_expected = len(self.sequence) - 1
        found = len(results['correct_sequential'])
        
        results['statistics'] = {
            'total_expected_sequential': total_expected,
            'found_sequential': found,
            'completeness': found / total_expected if total_expected > 0 else 0,
            'total_edges': len(self.edges),
            'long_range_edges': len(results['long_range'])
            # 'correct_edges' removido porque não faz sentido com o novo modelo
        }
        
        return results
    
    def _map_residues_to_systems(self) -> Dict[str, List[int]]:
        """Mapeia cada resíduo da sequência para os sistemas possíveis"""
        mapping = defaultdict(list)
        for i, (ss, ass) in enumerate(zip(self.spin_systems, self.assignments)):
            if ass and ass.residue_name != 'UNKNOWN':
                mapping[ass.residue_name].append(i)
        return mapping
    
    def build_expected_connections(self) -> List[Tuple[int, int]]:
        """Retorna lista de pares esperados (i, i+1) na sequência"""
        expected = []
        for i in range(len(self.sequence) - 1):
            expected.append((i, i+1))
        return expected
    
    def find_possible_system_pairs(self) -> List[Tuple[int, int, str, str]]:
        """Encontra todos os pares possíveis de sistemas baseado nos resíduos"""
        pairs = []
        for i, res_i in enumerate(self.sequence):
            for j, res_j in enumerate(self.sequence):
                if i == j:
                    continue
                
                # Busca sistemas que podem ser estes resíduos
                sys_i_list = self.residue_to_system.get(res_i, [])
                sys_j_list = self.residue_to_system.get(res_j, [])
                
                for sys_i in sys_i_list:
                    for sys_j in sys_j_list:
                        if sys_i != sys_j:
                            pairs.append((sys_i, sys_j, res_i, res_j, i, j))
        
        return pairs
    
    def _assign_positions_to_systems(self) -> Dict[int, int]:
        """Atribui cada sistema a uma posição única na sequência"""
        systems_with_score = []
        for i, (ss, ass) in enumerate(zip(self.spin_systems, self.assignments)):
            if ass and ass.residue_name != 'UNKNOWN':
                systems_with_score.append((i, ass.residue_name, ass.score))
        
        systems_with_score.sort(key=lambda x: x[2], reverse=True)
        
        position_to_system = {}
        system_to_position = {}
        
        for sys_id, res_name, score in systems_with_score:
            positions = [i for i, r in enumerate(self.sequence) if r == res_name]
            available = [p for p in positions if p not in position_to_system]
            
            if available:
                pos = available[0]
                position_to_system[pos] = sys_id
                system_to_position[sys_id] = pos
        
        return system_to_position
    
    def _build_connection_matrix(self) -> Dict[Tuple[int, int], Dict]:
        """Constrói matriz de conexões entre posições com scores específicos"""
        
        # Inicializa matriz
        matrix = {}
        for i in range(len(self.sequence)):
            for j in range(len(self.sequence)):
                if i != j:
                    matrix[(i, j)] = {'hn_hn': 0.0, 'ha_hn': 0.0, 'hb_hn': 0.0}
        
        # Preenche com dados dos NOEs
        for edge in self.edges:
            src = edge.source_system_id
            tgt = edge.target_system_id
            
            if src in self.system_to_position and tgt in self.system_to_position:
                pos_i = self.system_to_position[src]
                pos_j = self.system_to_position[tgt]
                
                key = (pos_i, pos_j)
                matrix[key]['hn_hn'] = max(matrix[key]['hn_hn'], edge.hn_hn_score)
                matrix[key]['ha_hn'] = max(matrix[key]['ha_hn'], edge.ha_hn_score)
                matrix[key]['hb_hn'] = max(matrix[key]['hb_hn'], edge.hb_hn_score)
        
        return matrix
    
    def get_sequential_table(self) -> List[Dict]:
        """Gera tabela de conexões sequenciais i,i+1"""
        
        table = []
        for i in range(len(self.sequence) - 1):
            key_fwd = (i, i+1)
            key_rev = (i+1, i)
            
            # Pega o melhor score nas duas direções
            hn_hn = max(self.connection_matrix.get(key_fwd, {}).get('hn_hn', 0),
                       self.connection_matrix.get(key_rev, {}).get('hn_hn', 0))
            ha_hn = max(self.connection_matrix.get(key_fwd, {}).get('ha_hn', 0),
                       self.connection_matrix.get(key_rev, {}).get('ha_hn', 0))
            hb_hn = max(self.connection_matrix.get(key_fwd, {}).get('hb_hn', 0),
                       self.connection_matrix.get(key_rev, {}).get('hb_hn', 0))
            
            has_connection = (hn_hn > 0 or ha_hn > 0 or hb_hn > 0)
            
            table.append({
                'pos': i+1,
                'res_i': self.sequence[i],
                'res_j': self.sequence[i+1],
                'hn_hn': hn_hn,
                'ha_hn': ha_hn,
                'hb_hn': hb_hn,
                'has_connection': has_connection
            })
        
        return table
    
    def print_sequential_table(self):
        """Imprime tabela de conexões sequenciais formatada"""
        
        table = self.get_sequential_table()
        
        print("\n" + "="*80)
        print("SEQUENTIAL CONNECTIONS TABLE (i → i+1)")
        print("="*80)
        print(f"\n{'Pos':<6} {'Residues':<14} {'HN-HN':<10} {'HA-HN':<10} {'HB-HN':<10} {'Status':<12}")
        print(f"{'-'*6} {'-'*14} {'-'*10} {'-'*10} {'-'*10} {'-'*12}")
        
        for row in table:
            status = "✓ FOUND" if row['has_connection'] else "✗ MISSING"
            hn_str = f"{row['hn_hn']:.2f}" if row['hn_hn'] > 0 else "—"
            ha_str = f"{row['ha_hn']:.2f}" if row['ha_hn'] > 0 else "—"
            hb_str = f"{row['hb_hn']:.2f}" if row['hb_hn'] > 0 else "—"
            
            print(f"{row['pos']}→{row['pos']+1:<3} {row['res_i']}→{row['res_j']:<8} "
                  f"{hn_str:<10} {ha_str:<10} {hb_str:<10} {status}")
        
        # Resumo
        found = sum(1 for row in table if row['has_connection'])
        total = len(table)
        print(f"\n📊 Summary: {found}/{total} sequential connections found ({found/total*100:.1f}%)")
    
    def print_long_range_connections(self, max_distance: int = 4):
        """Imprime conexões de longo alcance (limitadas a max_distance)"""
        
        print(f"\n" + "="*80)
        print(f"LONG-RANGE CONNECTIONS (|i-j| ≤ {max_distance}, i≠i+1)")
        print("="*80)
        
        connections = []
        for i in range(len(self.sequence)):
            for j in range(i+2, min(len(self.sequence), i + max_distance + 1)):
                key = (i, j)
                hn_hn = self.connection_matrix.get(key, {}).get('hn_hn', 0)
                ha_hn = self.connection_matrix.get(key, {}).get('ha_hn', 0)
                hb_hn = self.connection_matrix.get(key, {}).get('hb_hn', 0)
                
                if hn_hn > 0 or ha_hn > 0 or hb_hn > 0:
                    connections.append({
                        'i': i, 'j': j,
                        'res_i': self.sequence[i],
                        'res_j': self.sequence[j],
                        'hn_hn': hn_hn,
                        'ha_hn': ha_hn,
                        'hb_hn': hb_hn,
                        'distance': j - i
                    })
        
        if connections:
            print(f"\n{'Pos':<8} {'Residues':<14} {'Dist':<6} {'HN-HN':<10} {'HA-HN':<10} {'HB-HN':<10}")
            print(f"{'-'*8} {'-'*14} {'-'*6} {'-'*10} {'-'*10} {'-'*10}")
            
            for conn in connections:
                hn_str = f"{conn['hn_hn']:.2f}" if conn['hn_hn'] > 0 else "—"
                ha_str = f"{conn['ha_hn']:.2f}" if conn['ha_hn'] > 0 else "—"
                hb_str = f"{conn['hb_hn']:.2f}" if conn['hb_hn'] > 0 else "—"
                
                print(f"{conn['i']+1}→{conn['j']+1:<3} {conn['res_i']}→{conn['res_j']:<8} "
                      f"{conn['distance']:<6} {hn_str:<10} {ha_str:<10} {hb_str:<10}")
        else:
            print("\n   No long-range connections found within distance limit")
    
    def print_validation_report(self):
        """Imprime relatório completo de validação"""
        
        print("\n" + "="*80)
        print("SEQUENTIAL ASSIGNMENT VALIDATION REPORT")
        print("="*80)
        
        # Tabela sequencial
        self.print_sequential_table()
        
        # Conexões de longo alcance
        self.print_long_range_connections(max_distance=self.MAX_LONG_RANGE_DISTANCE)
        
        # Resumo final
        table = self.get_sequential_table()
        found = sum(1 for row in table if row['has_connection'])
        total = len(table)
        
        print(f"\n{'='*80}")
        print(f"FINAL SUMMARY")
        print(f"{'='*80}")
        print(f"   Sequential connections found: {found}/{total} ({found/total*100:.1f}%)")
        
        if found/total < 0.5:
            print(f"   ⚠️  Low completeness - consider adjusting NOE tolerances")
        elif found/total < 0.8:
            print(f"   📈 Good completeness - fine-tuning could improve")
        else:
            print(f"   🎉 Excellent completeness! Sequential assignment is reliable")
        
        # Lista de conexões faltantes
        missing = [row for row in table if not row['has_connection']]
        if missing:
            print(f"\n   Missing connections:")
            for row in missing:
                print(f"      Position {row['pos']}→{row['pos']+1}: {row['res_i']}→{row['res_j']}")