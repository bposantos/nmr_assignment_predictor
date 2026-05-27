# assignment/sequential_resolver.py

from collections import defaultdict

class SequentialResolver:
    """Resolve a ordem sequencial dos spin systems usando NOEs"""
    
    def __init__(self, sequence_composition, edges):
        self.sequence_composition = sequence_composition
        self.edges = edges
    
    def find_optimal_path(self, n_systems):
        """
        Encontra o caminho ótimo através do grafo de NOEs.
        """
        # Constrói grafo (considerando conexões com HN-HN)
        adj = defaultdict(list)
        for edge in self.edges:
            # Inclui qualquer conexão com HN-HN ou HA-HN significativo
            if edge.hn_hn_score > 0.3 or edge.ha_hn_score > 0.5:
                adj[edge.source_system_id].append((edge.target_system_id, edge.total_score))
        
        # Encontra todos os nós com conexões
        connected = set(adj.keys()) | set(sum(adj.values(), []))
        
        if not connected:
            return []
        
        # Busca o caminho mais longo
        best_path = []
        
        for start in connected:
            path = []
            used = set()
            current = start
            
            while current is not None:
                path.append(current)
                used.add(current)
                
                # Próximo não usado
                candidates = [(n, s) for n, s in adj.get(current, []) if n not in used]
                if not candidates:
                    break
                current = max(candidates, key=lambda x: x[1])[0]
            
            if len(path) > len(best_path):
                best_path = path
        
        return best_path
    
    def _find_start_node(self, adj, n_systems):
        """Encontra o nó mais provável como N-terminal"""
        # N-terminal geralmente tem menos conexões de entrada
        in_degree = defaultdict(int)
        for src, targets in adj.items():
            for tgt, _ in targets:
                in_degree[tgt] += 1
        
        # Sistema com menor in-degree é provável N-terminal
        min_in = min(in_degree.values()) if in_degree else 0
        for i in range(n_systems):
            if in_degree.get(i, 0) == min_in:
                return i
        return 0