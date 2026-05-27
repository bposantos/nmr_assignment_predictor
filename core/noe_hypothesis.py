# core/noe_hypothesis.py

from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class NoeHypothesis:
    """Hipótese de NOE entre dois sistemas de spin"""
    donor_system_id: int
    donor_atom: str  # 'HN', 'HA', 'HB', etc.
    acceptor_system_id: int
    acceptor_atom: str
    intensity: float
    peak_shift1: float
    peak_shift2: float
    confidence: float = 0.0
    
    def __post_init__(self):
        # Calcula confiança baseada na intensidade e tipo de átomo
        self.confidence = self._calculate_confidence()
    
    def _calculate_confidence(self) -> float:
        # Intensidade normalizada
        intensity_score = min(1.0, self.intensity / 1e10)
        
        # Bônus para pares específicos
        bonus = 0.0
        if self.donor_atom == 'HN' and self.acceptor_atom == 'HN':
            bonus = 0.4
        elif self.donor_atom == 'HA' and self.acceptor_atom == 'HN':
            bonus = 0.3
        elif self.donor_atom == 'HB' and self.acceptor_atom == 'HN':
            bonus = 0.1
        
        return min(1.0, intensity_score * 0.5 + bonus)


@dataclass
class SequentialEdge:
    """Aresta sequencial entre dois sistemas de spin"""
    source_system_id: int
    target_system_id: int
    
    hn_hn_score: float = 0.0
    ha_hn_score: float = 0.0
    hb_hn_score: float = 0.0
    
    supporting_peaks: List[NoeHypothesis] = field(default_factory=list)
    
    @property
    def total_score(self) -> float:
        """Score total baseado em todas as evidências"""
        return (self.hn_hn_score * 0.3 + 
                self.ha_hn_score * 0.5 + 
                self.hb_hn_score * 0.2)
    
    @property
    def is_strong(self) -> bool:
        """Verifica se a conexão é forte"""
        return self.total_score > 0.35