from database.chemical_shift_database import ChemicalShiftDatabase
from collections import defaultdict
import numpy as np
import os
from pathlib import Path
import sys

from graph.graph_builder import GraphBuilder
from assignment.assignment_engine_competitive import AssignmentEngineCompetitive

from parsers.tocsy_parser import parse_tocsy
from parsers.noesy_parser import parse_noesy
from parsers.hsqc_parser import parse_hsqc

from visualization.nmr_visualizer import NMRVisualizer

# NOVOS IMPORTS PARA NOESY SEQUENCIAL
from scoring.noe_sequential_scorer import NoeSequentialScorer
from assignment.sequential_resolver import SequentialResolver
from core.noe_hypothesis import NoeHypothesis, SequentialEdge
from assignment.sequential_validator import SequentialValidator

# Adiciona módulos ao path
sys.path.append(str(Path(__file__).parent))

from scoring.aromatic_handler import AromaticHandler
from graph.candidate_generator_corrected import CandidateGeneratorCorrected

# Optional ML imports
USE_ML = False

if USE_ML:
    from ml.feature_extractor import FeatureExtractor
    from ml.gnn_model import GNNModel
    from ml.post_processor import PostProcessor
    import torch


# =========================================================
# Input Files
# =========================================================

TOCSY_FILE = "/home/overwatch/Documents/2026_assignment/gnn_test_data/peptide_05/tocsy.list"
NOESY_FILE = "/home/overwatch/Documents/2026_assignment/gnn_test_data/peptide_05/noesy.list"
HSQC_FILE = "/home/overwatch/Documents/2026_assignment/gnn_test_data/peptide_05/hsqc.list"
SEQUENCE_FILE = "/home/overwatch/Documents/2026_assignment/gnn_test_data/peptide_05/sequence.txt"

# =========================================================
# Load Experimental Spectra
# =========================================================

print("Loading experimental spectra...")

tocsy_peaks = parse_tocsy(TOCSY_FILE)
noesy_peaks = parse_noesy(NOESY_FILE)
hsqc_peaks = parse_hsqc(HSQC_FILE)

print(f"Loaded {len(tocsy_peaks)} TOCSY peaks")
print(f"Loaded {len(noesy_peaks)} NOESY peaks")
print(f"Loaded {len(hsqc_peaks)} HSQC peaks")

# =========================================================
# Initialize Chemical Shift Database
# =========================================================

print("\nLoading chemical shift database...")
chemical_database = ChemicalShiftDatabase()

# =========================================================
# Build Experimental Graphs
# =========================================================

print("\nBuilding experimental graphs...")

builder = GraphBuilder(
    proton_tolerance=0.01,
    intensity_threshold=300
)

# Adiciona os peaks ao builder
builder.tocsy_peaks = tocsy_peaks
builder.noesy_peaks = noesy_peaks
builder.hsqc_peaks = hsqc_peaks  # <--- IMPORTANTE!

tocsy_graph = builder.build_tocsy_graph(tocsy_peaks)
noesy_graph = builder.build_noesy_graph(noesy_peaks)
hsqc_graph = builder.build_hsqc_graph(hsqc_peaks)

print(f"\nTOCSY graph: {len(tocsy_graph.nodes)} nodes / {len(tocsy_graph.edges)} edges")
print(f"NOESY graph: {len(noesy_graph.nodes)} nodes / {len(noesy_graph.edges)} edges")
print(f"HSQC graph: {len(hsqc_graph.nodes)} nodes / {len(hsqc_graph.edges)} edges")

# =========================================================
# Detect and Annotate Spin Systems
# =========================================================

print("\nDetecting spin systems...")
spin_systems = builder.extract_spin_systems(tocsy_graph, hsqc_graph, min_spins=2)
print(f"Detected {len(spin_systems)} spin systems")

print("\nAnnotating spin systems with HSQC data...")
for ss in spin_systems:
    builder.annotate_spin_system_with_hsqc(ss, tolerance=0.03)

print("\nBuilding spin system structures with H-C data...")
spin_structures = []

for ss in spin_systems:
    structure = builder.build_spin_system_structure(ss)
    if structure.hn_shift and structure.spins:
        spin_structures.append(structure)
        print(f"  System {len(spin_structures)}: HN={structure.hn_shift:.3f}, "
              f"{structure.n_spins} spins, annotated={structure.n_carbon_annotated}, "
              f"fingerprint={structure.fingerprint}")
    

# Limita ao número esperado de resíduos
try:
    with open(SEQUENCE_FILE, 'r') as f:
        expected_residue_count = sum(1 for line in f if line.strip() and line.strip() not in ['NH3', 'ACE'])
except:
    expected_residue_count = 13

if len(spin_structures) > expected_residue_count:
    print(f"Limiting from {len(spin_structures)} to {expected_residue_count} spin systems (by size)")
    spin_structures.sort(key=lambda x: x.n_spins, reverse=True)
    spin_structures = spin_structures[:expected_residue_count]


# =========================================================
# 5. Carrega composição da sequência (depois do database)
# =========================================================
engine = AssignmentEngineCompetitive(chemical_database, sequence_file=SEQUENCE_FILE)
sequence_composition = engine.expected_composition
print(f"Sequence composition for candidate filtering: {sequence_composition}")


# =========================================================
# Initialize Competitive Assignment Engine
# =========================================================

candidate_generator = CandidateGeneratorCorrected(
    database=chemical_database,  # Use chemical_database, não chemical_database (já está correto)
    sequence_composition=sequence_composition,
    max_candidates=5,
    min_score_threshold=0.05
)

engine = AssignmentEngineCompetitive(chemical_database, sequence_file=SEQUENCE_FILE)

# =========================================================
# Process All Spin Systems Competitively
# =========================================================

print("\n" + "="*60)
print("STARTING COMPETITIVE ASSIGNMENT")
print("="*60)

# Executa atribuição competitiva
results = engine.assign_all_systems(spin_structures, candidate_generator, max_iterations=3)

#############
# NOESY
#############

print("\n" + "="*60)
print("SEQUENTIAL ASSIGNMENT WITH NOESY")
print("="*60)

# Cria scorer NOESY
noe_scorer = NoeSequentialScorer(noesy_peaks)

# Constrói grafo sequencial
edges = noe_scorer.build_sequential_graph(spin_systems, results)

print(f"Found {len(edges)} sequential connections")

# Mostra conexões fortes
for edge in edges:
    if edge.is_strong:
        print(f"  S{edge.source_system_id+1} -> S{edge.target_system_id+1}: "
              f"score={edge.total_score:.2f} (HN-HN:{edge.hn_hn_score:.2f}, "
              f"HA-HN:{edge.ha_hn_score:.2f})")

# Resolve ordem sequencial
resolver = SequentialResolver(sequence_composition, edges)
sequential_order = resolver.find_optimal_path(len(spin_systems))

print(f"\nProposed sequential order: {[f'S{i+1}' for i in sequential_order]}")


print("\n" + "="*60)
print("SEQUENTIAL ASSIGNMENT VALIDATION")
print("="*60)

# Carrega a sequência real
sequence = []
with open(SEQUENCE_FILE, 'r') as f:
    for line in f:
        res = line.strip().upper()
        if res and res not in ['NH3', 'ACE']:
            sequence.append(res)

print(f"Real sequence: {sequence}")

# Cria validador
validator = SequentialValidator(sequence, spin_structures, results, edges)

# Gera relatório completo
validator.print_validation_report()

# =========================================================
# Display Results
# =========================================================

print("\n" + "="*60)
print("ASSIGNMENT RESULTS")
print("="*60)

for idx, (structure, result) in enumerate(zip(spin_structures, results)):
    engine.print_system_details(structure, result, idx + 1)

# Atribui os resultados aos spin_systems originais para visualização
for i, (structure, result) in enumerate(zip(spin_structures, results)):
    if i < len(spin_systems):
        spin_systems[i].assigned_residue = result.residue_name
        spin_systems[i].assignment_score = result.score

# =========================================================
# ANCHOR AROMATIC SPINS TO PHE SYSTEMS
# =========================================================

print("\n" + "="*60)
print("ANCHORING AROMATIC SPINS")
print("="*60)

from scoring.aromatic_anchor import AromaticAnchor

# Cria ancora
aromatic_anchor = AromaticAnchor(hsqc_peaks, noesy_peaks, spin_systems, results)

# Ancora sinais aromáticos aos sistemas PHE
anchored = aromatic_anchor.anchor_all_phe_systems()

# Adiciona os spins aromáticos aos sistemas
if anchored:
    print("\n  Adding aromatic spins to spin systems...")
    for phe_id, aromatic_spins in anchored.items():
        try:
            aromatic_anchor.add_aromatic_to_spin_system(phe_id, aromatic_spins)
        except Exception as e:
            print(f"    Warning: Could not add spins to S{phe_id+1}: {e}")
            # Mostra o que seria adicionado
            for a in aromatic_spins:
                proton_name = a.get('proton_name', a.get('atom_type', 'H?'))
                carbon_name = a.get('carbon_name', 'C?')
                print(f"      Would add: {proton_name} H={a['proton']:.3f}, C={a['carbon']:.1f}")

# Rebuild spin_structures
if anchored:
    print("\n  Rebuilding spin structures with aromatic spins...")
    spin_structures = []
    for ss in spin_systems:
        structure = builder.build_spin_system_structure(ss)
        if structure.hn_shift and structure.spins:
            spin_structures.append(structure)

# =========================================================
# Análise de Aromáticos Não Atribuídos
# =========================================================

print("\n" + "="*60)
print("AROMATIC ANALYSIS")
print("="*60)

unassigned_aromatics = []
for i, (structure, result) in enumerate(zip(spin_structures, results)):
    if result.residue_name == "UNKNOWN":
        is_aromatic, candidate = AromaticHandler.detect_aromatic_candidate(structure)
        if is_aromatic:
            unassigned_aromatics.append((structure, candidate))
            print(f"⚠️  System {i+1} (HN={structure.hn_shift:.3f}) appears to be {candidate} (not assigned)")

if unassigned_aromatics:
    print(f"\nFound {len(unassigned_aromatics)} potential aromatic systems not assigned")
    print("\n💡 These may be PHE/TYR residues where only HA/HB were detected")
    
    # Sugestão para reatribuição manual
    print("\nSuggested manual assignments:")
    for structure, candidate in unassigned_aromatics[:5]:
        print(f"  System HN={structure.hn_shift:.3f} -> {candidate}")

# =========================================================
# Análise de Sistemas Parciais
# =========================================================

print("\n" + "="*60)
print("PARTIAL SYSTEMS ANALYSIS")
print("="*60)

partial_systems = []
for i, (structure, result) in enumerate(zip(spin_structures, results)):
    if result.residue_name == "UNKNOWN":
        # Verifica padrão LEU (HG ~41 ppm)
        for spin in structure.spins:
            if hasattr(spin, 'c_shift') and spin.c_shift and 38 <= spin.c_shift <= 44:
                if spin.c_shift >= 40:  # HG de LEU
                    partial_systems.append((structure, "LEU"))
                    print(f"⚠️  System {i+1} (HN={structure.hn_shift:.3f}) - possible LEU (HG={spin.c_shift:.1f} ppm)")
                    break
                else:
                    partial_systems.append((structure, "ILE/VAL"))
                    print(f"⚠️  System {i+1} (HN={structure.hn_shift:.3f}) - possible ILE/VAL (CB={spin.c_shift:.1f} ppm)")
                    break

if partial_systems:
    print(f"\nFound {len(partial_systems)} partial systems with diagnostic patterns")        

# =========================================================
# Final Summary
# =========================================================

print("\n" + "="*50)
print("FINAL ASSIGNMENT SUMMARY")
print("="*50)

final_counts = defaultdict(int)
for result in results:
    final_counts[result.residue_name] += 1

# CORREÇÃO: Verificar se engine.expected_composition existe
if hasattr(engine, 'expected_composition'):
    expected_composition = engine.expected_composition
else:
    expected_composition = sequence_composition
    print("Note: Using sequence_composition from engine")

for residue, expected in sorted(expected_composition.items()):
    assigned = final_counts.get(residue, 0)
    status = "✓" if assigned == expected else f"({assigned}/{expected})"
    bar = "█" * assigned + "░" * (expected - assigned)
    print(f"{residue:4} : {bar:10} {status}")

total_assigned = sum(final_counts.values())
total_expected = sum(expected_composition.values())
print(f"\nTotal assigned: {total_assigned}/{total_expected}")

unknown_count = sum(1 for r in results if r.residue_name == "UNKNOWN")
if unknown_count > 0:
    print(f"\n⚠️  {unknown_count} spin systems could not be assigned")

# =========================================================
# Visualizations
# =========================================================

print("\nGenerating visualizations...")

visualizer = NMRVisualizer()


# Visualização detalhada de cada sistema (opcional, pode comentar se for muitos)
#for idx, ss in enumerate(spin_systems[:5]):  # Mostra só os 5 primeiros
#    if hasattr(ss, 'assigned_residue') and ss.assigned_residue != "UNKNOWN":
#        visualizer.plot_single_spin_system(ss, tocsy_peaks, hsqc_peaks, idx)

# Summary figure com todos os sistemas coloridos
visualizer.plot_spin_systems_summary(
    spin_systems[:len(spin_structures)],
    tocsy_peaks,
    hsqc_peaks,
    save_path="spin_system_summary.png"
)

print("\n✅ Assignment completed successfully!")