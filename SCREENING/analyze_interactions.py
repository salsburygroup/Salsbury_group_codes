from plip.structure.preparation import PDBComplex

def analyze_interactions(pdb_file):
    complex_structure = PDBComplex()
    complex_structure.load_pdb(pdb_file)
    complex_structure.analyze()

    interaction_sets = complex_structure.interaction_sets

    for ligand_key, interaction_obj in interaction_sets.items():
        ligand_id = str(ligand_key)
        print(f"\n--- Interactions for ligand {ligand_id} ---")

        # Hydrogen bonds: protein is donor
        if hasattr(interaction_obj, "hbonds_pdon") and interaction_obj.hbonds_pdon:
            print("\nHydrogen Bonds (Protein Donor → Ligand Acceptor):")
            for hbond in interaction_obj.hbonds_pdon:
                print(f"  {hbond.restype}{hbond.resnr} chain {hbond.reschain} → Ligand atom {hbond.atype} (distance: {hbond.distance_ad:.2f} Å)")

        # Hydrogen bonds: ligand is donor
        if hasattr(interaction_obj, "hbonds_ldon") and interaction_obj.hbonds_ldon:
            print("\nHydrogen Bonds (Ligand Donor → Protein Acceptor):")
            for hbond in interaction_obj.hbonds_ldon:
                print(f"  Ligand atom {hbond.atype} → {hbond.restype}{hbond.resnr} chain {hbond.reschain} (distance: {hbond.distance_ad:.2f} Å)")

        # Hydrophobic contacts
        if hasattr(interaction_obj, "hydrophobic_contacts") and interaction_obj.hydrophobic_contacts:
            print("\nHydrophobic Contacts:")
            for contact in interaction_obj.hydrophobic_contacts:
                print(f"  {contact.restype}{contact.resnr} chain {contact.reschain} ↔ Ligand atom (distance: {contact.distance:.2f} Å)")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python analyze_interactions.py docked_structure.pdb")
        sys.exit(1)

    pdb_path = sys.argv[1]
    analyze_interactions(pdb_path)

