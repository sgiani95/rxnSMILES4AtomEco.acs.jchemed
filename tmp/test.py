
reactions_smiles_pathway = {
    # Cumene hydroperoxide decomposition
    "CC(C)(OO)c1ccccc1>acid>CC(C)=O": "Cumene",
    # Isopropanol dehydrogenation
    "CC(C)O>Cu>CC(C)=O": "Isopropanol",
    # Propene oxidation
    "{2}C=CC.O=O>Pd/Cu>{2}CC(C)=O": "Propene"
}

print(reactions_smiles_pathway)
