import pytest
import mdtraj as md

@pytest.fixture(scope="session")
def trajectory():
    return md.load("trypsin_benzamidine.pdb")

@pytest.fixture(scope="session")
def non_water_atoms_indices(
        trajectory,
        atoms,
):
    return [
        atom.index for atom in atoms
        if "HOH" not in atom.residue.name
    ]

@pytest.fixture(scope="session")
def water_atoms_indices(
        trajectory,
        atoms,
):
    return trajectory.topology.select("water")

@pytest.fixture(scope="session")
def atoms(
        trajectory,
):
    return [atom for atom in trajectory.topology.atoms]

@pytest.fixture(scope="session")
def bonds(
        trajectory,
):
    return [bond for bond in trajectory.topology.bonds]

@pytest.fixture(scope="session")
def repartition_groups(
        bonds,
        atoms,
):
    non_water_bonds_bonds = [
        bond for bond in bonds if "HOH" not in bond[0].residue.name
    ]
    non_water_hydrogen_bonds = [
        bond for bond in non_water_bonds_bonds
        if bond[0].element.symbol == "H" or bond[1].element.symbol
           == "H"
    ]
    atom_masses = [atom.element.mass for atom in atoms]
    repartition_groups_dict = {}
    heavy_atom_indices = set([])
    for bond in non_water_hydrogen_bonds:
        if bond[0].element.symbol == "H":
            if bond[1].index not in heavy_atom_indices:
                heavy_atom_indices.add(bond[1].index)
                repartition_groups_dict[(bond[1].index, atom_masses[
                    bond[1].index])] = []
            repartition_groups_dict[(bond[1].index, atom_masses[
                bond[1].index])].append(
                (bond[0].index, atom_masses[bond[0].index])
            )
        else:
            if bond[0].index not in heavy_atom_indices:
                heavy_atom_indices.add(bond[0].index)
                repartition_groups_dict[(bond[0].index, atom_masses[
                    bond[0].index])] = []
            repartition_groups_dict[(bond[0].index, atom_masses[
                bond[0].index])].append(
                (bond[1].index, atom_masses[bond[1].index]))
    return repartition_groups_dict

