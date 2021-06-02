"""
check.py

Check SEEKR2 calculations for common problems both before and after
simulations.

check_pre_simulation_all : 
    These are checks that should be done before running any
    simulations (run stage). They are automatically run by prepare.py.
 - Check for bubbles in starting structures, or unusually low densities
   (high densities will probably fail on their own). Check that box
   boundaries aren't significantly larger than starting structure 
   boxes.
 - Check to ensure that the same (or similar) salt concentrations exist
   between the MD and BD stages.
 - Check that the system exists within the expected Voronoi cell and
   give detailed information of what is wrong if it is not in there.
 - Check that atom selections exist on the same molecule, and that
   atoms from water molecules are not included accidentally.
 - Check that BD milestone's atom selections are the same between
   MD and BD.

check_post_simulation_all : 
    These are checks that should be performed after any or all
    simulations are performed, and before the analysis stage.
    They are automatically run by analyze.py
 - Check that all simulation trajectories keep the system where
   expected: 
    - Elber umbrella stage stays close to milestone.
    - Elber reversal/forward stage ends close to the milestone it was
      recorded to have crossed.
    - MMVT trajectories stay within the Voronoi cell, and they end
      close to the boundary recorded by the crossing.
  - Re-check simulation density? ***
  - Check that BD simulations end on the correct milestone
  
This module maybe also be run from the command line. Example:

$ python check.py /path/to/model.xml
"""

import os
import argparse
import collections
import glob
import warnings
import bubblebuster
warnings.filterwarnings("ignore")

import numpy as np
import parmed
import mdtraj

import seekr2.modules.common_base as base
import seekr2.modules.elber_base as elber_base
import seekr2.modules.mmvt_base as mmvt_base
from seekr2.modules.common_base import Model 
# The charges of common ions
ION_CHARGE_DICT = {"li":1.0, "na":1.0, "k":1.0, "rb":1.0, "cs":1.0, "fr":1.0,
                   "be":2.0, "mg":2.0, "ca":2.0, "sr":2.0, "ba":2.0, "ra":2.0,
                   "f":-1.0, "cl":-1.0, "br":-1.0, "i":-1.0, "at":-1.0,
                   "he":0.0, "ne":0.0, "ar":0.0, "kr":0.0, "xe":0.0, "rn":0.0}

AVOGADROS_NUMBER = 6.022e23
SAVE_STATE_DIRECTORY = "states/"

def load_structure_with_parmed(model, anchor):
    """
    Given the simulation inputs, load an anchor's structure for one of
    the checks and return the parmed structure object.
    """
    
    if anchor.amber_params is not None:
        #inpcrd_structure = None
        building_directory = os.path.join(
            model.anchor_rootdir, anchor.directory, 
            anchor.building_directory)
        if anchor.amber_params.prmtop_filename is not None:
            prmtop_filename = os.path.join(
                building_directory, anchor.amber_params.prmtop_filename)
        else:
            return None
        #if anchor.amber_params.inpcrd_filename is not None \
        #        and anchor.amber_params.inpcrd_filename != "":
        #    inpcrd_filename = os.path.join(
        #        building_directory, anchor.amber_params.inpcrd_filename)
        #    inpcrd_structure = parmed.load_file(prmtop_filename, 
        #                                        xyz=inpcrd_filename)
        #else:
        #    inpcrd_filename = None
        
        if anchor.amber_params.pdb_coordinates_filename is not None \
                and anchor.amber_params.pdb_coordinates_filename != "":
            pdb_filename = os.path.join(
                building_directory, 
                anchor.amber_params.pdb_coordinates_filename)
            structure = parmed.load_file(pdb_filename)
            if anchor.amber_params.box_vectors is not None:
                structure.box_vectors = anchor.amber_params.box_vectors.to_quantity()
        else:
            return None

        #elif inpcrd_filename is not None and inpcrd_filename != "":
        #    #structure = parmed.load_file(prmtop_filename, xyz=inpcrd_filename)
        #    return None
        #else:
        #    # anchor has no structure files
        #    return None
        return structure
        
    elif anchor.forcefield_params is not None:
        forcefield_filenames = []
        for forcefield_filename in \
                anchor.forcefield_params.built_in_forcefield_filenames:
            forcefield_filenames.append(forcefield_filename)
        for forcefield_filename in \
                anchor.forcefield_params.custom_forcefield_filenames:
            forcefield_filenames.append(os.path.join(
                building_directory, forcefield_filename))
        parameter_set = parmed.openmm.parameters.OpenMMParameterSet(
            forcefield_filenames)
        pdb_filename = os.path.join(building_directory, 
                               anchor.forcefield_params.pdb_filename)
        structure = parmed.load_file(parameter_set)
        pdb_structure = parmed.load_file(pdb_filename)
        structure.coordinates = pdb_structure.coordinates
        if anchor.forcefield_params.box_vectors is not None:
                structure.box_vectors \
                    = anchor.forcefield_params.box_vectors.to_quantity()
        return structure
    
    elif anchor.charmm_params is not None:
        raise Exception("Charmm systems not yet implemented")
        
    else:
        return None

def load_structure_with_mdtraj(model, anchor, mode="pdb", coords_filename=None):
    """
    Given the simulation inputs, load an anchor's structure for one of
    the checks and return the mdtraj Trajectory() object.
    """
    
    building_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.building_directory)
    prod_directory = os.path.join(
        model.anchor_rootdir, anchor.directory, anchor.production_directory)
    if mode == "pdb":
        pass
    elif mode == "elber_umbrella":
        umbrella_dcd_path = os.path.join(prod_directory)
        umbrella_basename = elber_base.ELBER_UMBRELLA_BASENAME+"*.dcd"
        umbrella_traj_glob = os.path.join(umbrella_dcd_path, umbrella_basename)
        umbrella_traj_filenames = glob.glob(umbrella_traj_glob)
        if len(umbrella_traj_filenames) == 0:
            return None
        
        # search for and remove any empty trajectory files
        indices_to_pop = []
        for i, umbrella_traj_filename in enumerate(umbrella_traj_filenames):
            if os.path.getsize(umbrella_traj_filename) == 0:
                indices_to_pop.append(i)
        
        for i in indices_to_pop[::-1]:
            umbrella_traj_filenames.pop(i)
            
        assert len(umbrella_traj_filenames) > 0, "Only empty umbrella " \
        "trajectories were found. You can force SEEKR to skip these checks "\
        "by using the --skip_checks (-s) argument"
    
    elif mode == "state_xml":
        assert coords_filename is not None
        
    elif mode == "mmvt_traj":
        mmvt_traj_dcd_path = os.path.join(prod_directory)
        mmvt_traj_basename = mmvt_base.OPENMMVT_BASENAME+"*.dcd"
        mmvt_traj_glob = os.path.join(mmvt_traj_dcd_path, mmvt_traj_basename)
        mmvt_traj_filenames = glob.glob(mmvt_traj_glob)
        if len(mmvt_traj_filenames) == 0:
            return None
        
        # search for and remove any empty trajectory files
        indices_to_pop = []
        for i, mmvt_traj_filename in enumerate(mmvt_traj_filenames):
            if os.path.getsize(mmvt_traj_filename) == 0:
                indices_to_pop.append(i)
        
        for i in indices_to_pop[::-1]:
            mmvt_traj_filenames.pop(i)
    else:
        raise Exception("Check mode not implemented: {}".format(mode))
        
    if anchor.amber_params is not None:
        
        if anchor.amber_params.prmtop_filename is not None:
            prmtop_filename = os.path.join(
                building_directory, anchor.amber_params.prmtop_filename)
        else:
            return None
        #if anchor.amber_params.inpcrd_filename is not None:
        #    inpcrd_filename = os.path.join(
        #        building_directory, anchor.amber_params.inpcrd_filename)
        #else:
        #    inpcrd_filename = None
        
        if mode == "pdb":
            if anchor.amber_params.pdb_coordinates_filename is not None \
                    and anchor.amber_params.pdb_coordinates_filename != "":
                pdb_filename = os.path.join(
                    building_directory, 
                    anchor.amber_params.pdb_coordinates_filename)
                traj = mdtraj.load(pdb_filename) #, top=prmtop_filename)
            #elif inpcrd_filename is not None and inpcrd_filename != "":
            #    #structure = parmed.load_file(prmtop_filename, xyz=inpcrd_filename)
            #    return None
            else:
                # anchor has no structure files
                return None
        
        elif mode == "elber_umbrella":
            assert len(umbrella_traj_filenames) > 0, "Only empty umbrella " \
                "trajectories were found. You can force SEEKR to skip these "\
                "checks by using the --skip_checks (-s) argument"
            traj = mdtraj.load(umbrella_traj_filenames, top=prmtop_filename)
        
        elif mode == "state_xml":
            traj = mdtraj.load_xml(coords_filename, top=prmtop_filename)
        
        elif mode == "mmvt_traj":
            assert len(mmvt_traj_filenames) > 0, "Only empty mmvt " \
                "trajectories were found. You can force SEEKR to skip these "\
                "checks by using the --skip_checks (-s) argument"
            traj = mdtraj.load(mmvt_traj_filenames, top=prmtop_filename)
            
        return traj
        
    elif anchor.forcefield_params is not None:
        forcefield_filenames = []
        for forcefield_filename in \
                anchor.forcefield_params.built_in_forcefield_filenames:
            forcefield_filenames.append(forcefield_filename)
        for forcefield_filename in \
                anchor.forcefield_params.custom_forcefield_filenames:
            forcefield_filenames.append(os.path.join(
                building_directory, forcefield_filename))
        parameter_set = parmed.openmm.parameters.OpenMMParameterSet(
            forcefield_filenames)
        pdb_filename = os.path.join(building_directory, 
                               anchor.forcefield_params.pdb_filename)
        if mode == "pdb":
            traj = mdtraj.load(pdb_filename)
        elif mode == "elber_umbrella":
            assert len(umbrella_traj_filenames) > 0, "Only empty umbrella " \
                "trajectories were found. You can force SEEKR to skip these "\
                "checks by using the --skip_checks (-s) argument"
            traj = mdtraj.load(umbrella_traj_filenames, top=pdb_filename)
        elif mode == "state_xml":
            traj = mdtraj.load_xml(coords_filename, top=prmtop_filename)
        elif mode == "mmvt_traj":
            assert len(mmvt_traj_filenames) > 0, "Only empty mmvt " \
                "trajectories were found. You can force SEEKR to skip these "\
                "checks by using the --skip_checks (-s) argument"
            traj = mdtraj.load(mmvt_traj_filenames, top=pdb_filename)
        
        return traj
    
    elif anchor.charmm_params is not None:
        raise Exception("Charmm systems not yet implemented")
        
    else:
        return None

def is_ion(atom):
    """If a lone atom has no bonds, assume it's an ion."""
    
    if len(atom.bond_partners) == 0:
        return True
    else:
        return False

def check_pre_sim_md_anchor_dir_per_md_anchor(model):
    """
    Makes sure that there is an anchor directory in the anchor root
    directory for every molecular dynamics anchor defined in the
    Model object.
    

    Parameters
    ----------
    model : Model
        The SEEKR2 model object containing all calculation information.

    Returns
    -------
    : bool
        False if the number of md anchor directories does not equal the
        number of md anchors.

    """
    num_md_anchor_directories = len([
        md_anchor_dir for md_anchor_dir 
        in os.listdir(model.anchor_rootdir) 
            if "anchor_" in md_anchor_dir
    ])
    num_md_anchors = len([
        anchor for anchor in model.anchors 
            if anchor.md == True
    ])
    if num_md_anchor_directories == num_md_anchors:
        return True
    return False

def check_pre_sim_bubbles(model):
    """
    Checks starting pdb structures for water box bubbles.

    Parameters
    ----------
    model : Model
        The SEEKR2 model object containing all calculation information.

    Returns
    -------
    : bool
        False if a bubble is detected and True otherwise.

    """
    
    for anchor in model.anchors:
        building_directory = os.path.join(
            model.anchor_rootdir, 
            anchor.directory, 
            anchor.building_directory
        )
        if anchor.amber_params is not None:
            if (anchor.amber_params.pdb_coordinates_filename is not None \
                and anchor.amber_params.pdb_coordinates_filename != ""
            ):
                pdb_filename = os.path.join(
                    building_directory,
                    anchor.amber_params.pdb_coordinates_filename
                )
                bvecs = anchor.amber_params.box_vectors
                bvecs = np.array([
                    [bvecs.ax, bvecs.ay, bvecs.az],
                    [bvecs.bx, bvecs.by, bvecs.bz],
                    [bvecs.cx, bvecs.cy, bvecs.cz],
                ])
                pdb_periodic_box_properties = \
                    bubblebuster.periodic_box_properties(
                        pdb_filename,
                        box_vectors=np.array(bvecs, dtype=np.float32)
                    )
                if pdb_periodic_box_properties.has_bubble:
                    warnstr = "CHECK FAILURE: Water box bubble detected "\
                    "in one or more starting structures." 
                    print(warnstr)
                else:
                    box_vectors = np.array([
                        *bvecs
                    ], dtype=np.float32)
                    pdb_periodic_box_properties = \
                        bubblebuster.periodic_box_properties(
                            pdb_filename,
                            box_vectors=box_vectors
                        )
                    if pdb_periodic_box_properties.has_bubble:
                        warnstr = "CHECK FAILURE: Water box bubble detected "\
                        "in one or more starting structures." 
                        print(warnstr)
                        return False
    return True

def check_hmr_timestep(model):
    """
    Checks the simulation time step size if HMR is implemented.
    If HMR is implemented and the time step equals 0.002 picoseconds,
    the user likely forgot to increase the time step in the input 
    XML file. 

    Parameters
    ----------
    model : Model
        The SEEKR2 model object containing all calculation information.

    Returns
    -------
    : bool
        False if HMR is implemented and the time step equals 0.002 
        picoseconds, otherwise True.

    """
    if model.openmm_settings:
        if model.openmm_settings.hydrogenMass is not None:
            if (model.openmm_settings.langevin_integrator.timestep \
                <= 0.002
            ):
                warnstr = "CHECK FAILURE: HMR is implemented, but a "\
                "0.002 ps time step is being used. Increase the time step "\
                "in the input XML file to take advantage of the "\
                "performance benefits HMR offers." 
                print(warnstr)
                return False
    return True
    

def check_pre_sim_MD_and_BD_salt_concentration(model):
    """
    Users might inadvertently define different salt concentrations 
    between the MD and BD stages. 
    Examine BD settings and count the numbers of ions in MD starting
    structures to ensure that ion concentrations are relatively 
    consistent to avoid situations where different salt concentrations
    exist between the various anchors and scales.
    """
    
    RELATIVE_TOLERANCE = 0.25
    ABSOLUTE_TOLERANCE = 0.1
    if model.k_on_info:
        bd_ionic_strength = 0.0
        for ion in model.k_on_info.ions:
            bd_ionic_strength += 0.5 * ion.conc * ion.charge**2
    else:
        # No BD to check
        return True
    for anchor in model.anchors:
        md_ionic_strength = 0.0
        structure = load_structure_with_parmed(model, anchor)
        if structure is None:
            continue
        box_6_vector = structure.get_box()
        assert box_6_vector is not None, "Unable to compute box volume for "\
            "structures in anchor {}".format(anchor.index)
        box_vectors = base.Box_vectors()
        box_vectors.from_6_vector(box_6_vector[0])
        box_volume = box_vectors.get_volume()
        particle_concentration = 1.0e24 / (AVOGADROS_NUMBER * box_volume)
        for index, atom in enumerate(structure.atoms):
            if is_ion(atom):
                found_ion = False
                for key in ION_CHARGE_DICT:
                    if atom.name.lower().startswith(key):
                        charge = ION_CHARGE_DICT[key]
                        md_ionic_strength += particle_concentration * charge**2
                        found_ion = True
                    
                if not found_ion:
                    print("Found unbonded atom with index: "\
                          "{} and name: {}.".format(index, atom.name)\
                          +"Charge is uncertain, assuming zero.")
        if not np.isclose(md_ionic_strength, bd_ionic_strength, 
                          rtol=RELATIVE_TOLERANCE, atol=ABSOLUTE_TOLERANCE):
            print("""CHECK FAILURE: BD simulation has significantly different
    ionic strength of {} M*e^2 than MD simulation ionic strength of 
    {} M*e^2 for anchor {}. Please check the ion concentrations in the 
    BD simulation settings, and also count the number of ions 
    in the MD simulations.""".format(bd_ionic_strength, 
                                 md_ionic_strength,
                                 anchor.index))
            return False
        
    return True

def check_systems_within_Voronoi_cells(model):
    """
    Users might provide the wrong starting structure for a given 
    anchor, and the structure may actually belong in a different one.
    When the model defines Voronoi cells (such as in MMVT), check that
    the starting structures lie within the expected Voronoi cells, and
    suggest corrections if the check fails. Otherwise, the SEEKR2
    backend would fail with a non-helpful error message.
    """
    
    if model.get_type() != "mmvt":
        # only apply to MMVT systems
        return True
    for anchor in model.anchors:
        traj = load_structure_with_mdtraj(model, anchor)
        if traj is None:
            continue
        for milestone in anchor.milestones:
            cv = model.collective_variables[milestone.cv_index]
            result = cv.check_mdtraj_within_boundary(traj, milestone.variables, 
                                                     verbose=True)
            if result == False:
                correct_anchor = None
                for anchor2 in model.anchors:
                    within_milestones = True
                    for milestone in anchor2.milestones:
                        result2 = cv.check_mdtraj_within_boundary(
                            traj, milestone.variables)
                        if not result2:
                            within_milestones = False
                            break
                    if within_milestones:
                        correct_anchor = anchor2
                        break
                
                warnstr = """CHECK FAILURE: The provided initial starting 
    structure for anchor {} does not lie within the 
    anchor boundaries. The simulation would fail. Please 
    carefully check this anchor's starting structure, as 
    well as the collective variable's (CV) atom selections, 
    and anchor/milestone variables.""".format(anchor.index)
                print(warnstr)
                if correct_anchor is not None:
                    print("It looks like the failed structure might belong in "\
                          "anchor {}.".format(correct_anchor.index))
                return False
    
    return True

def recurse_atoms(atom, _visited_indices=set()):
    """
    Recursively visit all atoms within a molecule for the purposes
    of determining the molecules (all sets of atoms connected by bonds)
    in the system.
    """
    
    _visited_indices.add(atom.idx)
    for bonded_atom in atom.bond_partners:
        if not bonded_atom.idx in _visited_indices:
            branch_indices = recurse_atoms(bonded_atom, _visited_indices)
            _visited_indices.update(branch_indices)
    return _visited_indices

def check_atom_selections_on_same_molecule(model):
    """
    The user might accidentally define atom selections that span
    multiple molecules. Check this possibility by finding all molecules
    in the system and ensure that atom selections only exist on one
    molecule.
    """
    
    warnstr1 = """CHECK FAILURE: the atom selection for collective variable 
    (CV) number {} is split over multiple molecules. Atom index {}, 
    which has the name {} and serial id {}, was the first atom to 
    be detected in a different molecule than the others. Please
    check the structure for anchor {} to ensure that atom indexing 
    is correct. Keep in mind: SEEKR2 atom indexing starts at 0 (but 
    PDB files often start with a different atom serial index, such
    as 1 or another number)."""
    warnstr2 = """CHECK FAILURE: the atom index {} for collective variable
    (CV) number {} does not exist in the structure for anchor {}."""
    for anchor in model.anchors:
        structure = load_structure_with_parmed(model, anchor)
        if structure is None:
            continue
        molecules = []
        traversed_atoms = set()
        for index, atom in enumerate(structure.atoms):
            if index in traversed_atoms:
                continue
            molecule = recurse_atoms(atom, set())
            traversed_atoms.update(molecule)
            molecules.append(molecule)
        
        for cv in model.collective_variables:
            cv_in_anchor = False
            for milestone in anchor.milestones:
                if milestone.cv_index == cv.index:
                    cv_in_anchor = True
                    
            if not cv_in_anchor:
                continue
            
            atom_groups = cv.get_atom_groups()
            for atom_group in atom_groups:
                in_molecule = None
                for atom_index in atom_group:
                    index_found = False
                    for mol_id, molecule in enumerate(molecules):
                        if atom_index in molecule :
                            if in_molecule is None:
                                in_molecule = mol_id
                                index_found = True
                            elif in_molecule == mol_id:
                                index_found = True
                            else:
                                print(warnstr1.format(
                                    cv.index, atom_index, 
                                    structure.atoms[atom_index].name, 
                                    structure.atoms[atom_index].number, 
                                    anchor.index))
                                return False
                    if not index_found:
                        print(warnstr2.format(atom_index, cv.index, 
                                              anchor.index))
                        return False
    return True

def check_atom_selections_MD_BD(model):
    """
    Users might accidentally select atoms that are not equivalent
    between the MD and BD stages. Detect whether this is the case.
    """
    
    warnstr1 = """CHECK FAILURE: The atom selection as defined for BD
    milestone {} includes atom index {} that does not exist 
    in the file {}. Please check the atom indices in Browndye
    inputs and also check the atom indices in your input PQR
    files."""
    warnstr2 = """CHECK FAILURE: Different atoms are selected between
    the MD and BD stages. The BD selection contains {} while the
    MD selection contains {}. Please check the input molecular
    structure files, as well as the input atom indices (keep in
    mind that SEEKR2 requires atom indexing to start from 0 in 
    each molecular structure input file, but PDB and PQR files
    might start their numbering at 1 or another number."""
    if model.browndye_settings is not None:
        b_surface_dir = os.path.join(
            model.anchor_rootdir, model.k_on_info.b_surface_directory)
        rec_pqr_path = os.path.join(
            b_surface_dir, model.browndye_settings.receptor_pqr_filename)
        lig_pqr_path = os.path.join(
            b_surface_dir, model.browndye_settings.ligand_pqr_filename)
        rec_pqr_structure = parmed.load_file(rec_pqr_path)
        lig_pqr_structure = parmed.load_file(lig_pqr_path)
        for bd_index, bd_milestone in enumerate(model.k_on_info.bd_milestones):
            receptor_indices = bd_milestone.receptor_indices
            for receptor_index in receptor_indices:
                if receptor_index < 0 \
                        or receptor_index > len(rec_pqr_structure.atoms):
                    print(warnstr1.format(bd_index, receptor_index, 
                                          rec_pqr_path))
                    return False
            ligand_indices = bd_milestone.ligand_indices
            for ligand_index in ligand_indices:
                if ligand_index < 0 \
                        or ligand_index > len(lig_pqr_structure.atoms):
                    print(warnstr1.format(bd_index, ligand_index, 
                                          lig_pqr_path))
                    return False
    else:
        # No BD to check
        return True
    
    for anchor in model.anchors:
        md_structure = load_structure_with_parmed(model, anchor)
        if md_structure is None:
            continue
        md_bd_cv_pairs = []
        for md_cv in model.collective_variables:
            for bd_index, bd_milestone in enumerate(
                    model.k_on_info.bd_milestones):
                if bd_milestone.outer_milestone.cv_index == md_cv.index:
                    md_bd_cv_pairs.append([md_cv.index, bd_index])
        
        for md_bd_pair in md_bd_cv_pairs:
            md_cv = model.collective_variables[md_bd_pair[0]]
            bd_milestone = model.k_on_info.bd_milestones[md_bd_pair[1]]
            md_atom_groups = sorted(md_cv.get_atom_groups(), 
                                    key=lambda L: (len(L), L))
            bd_atom_groups = sorted([bd_milestone.receptor_indices, 
                                    bd_milestone.ligand_indices],
                                    key=lambda L: (len(L), L))
            # check atom selection lengths
            for md_atom_group, bd_atom_group in zip(
                    md_atom_groups, bd_atom_groups):
                if len(md_atom_group) != len(bd_atom_group):
                    bd_error_str = "{} atoms".format(len(bd_atom_group))
                    md_error_str = "{} atoms".format(len(md_atom_group))
                    print(warnstr2.format(bd_error_str, md_error_str))
                    return False
                md_atomic_number_dict = collections.defaultdict(int)
                bd_atomic_number_dict = collections.defaultdict(int)
                for md_atom_idx in md_atom_group:
                    md_atom = md_structure.atoms[md_atom_idx]
                    md_atomic_number_dict[md_atom.atomic_number] += 1
                rec_pqr_indices = bd_milestone.receptor_indices
                lig_pqr_indices = bd_milestone.ligand_indices
                if bd_atom_group == rec_pqr_indices:
                    for bd_atom_idx in rec_pqr_indices:
                        bd_atom = rec_pqr_structure.atoms[bd_atom_idx]
                        bd_atomic_number_dict[bd_atom.atomic_number] += 1
                elif bd_atom_group == lig_pqr_indices:
                    for bd_atom_idx in lig_pqr_indices:
                        bd_atom = lig_pqr_structure.atoms[bd_atom_idx]
                        bd_atomic_number_dict[bd_atom.atomic_number] += 1
                else:
                    raise Exception("BD/MD mismatch in atom groups")
                
                if md_atomic_number_dict != bd_atomic_number_dict:
                    bd_err_str = ""
                    for key in bd_atomic_number_dict:
                        this_str = "{} atoms with atomic number {}, ".format(
                            bd_atomic_number_dict[key], key)
                        bd_err_str += this_str
                    md_err_str = ""
                    for key in md_atomic_number_dict:
                        this_str = "{} atoms with atomic number {}, ".format(
                            md_atomic_number_dict[key], key)
                        md_err_str += this_str
                    print(warnstr2.format(bd_err_str, md_err_str))
                    return False
    return True

def check_pre_simulation_all(model):
    """
    After the completion of the prepare stage, check inputs for some
    of the most common problems and mistakes a user is likely to 
    encounter. If a check fails, raise an Exception.
    
    Parameters:
    -----------
    model : Model()
        The SEEKR2 model object containing all calculation information.
        
    Returns:
    None
    """
    
    check_passed_list = []
    check_passed_list.append(
        check_pre_sim_md_anchor_dir_per_md_anchor(model)
    )
    check_passed_list.append(check_pre_sim_bubbles(model))
    check_passed_list.append(check_hmr_timestep(model))
    check_passed_list.append(check_pre_sim_MD_and_BD_salt_concentration(model))
    check_passed_list.append(check_atom_selections_on_same_molecule(model))
    check_passed_list.append(check_systems_within_Voronoi_cells(model))
    check_passed_list.append(check_atom_selections_MD_BD(model))
    
    no_failures = True
    for check_passed in check_passed_list:
        if not check_passed:
            no_failures = False
            
    if no_failures:
        print("All pre-simulation checks passed.")
        return
    else:
        check_fail_str = "One or more fatal pre-simulation checks failed. It "\
        "is highly recommended that you address and correct each of these "\
        "problems. However, you can force SEEKR to skip these checks by using "\
        "the --skip_checks (-s) argument on prepare.py."
        print(check_fail_str)
        raise Exception("The SEEKR2 calculation can not proceed due "\
                        "to failed pre-simulation checks.")
    return

def check_elber_umbrella_stage(model):
    """
    If an umbrella force is improperly constructed or because of a bug,
    a system may deviate too far from the milestone during an umbrella
    stage. Detect this situation to alert the user.
    """
    
    warnstr = """CHECK FAILURE: Elber umbrella stage trajectory for anchor {} 
    deviates significantly from the location of the central milestone.
    This could be caused by a bad umbrella force constant (too large or
    too small) umbrella_force_constant value.
    If may also happen if the model.xml file was improperly modified,
    or perhaps due to a bug."""
    if model.get_type() != "elber":
        # this test would be irrelevant
        return True
    for anchor in model.anchors:
        traj = load_structure_with_mdtraj(model, anchor, mode="elber_umbrella")
        if traj is None:
            continue
        for milestone in anchor.milestones:
            if milestone.alias_index == 2:
                cv = model.collective_variables[milestone.cv_index]
                result = cv.check_mdtraj_close_to_boundary(
                    traj, milestone.variables, verbose=True)
                if result == False:
                    print(warnstr.format(anchor.index))
                    return False
    return True

def check_xml_boundary_states(model):
    """
    SEEKR2 calculations generate OpenMM XML state files when 
    boundaries are encountered. Ensure that these boundary encounters
    are properly accounted and that they exist where expected.
    """
    
    warnstr = """CHECK FAILURE: Saved boundary states in anchor {} were not
    saved close to the expected milestone(s). File name: {}. 
    This could be caused if the model.xml file was improperly 
    modified, or perhaps due to a bug."""
    for anchor in model.anchors:
        states_dir = os.path.join(
            model.anchor_rootdir, anchor.directory, 
            anchor.production_directory, SAVE_STATE_DIRECTORY)
        for milestone in anchor.milestones:
            if model.get_type() == "elber" and milestone.alias_index == 2:
                continue
            state_file_glob = os.path.join(
                states_dir, "*_*_%d" % milestone.alias_index)
            state_files = glob.glob(state_file_glob)
            if len(state_files) == 0:
                continue
            for state_file in state_files:
                traj = load_structure_with_mdtraj(
                    model, anchor, mode="state_xml", coords_filename=state_file)
                cv = model.collective_variables[milestone.cv_index]
                result = cv.check_mdtraj_close_to_boundary(
                    traj, milestone.variables, verbose=True)
                save_pdb_filename = state_file+".pdb"
                traj.save_pdb(save_pdb_filename)
                if result == False:
                    save_pdb_filename = os.path.join(states_dir, "problem.pdb")
                    traj.save_pdb(save_pdb_filename)
                    print(warnstr.format(anchor.index, state_file))
                    print("Saved problematic structure to:", save_pdb_filename)
                    return False
    
    return True

def check_mmvt_in_Voronoi_cell(model):
    """
    For some reason, a MMVT system may drift out of a Voronoi cell.
    This would most likely indicate a bug, so detect this situation.
    """
    
    if model.get_type() != "mmvt":
        # irrelevant test
        return True
    for anchor in model.anchors:
        traj = load_structure_with_mdtraj(model, anchor, mode="pdb")
        if traj is None:
            continue
        for milestone in anchor.milestones:
            cv = model.collective_variables[milestone.cv_index]
            result = cv.check_mdtraj_within_boundary(traj, milestone.variables, 
                                                verbose=True)
            if result == False:
                warnstr = """CHECK FAILURE: The MMVT trajectory(ies)
    for anchor {} do not lie within the 
    anchor boundaries. This could be caused by
    careless tampering with files or the model.xml
    file, or possibly a bug.""".format(anchor.index)
                print(warnstr)
                return False
    
    return True

def find_parmed_structure_com(structure, indices):
    """
    For a parmed structure, find the center of mass (COM) of a set
    of atoms defined by "indices".
    """
    
    total_mass = 0.0
    com_vector = np.zeros(3)
    for i, atom in enumerate(structure.atoms):
        total_mass += atom.mass
        com_vector += atom.mass * structure.coordinates[i]
    
    assert total_mass > 0.0
    return com_vector / total_mass

def check_bd_simulation_end_state(model):
    """
    The BD stage simulation program, Browndye2, can save encounter
    complexes when the system encounters a milestone. Check these
    encounter complex structures to ensure that they exist close to
    the expected milestone.
    """
    
    ATOL = 0.1
    if model.k_on_info:
        bd_ionic_strength = 0.0
        for ion in model.k_on_info.ions:
            bd_ionic_strength += 0.5 * ion.conc * ion.charge**2
    else:
        # No BD to check
        return True
    
    for bd_index, bd_milestone in enumerate(model.k_on_info.bd_milestones):
        outer_radius = bd_milestone.outer_milestone.variables["radius"]
        fhpd_directory = os.path.join(
            model.anchor_rootdir, bd_milestone.directory, 
            bd_milestone.fhpd_directory)
        directory_glob = os.path.join(fhpd_directory, "lig*")
        directories = glob.glob(directory_glob)
        for directory in directories:
            lig_pqr_filename = os.path.join(directory, "ligand.pqr")
            rec_pqr_filename = os.path.join(directory, "receptor.pqr")
            lig_struct = parmed.load_file(lig_pqr_filename)
            rec_struct = parmed.load_file(rec_pqr_filename)
            lig_com = find_parmed_structure_com(
                lig_struct, bd_milestone.ligand_indices)
            rec_com = find_parmed_structure_com(
                rec_struct, bd_milestone.receptor_indices)
            distance = np.linalg.norm(rec_com - lig_com) / 10.0 # A_to_nm
            if not np.isclose(distance, outer_radius, atol=ATOL):
                warnstr = """CHECK FAILURE: The BD milestone number {}
    has saved FHPD structures from the b-surface simulations
    that were significantly far away from the outermost 
    milestone of this binding site. Observed distance: {}. 
    Expected distance: {}.""".format(bd_milestone.index, distance, outer_radius)
                print(warnstr)
                return False
    return True

def check_post_simulation_all(model, long_check=False):
    """
    After the completion of the run stage, check simulation files
    for some of the most common problems and mistakes a user is likely 
    to encounter. If a check fails, raise an Exception.
    
    Parameters:
    -----------
    model : Model()
        The SEEKR2 model object containing all calculation information.
    
    long_check : bool, Default False
        Whether to conduct a detailed check of post-simulation outputs.
        If set to True, the checks are likely to take a significantly
        longer amount of time, but is more likely to detect problems.
    
    Returns:
    None
    """
    
    check_passed_list = []
    check_passed_list.append(check_elber_umbrella_stage(model))
    check_passed_list.append(check_mmvt_in_Voronoi_cell(model))
    check_passed_list.append(check_bd_simulation_end_state(model))
    if long_check:
        check_passed_list.append(check_xml_boundary_states(model))
    
    no_failures = True
    for check_passed in check_passed_list:
        if not check_passed:
            no_failures = False
            
    if no_failures:
        print("All post-simulation checks passed.")
        return
    else:
        check_fail_str = "One or more fatal post-simulation checks failed. It "\
        "is highly recommended that you address and correct each of these "\
        "problems. However, you can force SEEKR to skip these checks by using "\
        "the --skip_checks (-s) argument on analyze.py."
        print(check_fail_str)
        raise Exception("The SEEKR2 calculation can not proceed due "\
                        "to failed post-simulation checks.")
    return

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description="Check SEEKR2 calculations for common problems both "\
        "before and after simulations.")
    
    argparser.add_argument(
        "model_file", metavar="MODEL_FILE", type=str, 
        help="name of model file for OpenMMVT calculation. This would be the "\
        "XML file generated in the prepare stage.")
    
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    xmlfile = args["model_file"]
    model = base.Model()
    model.deserialize(xmlfile)
    if model.anchor_rootdir == ".":
        model_dir = os.path.dirname(xmlfile)
        model.anchor_rootdir = os.path.abspath(model_dir)
    check_pre_simulation_all(model)
    check_post_simulation_all(model, long_check=True)
    
