import tequila as tq
import sunrise as sun
from sunrise.measurement import get_hcb_part
import numpy
import csv
import time
import sys
# print(sys.getrecursionlimit())
sys.setrecursionlimit(30000)

def generate_geometry(n, iter, max_iter, d_min=0.5, d_max=4.0):
    if max_iter == 1:
        R = d_min
    else:
        R = d_min + (d_max - d_min) * iter / (max_iter - 1)

    geom = ""
    for i in range(n):
        geom += f"H 0.0 0.0 {i * R}\n"
    return geom, R

def get_edges_and_guess(n):
    edges = [(i, i + 1) for i in range(0, n, 2)]
    guess = numpy.eye(n)
    for edge in edges:
        guess[edge[0]][edge[1]] = 1
        guess[edge[1]][edge[0]] = -1
    return edges, guess.T

def generate_data_point(n, iter, max_iter, d_min=0.5, d_max=4.0, nroots=1, get_fci=False, get_var=False, verbose=False, get_coupling=False):
    start=time.time()
    print(f"\nIteration : {iter + 1} / {max_iter}")

    # Build tequila/sunrise molecule
    geometry, distance = generate_geometry(n, iter, max_iter, d_min=d_min, d_max=d_max)
    print(f"Interatomic distance = {distance:.5f}")
    mol = sun.Molecule(geometry=geometry,basis_set='sto-3g',transformation='reordered-jordan-wigner').use_native_orbitals()

    results = {
        "n": n,
        "distance": f"{distance:.3f}"
        }
    
    # Build SPA circuit
    edges, guess = get_edges_and_guess(n)
    U = mol.make_spa_ansatz(edges=edges, hcb=True)

    # Optimise orbitals
    oo_start=time.time()
    opt = sun.SPAFP.run_spa(mol=mol, edges=edges, initial_guess=guess)
    mol = opt.molecule
    results["oo_t"] = f"{time.time() - oo_start:.6f}"
    if verbose:
        print(f"Orb opt took {time.time()-oo_start:.6f}s") 

    # Fast HCB-SPA VQE
    spa_start=time.time()
    H = mol.make_hardcore_boson_hamiltonian()
    grouping = sun.SPAFP.make_decomposed_clusters(U)
    vqe_solver = sun.SPAFP.SPASolver(decompose=True,grouping=grouping)
    res = vqe_solver(H=H, circuit=U, molecule=mol)
    results["spa"] = f"{res.energy:.10f}"
    results["spa_t"] = f"{time.time() - oo_start:.6f}"
    if verbose:
        print(f"VQE SPA  : {res.energy:+2.10f}")
        print(f"SPA energy took {time.time()-spa_start:.6f}s") 

    if get_fci:
        # Get HCB-SPA wavefunction
        spa_start = time.time()
        wfn_spa_hcb = tq.simulate(U, variables=res.variables)
        print(wfn_spa_hcb)
        results["spa_wfn_t"] = f"{time.time() - spa_start:.6f}"
        if verbose:
            print(f"SPA wfn took {time.time() - spa_start:.6f}s") 
    
        # Compute FCI energy and wavefunction
        fci_start = time.time()
        nroots_map = {4: 6, 6: 20, 8: 70}
        if distance >= 2.5 and n > 8:
            print(f"\n!!! Warning !!! \nFCI for H{n} at distance {distance:.2f} Å is degenerate. \n"
                "Fidelity requires many FCI roots and is not computed due to high cost.\n")
        elif distance >= 2.7 and n in nroots_map:
            nroots = nroots_map[n]
        ci0 = wfn_spa_hcb if n > 6 else None
        fci, wfn_fci = mol.compute_energy("fci", get_wfn=True, nroots=nroots, ci0=ci0, use_hcb=True)
        fci0 = fci if nroots == 1 else fci[0]     
        results["fci"] = f"{fci0:.10f}"
        results["fci_t"] = f"{time.time() - fci_start:.6f}"
        if verbose:
            print(f"FCI      : {fci0:.10f}")
            print(f"Error    : {res.energy-fci0:.10f}")
            print(f"FCI energy+wfn took {time.time()-fci_start:.6f}s")
    
        # Compute fidelity (SPA/FCI overlap)
        fid_start = time.time()
        if nroots == 1:
            fidelity = abs(wfn_spa_hcb.inner(wfn_fci))**2
        else:
            fidelity = 0.0
            for i in range(nroots):
                if abs(fci[i]-fci[0]) < 0.0016: 
                    fidelity += abs(wfn_spa_hcb.inner(wfn_fci[i]))**2
        results["fid"] = f"{fidelity:.6f}"
        results["fid_t"] = f"{time.time() - fid_start:.6f}"
        if verbose:
            print(f"fidelity : {fidelity:.6f}")
            print(f"Fidelity took {time.time()-fid_start:.6f}s")

    # Compute variance
    if get_var:
        var_start=time.time()
        H = mol.make_hamiltonian()
        H2 = tq.ExpectationValue(H=H**2, U=U+mol.hcb_to_me())
        results["h2_t"] = f"{time.time() - var_start:.6f}"
        if verbose:
            print(f"Building H2 took {time.time()-var_start:.6f}s")
        E = tq.ExpectationValue(H=H, U=U+mol.hcb_to_me())
        var = sun.simulate(H2-E**2, variables=res.variables)
        results["var"] = f"{var:.6f}"
        results["var_t"] = f"{time.time() - var_start:.6f}"
        if verbose:
            print(f"Variance : {var:.6f}")
            print(f"Variance took {time.time()-var_start:.6f}s")

    # Compute the coupling between SPA/non-seniority zero
    if get_coupling:
        coupling_start = time.time()
        hcb_mol, res_mol = get_hcb_part(mol)
        H_res = res_mol.make_hamiltonian()
        H2_res = tq.ExpectationValue(H=H_res**2, U=U+mol.hcb_to_me())
        results["h2_res_t"] = f"{time.time() - coupling_start:.6f}"
        if verbose:
            print(f"Building H_res + H2_res took {time.time()-coupling_start:.6f}s")
        coupling_sq = sun.simulate(H2_res, variables=res.variables)
        coupling = numpy.sqrt(coupling_sq)
        results["coupling"] = f"{coupling:.6f}"
        results["coupling_t"] = f"{time.time() - coupling_start:.6f}"
        if verbose:
            print(f"Coupling : {coupling:.6f}")
            print(f"Coupling took {time.time()-coupling_start:.6f}s")

    results["total_t"] = f"{time.time() - start:.6f}"
    if verbose:
            print(f"Data point took {time.time()-start:.6f}s")
    return results


def run_dissociation(n, max_iter, d_min=0.5, d_max=4.0, nroots=1, filename=None, get_fci=False, get_var=False, verbose=False, get_coupling=False):
    """
    :param n: Number of H atoms in H_n chain.
    :param max_iter: Total number of dissociation points.
    :param d_min: Minimum interatomic distance (in angstroms).
    :param d_max: Maximum interatomic distance (in angstroms).
    :param nroots: Number of FCI eigenstates used as reference.
            Use nroots > 1 in the presence of near-degeneracies.
    :param filename: Name of the output CSV file.
             If None, default is f"results_h{n}.csv".
    :param get_fci: If True, compute and store FCI energies and fidelities.
    :param get_var: If True, compute and store the energy variance.
    :param get_coupling: If True, compute and store coupling between SPA and non-seniority zero
    :param verbose: If True, prints results
    
    :return: Generates a CSV file with one row per dissociation point.
    """
    if filename is None:
        filename = f"results_h{n}.csv"
    
    with open(filename, "w", newline="") as file:
        columns = ["distance", "spa", "fci", "fid", "var", "coupling"]
        writer = csv.DictWriter(file, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        for iter in range(max_iter):
            data_dict = generate_data_point(n, iter, max_iter, d_min=d_min, d_max=d_max, nroots=nroots, get_fci=get_fci, get_var=get_var, verbose=verbose, get_coupling=get_coupling)
            writer.writerow(data_dict)

def run_single_point(n, distance, nroots=1, get_fci=False, get_var=False, verbose=False, get_coupling=False):
    """
    :param n: Number of H atoms in H_n chain.
    :param distance: Interatomic distance (in angstroms).
    :param nroots: Number of FCI eigenstates used as reference.
            Use nroots > 1 in the presence of near-degeneracies. 

    :return: dictionary containing all the computed data
    """
    return generate_data_point(n=n, iter=0, max_iter=1, d_min=distance, nroots=nroots, get_fci=get_fci, get_var=get_var, verbose=verbose, get_coupling=get_coupling)

def run_scaling(n_min=2, n_max=10, distance=1.0, nroots=1, filename_t="timing_vs_n.csv", filename="results_vs_n.csv", verbose=False):
    """
    Arguments
    ----------
    n_min, n_max: Minimum and maximum number of H atoms in H_n chain.
    distance: Interatomic distance (in angstroms).
    nroots: Number of FCI eigenstates used as reference.
            Use nroots > 1 in the presence of near-degeneracies.
    filename_t: Name of the output CSV file for timing data (default "timing_vs_n.csv").
    filename: Name of the output CSV file for SPA/FCI/fidelity/variance results (default "results_vs_n.csv"). 
    
    Returns
    -------
    Generates two CSV files, one row per H_n chain: one for results and one for timings.
    """
    columns = ["n", "spa", "fci", "fid", "var", "coupling"]
    columns_t = ["n", "total_t", "oo_t", "spa_t", "spa_wfn_t", "fci_t", "fid_t", "h2_t", "var_t", "h2_res_t", "coupling_t"]

    with open(filename, "w", newline="") as file, \
         open(filename_t, "w", newline="") as file_t:

        writer = csv.DictWriter(file, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer_t = csv.DictWriter(file_t, fieldnames=columns_t, extrasaction='ignore')
        writer_t.writeheader()

        for n in range(n_min, n_max + 1, 2):
            print(f"\n===== Running H{n} =====")
            get_fci = True if n <= 14 else False
            get_var = True if n <= 10 else False
            get_coupling = True if n <= 10 else False
            data_dict = run_single_point(n=n, distance=distance, nroots=nroots, get_fci=get_fci, get_var=get_var, verbose=verbose, get_coupling=get_coupling)
            writer.writerow(data_dict)
            writer_t.writerow(data_dict)
