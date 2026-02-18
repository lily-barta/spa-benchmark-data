import tequila as tq
import sunrise as sun
import csv
import time
import spafastprototype as fspa
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

def generate_data_point(n, iter, max_iter, d_min=0.5, d_max=4.0, nroots=1, get_fci=False, get_var=False):
    start=time.time()
    print(f"\nIteration : {iter + 1} / {max_iter}")

    ## PREPARE MOLECULE AND CIRCUIT
    geometry, distance = generate_geometry(n, iter, max_iter, d_min=d_min, d_max=d_max)
    print(f"Interatomic distance = {distance:.5f}")
    mol = sun.Molecule(geometry=geometry,basis_set='sto-3g',nature='h',transformation='reordered-jordan-wigner').use_native_orbitals()
    
    edges, guess = mol.get_spa_edges_and_guess()
    U = mol.make_spa_ansatz(edges=edges, hcb=True)

    ## OPTIMISE ORBITALS
    oo_start=time.time()
    opt = fspa.fast_spa_orb_opt(mol=mol, edges=edges, initial_guess=guess.T)
    mol = opt.molecule
    oo_end = time.time()
    print(f"Orbital optimisation took {oo_end-oo_start}s") 

    ## COMPUTE SPA ENERGY 
    spa_start=time.time()
    res = fspa.fast_spa_vqe(mol, U)
    spa_end = time.time()
    print(f"VQE SPA  : {res.energy:+2.10f}")
    print(f"SPA energy took {spa_end-spa_start}s") 

    results = {
        "n": n,
        "distance": distance,
        "spa": res.energy,
        "oo_t": oo_end - oo_start,
        "spa_t": spa_end - spa_start
    }

    if get_fci:
        ## COMPUTE SPA WFN
        spa_start = time.time()
        wfn_spa_hcb = sun.simulate(U, variables=res.variables)
        results["spa_wfn_t"] = time.time() - spa_start
        print(f"SPA wfn took {time.time() - spa_start}s") 
    
        ## COMPARE TO FCI ENERGY
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
        results["fci"] = fci0
        results["fci_t"] = time.time() - fci_start       
        print(f"FCI      : {fci0:.10f}")
        print(f"Error    : {res.energy-fci0:.10f}")
        print(f"FCI took {time.time() - fci_start}s") 
    
        ## COMPUTE FIDELITY
        fid_start = time.time()
        if nroots == 1:
            fidelity = abs(wfn_spa_hcb.inner(wfn_fci))**2
        else:
            fidelity = 0.0
            for i in range(nroots):
                if abs(fci[i]-fci[0]) < 0.0016: 
                    fidelity += abs(wfn_spa_hcb.inner(wfn_fci[i]))**2
        results["fid"] = fidelity
        results["fid_t"] = time.time() - fid_start
        print(f"fidelity : {fidelity:.6f}")
        print(f"Fidelity took {time.time()-fid_start}s")

    ## COMPUTE VARIANCE
    if get_var:
        var_start=time.time()
        H = mol.make_hamiltonian()
        H2 = tq.ExpectationValue(H=H**2, U=U+mol.hcb_to_me())
        results["h2_t"] = time.time() - var_start
        print(f"Building H2 took {time.time()-var_start}s")
        E = tq.ExpectationValue(H=H, U=U+mol.hcb_to_me())
        var = sun.simulate(H2-E**2, variables=res.variables)
        results["var"] = var
        results["var_t"] = time.time() - var_start
        print(f"Variance : {var:.6f}")
        print(f"Variance took {time.time()-var_start}s")

    results["total_t"] = time.time() - start
    print(f"Data point took {time.time()-start}s")
    return results

def run_dissociation(n, max_iter, d_min=0.5, d_max=4.0, nroots=1, filename=None, get_fci=False, get_var=False):
    """
    Arguments
    ----------
    n: Number of H atoms in H_n chain.
    max_iter: Total number of dissociation points.
    d_min, d_max: Minimum and maximum interatomic distance (in angstroms).
    nroots: Number of FCI eigenstates used as reference.
            Use nroots > 1 in the presence of near-degeneracies.
    filename: Name of the output CSV file.
             If None, default is f"results_h{n}.csv".
    get_fci: If True, compute and store FCI energies and fidelities.
    get_var: If True, compute and store the energy variance.
    
    Returns
    -------
    Generates a CSV file with one row per dissociation point.
    """
    if filename is None:
        filename = f"results_h{n}.csv"
    
    with open(filename, "w", newline="") as file:
        columns = ["distance", "spa", "fci", "fid", "var"]
        writer = csv.DictWriter(file, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        for iter in range(max_iter):
            data_dict = generate_data_point(n, iter, max_iter, d_min=d_min, d_max=d_max, nroots=nroots, get_fci=get_fci, get_var=get_var)
            writer.writerow(data_dict)

def run_single_point(n, distance, nroots=1, get_fci=False, get_var=False):
    """
    Arguments
    ----------
    n: Number of H atoms in H_n chain.
    distance: Interatomic distance (in angstroms).
    nroots: Number of FCI eigenstates used as reference.
            Use nroots > 1 in the presence of near-degeneracies. 

    Returns
    -------
    dictionary containing all the computed data
    """
    return generate_data_point(n=n, iter=0, max_iter=1, d_min=distance, nroots=nroots, get_fci=get_fci, get_var=get_var)

def run_scaling(n_min=2, n_max=10, distance=1.0, nroots=1, filename_t="timing_vs_n.csv", filename="results_vs_n.csv"):
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
    columns = ["n", "spa", "fci", "fid", "var"]
    columns_t = ["n", "total_t", "oo_t", "spa_t", "spa_wfn_t", "fci_t", "fid_t", "h2_t", "var_t"]

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
            data_dict = run_single_point(n=n, distance=distance, nroots=nroots, get_fci=get_fci, get_var=get_var)
            writer.writerow(data_dict)
            writer_t.writerow(data_dict)
