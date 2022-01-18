pdb_dir = './datafiles/raw_pdb/'
proc_dir = './datafiles/processed_pdb/'
EGFR_KLIFS_IDENT = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 33, 34, 35, 43, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 61, 62, 63, 64, 65, 66, 67, 68, 69, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 114, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 139, 142, 143, 144, 145, 146, 147, 148, 149, 163, 165, 167, 168, 169, 171, 174, 185, 187, 190, 192, 195, 220, 227, 262, 274]

proc_pdb_list = []
NCAA = ['TPO', 'CSD', 'KCX', 'CSD', 'PTR', 'OCS', 'ALY']
pdb_proc = dist_analy.import_pdb.PDB_Processer(NCAA)
for pdb_fn in found_pdbs:
    pdb_list = pdb_proc.process_pdb(pdb_fn+'.pdb', pdb_dir, proc_dir, egfr_uniprot)
    for proc in pdb_list:
        print(proc)
        proc_pdb_list.append(proc)
