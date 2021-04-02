from hatchet.utils.solve import solve
from hatchet.utils.solve.utils import segmentation

if __name__ == '__main__':

    _, cA, cB, u = solve(
        clonal='28:1:1',                      # -c, examples: 28:1:1/28:2:2,24:2:0,15:1:2/None
        seg_file='bbc/bulk.seg',
        n=3,                                  # -n, Number of distinct clones, 2-6 sweep by default
        solve_mode='cd',                      # 'cd' or 'ilp'; corresponds closely to
                                              # -M, 0=ILP+COORD_DESC, 1=ILP, 2=COORD_DESC (default 2)
        d=-1,                                 # -d, Maximum number of distinct copy-number states, -1 = no limit
        cn_max=6,                             # -e, Maximum copy number
        mu=0.03,                              # -u, Minimum tumor-clone threshold
        diploid_threshold=0.1,                # -t, hardcoded; argument exposed but not used in C++
        ampdel=True,                          # -f, mutated allele amplification/deletion across tumor clones
        n_seed=400,                           # -p, no. of seeds for coordinate descent
        n_worker=8,                           # -j, no. of parallel workers for coordinate descent
    )

    segmentation_df = segmentation(cA, cB, u, bbc_file='bbc/bulk.bbc', bbc_out_file='out.bbc.ucn.tsv',
                                   seg_out_file='out.seg.ucn.tsv')
