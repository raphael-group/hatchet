import os
import os.path
import subprocess as pr
import gzip

import hatchet.utils.ArgParsing as ap
from hatchet.utils.Supporting import *
import hatchet.utils.Supporting as sp
from hatchet import config


def main(args=None):
    log(msg='# log notes\n', level='STEP')
    args = ap.parse_download_panel_arguments(args)
    logArgs(args)

    os.makedirs(args['refpaneldir'], exist_ok=True)

    # download reference panel, prepare files for liftover
    if args['refpanel'] == '1000GP_Phase3':
        # download 1000GP ref panel
        sp.download(
            url=config.urls.onekgp,
            dirpath=args['refpaneldir'],
            sentinel_file='1000GP_Phase3.sample',
        )
    else:
        error(
            'Currently, only the 1000 genome panel aligned to GRCh37 without "chr" prefix is supported\n',
            raise_exception=True,
        )

    # download necessary liftover files; 1000GP in hg19 coordinates
    # if users aligned reads to the hg38 build, we need to liftover coordinates to the reference panel (hg38 -> hg19)
    # since the 1000GP panel is in hg19 coordinates, we need to download (1) hg19  genome and (2) chain files
    # for liftover via picard
    dwnld_refpanel_genome(path=args['refpaneldir'])
    dwnld_chains(dirpath=args['refpaneldir'])

    # if users aligned reads to the same reference genome as used in the reference panel, liftover isn't required, but
    # there could be different naming conventions of chromosomes, with or without the 'chr' prefix. The 1000GP reference
    # panel does NOT use 'chr' prefix, so input into shapeit also should not have this
    mk_rename_file(path=args['refpaneldir'])


def dwnld_chains(dirpath):
    def mod_chain(infile, sample_chr, refpanel_index, sample_index):
        if sample_chr:
            name = infile.strip('.gz').replace('over', 'chr')
        else:
            name = infile.strip('.gz').replace('over', 'no_chr')

        with open(name, 'w') as new:
            with gzip.open(infile, 'rt') as f:
                for l in f:
                    if l.startswith('chain'):
                        l = l.split()
                        if not config.urls.refpanel_genome_chr_notation:
                            l[refpanel_index] = l[refpanel_index].replace(
                                'chr', ''
                            )
                        if not sample_chr:
                            l[sample_index] = l[sample_index].replace(
                                'chr', ''
                            )
                        new.write(' '.join(l) + '\n')
                    else:
                        new.write(l)
        return name

    hg38tohg19 = download(
        url=config.urls.refpanel_hg38tohg19,
        dirpath=dirpath,
        overwrite=False,
        extract=False,
    )
    hg19tohg38 = download(
        url=config.urls.refpanel_hg19tohg38,
        dirpath=dirpath,
        overwrite=False,
        extract=False,
    )

    # make all necessary chain files to convert from hg38 (w/ or w/out chr notation) to hg19 (no chr notation),
    # and also to lift back over from hg19 (no chr notation) to hg38 (w/ or w/out chr notation).

    # modify chr notation of hg38ToHg19, ref panel chr in 7th field, sample chr in 2nd field
    mod_chain(hg38tohg19, sample_chr=True, refpanel_index=7, sample_index=2)
    mod_chain(hg38tohg19, sample_chr=False, refpanel_index=7, sample_index=2)

    # modify chr notation of hg19ToHg38, ref panel chr in 2nd field, sample chr in 7th field
    mod_chain(hg19tohg38, sample_chr=True, refpanel_index=2, sample_index=7)
    mod_chain(hg19tohg38, sample_chr=False, refpanel_index=2, sample_index=7)


def dwnld_refpanel_genome(path):
    newref = os.path.join(path, 'hg19_no_chr.fa')
    if not os.path.isfile(newref):

        # If the genome reference file used in other parts of HATCHet matches the one we want, use it
        reference_file = config.paths.reference
        if (
            os.path.isfile(reference_file)
            and checksum(reference_file)
            == config.urls.refpanel_genome_checksum
        ):
            out = reference_file
        else:
            out = download(
                config.urls.refpanel_genome, dirpath=path, extract=False
            )

        _open, _mode = open, 'r'
        if out.endswith('gz'):
            _open, _mode = gzip.open, 'rt'

        # change chr notation
        with open(newref, 'w') as new:
            with _open(out, _mode) as f:
                for line in f:
                    if line.startswith('>'):
                        new.write(line.replace('chr', ''))
                    else:
                        new.write(line)

        # make dict file
        dict_file = newref.replace('.fa', '.dict')
        samtools = os.path.join(config.paths.samtools, 'samtools')
        cmd = f'{samtools} dict {newref} > {dict_file}'
        errname = os.path.join(path, 'samtools.log')
        with open(errname, 'w') as err:
            run = pr.run(
                cmd,
                stdout=err,
                stderr=err,
                shell=True,
                universal_newlines=True,
            )
        if run.returncode != 0:
            raise ValueError(
                sp.error(
                    f'Samtools dict creation failed, please check errors in {errname}!'
                )
            )
        else:
            os.remove(errname)
    return newref


def mk_rename_file(path):
    # makes rename_chrs1.txt for removing "chr", rename_chrs2.txt for adding "chr"
    names = [os.path.join(path, f'rename_chrs{i}.txt') for i in range(1, 3)]
    for i, n in enumerate(names):
        with open(n, 'w') as f:
            for j in range(1, 23):
                f.write(f'chr{j} {j}\n') if i == 0 else f.write(
                    f'{j} chr{j}\n'
                )
    return names


if __name__ == '__main__':
    main()
