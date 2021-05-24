import os
import sys
import argparse
import sourmash
from sourmash import sourmash_args
from collections import Counter




def drop_below_mincount(hashCounter, min_count):
    for hashes, cnts in hashCounter.copy().items():
        if cnts < min_count:
            hashCounter.pop(hashes)
    return hashCounter

def build_hashCounter(sigs, hashCounter):
    for sig in sigs:
        hashCounter.update(sig.minhash.hashes) # add hashes to Counter
    return hashCounter


def main(args):
    min_count = args.minimum_count
    # select sigs that match specified ksize, moltype.
    sigs = sourmash_args.load_file_as_signatures(args.signatures_file, ksize=args.ksize, select_moltype=args.alphabet)

    # count hashes
    counts = Counter()
    counts = build_hashCounter(sigs, counts)
    print(f"dropping unique hashes for molecule: {args.alphabet}, ksize: {args.ksize}")
    # drop hashes below min_count threshold
    counts = drop_hashes(counts, min_count)

    # write out sig
    if args.output_sig:
        fresh_mh = sigs[0].minhash.copy_and_clear()
        # add hashes to fresh_mh
        fresh_mh.add_many(counts.keys())
        # build sig
        new_sig= sourmash.signature.SourmashSignature(fresh_mh, name=f"aggregated_hashvals_above_{min_count}", filename=args.signatures_file)
        with open(args.output_sig, "w") as out:
            sourmash.signature.save_signatures([new_sig], sigout)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--signatures-file', help='file containing list of signatures', required=True)
    p.add_argument('--minimum-count', help='keep hashes >= this minimum abundance', default=2)
    p.add_argument('--alphabet', help='alphabet', default="dna")
    p.add_argument('--ksize', help='ksize', default=31)

    # output options:
    p.add_argument('--output', help='store abundance filtered hashes as signature')
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
