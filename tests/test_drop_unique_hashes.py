""" test drop_unique_hashes.py"""

from collections import Counter
from sourmash import MinHash, SourmashSignature
from distillerycats.drop_unique_hashes import build_hashCounter, drop_below_mincount


def test_build_hashCounter():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=True)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=True)
    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 2, 5))
    true_res=Counter({1: 2, 2: 2, 3: 1, 4: 1, 5: 1})

    ss1 = SourmashSignature(mh1)
    ss2 = SourmashSignature(mh2)

    counts = Counter()
    hc = build_hashCounter([ss1,ss2], counts)
    print("Hash Counter: ", hc)
    assert hc == true_res

def test_drop_below_mincount():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=True)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=True)
    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 2, 5))

    ss1 = SourmashSignature(mh1)
    ss2 = SourmashSignature(mh2)

    counts = Counter()
    hc = build_hashCounter([ss1,ss2], counts)
    kept_hashes = drop_below_mincount(hc, 2)
    true_kept = Counter({1: 2, 2: 2})
    print("kept hashes: ", kept_hashes)
    assert kept_hashes == true_kept

def test_drop_below_mincount_threshold():
    mh1 = MinHash(0, 21, scaled=1, track_abundance=True)
    mh2 = MinHash(0, 21, scaled=1, track_abundance=True)
    mh1.add_many((1, 2, 3, 4))
    mh2.add_many((1, 1, 2, 5))

    ss1 = SourmashSignature(mh1)
    ss2 = SourmashSignature(mh2)

    counts = Counter()
    hc = build_hashCounter([ss1,ss2], counts)
    kept_hashes = drop_below_mincount(hc, 3)
    true_kept = Counter({1: 3})
    print("kept hashes: ", kept_hashes)
    assert kept_hashes == true_kept
