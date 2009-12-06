#!/usr/bin/awk -f
#strips off extended version number after dash from qName

function stripVersion(uacc) {
    split(uacc, a, "-");
    return a[1];
}

BEGIN {
    FS = "\t";
    OFS = "\t";
    prevqs = 0; prevqe = 0; prevts = 0 ; prevte = 0;
    prevkey = "";
}

{$10 = stripVersion($10); print $0}
