#!/bin/bash
# Work with SINE_orth_loc alignment bundles (aln_<sp1>-<sp2>_<PM|MP|SINE|rejected>.aln[.gz]).
# A bundle is plain text: "##FILE <name>" followed by that alignment in FASTA,
# repeated for every alignment; ".gz" bundles are read transparently.

usage() {
    cat << EOF
Usage:
  $0 list    BUNDLE...                  names of alignments in the bundles
  $0 count   BUNDLE...                  number of alignments per bundle
  $0 get     NAME BUNDLE...             print one alignment (FASTA) to stdout
  $0 extract [-p REGEX] DIR BUNDLE...   write alignments (names matching REGEX) as files into DIR
  $0 pack    [--remove] DIR PREFIX      bundle loose .PM/.MP/.SINE files of an older run in DIR
                                        into PREFIX_PM.aln.gz, PREFIX_MP.aln.gz, PREFIX_SINE.aln.gz
                                        (--remove deletes the packed files afterwards)
EOF
    exit 1
}

[[ $# -lt 2 ]] && usage
cmd=$1; shift

case $cmd in
    list)
        for b in "$@"; do zcat -f "$b" | awk '/^##FILE /{print substr($0,8)}'; done
        ;;
    count)
        for b in "$@"; do printf '%s\t%s\n' "$(zcat -f "$b" | grep -c '^##FILE ')" "$b"; done
        ;;
    get)
        name=$1; shift
        [[ $# -lt 1 ]] && usage
        for b in "$@"; do zcat -f "$b" | awk -v n="$name" '/^##FILE /{p=(substr($0,8)==n); next} p'; done
        ;;
    extract)
        regex=""
        if [[ $1 == "-p" ]]; then regex=$2; shift 2; fi
        [[ $# -lt 2 ]] && usage
        dir=$1; shift
        mkdir -p "$dir" || exit 1
        for b in "$@"; do
            # regex passed via the environment: awk -v would interpret its backslashes
            zcat -f "$b" | RE="$regex" awk -v d="$dir" -v src="$b" 'BEGIN {re=ENVIRON["RE"]}
                /^##FILE / {if (o) close(o); n=substr($0,8); o=""
                            if (n ~ /\// || n == "" || (re != "" && n !~ re)) next
                            o=d"/"n; printf "" > o; k++; next}
                o {print > o}
                END {print k+0" alignments extracted from "src > "/dev/stderr"}'
        done
        ;;
    pack)
        remove=false
        if [[ $1 == "--remove" ]]; then remove=true; shift; fi
        [[ $# -ne 2 ]] && usage
        dir=$1; prefix=$2
        for t in PM MP SINE; do
            out="${prefix}_$t.aln"
            echo "##SINE_orth_loc bundle v1" > "$out" || exit 1
            find "$dir" -maxdepth 1 -type f -name "*.$t" -print0 | sort -z |
             xargs -0 -r awk 'FNR==1 {n=FILENAME; sub(/.*\//,"",n); print "##FILE " n} {print}' >> "$out" || exit 1
            gzip -f "$out" || exit 1
            k=$(zcat "$out.gz" | grep -c '^##FILE ')
            echo "$k alignments -> $out.gz"
            if $remove; then
                # delete only files whose names are recorded in the bundle just written
                zcat "$out.gz" | awk -v d="$dir" '/^##FILE /{printf "%s/%s\0", d, substr($0,8)}' | xargs -0 -r rm -f
            fi
        done
        ;;
    *)
        usage
        ;;
esac
