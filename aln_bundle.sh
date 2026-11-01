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
                                        into PREFIX_PM.aln.gz, PREFIX_MP.aln.gz, PREFIX_SINE.aln.gz,
                                        and its rejected ones (*.cl.*.badLF, badRF, shortRF, lfSINE,
                                        rfSINE, tiePM, contig_end, contig_end_lf, overlong) into PREFIX_rejected.aln.gz
                                        (--remove deletes the packed files afterwards)
EOF
    exit 1
}

[[ $# -lt 2 ]] && usage
cmd=$1; shift
set -o pipefail
# every bundle must be a readable file (a wrong path must not look like an empty bundle)
check_bundles() { local b; for b in "$@"; do [[ -f "$b" && -r "$b" ]] || { echo "$0: bundle not found or not readable: $b" >&2; exit 1; }; done; }

case $cmd in
    list)
        check_bundles "$@"
        for b in "$@"; do zcat -f "$b" | awk '/^##FILE /{print substr($0,8)}' || { echo "$0: cannot read $b" >&2; exit 1; }; done
        ;;
    count)
        check_bundles "$@"
        for b in "$@"; do
            k=$(zcat -f "$b" | awk '/^##FILE /{k++} END {print k+0}') || { echo "$0: cannot read $b" >&2; exit 1; }
            printf '%s\t%s\n' "$k" "$b"
        done
        ;;
    get)
        name=$1; shift
        [[ $# -lt 1 ]] && usage
        check_bundles "$@"
        found=0
        for b in "$@"; do
            out=$(zcat -f "$b" | awk -v n="$name" '/^##FILE /{p=(substr($0,8)==n); if (p) f=1; next} p; END {exit !f}')
            rc=$?; [[ -n "$out" ]] && printf '%s\n' "$out"
            if (( rc == 0 )); then found=1; elif ! zcat -f "$b" > /dev/null; then echo "$0: cannot read $b" >&2; exit 1; fi
        done
        (( found )) || { echo "$0: $name not in the given bundles" >&2; exit 1; }
        ;;
    extract)
        regex=""
        if [[ $1 == "-p" ]]; then regex=$2; shift 2; fi
        [[ $# -lt 2 ]] && usage
        dir=$1; shift
        check_bundles "$@"
        mkdir -p "$dir" || exit 1
        # names already written by this call (a name met again, in the same or a later bundle, e.g. a run's bundle
        # and its by_reference copy, or C<n>R ids of two runs, is written as <name>~2, <name>~3, ... not over it)
        seen=$(mktemp) || exit 1
        for b in "$@"; do
            # regex passed via the environment: awk -v would interpret its backslashes
            zcat -f "$b" | RE="$regex" SEEN="$seen" awk -v d="$dir" -v src="$b" 'BEGIN {re=ENVIRON["RE"]; S=ENVIRON["SEEN"]
                    while ((getline l < S) > 0) seen[l]++; close(S)}
                /^##FILE / {if (o) close(o); n=substr($0,8); o=""
                            if (n ~ /\// || n == "" || (re != "" && n !~ re)) next
                            f=n; if (++seen[n] > 1) {f=n "~" seen[n]
                             print "duplicate name " n " (occurrence " seen[n] ") in " src ": written as " f > "/dev/stderr"}
                            print n >> S; o=d"/"f; printf "" > o; k++; next}
                o {print > o}
                END {print k+0" alignments extracted from "src > "/dev/stderr"}'
        done
        rm -f "$seen"
        ;;
    pack)
        remove=false
        if [[ $1 == "--remove" ]]; then remove=true; shift; fi
        [[ $# -ne 2 ]] && usage
        dir=$1; prefix=$2
        for t in PM MP SINE rejected; do
            out="${prefix}_$t.aln"
            echo "##SINE_orth_loc bundle v1" > "$out" || exit 1
            names=(-name "*.$t")
            [[ $t == rejected ]] && names=(-name "*.cl.*.badLF" -o -name "*.cl.*.badRF" -o -name "*.cl.*.shortRF" -o -name "*.cl.*.lfSINE"
                -o -name "*.cl.*.rfSINE" -o -name "*.cl.*.tiePM" -o -name "*.cl.*.contig_end" -o -name "*.cl.*.contig_end_lf" -o -name "*.cl.*.overlong")
            find "$dir" -maxdepth 1 -type f \( "${names[@]}" \) -print0 | sort -z |
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
