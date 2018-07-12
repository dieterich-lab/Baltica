ORG=$1
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D $ORG -P 3306   -e "select chrom,txStart,txEnd,K.name,X.geneSymbol,strand,exonStarts,exonEnds from knownGene as K,kgXref as X where  X.kgId=K.name;" > tmp.notbed

grep -v txStart tmp.notbed | awk '
        BEGIN { OFS = "\t"; FS = "\t"} ;
            {
                delete astarts;
                delete aends;
                split($7, astarts, /,/);
                split($8, aends, /,/);
                starts=""
                sizes=""
                exonCount=0
                for(i=1; i <= length(astarts); i++){
                    if (! astarts[i]) continue
                    sizes=sizes""(aends[i] - astarts[i])","
                    starts=starts""(astarts[i] = astarts[i] - $2)","
                    exonCount=exonCount + 1
                }
                print $1,$2,$3,$5","$4,1,$6,$2,$3,"0",exonCount,sizes,starts
            }' | sort -k1,1 -k2,2n > knownGene.$ORG.bed12
