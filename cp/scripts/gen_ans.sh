cp correct tests/correct
cd tests
for (( t = 1; t <= $1; t++ ));
do
	a=`cat $t.in`
	if [[ ("$a" = "") ]];
	then
		continue
	fi
	./correct < "$t.in" > "$t.out"
    # echo "$t " `./solution < "$t.in"`
done
