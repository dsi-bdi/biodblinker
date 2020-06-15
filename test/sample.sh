#!/bin/bash

IFS=$'\n'
lines=$(gunzip -c $1 | shuf -n $2)
#readarray -d '\n' -t lines <<< $(gunzip -c $1 | shuf -n $2)
for l in $lines
do
	delim=$'\t'
	string="${l%%*( )}"$delim
	parts=()
	while [[ $string ]]; do
		parts+=( "${string%%"$delim"*}" )
  		string=${string#*"$delim"}
	done
	id="'${parts[0]%%*()}'"
	delim=$';'
	string="${parts[1]%%*( )}"$delim
	
	mapped=()
	while [[ $string ]]; do
		mapped+=( "${string%%"$delim"*}" )
  		string=${string#*"$delim"}
	done
	
	values="'${mapped[0]%%*( )}'"
	len=${#mapped[@]}

	for (( i=1; i<${len}; i++))
	do
		values+=", '${mapped[$i]%%*( )}'"
	done

	echo "([${id%%*( )}], [[${values%%*( )}]]),"
done

#echo $lines