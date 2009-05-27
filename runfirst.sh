#!/bin/bash

echo "This script will create a directory called 'dune' containing "
echo "symbolic links to include directories of the modules "
echo "dune-common, dune-grid and dune-istl that are required by"
echo "dune-grid-glue."

if [ -f "dune" ]
then
	echo "Error! Cannot create directory 'dune' because of file of the same name."
	exit 1
elif [ ! -d "dune" ]
then
	echo "Creating directory 'dune'..."
	mkdir dune
fi

cd dune
for i in common grid istl
do
	if [ -d "../../dune-${i}/${i}" ] && [ ! -e "${i}" ]
	then
		echo "Creating link for '${i}'..."
		ln -s ../../dune-${i}/${i} ${i}
	else
		echo "Cannot create symbolic link '${i}'. Check includes manually!"
	fi
done
cd ..
