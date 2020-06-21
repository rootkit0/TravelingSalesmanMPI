#!/bin/bash

for i in script1.sh script2.sh script3.sh script4.sh script5.sh
do
	qsub $i
done