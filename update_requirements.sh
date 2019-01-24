#!/bin/bash

# update requirements for packaging
pip freeze > requirements1.txt

# find erroneous "pkg-resources==0.0.0" and remove - Linux bug
grep -vwE "(pkg-resources==0.0.0)" requirements1.txt > requirements.txt
rm requirements1.txt
