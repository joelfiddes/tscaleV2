#!/bin/bash
# ./commit.sh "my commit message"

if [ $# -eq 0 ]
  then
    echo "No commit message given"; exit
fi

# need to be in venv to make sure package list not entire system
source ./venv/bin/activate

./update_requirements.sh
echo "updated requirements"
git add -A
git commit -m $1
git push

