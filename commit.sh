#!/bin/bash
# ./commit.sh "my commit message"

if [ $# -eq 0 ]
  then
    echo "No commit message given"; exit
fi

source ./venv/bin/activate
./update_requirements.sh
echo "updated requirements"
git add -A
git commit -m $1
git push

