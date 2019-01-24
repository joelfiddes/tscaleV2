./update_requirements.sh
echo "updated requirements"
git add -A

git commit -m $1

git push

