./update_requiremnts.sh
echo "updated requirements"
git add -A
echo "added to local repo"
git commit -m $1
echo "commit with message" $1
git push
echo "pushed to github"
