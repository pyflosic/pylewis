cd ..
pylewis="$(dirname "$(pwd)")"
cd C60
export PYTHONPATH=$pylewis/src/:$PYTHONPATH

python $1
