cd ..
pylewis="$(dirname "$(pwd)")"
cd N2
export PYTHONPATH=$pylewis/src/:$PYTHONPATH

python $1
