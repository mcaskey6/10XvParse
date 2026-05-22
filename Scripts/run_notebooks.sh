#! /bin/bash
set -e

JUST_COMBO="true"

TECH_NOTEBOOKS=(
    "/home/mcaskey/10XvParse/Notebooks/Analysis_2/10x.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_2/parse.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_2/parse_mini.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_3/parse.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_4/10x.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_4/parse.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_5/parse.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_6/10x.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_6/parse.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_7/parse.ipynb"
)

COMBO_NOTEBOOKS=(
    "/home/mcaskey/10XvParse/Notebooks/Analysis_2/combo.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_2/combo_mini.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_3/combo_H1.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_3/combo_H2.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_4/combo.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_5/combo.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_6/combo.ipynb"
    "/home/mcaskey/10XvParse/Notebooks/Analysis_7/combo.ipynb"
)

if [ "$JUST_COMBO" = "true" ]; then
    NOTEBOOKS=("${COMBO_NOTEBOOKS[@]}")
else
    NOTEBOOKS=("${TECH_NOTEBOOKS[@]}" "${COMBO_NOTEBOOKS[@]}")
fi

for nb in "${NOTEBOOKS[@]}"; do
    echo "Running: $nb"
    jupyter nbconvert --to notebook --execute --inplace \
        --ExecutePreprocessor.timeout=3600 "$nb" 2> >(tr '\r' '\n' | grep '\[NbConvertApp\]' >&2)
    echo "Done: $nb"
done

echo "All notebooks completed."