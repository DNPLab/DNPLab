# clean up imports
isort -rc .

# clean up formats (must run after isort)
black .

# check syntax errors
git diff upstream/develop -u -- "*.py" | flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics


# pytest
pytest .
