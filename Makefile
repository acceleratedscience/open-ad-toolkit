.PHONY: check-lint test setup

setup:
	poetry install
	poetry run python -m ipykernel install --user --name=ad-kernel

check-lint:
	poetry run black --check .

test:
	poetry run coverage run --branch --source=./ad4e_opentoolkit/app/ -m pytest --durations=10 --color=yes tests/unit
	poetry run coverage report
