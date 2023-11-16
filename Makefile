.PHONY: check-lint test setup lint

setup:
	poetry install
	poetry run python -m ipykernel install --user --name=ad-kernel

check-lint:
	poetry run black --check .

lint:
	poetry run black .

test:
	poetry run coverage run --branch --source=./openad_opentoolkit/app/ -m pytest --durations=10 --color=yes tests/unit
	poetry run coverage report
