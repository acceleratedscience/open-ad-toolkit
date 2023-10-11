.PHONY: check-lint test

check-lint:
	poetry run black --check .

test:
	poetry run coverage run --branch --source=./ad4e_opentoolkit/app/ -m pytest --durations=10 --color=yes tests/unit
