.PHONY: check-lint test setupt setup lint

setup_file_path := $(PWD)/setup.sh

setupt:
	@if ! [ -x "$(setup_file_path)" ]; then \
		chmod +x ./setup.sh; \
	fi
	@./setup.sh

setup:
	poetry install
	poetry run python -m ipykernel install --user --name=ad-kernel

check-lint:
	poetry run black --check .

lint:
	poetry run black .

test:
	poetry run coverage run --branch --source=./ad4e_opentoolkit/app/ -m pytest --durations=10 --color=yes tests/unit
	poetry run coverage report
