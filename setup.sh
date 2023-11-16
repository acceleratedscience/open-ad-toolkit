#!/bin/bash
set -e

# helped functions
to_lowercase() {
    local input="$1"
    echo "$input" | tr '[:upper:]' '[:lower:]'
}

check_string_exist() {
    local str1="$1"
    local arr2=("${@:2}")

    for str2 in "${arr2[@]}"; do
        if [ "$str1" == "$str2" ]; then
            return 0  # Return true if the string exists in the second array
        fi
    done

    return 1  # Return false if the string does not exist in the second array
}

# start setup
clear

if command -v python3 &> /dev/null; then
    echo "Python detected, moving forward with setup"
else
    echo "Python not detected"
    echo "Please install python and rerun setup"
	exit 1
fi

package_name="poetry"
array2=("String1" "String2" "String3")
echo "Please choose your setup option: "
echo "1) Use virtual environment using Poetry"
echo "2) Use your global space"
read -p "# Input: " option

clear

if [ "$option" -gt 0 ] && [ "$option" -le 2 ]; then
	echo "You have chosen $option."
	read -p "Confirm (Y/N):" confirm
else
	echo "$option is not a valid option."
	exit 1
fi

clear


confirmlower=$(to_lowercase "$confirm")
possible_confirm=("y" "n")

if ! check_string_exist "$confirmlower" "${possible_confirm[@]}"; then
	echo "$confirmlower is not a valid option."
	exit 1
fi

if [ "$confirmlower" == "n" ]; then
	echo "Aborted"
	exit 0
fi

if [ "$option" -eq 1 ]; then
	if pip show "$package_name" &> /dev/null; then
		echo "$package_name detected."
	else
		read -p "$package_name is not installed, and this is needed to proceed. Should we install it? (Y/N):" poetry_confirm
		poetry_confirm_lower=$(to_lowercase "$poetry_confirm")

		if ! check_string_exist "$poetry_confirm_lower" "${possible_confirm[@]}"; then
			echo "$poetry_confirm_lower is not a valid option."
			exit 1
		fi

		if [ "$poetry_confirm_lower" == "y" ]; then
			pip install poetry
		else 
			echo "Aborted"
			exit 0
		fi
	fi
	poetry config virtualenvs.create true
	# install opentoolkit
	poetry install
	poetry run python -m ipykernel install --user --name=ad-kernel
	# Virtual ENV
	clear
	echo "Success!"
	echo "Enter the virtual environment by entering 'poetry shell' and start using Opentoolkey!"
fi

if [ "$option" -eq 2 ]; then
	# Global ENV
	pip install .
	python -m ipykernel install --user --name=ad-kernel
	clear
	echo "Success!"
	echo "Start using Opentoolkit!"
fi
