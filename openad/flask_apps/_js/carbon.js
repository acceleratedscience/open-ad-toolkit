class Carbon {
	constructor() {
		this.init()
	}

	init() {
		this.initDropdowns()
	}

	/**
	 * Dropdowns
	 */

	initDropdowns() {
		document.querySelectorAll('.ibm-dropdown').forEach($dropdown => {
			const $select = $dropdown.querySelector('select')
			$select.removeEventListener('change', this.onDropdownChange)
			$select.addEventListener('change', this.onDropdownChange)
			this.onDropdownChange({ target: $select })
		})
	}

	// Update display when dropdown value changes.
	onDropdownChange(e) {
		const $display = e.target.parentNode.querySelector('.ibm-dd-display')
		if (!$display) return
		const val = e.target.value
		const displayVal = e.target.querySelector(`option[value="${val}"]`).innerText
		$display.innerText = displayVal
	}

	// Programatically update a dropdown value.
	updateDropdown(querySelector, val) {
		const $dropdown = document.querySelector(querySelector)

		// Check if the dropdown exists.
		if (!$dropdown) {
			console.error(`Invalid dropdown selector '${querySelector}'`)
			return
		}

		// Check if the dropdoen has the value.
		const $option = $dropdown.querySelector(`option[value="${val}"]`)
		if (!$option) {
			console.error(`Invalid dropdown value '${val}'`)
			return
		}

		const event = new Event('change', { bubbles: true, cancelable: true })
		$dropdown.value = val
		$dropdown.dispatchEvent(event)
	}
}

let carbonUi
document.addEventListener('DOMContentLoaded', () => {
	carbonUi = new Carbon()
})
