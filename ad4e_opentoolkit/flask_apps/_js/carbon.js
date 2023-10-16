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

	onDropdownChange(e) {
		const val = e.target.value
		const displayVal = e.target.querySelector(`option[value="${val}"]`).innerText
		const $display = e.target.parentNode.querySelector('.ibm-dd-display')
		$display.innerText = displayVal
	}
}

document.addEventListener('DOMContentLoaded', () => {
	new Carbon()
})
