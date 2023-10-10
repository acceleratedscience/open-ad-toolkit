// This is a custom module that extends
// Tabulator with some extra functionality.
class Table extends Tabulator {
	constructor(element, options) {
		// Modifications
		super(element, options)
		super.on('renderComplete', () => {
			this.storeCols()
		})
		super.on('columnResized', () => {
			this.element.classList.add('resized')
		})

		// Variables
		this.editMode = false

		// Parse custom options.
		this.init(options)

		// super.on('cellEditing', cell => {
		// 	console.log(1)
		// 	return false
		// })
	}

	init(options) {
		if (options.editMode) {
			this.toggleEditMode(true)
		}
	}

	// on(eventName, callback) {
	// 	console.log(333, eventName)
	// 	if (eventName != 'cellEditing') {
	// 		super.on(eventName, callback)
	// 	} else {
	// 		console.log(1234)
	// 	}
	// }

	// isEditMode(cell) {
	// 	console.log(111)
	// 	// cell = { ...cell }
	// 	// console.log(cell._cell.row.position, cell._cell.column)
	// 	// console.log(cell)
	// 	return false
	// }

	// Store column default widths.
	storeCols() {
		this.defaultColumnWidths = {}
		this.getColumns().forEach(col => {
			this.defaultColumnWidths[col.getField()] = col.getWidth()
		})
	}

	// Reset column widths to default.
	resetCols() {
		this.getColumns().forEach(col => {
			col.setWidth(this.defaultColumnWidths[col.getField()])
			this.element.classList.remove('resized')
		})
	}

	// Toggle edit mode.
	toggleEditMode(bool) {
		this.editMode = bool == undefined ? !this.editMode : bool
		if (this.editMode) {
			this.element.classList.add('edit-mode')
			history.pushState('', document.title, window.location.pathname + window.location.search + '#edit') // Add hash
		} else {
			this.element.classList.remove('edit-mode')
			history.pushState('', document.title, window.location.pathname + window.location.search) // Remove hash
		}
	}
}

Table.moduleName = 'custom'
