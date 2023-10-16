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
		this.addedIndexCol = false

		// Parse custom options
		this.init(options)
	}

	// Initialization
	init(options) {
		this.on('dataProcessed', () => {
			// Turn on edit mode if hash is present
			if (options.editMode) {
				this.toggleEditMode(true)
			}
		})
	}

	// Trash
	// on(eventName, callback) {
	// 	console.log(333, eventName)
	// 	if (eventName != 'cellEditing') {
	// 		super.on(eventName, callback)
	// 	} else {
	// 		console.log(1234)
	// 	}
	// }

	// Store column default widths
	storeCols() {
		this.defaultColumnWidths = {}
		this.getColumns().forEach(col => {
			this.defaultColumnWidths[col.getField()] = col.getWidth()
		})
	}

	// Reset column widths to default
	resetCols() {
		this.getColumns().forEach(col => {
			col.setWidth(this.defaultColumnWidths[col.getField()])
			this.element.classList.remove('resized')
		})
	}

	// Toggle edit mode
	toggleEditMode(bool, revertChanges) {
		this.editMode = bool == undefined ? !this.editMode : bool
		if (this.editMode) {
			this.storeData()
			this.element.classList.add('edit-mode')
			history.pushState('', document.title, window.location.pathname + window.location.search + '#edit') // Add hash
		} else {
			if (revertChanges) this.revertData()
			this.element.classList.remove('edit-mode')
			history.pushState('', document.title, window.location.pathname + window.location.search) // Remove hash
		}
	}

	// Store copy of data so we can revert changes
	storeData() {
		this.dataBeforeEdit = this.getData()
		console.log('Store:', this.dataBeforeEdit)
	}

	// Revert changes on cancel
	revertData() {
		this.setData(this.dataBeforeEdit)
		console.log('Revert:', this.data, this.dataBeforeEdit)
		this.redraw(true)
	}
}

Table.moduleName = 'custom'
