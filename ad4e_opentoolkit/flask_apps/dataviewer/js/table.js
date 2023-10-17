// This is a custom module that extends
// Tabulator with some extra functionality.
class Table extends Tabulator {
	constructor(element, options) {
		// Modifications
		super(element, options)
		super.on('tableBuilt', () => {
			this.addIndexCol()
		})
		super.on('renderComplete', () => {
			this.storeColWidths()
		})
		super.on('columnResized', () => {
			this.element.classList.add('resized')
		})

		// Variables
		this.editMode = false
		this.colDefaultWidths = {}
		this.addedIndex = false

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

	// Store column default widths
	storeColWidths() {
		this.getColumns().forEach(col => {
			this.colDefaultWidths[col.getField()] = col.getWidth()
		})
	}

	// Add index column
	addIndexCol(force) {
		if (this.addedIndex) return
		const data = this.getData()

		let addCol = false
		if (force) {
			// Column is added manually via options panel.
			addCol = true
		} else if (!force) {
			// Column is added on init, but only if one doesn't already exist.
			const rowSample = data[1]
			let hasIndex = !!rowSample['#']
			if (!hasIndex) {
				for (const key in rowSample) {
					if (key.toLowerCase() == 'index') {
						hasIndex = true
						break
					}
				}
				addCol = !hasIndex
			}
		}

		// Add index column if missing
		if (addCol) {
			this.addedIndex = true
			data.forEach((row, i) => {
				row['#'] = i + 1
			})
			this.setData(data)

			const colSample = this.getColumns()[0]
			const newCol = {
				...colSample,
				sorter: 'number',
				title: '#',
				field: '#',
			}
			this.addColumn(newCol, true)
		}
	}

	// Remove index column
	removeIndexCol() {
		const indexCol = this.getColumns()[0]
		if (indexCol.getField() == '#') {
			this.addedIndex = false
			indexCol.delete()
			const data = this.getData()
			data.forEach((row, i) => {
				delete row['#']
			})
			this.setData(data)
		}
	}

	// Reset column widths to default
	resetCols() {
		this.getColumns().forEach(col => {
			col.setWidth(this.colDefaultWidths[col.getField()])
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
