TEST = null
/**
 * Table is a subclass that extends
 * Tabulator with some extra functionality.
 */

class Table extends Tabulator {
	constructor(element, options) {
		// Initiation
		super(element, options)
		super.on('tableBuilt', () => {
			// Add index column is missing
			this.addIndexCol()

			// Turn on edit mode if hash is present
			if (options.editMode) {
				console.log('$')
				this.toggleEditMode(true)
			}

			// Store column widths so we can reset them
			this.storeColWidths()
		})

		// Modofation
		super.on('columnResized', () => {
			this.element.classList.add('resized-col')
		})
		super.on('rowResized', () => {
			this.element.classList.add('resized-row')
		})

		// Variables
		this.editMode = false // Used to determine if we're in edit mode - see isEditMode()
		this.colDefaultWidths = {} // Used to store column widths so we can reset them
		// this.rowDefaultHeights = {} // Used to store row heights so we can reset them
		this.addedIndex = false // Used to keep track if we added an index column
		this.lastSelectedRowIndex = null // Number used to determine where to start from when shift-selecting
		this.lastSelectedRowSelState = null // Boolean used to determine if shift-select row should batch-select or batch-deselect
		this.selectMode = false // Boolean used to ignore row resize handles while selecting rows
	}

	// Store column default widths
	storeColWidths() {
		this.getColumns().forEach(col => {
			this.colDefaultWidths[col.getField()] = col.getWidth()
		})
	}

	// // Store row default heights
	// storeRowHeights() {
	// 	this.getRows().forEach(row => {
	// 		this.rowDefaultHeights[row.getPosition()] = row.getHeight()
	// 	})
	// }

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
			this.element.classList.remove('resized-col')
		})
	}

	// // Reset row heights to default
	// resetRows() {
	// 	this.getRows().forEach(row => {
	// 		row.setHeight(this.rowDefaultHeights[row.getField()])
	// 		this.element.classList.remove('resized-row')
	// 	})
	// }

	// Check if table is in edit mode
	isEditMode() {
		return this.editMode
	}

	// Toggle edit mode
	// - - -
	// Note: Tabulator doesn't support toggling edit mode, so we use a little hack.
	// The actual toggling of edit mode is done via isEditMode when creating the table.
	toggleEditMode(bool, revertChanges) {
		console.log('>', bool)
		this.editMode = bool == undefined ? !this.editMode : bool
		if (this.editMode) {
			this.storeData() // Store data so we can revert on cancel
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
	}

	// Revert changes on cancel
	async revertData() {
		if (this.hasEdits(this.dataBeforeEdit, this.getData())) {
			const selectedRows = await this.getSelectedRows()
			const selectedRowsIndexes = selectedRows.map(row => row.getPosition())
			this.setData(this.dataBeforeEdit)
			this.redraw()

			// Re-select rows that were selected before
			const newRows = table.getRows()
			const newSelectedRows = selectedRowsIndexes.map(index => newRows[index - 1])
			this.selectRow(newSelectedRows)
		}
	}

	// Check if table data was edited
	hasEdits(data1, data2) {
		if (data1.length != data2.length) return true
		for (let i = 0; i < data1.length; i++) {
			const row1 = data1[i]
			const row2 = data2[i]
			for (const key in row1) {
				if (row1[key] != row2[key]) return true
			}
		}
		return false
	}

	// Recreate the built in sort
	sort(col) {
		if (this.editMode) {
			alert('Please exit edit mode before sorting.')
			return
		}
		const colName = col.getField()
		const $col = col.getElement()
		const sortOrder = $col.getAttribute('aria-sort') == 'ascending' ? 'desc' : 'asc'
		$col.classList.add('sortable', `sort-${sortOrder}`)
		// tabulator-col-sorter
		this.setSort(colName, sortOrder)
	}

	// // Tag on to the Tabulator redraw function
	// redraw(bool) {
	// 	console.log('redrawww')
	// 	// Reset the edit parameters to their default values
	// 	this.lastSelectedRowIndex = null
	// 	this.lastSelectedRowSelState = null
	// 	this.selectMode = false
	// 	super.redraw(bool)
	// }
}

Table.moduleName = 'custom'
