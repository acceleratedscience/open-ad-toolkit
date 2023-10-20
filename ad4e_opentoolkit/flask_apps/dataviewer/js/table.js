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

	// Homebrewed version of movableRows.
	// - - -
	// We had to build our own because the build-in movableRows
	// can't be made consitional, and it creates undesrieable
	// side effects in edit mode. Specifically, it will interfere
	// with the textarea resize handle.
	moveRowStart(e, row) {
		if (this.editMode) return

		// Store elements
		const $row = row.getElement()
		const $rowClone = $row.cloneNode(true)
		const $parent = this.element.querySelector('.tabulator-tableholder')

		// Calculate values
		const rowY = $row.getBoundingClientRect().top
		let parentY = $parent.getBoundingClientRect().top
		let y = rowY - parentY - 1
		const mouseRowOffset = e.clientY - rowY
		let jump = 0 // How many positions the row has jumped
		const table = this

		// Update DOM
		$row.classList.add('while-dragging')
		$rowClone.classList.remove('tabulator-selectable')
		$rowClone.classList.add('dragging')
		$rowClone.style.top = y + 'px'
		$parent.appendChild($rowClone)

		// Make interactive
		document.addEventListener('mousemove', _moveRowDrag)
		document.addEventListener('mouseup', _moveRowEnd)

		// // For testing
		// const $indicator = document.createElement('div')
		// $indicator.style.height = '1px'
		// $indicator.style.width = '100px'
		// $indicator.style.background = 'red'
		// $indicator.style.position = 'absolute'
		// $indicator.style.zIndex = 1000
		// $parent.appendChild($indicator)
		// const $reference1 = document.createElement('div')
		// $reference1.style.height = '1px'
		// $reference1.style.width = '100px'
		// $reference1.style.background = 'blue'
		// $reference1.style.position = 'absolute'
		// $reference1.style.left = 0
		// $reference1.style.zIndex = 1000
		// $parent.appendChild($reference1)
		// const $reference2 = document.createElement('div')
		// $reference2.style.height = '1px'
		// $reference2.style.width = '100px'
		// $reference2.style.background = 'green'
		// $reference2.style.position = 'absolute'
		// $reference2.style.left = 0
		// $reference2.style.zIndex = 1000
		// $parent.appendChild($reference2)

		// On drag listener
		function _moveRowDrag(e) {
			parentY = $parent.getBoundingClientRect().top
			y = e.clientY - parentY - mouseRowOffset
			const minY = -1
			const maxY = $parent.clientHeight - $rowClone.clientHeight - 2
			y = Math.max(Math.min(y, maxY), minY)
			$rowClone.style.top = y + 'px'
			_moveRow(y)
		}

		// On drag - move row to closest position
		function _moveRow(y) {
			const yRowMiddle = y + $rowClone.clientHeight / 2
			const moveRow = _findTargetIndex(yRowMiddle)
			jump += moveRow

			if (moveRow == 1) {
				const $nextRow = $row.nextElementSibling
				$nextRow.after($row)
			} else if (moveRow == -1) {
				const $prevRow = $row.previousElementSibling
				$prevRow.before($row)
			}
		}

		function _findTargetIndex(yRowMiddle) {
			const targetRow = document.querySelector('.tabulator-row.while-dragging')
			const yTargetRowTop = targetRow.offsetTop
			const yTargetRowBottom = targetRow.offsetTop + targetRow.clientHeight

			// // For testing
			// $indicator.style.top = yRowMiddle + 'px'
			// $reference1.style.top = yTargetRowTop + 'px'
			// $reference2.style.top = yTargetRowBottom + 'px'

			if (yRowMiddle > yTargetRowBottom) {
				return 1
			} else if (yRowMiddle < yTargetRowTop) {
				return -1
			} else {
				return 0
			}
		}

		// On drop
		function _moveRowEnd() {
			document.removeEventListener('mousemove', _moveRowDrag)
			document.removeEventListener('mouseup', _moveRowEnd)
			$rowClone.remove()
			$row.classList.remove('while-dragging')

			const data = table.getData()
			// console.log('')
			// console.log(data.map(itm => itm.Index))
			const currentPos = row.getPosition() - 1
			const newPos = currentPos + jump
			data.splice(newPos, 0, data.splice(currentPos, 1)[0])
			table.setData(data)
			// console.log(currentPos, '-->', newPos)
			// console.log(data.map(itm => itm.Index))

			// // For testing
			// $indicator.remove()
			// $reference1.remove()
		}
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
